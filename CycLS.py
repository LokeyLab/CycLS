"""A program for identifying a cyclic peptides from a known list of compounds by MSMS.
Not meant for identification of natural products, but useful to identify members of a library.
The program requires an MZML file from which to draw spectra, a string containing instructions to generate the library,
and a database containing the smiles strings of the amino acids used in the library configured in an N->C terminal manner
as L-Ala  N[C@@H](C)C(=O)O.

Output is given in full in "*_Out.xlsx" and summarized in "*_Results.xlsx", where "*" is "Sequencing" by default.
As the m/z precision of the MS1 and MS2 spectra will impact sequencing accuracy and processing time, they should be optimized for your mass spectrometry system.
If the MS2 precision is too wide, then incorrect peptides will be rated more highly; too stringent and key evidence will be ignored.
Furthermore, not all relevant neutral loss types may be implemented if the residues used differ significantly from those included in the default residue definition file (or through oversight); 
this may be especially true for mass spectrometry systems not using ESI ionization and CID fragmentation for tandem mass spectrometry.
"""
#Copywrite Chad Townsend, 2018
#Reachable at cetownse@ucsc.edu

import collections, math, itertools, sys, argparse, multiprocessing, time, gzip, signal, copy, operator, bisect, random, statistics, openpyxl, os, random
import pymzml
import numpy as np
import statsmodels.api as sm
import matplotlib.pyplot as plt
import seaborn
from rdkit import Chem
from rdkit.Chem import Crippen

def openfiles (filename):
    """
    Takes a file name and opens it for reading, choosing
    the appropriate method depending on whether it was
    gzipped or not. Will error for files not found.
    Only looks in the current directory.
    """
    try:
        if filename.endswith('.gz'):
            filename = gzip.GzipFile(filename, 'r')
        else:
            filename = open(filename, 'r')
    except IOError:
        print ("Error: File \'{}\' does not exist.".format(filename), file=sys.stderr)
        sys.exit(1)
    return filename
    
def parse_args ():
    """Takes arguments from the command line and checks their validity. 
    Handles help documentation.
    """
    parser = argparse.ArgumentParser (description=__doc__)
    incompatible = parser.add_mutually_exclusive_group()
    parser.add_argument ('targets',
            help='Input the exact mass of the parent ion ([M+H]) of interest as the target(s) for analysis. \
            Comma separated values will be treated as a list of targets. \
            An asterisk will ask the program to search for targets itself. \
            Only one scan range will be accepted in such cases and will serve to limit the search space to that set of scans. \
            Ex: \"647.523,831.978,745.069\" would search for those three masses alone.')
    parser.add_argument ('mzml',
            help='The mzml file from which spectra are to be drawn.')
    parser.add_argument ('constraint',
            help='The constraint string governing library generation. Building blocks at a single location are separated by commas. \
            Positions are separated by semicolons. There should be no whitespace. \
            Ex: \"L,A,D;E,Q,K;P,G,R\" is a tripeptide with three possibilities for building blocks at each position.')
    incompatible.add_argument ('-r', '--rules', dest='rules', 
            help='Input an additional string with position-independent constraints on library composition. \
            Rules are input as a single building block, a comparison operator (=,<,>,<=,>=),\
            and an integer value. Rules for AlogP and mass can be accessed through the names \'AlogP\' and \
            \'MolWeight\' respectively and can be compared to float values. Multiple rules should be separated \
            by a semicolon. Because there will be spaces between characters, the set of rules should be \
            encompassed by quotes. Ex: "G < 3;R > 2;AlogP > 2.27" will generate a library from the \
            constraint string and filter down to members with less than 3 glycines, more than 2 arginines \
            and with AlogP greater than 2.27.')
    parser.add_argument ('-s', '--scanranges', dest='scanranges', default=None,
            help='Sets the minimum and maximum scan id\'s to be considered. Ex: 175-354.\
            Comma separated values will be treated as a list of scan ranges.\
            Scan ranges will be distributed to targets in the order they are entered and \
            there must be an equal number of targets and scan ranges.')
    parser.add_argument ('-o', '--out', dest='outfile', default='Sequencing',
            help='Sets the prefix of the output file. Defaults to \'Sequencing\', resulting in the files \'Sequencing_Results.xlsx\' and \'Sequencing_Out.xlsx\'.')
    parser.add_argument ('-d', '--database', dest='aadatabase', type=openfiles,
            default='aadatabase.txt', help='Sets the name of the amino acid name to smiles string \
            database to be accessed for library generation. Defaults to \'aadatabase.txt\'.')
    parser.add_argument ('-l', '--linear', action='store_true',dest='linear', 
            help='The library is linear. Otherwise assumed to be cyclic.')
    parser.add_argument ('-v', '--verbose',dest='verbose',default=0,type=int,
            help='Verbosity level. Prints out general status anouncments at 1. Higher levels give more detail.')
    parser.add_argument ('-q', '--query', action='store_true',dest='query', 
            help='Activate query mode, in which an interactive state occurs after normal program operation. Compound names entered will \
            yield fragment names, masses, and which spectra they hit in.')
    parser.add_argument ('-u', dest='workernum', type=int,
            help='Sets the number of worker processes. Defaults to the number of CPU\'s minus one.')
    parser.add_argument ('-n', dest='noiselevel', type=float, default=100.0,
            help='Sets the threshold intensity below which peaks are thrown out if above 1.0. If below 1.0, acts as the \
            maximum probability of a peak which survives filtration being due to noise. 0.1-0.5 tends to give good results.')
    parser.add_argument ('-p', dest='precision', default=None,
            help='Precision of MS2 and MS1 spectra m/z values respectively; this will vary with instrument and protocol. Defaults to \'0.3,.02\'.')
    incompatible.add_argument ('-t', '--truncate', action='store_true',dest='trunc', 
            help='Searches for truncations of the library in addition to full length library members. May significantly increase run time for large libraries. Incompatible with rule usage.')
    parser.add_argument('-e', dest='evalues', action='store_true',
            help='Given that the search database for a spectrum is of sufficient size, attempts to calculate E-values for the scores of each candidate molecule. Decoy database sizes under 1000 do not generate useful E-values. Experimental feature.')
    options = parser.parse_args()
    #-n
    if not options.workernum:
        try:
            options.workernum = (multiprocessing.cpu_count()-1)
        except NotImplementedError:
            print("Warning: Could not detect number of CPU's. Setting the number of worker processes to 1.")
            options.workernum = 1
    #-p
    if not options.precision:
        options.precision = [0.3,0.02]
    else:
        options.precision = options.precision.split(',')
        options.precision = [float(x) for x in options.precision]
        if len(options.precision) != 2:
            print('Error: -p argument formatted incorrectly.')
            parser.print_help()
            sys.exit(1)
    #Targets and -s
    #Target search mode
    if options.targets == '*':
        #Search scan limit specified
        if options.scanranges:
            options.scanranges = options.scanranges.split(',')
            if len(options.scanranges) != 1:
                print('Error: -s argument contains more than 1 scan range when target search mode is engaged.')
                parser.print_help()
                sys.exit(1)
            else:
                options.scanranges = options.scanranges[0].split('-')
                try:
                    options.scanranges = [int(thing) for thing in options.scanranges]
                except ValueError:
                    print('Error: -s argument formatted incorrectly.')
                    parser.print_help()
                    sys.exit(1)
                if len(options.scanranges) != 2:
                    print('Error: -s argument formatted incorrectly.')
                    parser.print_help()
                    sys.exit(1)
        else:
            options.scanranges = '*'
    #Targets specified
    else:
        options.targets = options.targets.split(',')
        try:
            for i,target in zip(itertools.count(),options.targets):
                options.targets[i] = float(target)
        except ValueError:
            print('Error: Targets formatted incorrectly.')
            parser.print_help()
            sys.exit(1)
        #Scan ranges specified
        if options.scanranges:
            options.scanranges = options.scanranges.split(',')
            if len(options.scanranges) != len(options.targets):
                print('Error: The number of scan ranges do not match the number of targets.')
                parser.print_help()
                sys.exit(1)
            if len(options.scanranges) == 1:
                options.scanranges = [options.scanranges]
            try:
                for i,scanrange in zip(itertools.count(),options.scanranges):
                    options.scanranges[i] = scanrange[0].split('-')
                    options.scanranges[i] = [int(thing) for thing in options.scanranges[i]]
            except ValueError:
                print('Error: -s argument formatted incorrectly.')
                parser.print_help()
                sys.exit(1)
            for range in options.scanranges:
                if len(range) != 2:
                    print('Error: -s argument formatted incorrectly.')
                    parser.print_help()
                    sys.exit(1)
        #No scan ranges specified
        else:
            options.scanranges=['*']*len(options.targets)
    return options

def getfrags (name, aatomass, linear):
    """Create all possible fragments of the input molecule
    and organize them by name and mass. Ions are +1 with H+.
    """
    byname = {}#name:getmass(name,aatomass,linear)+1.0078250321 #+H
    positionlist = name.split(',')
    l = len(positionlist)
    if linear:
        for i in range(1,l):
            if len(positionlist[:i]) > 1:
                beta = ','.join(positionlist[:i])
                byname[beta] = getmass(beta,aatomass,True)-15.9949146221-1.0078250321 #-OH
            if len(positionlist[i:]) > 1:
                gamma = ','.join(positionlist[i:])
                byname[gamma+'_y'] = getmass(gamma,aatomass,True)+1.0078250321 #+H
        #Include full length gamma fragment here for neutral loss. Otherwise we miss (full length gamma C-terminal loss) fragments.
        byname[name+'_y'] = getmass(name,aatomass,True)+1.0078250321 #+H
    else:
        loop = positionlist*2
        alreadydid = set()
        for fragsize in range(2, l):
            for subfragstart in range(l):
                subfrag = loop[subfragstart:subfragstart+fragsize]
                flat = ','.join(subfrag) #could be any unique key generated
                if flat not in alreadydid:
                    byname[flat] = getmass(flat, aatomass,False)+1.0078250321 #H+ ions
                    alreadydid.add(flat)
    return byname

def getmass (name, aatomass, linear):
    """Derive the mass of the input molecule from its name and separate residue masses.
    """
    splitname = name.split(',')
    namelen = len(splitname)
    amidebonds = namelen-1 if linear else namelen
    return sum([aatomass[aa] for aa in splitname])-amidebonds*18.0105646863
        
def merge_dicts (*dict_args):
    """merges any number of dictionaries. Where two or more dictionaries 
    contain the same key, dictionaries later in the arguments overwrite 
    values in dictionaries earlier in the arguments.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result

def flatten (x):
    """Flattens unevenly layered shallow lists. Does not conserve order.
    """
    result = []
    for el in x:
        if hasattr(el, "__getitem__") and not isinstance(el, str):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result

def uniquemass (positionlist, aatomass, linear):
    """generates the smallest set of unique masses by locating and combining positions with identical composition, then pruning duplicates that remain.
    """
    amidebonds = len(positionlist)-1 if linear else len(positionlist)
    #Convert amino acid names from each position in positionlist to masses, then sort them.
    mposlist = [sorted([str(aatomass[aa]) for aa in pos], key=lambda x: float(x)) for pos in positionlist]
    #Then concatenate them into strings and count unique strings to collapse redundant positions into combinations instead of permutations.
    flatcnt = collections.Counter([','.join(sub) for sub in mposlist])
    #For redundant positions, generate a tuple of combinations of the positions as floats, else just convert to a list of floats.
    comblist = [[[float(b) for b in c] for c in itertools.combinations_with_replacement(sub.split(','),flatcnt[sub])] if flatcnt[sub] > 1 else [float(b) for b in sub.split(',')] for sub in flatcnt]
    #Combinations are inside tuples, so we need to flatten the list before we can sum it for each product. Format fixes float math imprecision for otherwise unrecognized duplicates.
    permlist = [float('{:.{digits}f}'.format((sum(flatten(p))-amidebonds*18.0105646863),digits=5)) for p in itertools.product(*comblist)]
    #If there are duplicate masses in different positions then there will still be duplicate masses. Remove them.
    return list(set(permlist))

def stringid (target,spectra):
    """Generate a string from target, time, and scans.
    """
    return '{},{},{}'.format(target,spectra[0]['scan start time']*60,sorted([s['id'] for s in spectra]))
    
def copyspectrum (spec):
    """Copies all the bits of a spectrum that I care about.
    Used to avoid a slow deepcopy triggered by spectrum.deRef()
    """
    newspectrum = pymzml.spec.Spectrum(spec.measuredPrecision)
    for k in ['ms level','filter string','id','precursors','scan start time','selected ion m/z','total ion current','charge state','MS:1000512']:
        try:
            newspectrum[k] = spec[k]
        except KeyError:
            pass
    newspectrum._peaks = spec.peaks
    newspectrum._i = spec.i
    newspectrum._mz = spec.mz
    return newspectrum
    
def processSpectra (mzml, precision, scanranges, targets, positionlist, aatomass, linear, noiselevel, verbose):
    """Produces a list of processed spectra data structures from the input mzml file, filtering out MSMS spectra which could not be library members.
            Determine mass range of interest.
            Determine presence of isotopes (currently not used).
            Read in spectral data (dumping lots of info we don't need).
            Resolve multiply-charged MS2 spectra (only to +2 currently), then infer an accurate parent mass to the spectrum by using the MS1 scan which was used to determine it as a target.
            If scan ranges were specified, then prune away scans outside them.
            Approved MS2 spectra originating from a single MS1 spectra are grouped into time-neighborhoods.
            Generate all library peptide masses at high resolution, then filter away MS2 spectra that don't match any.
            Group spectra of the same high res mass.
            Subdivide groups by time-adjacency (acceptable adjacency rules learned from the input file automatically, defined as no more than 1 FULL scan cycle apart) and by intensity-slope (to separate nearby peaks).
            Combine spectra which are still grouped by spectral addition.
    """
    MS1precision=precision[1]
    MS2precision=precision[0]
    #Generate search bounds and constants
    #What masses should we care about?
    amidebonds = len(positionlist)-1 if linear else len(positionlist)
    positionalmasslist = [[aatomass[p] for p in pos] for pos in positionlist]
    mintarget = sum([min(pos) for pos in positionalmasslist])-amidebonds*18.0105646863+1.0078250321-MS1precision #smallest [M+H]
    maxtarget = sum([max(pos) for pos in positionalmasslist])-amidebonds*18.0105646863+22.98976967+MS1precision #largest [M+Na]
    #Isotopes?
    isotopes=False
    for pos in positionlist:
        for aa in pos:
            if '*' in aa:
                isotopes=True
                break
        if isotopes:
            break
    # if isotopes and verbose:
        # print('Isotopes are present in the library. MS2 deisotoping aborted.')
    # if MS2precision > 0.25 and verbose:
        # print('MS2 precision is insufficient to reliably deisotope. MS2 deisotoping aborted.')
    #gobble up the spectra; there really shouldn't be an mzml file larger than your RAM here. That would be quite impressive.
    reader = pymzml.run.Reader(mzml, MSn_Precision=MS2precision/maxtarget)
    spectralist = [copyspectrum(spec) for spec in reader]
    processedspectra = {}
    allpeaks=[]
    if verbose:
        print('Processing all MSMS spectra within extended library mass bounds {:.2f} and {:.2f}.'.format(mintarget,maxtarget))
    #Inspect the MS2 spectra and infer their properties from the corresponding MS1 spectrum.
    neighborlist = []
    nlenlist = []
    correctedMS2list = []
    parentspectrum = None
    #Correct MS2 precursors.
    for spectrum in spectralist:
        if spectrum['ms level'] == 1:
            if neighborlist == []:
                parentspectrum = copy.copy(spectrum)
            elif neighborlist != []:
                precursorCorrect(parentspectrum, neighborlist, mintarget, maxtarget, MS2precision, MS1precision, isotopes)
                correctedMS2list.extend(filterbyscan(targets, scanranges, neighborlist, mintarget, maxtarget, MS2precision, MS1precision))
                nlenlist.append(len(neighborlist))
                parentspectrum = copy.copy(spectrum)
                neighborlist = []
        if spectrum['ms level'] == 2:
            neighborlist.append(spectrum)
            allpeaks.extend([x[1] for x in spectrum.peaks])
    #Clean up after the loop since we never hit a capping MS1 event.
    if neighborlist != []:
        precursorCorrect(parentspectrum, neighborlist, mintarget, maxtarget, MS2precision, MS1precision, isotopes)
        correctedMS2list.extend(filterbyscan(targets, scanranges, neighborlist, mintarget, maxtarget, MS2precision, MS1precision))
        nlenlist.append(len(neighborlist))
    #Needed for grouping on time-adjacency
    neighborhoodlen = max(nlenlist)*2    
    #Sort into high res mass bins
    highresmasses = collections.defaultdict(list)
    for MSMS in correctedMS2list:
        #k: corrected precursor mass, v:spectrum
        highresmasses[MSMS['selected ion m/z']].append(MSMS)
    #Merge high res masses that are close.
    alreadymatched=set()
    hrmasslist = sorted(highresmasses.keys())
    for i in range(len(hrmasslist)):
        if i in alreadymatched:
            continue
        else:
            for j in range(i+1,len(hrmasslist)):
                if hrmasslist[j] < hrmasslist[i]+MS1precision:
                    highresmasses[hrmasslist[i]].extend(highresmasses[hrmasslist[j]])
                    alreadymatched.add(j)
                else:
                    break
    for number in alreadymatched:
        del highresmasses[hrmasslist[number]]     
    #Determine possible targets if we don't know them already.
    if targets == '*':
        #See if high res masses are in library.
        uniquemasslist = sorted(uniquemass(positionlist, aatomass, linear))
        for i,j in zip(itertools.count(),uniquemasslist):
            uniquemasslist[i]+=1.0078250321 #Assuming H+ Ions
        hrmasslist = sorted(highresmasses.keys())
        m = 0
        mlen = len(uniquemasslist)
        targets = []
        for hrm in hrmasslist:
            while uniquemasslist[m]+MS1precision < hrm:
                if m+1 < mlen:
                    m+=1
                else:
                    break
            if uniquemasslist[m]-MS1precision > hrm:
                continue
            if uniquemasslist[m]+MS1precision > hrm:
                targets.append(hrm)
    #Now group them into target:spectra pairings
    if targets == '*' or '*' in scanranges: #Either all scanranges are '*' or none are.
        groupedspectra = groupSpectra(targets, neighborhoodlen, highresmasses)
    else: #Already grouped! Reogranize it.
        groupedspectra = {'targets':[],'spectra':[]}
        for t in highresmasses:
            groupedspectra['targets'].append(t)
        groupedspectra['targets'] = sorted(groupedspectra['targets'])
        for t in groupedspectra['targets']:
            groupedspectra['spectra'].append(highresmasses[t])
    #Report findings if verbose.
    if verbose:
        print('Found {} targets:'.format(len(groupedspectra['targets'])))
        print(groupedspectra['targets'])
        print('At scan ranges:')
        print([[s['id'] for s in g] for g in groupedspectra['spectra']])
    #Error-check
    if groupedspectra['targets'] == []:
        print('No targets were found.')
        sys.exit(0)
    #Combine spectra.
    groupedspectra['combined']=[]
    for target,spectra in zip(groupedspectra['targets'],groupedspectra['spectra']):
        groupedspectra['combined'].append(spectralcombine(target,spectra,precision,verbose))
    #Filter out noise peaks
    groupedspectra['combined'] = filterpeaks(groupedspectra['combined'],allpeaks,verbose,noiselevel)
    prettyspectra = {} #k:uid v:spectrum
    for target,spectra,combined in zip(groupedspectra['targets'],groupedspectra['spectra'],groupedspectra['combined']):
        if combined != []:
            if len(combined.peaks) > 1:
                prettyspectra[stringid(target,spectra)]=combined
            elif verbose:
                print('No remaining peaks in target {} after noise filtration.'.format(combined['uid']))
    #Graphing
    if verbose > 1:
        for target,spectra,combined in zip(groupedspectra['targets'],groupedspectra['spectra'],groupedspectra['combined']):
            graphspectra(allpeaks,target,spectra,combined)
    return prettyspectra

def getPrecursorFromSpectrum (spectrum):
    """Grabs the precursor value recorded in the Thermo filter string attribute when possible.
    This is not the same as the number stored in the 'precursor' dictionary value in certain circumstances, but is the accurate one of the two in those cases.
    """
    try:
        return float(spectrum['MS:1000512'].split()[7].split('@')[0])
    except KeyError:
        return float(spectrum["selected ion m/z"])

def getMassrangeFromSpectrum (spectrum):
    """Grabs the mass range from the filter string when possible. It's the only place I've found this value in a spectrum object (inexplicably).
    """
    try:
        r = spectrum['MS:1000512'].split()[-1].split('-')
        r0 = float(r[0][1:])
        r1 = float(r[1][:-1])
        return [r0,r1]
    except KeyError:
        return [spectrum.mz[0]-1,spectrum.mz[-1]+1]
    
def precursorCorrect (MS1, MS2list, mintarget, maxtarget, MS2precision, MS1precision, isotopes):
    """Mutator function of MS2 that corrects the precursor ion to a more accurate value by inferring it from the parental MS1 spectrum.
    Also corrects for doubly charged precursor ions (but not higher for now) and isotopic variations.
    May also include de-isotoping in the future.
    """
    precdigits = len(str(MS1precision).split('.')[1])
    mz = MS1.mz
    intensity = MS1.i
    for MS2 in MS2list:
        initialprecursor = getPrecursorFromSpectrum(MS2)
        lower = bisect.bisect_right(mz, initialprecursor-MS2precision)
        upper = bisect.bisect_left(mz, initialprecursor+MS2precision)
        if upper == len(mz)-1:
            index = intensity.index(max(intensity[lower:]))
        elif upper == lower:
            index = upper
        else:
            index = intensity.index(max(intensity[lower:upper]))
        corrected = mz[index]
        #1. Determine charge state for the precursor
        split = bisect.bisect_right(MS2.mz,corrected+MS1precision+corrected/10) #Make sure to not catch isotopes by adding mz/10
        if MS2['total ion current'] == 0:
            intensityfrac = 0
        else:
            intensityfrac = sum(MS2.i[split:])/MS2['total ion current']
        possiblecharges = []
        MS2['charge state'] = 1
        if intensityfrac > 0.25: #Multiply charged? yes/no
            chargestosearch = [2,3,4,5]
            for c in chargestosearch: #How many charges?
                proposedmz = (corrected*c)-(1.0078250321*c-1)
                if mintarget < proposedmz < maxtarget:
                    split = bisect.bisect_right(MS2.mz,proposedmz+MS1precision+proposedmz/10)
                    intensityfrac = sum(MS2.i[split:])/MS2['total ion current']
                    #If intensityfrac is significant, then the true precursor is larger (more highly charged).
                    #If not, we have found the charge.
                    if intensityfrac < 0.25: 
                        MS2['charge state'] = c
                        break
        #ON HOLD; Possibly of low importance if data is acquired reasonably. Only possibly useful when we get up to 11mers or so.
        # #1. Look for isotopic envelopes in the MS2 spectrum.
        # if (not isotopes) and MS2precision < 0.25:
            # envlist = []
            # firstpeak = None
            # for p in MS2.peaks:
                # if p[1] > 50:
                    
        # #2. If they're found, what size are they?
        # #3. Derive monoisotopic mass from that!
        # #4. Deisotope the MS2 spectrum
        # MS2 = deisotopeSpectrum(MS2,MS2precision)
        #7. Correct for charge
        corrected = (corrected*MS2['charge state'])-(1.0078250321*(MS2['charge state']-1))
        MS2['selected ion m/z'] = float('{:.{digits}f}'.format(corrected,digits=precdigits))

def filterbyscan (targets, scanranges, MS2list, mintarget, maxtarget, MS2precision, MS1precision):
    """Takes into account the various target:scanrange correspondances and filters for MS2 events which satisfy them.
    """
    approved = []
    for spectrum in MS2list:
        if targets == '*': #If we don't know what we're looking for yet.
            if scanranges != '*': #But we know where to look.
                if scanranges[0] <= spectrum['id'] <= scanranges[1]:
                    if mintarget < spectrum['selected ion m/z'] < maxtarget:
                        approved.append(spectrum)
            elif mintarget < spectrum['selected ion m/z'] < maxtarget: #And we don't know where to look.
                approved.append(spectrum)
        else: #We do know what we're looking for.
            for target,sr in zip(targets,scanranges):
                if target-MS1precision < spectrum['selected ion m/z'] < target+MS1precision:
                    if sr != '*':
                        if sr[0] <= spectrum['id'] <= sr[1]:                
                            approved.append(spectrum)
                            break
                    else:
                        approved.append(spectrum)
                        break
    return approved
        
def groupSpectra (targets, neighborhoodlen, highresmasses):
    """Splits spectra into different groupings based on adjacency.
    """
    #split into target, spectra paired lists based on scan adjacencies.
    paired= {'targets':[],'spectra':[]}
    finalspectra=[]
    for target in targets:
        neighbors = []
        lastid = None
        lasttotI = None
        lastslopedown = None
        for spectrum in sorted(highresmasses[target], key=lambda x: x['id']):
            scan_id = spectrum['id']
            #If we have nothing to compare to or we are close enough to group, append it.
            if neighbors == []: 
                neighbors.append(spectrum)
            elif lastid:
                if scan_id-lastid <= neighborhoodlen:
                    slopedown = True if lasttotI*1.05 > spectrum['total ion current'] else False
                    #If within proper distance and slope goes from down to up, start a new group.
                    if lastslopedown and not slopedown:
                        paired['targets'].append(target)
                        paired['spectra'].append(neighbors)
                        neighbors = [spectrum]
                        lastslopedown = None
                    #Slope change does not indicate a new peak
                    else:
                        neighbors.append(spectrum)
                        lastslopedown = slopedown
                #Otherwise, start a new group with this spectrum
                else:
                    paired['targets'].append(target)
                    paired['spectra'].append(neighbors)
                    neighbors = [spectrum]
                    lastslopedown = None
            lastid = scan_id
            lasttotI = spectrum['total ion current']
        if neighbors != []:
            paired['targets'].append(target)
            paired['spectra'].append(neighbors)
    return paired
    
def graphspectra (modelinput,target,spectra,combined):
    """Graphs spectra and shows where various noise thresholds fall on them.
       Useful for troubleshooting.
    """
    import matplotlib
    from matplotlib import pyplot
    if spectra == [] or len(combined.peaks) <= 1:
        print('Nothing to plot for target {}'.format(combined['uid']))
        return
    usecombined = True
    if len(combined.peaks)==0:
        usecombined = False
    print("Generating Plots.")
    numplots = len(spectra)
    if usecombined:
        numplots +=1
    uberplot = pyplot.figure(figsize=(15,3*numplots))
    pyplot.suptitle('MSMS spectra for target {}'.format(stringid(target,spectra)),size='x-large')
    pyplot.axis('off')
    plotnum=1
    thresholdsy=[]
    thresholdsy.append(noiseleveltointensitythreshold(modelinput,0,.1))
    thresholdsy.append(noiseleveltointensitythreshold(modelinput,0,.05))
    thresholdsy.append(noiseleveltointensitythreshold(modelinput,0,.01))
    massrangesmin = []
    massrangesmax = []
    for spectrum in spectra:
        mssrng= getMassrangeFromSpectrum(spectrum)
        massrangesmin.append(mssrng[0])
        massrangesmax.append(mssrng[1])
    massrange = [min(massrangesmin),max(massrangesmax)]
    ax = uberplot.add_subplot(numplots,1,plotnum)
    pyplot.title('Combined spectrum for target {}'.format(target))
    ax.stem(combined.mz,combined.i,color='b',markerfmt=' ')
    ax.hlines(thresholdsy,massrange[0],massrange[1],color='r',label='90%,95%,99% confidence in peaks not being noise.')
    ax.legend()
    ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(10))
    ax.set_xlim(massrange)
    plotnum +=1
    for spectrum in sorted(spectra,key=lambda x: x['id']):
        ax = uberplot.add_subplot(numplots,1,plotnum)
        pyplot.title('Scan {} of parent ion {}'.format(spectrum['id'],spectrum['selected ion m/z']),size='large')
        ax.stem(spectrum.mz,spectrum.i,color='b',markerfmt=' ')
        ax.hlines(thresholdsy,massrange[0],massrange[1],color='r',label='90%,95%,99% confidence in peaks not being noise.')
        if plotnum == numplots:
            pyplot.xlabel('m/z',size='large')
        ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(10))
        ax.ticklabel_format(style='plain')
        ax.set_xlim(massrange)
        plotnum +=1
    uberplot.subplots_adjust(left=0.15, bottom=0.05, right=0.97, top=0.92, wspace=0.0, hspace=0.35)
    uberplot.text(0.02, 0.5, 'Intensity', va='center', rotation='vertical', size='large')
    uberplot.savefig('MSMS spectra for target {}.png'.format(stringid(target,spectra)))
    pyplot.close(uberplot)

def filterpeaks (spectra,modelinput,verbose,noiselevel):
    """Mutator of spectra that filters out noise peaks by comparing to a specified or model-based threshold value.
    The fifty peaks of highest intensity are then accepted.
    """
    newspectra = []
    maxpeaks = 100 #Arbitrary
    if noiselevel > 1.0:
        minintensity=noiselevel
    else:
        minintensity = noiseleveltointensitythreshold(modelinput,verbose,noiselevel)
    if verbose:
        print('Filtering out peaks below intensity {}'.format(minintensity))
    for spectrum in spectra:
        newlist = []
        for peak in spectrum.peaks:
            if peak[1] > minintensity:
                newlist.append(peak)
        if len(newlist) > maxpeaks:
            newlist.sort(key=lambda x: x[1],reverse=True)
            newlist[:] = newlist[:maxpeaks]
            newlist.sort(key=lambda x: x[0])
        elif newlist == []:
            spectrum.peaks = [(0,0)]
        spectrum.peaks = newlist #This should update mz and i properly through the setter function.
        spectrum['total ion current']=sum(spectrum.i)
    return spectra
    
def noiseleveltointensitythreshold (modelinput,verbose,noiselevel):
    """Generates a noise model and locates the intensity at which the p-value of a peak being noise is equal to noiselevel.
    The assumption is that the top 10% includes all signal peaks, and is therefore not relevant to the noise model we are building.
    """
    if verbose:
        print('Generating noise model from MSMS data.')
    modelinput.sort()
    ninetenthslen = (len(modelinput)*9)//10
    x = np.asarray(modelinput[:ninetenthslen])
    if len(x) > 30000: #Arbitrary; adjust based on speed of cdf evaluation.
        x = np.random.choice(x,size=30000,replace=False)
    logx = np.log(x) #Transform from [0,inf] to [-inf,inf], which is the space the kernel is in
    kde = sm.nonparametric.KDEUnivariate(logx)
    kde.fit(gridsize=1000)
    truesupport = math.e**kde.support #Convert back to bounded [0,inf] space
    if verbose:
        print('Finding intensity threshold for peaks for which the likelihood of a peak at that intensity or greater being noise is less than {}%'.format(noiselevel*100))
    idx = (np.abs(kde.cdf-(1-noiselevel))).argmin()
    minintensity = truesupport[idx]
    if verbose:
        print('Intensity threshold found to be {}'.format(minintensity))
    return minintensity

def spectralcombine (target,spectra,precision,verbose):
    """Combines spectra for a single target in order to amplify signal to noise ratio.
    Makes the assumption that all signal peaks will be present in the most intense MSMS spectrum of the target mass.
    """
    uniqueidentifier = stringid(target,spectra)
    MS2precision = precision[0]
    MS1precision = precision[1]
    if verbose:
        print('Combining spectra for {}'.format(uniqueidentifier))
    intensitythreshold=1
    intensityweight=0.1
    combinedthreshold=MS2precision+intensitythreshold*intensityweight
    if len(spectra) == 1:
        spectra[0]['uid'] = uniqueidentifier
        return spectra[0]
    elif len(spectra) == 0:
        return spectra
    else:
        for spectrum in spectra:
            #Do normalization here
            spectrum['normalizedIntensities']=[math.log(x)-math.log(spectrum['total ion current']) for x in spectrum.i]
        spectra.sort(key= lambda x: x['total ion current'],reverse=True)
        #Peak-by-peak pairing method. 
        #Distance and intensity similarity thresholds used to find most similar peak and combine.
        masterspectrum = spectra[0]
        correspondances = [[] for thing in spectra[1:]]
        for i,slavespectrum in zip(itertools.count(),spectra[1:]):
            masterdistmat = np.tile(np.array(masterspectrum.mz),(len(slavespectrum.mz),1))#Horizontal copies out to len(peaklist[1])
            slavedistvector = np.array(slavespectrum.mz).transpose()#Vertical copies by numpy broadcast
            distmatrix = np.absolute(masterdistmat-slavedistvector[:,np.newaxis])#Reshape may be faster
            masterintmat = np.tile(np.array(masterspectrum['normalizedIntensities']),(len(slavespectrum['normalizedIntensities']),1))
            slaveintvector = np.array(slavespectrum['normalizedIntensities']).transpose()
            intmatrix = np.absolute(masterintmat-slaveintvector[:,np.newaxis])#Reshape may be faster
            scorematrix = distmatrix+intmatrix*intensityweight
            for j in range(len(masterspectrum.mz)):
                col=list(scorematrix[:,j])
                minscore = 9999
                minloc = 0
                for k,score in zip(itertools.count(0),col):
                    if score < minscore:
                        minscore = score
                        minloc=k
                #delta dist threshold = precison.
                if minscore < combinedthreshold:
                    correspondances[i].append([j,minloc])
        combined=copy.copy(spectra[0])
        plist = [list(p) for p in combined.peaks]
        for i,matchlist in zip(itertools.count(1),correspondances):
            for match in matchlist:
                plist[match[0]][1] += spectra[i].i[match[1]]
        combined.peaks = plist
        combined['total ion current']=sum(combined.i)
        combined['uid']=uniqueidentifier
        combined.pop('normalizedIntensities',None)
        return combined
    
class neutralloss ():
    """Generates ion variants based upon neutral loss and fragment composition.
    This includes alpha ions, counterion-less beta ion bx-yz fragments, ammonia and water loss, and water gain.
    """
    def __init__ (self, aadb, aatomass):
        self.H2O = 18.0105646863
        self.H = 1.0078250321
        self.Na = 22.98976967
        self.CO = 27.9949146221
        self.NH3 = 17.0265491015
        self.ammonialoss = set(['R','K','Q','N','r','k','q','n','MR','MK','MQ','MN','Mr','Mk','Mq','Mn'])
        self.waterloss = set(['S','T','s','t','MS','MT','Ms','MT'])
        self.Eset = set(['E','e','ME','Me'])
        self.aadb = aadb
        self.ctermloss = set(['k','K','Mk','MK','r','R','Mr','MR'])
        self.aatomass = aatomass
    def get (self, parent, name, mass):
        """Takes parent name, fragment name, fragment mass.
        Returns dictionary of name:mass for neutral losses of that fragment including the original fragment.
        May require modification to include/disinclude neutral losses commonly observed using your mass spectrometry system and residues.
        """
        if name.endswith('_y'):
            beta = False
            namelist = name[:-2].split(',')
            namelen = len(namelist)
        else:
            beta = True
        parentlist = parent.split(',')
        parentlen= len(parentlist)
        fragdict = {name:mass}
        #Do C-terminus loss such that it points back to the original fragment.
        if not beta:
            parentset = set(parentlist)
            if parentset & self.ctermloss: 
                if namelen-1 > 1:
                    nocterm = ','.join(namelist[:-1])
                    fragdict[nocterm+'_y'] = getmass(nocterm,self.aatomass,True)+1.0078250321 #+H
        #Neutral loss must be treated separately for each fragment if there is a c-terminal loss because the new fragment has a new sequence.
        finaldict = {}
        for fragname in fragdict:
            nldict = {fragname:fragdict[fragname]}
            if beta:
                fragnamelist = name.split(',')
            else:
                fragnamelist = fragname[:-2].split(',')
            fragnamelen = len(fragnamelist)
            #Count residues which commonly have neutral losses.
            ammonialosscount = 0
            waterlosscount = 0
            for residue in fragnamelist[:-1]:
                if residue in self.ammonialoss:
                    ammonialosscount +=1
                if residue in self.waterloss:
                    waterlosscount +=1
            #Positionally variant neutral loss
            if fragnamelist[-1] in self.Eset:
                waterlosscount+=1
            if not beta: #Pathways only relevant to gamma ions.
                tempdict = {}
                if fragnamelist[-1] in self.ammonialoss:
                    tempdict[fragname+'-OH'] = nldict[fragname]-self.H2O #Caprolactam and other side-chain OH removals on C-terminus of Y fragments.
                #If beta, then the caprolactamization and diketopiperizine pathways will not change the mass.
                if 5 > parentlen-fragnamelen > 1:
                    # The below is to make this cis-amide dependent.
                    # cis=False
                    # for residue in range(0,(parentlen-namelen)):
                        # if self.aadb[parentlist[parentlen-namelen-2]][1]=='(':
                            # cis=True
                            # break
                    #if cis:
                    for mol in nldict:
                        newname = mol+'+H' #Diketopiperizine pathways
                        tempdict[newname] = nldict[mol]+self.H
                    nldict = merge_dicts(nldict, tempdict)
            for mol in nldict: #Follow through with the neutral losses detected above on the original fragment and gamma variants.
                tempdict = {}
                for i in range(1,ammonialosscount+1):
                    newname = mol+"-NH3"*i
                    tempdict[newname] = nldict[mol]-self.NH3*i
                for i in range(1, waterlosscount+1):
                    newname = mol+"-H2O"*i
                    tempdict[newname] = nldict[mol]-self.H2O*i
            nldict = merge_dicts(nldict, tempdict)
            finaldict = merge_dicts(finaldict, nldict)
        return finaldict

def peakmatch (compoundlist,aatomass,linear,precision,spectra,aadb): 
    """Find all peaks which match for a list of molecules across all spectra.
    returns a dictionary of uids:dict of compound names:hits records.
    """
    uids = sorted(spectra.keys())
    targets = [float(u.split(',')[0]) for u in uids]
    nloss = neutralloss(aadb,aatomass)
    #One dictionary of fragment masses per target, with each mass having a dictionary of compound names:Counter for subfrag times encountered. 
    #This catches multiple fragments with the same exact mass and traces neutral losses back to their main fragment
    uberfragdictpertarget = [collections.defaultdict(lambda : collections.defaultdict(collections.Counter)) for i in range(len(targets))]
    compoundstoanalyze = collections.defaultdict(set)
    for compound in compoundlist:
        pmass = getmass(compound.name,aatomass,linear)
        MH = pmass+1.0078250321
        MHlow,MHhigh = MH-precision[1],MH+precision[1]
        for i,target in zip(itertools.count(),targets):
            if MHlow <= target <= MHhigh:
                compoundstoanalyze[compound].add(i)
    for compound in compoundstoanalyze:
        subfrags = getfrags(compound.name, aatomass, linear)
        for sub in subfrags:
            for neutral in nloss.get(compound.name,sub,subfrags[sub]).values():
                for i in compoundstoanalyze[compound]:
                    uberfragdictpertarget[i][round(neutral,10)][compound][sub]+=1
    hitsdict = {u:{} for u in spectra.keys()}
    for i,uid in zip(itertools.count(),uids):
        spectrum = spectra[uid]
        hitrecord = hitsdict[uid]
        mssrng = getMassrangeFromSpectrum(spectrum)
        allmasses = sorted(uberfragdictpertarget[i].keys())
        if allmasses == []:
            continue
        lower = bisect.bisect_right(allmasses,mssrng[0])
        upper = bisect.bisect_left(allmasses,mssrng[1])
        if upper == len(allmasses)-1:
            masses = allmasses[lower:]
        else:
            masses = allmasses[lower:upper]
        match(uberfragdictpertarget[i],spectrum,masses,precision,hitrecord)
    return hitsdict

def match (uberfragdict,spectrum,masses,precision,hitrecord):
    """Runs through two sorted lists without backtracking to compare one list of points to another list of intervals.
    Mutates hitrecord.
    """
    p = 0
    plen = len(spectrum.mz)
    peaksmatched = collections.defaultdict(set)
    for i in masses:
        #If we are at a mass higher than the upper bound, we need to advance the bounds.
        while spectrum.mz[p]+precision[0] < i:
            if p+1 < plen:
                p+=1
            else:
                break
        #If we are at a mass lower than the lower bound, we need to advance in mass.
        if spectrum.mz[p]-precision[0] > i:
            continue
        #Are we between the bounds?
        if spectrum.mz[p]+precision[0] > i:
            #add to list of hits
            fraghit(uberfragdict,spectrum,i,p,peaksmatched,hitrecord)

def fraghit (uberfragdict,spectrum,mass,peaknum,peaksmatched,hitrecord):
    """Records the subfragments or neutral loss fragments of a compound which are observed.
    Records the total the intensity matched (across all peaks matched to for a compound) without allowing intensity duplications.
    The case where multiple of the same fragment generate from a single compound is counted correctly.
    """
    for compound in uberfragdict[mass]:
        if compound not in hitrecord.keys():
            hitrecord[compound] = [uberfragdict[mass][compound],[spectrum.i[peaknum]]]
        else:
            hitrecord[compound][0] = hitrecord[compound][0]+uberfragdict[mass][compound]
            if peaknum not in peaksmatched[compound]:
                hitrecord[compound][1].append(spectrum.i[peaknum])
        peaksmatched[compound].add(peaknum)

class Worker (multiprocessing.Process):
    def __init__(self, inqueue, outqueue, params):
        super(Worker, self).__init__()
        self.inqueue=inqueue
        self.outqueue=outqueue
        self.params=params
    def run(self):
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        aatomass,linear,precision,spectra,aadb = self.params
        while True:
            compoundlist = self.inqueue.get()
            if compoundlist == None:
                break
            else:
                result = peakmatch(compoundlist,aatomass,linear,precision,spectra,aadb)
                self.outqueue.put(result) #Avoid overloading.
        self.outqueue.put(None)

class exactmass ():
    def __init__ (self, aadict):
        self.aadict = aadict
        
        self.symboltomass = {
        'Ac' : {227: 227.027747},
        'Ag' : {107: 106.905093, 109: 108.904756},
        'Al' : {27: 26.98153844},
        'Am' : {243: 243.0613727},
        'Ar' : {40: 39.962383123, 36: 35.96754628, 38: 37.9627322},
        'As' : {75: 74.9215964},
        'At' : {210: 209.987131},
        'Au' : {197: 196.966552},
        'B' : {10: 10.012937, 11: 11.0093055},
        'Ba' : {130: 129.90631, 132: 131.905056, 134: 133.904503, 135: 134.905683, 136: 135.90457, 137: 136.905821, 138: 137.905241},
        'Be' : {9: 9.0121821},
        'Bh' : {264: 264.12473},
        'Bi' : {209: 208.980383},
        'Bk' : {247: 247.070299},
        'Br' : {81: 80.916291, 79: 78.9183376},
        'C' : {12: 12.0, 13: 13.0033548378},
        'Ca' : {40: 39.9625912, 42: 41.9586183, 43: 42.9587668, 44: 43.9554811, 46: 45.9536928, 48: 47.952534},
        'Cd' : {106: 105.906458, 108: 107.904183, 110: 109.903006, 111: 110.904182, 112: 111.9027572, 113: 112.9044009, 114: 113.9033581, 116: 115.904755},
        'Ce' : {136: 135.90714, 138: 137.905986, 140: 139.905434, 142: 141.90924},
        'Cf' : {251: 251.07958},
        'Cl' : {35: 34.96885271, 37: 36.9659026},
        'Cm' : {247: 247.070347},
        'Co' : {59: 58.9332002},
        'Cr' : {50: 49.9460496, 52: 51.9405119, 53: 52.9406538, 54: 53.9388849},
        'Cs' : {133: 132.905447},
        'Cu' : {65: 64.9277937, 63: 62.9296011},
        'Db' : {262: 262.11415},
        'Dy' : {160: 159.925194, 161: 160.92693, 162: 161.926795, 163: 162.928728, 164: 163.929171, 156: 155.924278, 158: 157.924405},
        'Er' : {162: 161.928775, 164: 163.929197, 166: 165.93029, 167: 166.932045, 168: 167.932368, 170: 169.93546},
        'Es' : {252: 252.08297},
        'Eu' : {153: 152.921226, 151: 150.919846},
        'F' : {19: 18.9984032},
        'Fe' : {56: 55.9349421, 57: 56.9353987, 58: 57.9332805, 54: 53.9396148},
        'Fm' : {257: 257.095099},
        'Fr' : {223: 223.0197307},
        'Ga' : {69: 68.925581, 71: 70.924705},
        'Gd' : {160: 159.927051, 152: 151.919788, 154: 153.920862, 155: 154.922619, 156: 155.92212, 157: 156.923957, 158: 157.924101},
        'Ge' : {72: 71.9220762, 73: 72.9234594, 74: 73.9211782, 76: 75.9214027, 70: 69.9242504},
        'H' : {1: 1.0078250321, 2: 2.014101778},
        'He' : {3: 3.0160293097, 4: 4.0026032497},
        'Hf' : {174: 173.94004, 176: 175.9414018, 177: 176.94322, 178: 177.9436977, 179: 178.9458151, 180: 179.9465488},
        'Hg' : {196: 195.965815, 198: 197.966752, 199: 198.968262, 200: 199.968309, 201: 200.970285, 202: 201.970626, 204: 203.973476},
        'Ho' : {165: 164.930319},
        'Hs' : {269: 269.13411},
        'I' : {127: 126.904468},
        'In' : {113: 112.904061, 115: 114.903878},
        'Ir' : {193: 192.962924, 191: 190.960591},
        'K' : {40: 39.96399867, 41: 40.96182597, 39: 38.9637069},
        'Kr' : {78: 77.920386, 80: 79.916378, 82: 81.9134846, 83: 82.914136, 84: 83.911507, 86: 85.9106103},
        'La' : {138: 137.907107, 139: 138.906348},
        'Li' : {6: 6.0151223, 7: 7.016004},
        'Lr' : {262: 262.10969},
        'Lu' : {176: 175.9426824, 175: 174.9407679},
        'Md' : {258: 258.098425},
        'Mg' : {24: 23.9850419, 25: 24.98583702, 26: 25.98259304},
        'Mn' : {55: 54.9380496},
        'Mo' : {96: 95.9046789, 97: 96.906021, 98: 97.9054078, 100: 99.907477, 92: 91.90681, 94: 93.9050876, 95: 94.9058415},
        'Mt' : {268: 268.13882},
        'N' : {14: 14.0030740052, 15: 15.0001088984},
        'Na' : {23: 22.98976967},
        'Nb' : {93: 92.9063775},
        'Nd' : {142: 141.907719, 143: 142.90981, 144: 143.910083, 145: 144.912569, 146: 145.913112, 148: 147.916889, 150: 149.920887},
        'Ne' : {20: 19.9924401759, 21: 20.99384674, 22: 21.99138551},
        'Ni' : {64: 63.9279696, 58: 57.9353479, 60: 59.9307906, 61: 60.9310604, 62: 61.9283488},
        'No' : {259: 259.10102},
        'Np' : {237: 237.0481673},
        'O' : {16: 15.9949146221, 17: 16.9991315, 18: 17.9991604},
        'Os' : {192: 191.961479, 184: 183.952491, 186: 185.953838, 187: 186.9557479, 188: 187.955836, 189: 188.9581449, 190: 189.958445},
        'P' : {31: 30.97376151},
        'Pa' : {231: 231.0358789},
        'Pb' : {208: 207.976636, 204: 203.973029, 206: 205.974449, 207: 206.975881},
        'Pd' : {102: 101.905608, 104: 103.904035, 105: 104.905084, 106: 105.903483, 108: 107.903894, 110: 109.905152},
        'Pm' : {145: 144.912744},
        'Po' : {209: 208.982416},
        'Pr' : {141: 140.907648},
        'Pt' : {192: 191.961035, 194: 193.962664, 195: 194.964774, 196: 195.964935, 198: 197.967876, 190: 189.95993},
        'Pu' : {244: 244.064198},
        'Ra' : {226: 226.0254026},
        'Rb' : {85: 84.9117893, 87: 86.9091835},
        'Re' : {185: 184.9529557, 187: 186.9557508},
        'Rf' : {261: 261.10875},
        'Rh' : {103: 102.905504},
        'Rn' : {222: 222.0175705},
        'Ru' : {96: 95.907598, 98: 97.905287, 99: 98.9059393, 100: 99.9042197, 101: 100.9055822, 102: 101.9043495, 104: 103.90543},
        'S' : {32: 31.97207069, 33: 32.9714585, 34: 33.96786683, 36: 35.96708088},
        'Sb' : {121: 120.903818, 123: 122.9042157},
        'Sc' : {45: 44.9559102},
        'Se' : {74: 73.9224766, 76: 75.9192141, 77: 76.9199146, 78: 77.9173095, 80: 79.9165218, 82: 81.9167},
        'Sg' : {266: 266.12193},
        'Si' : {28: 27.9769265327, 29: 28.97649472, 30: 29.97377022},
        'Sm' : {144: 143.911995, 147: 146.914893, 148: 147.914818, 149: 148.91718, 150: 149.917271, 152: 151.919728, 154: 153.922205},
        'Sn' : {112: 111.904821, 114: 113.902782, 115: 114.903346, 116: 115.901744, 117: 116.902954, 118: 117.901606, 119: 118.903309, 120: 119.9021966, 122: 121.9034401, 124: 123.9052746},
        'Sr' : {88: 87.9056143, 84: 83.913425, 86: 85.9092624, 87: 86.9088793},
        'Ta' : {180: 179.947466, 181: 180.947996},
        'Tb' : {159: 158.925343},
        'Tc' : {98: 97.907216},
        'Te' : {128: 127.9044614, 130: 129.9062228, 120: 119.90402, 122: 121.9030471, 123: 122.904273, 124: 123.9028195, 125: 124.9044247, 126: 125.9033055},
        'Th' : {232: 232.0380504},
        'Ti' : {48: 47.9479471, 49: 48.9478708, 50: 49.9447921, 46: 45.9526295, 47: 46.9517638},
        'Tl' : {203: 202.972329, 205: 204.974412},
        'Tm' : {169: 168.934211},
        'U' : {234: 234.0409456, 235: 235.0439231, 238: 238.0507826},
        'V' : {50: 49.9471628, 51: 50.9439637},
        'W' : {184: 183.9509326, 186: 185.954362, 180: 179.946706, 182: 181.948206, 183: 182.9502245},
        'Xe' : {128: 127.9035304, 129: 128.9047795, 130: 129.9035079, 131: 130.9050819, 132: 131.9041545, 134: 133.9053945, 136: 135.90722, 124: 123.9058958, 126: 125.904269},
        'Y' : {89: 88.9058479},
        'Yb' : {168: 167.933894, 170: 169.934759, 171: 170.936322, 172: 171.9363777, 173: 172.9382068, 174: 173.9388581, 176: 175.942568},
        'Zn' : {64: 63.9291466, 66: 65.9260368, 67: 66.9271309, 68: 67.9248476, 70: 69.925325},
        'Zr' : {96: 95.908276, 90: 89.9047037, 91: 90.905645, 92: 91.9050401, 94: 93.9063158}}
        
        self.numbertomass = {
        1 : {1: 1.0078250321, 2: 2.014101778},
        2 : {3: 3.0160293097, 4: 4.0026032497},
        3 : {6: 6.0151223, 7: 7.016004},
        4 : {9: 9.0121821},
        5 : {10: 10.012937, 11: 11.0093055},
        6 : {12: 12.0, 13: 13.0033548378},
        7 : {14: 14.0030740052, 15: 15.0001088984},
        8 : {16: 15.9949146221, 17: 16.9991315, 18: 17.9991604},
        9 : {19: 18.9984032},
        10 : {20: 19.9924401759, 21: 20.99384674, 22: 21.99138551},
        11 : {23: 22.98976967},
        12 : {24: 23.9850419, 25: 24.98583702, 26: 25.98259304},
        13 : {27: 26.98153844},
        14 : {28: 27.9769265327, 29: 28.97649472, 30: 29.97377022},
        15 : {31: 30.97376151},
        16 : {32: 31.97207069, 33: 32.9714585, 34: 33.96786683, 36: 35.96708088},
        17 : {35: 34.96885271, 37: 36.9659026},
        18 : {40: 39.962383123, 36: 35.96754628, 38: 37.9627322},
        19 : {40: 39.96399867, 41: 40.96182597, 39: 38.9637069},
        20 : {40: 39.9625912, 42: 41.9586183, 43: 42.9587668, 44: 43.9554811, 46: 45.9536928, 48: 47.952534},
        21 : {45: 44.9559102},
        22 : {48: 47.9479471, 49: 48.9478708, 50: 49.9447921, 46: 45.9526295, 47: 46.9517638},
        23 : {50: 49.9471628, 51: 50.9439637},
        24 : {50: 49.9460496, 52: 51.9405119, 53: 52.9406538, 54: 53.9388849},
        25 : {55: 54.9380496},
        26 : {56: 55.9349421, 57: 56.9353987, 58: 57.9332805, 54: 53.9396148},
        27 : {59: 58.9332002},
        28 : {64: 63.9279696, 58: 57.9353479, 60: 59.9307906, 61: 60.9310604, 62: 61.9283488},
        29 : {65: 64.9277937, 63: 62.9296011},
        30 : {64: 63.9291466, 66: 65.9260368, 67: 66.9271309, 68: 67.9248476, 70: 69.925325},
        31 : {69: 68.925581, 71: 70.924705},
        32 : {72: 71.9220762, 73: 72.9234594, 74: 73.9211782, 76: 75.9214027, 70: 69.9242504},
        33 : {75: 74.9215964},
        34 : {74: 73.9224766, 76: 75.9192141, 77: 76.9199146, 78: 77.9173095, 80: 79.9165218, 82: 81.9167},
        35 : {81: 80.916291, 79: 78.9183376},
        36 : {78: 77.920386, 80: 79.916378, 82: 81.9134846, 83: 82.914136, 84: 83.911507, 86: 85.9106103},
        37 : {85: 84.9117893, 87: 86.9091835},
        38 : {88: 87.9056143, 84: 83.913425, 86: 85.9092624, 87: 86.9088793},
        39 : {89: 88.9058479},
        40 : {96: 95.908276, 90: 89.9047037, 91: 90.905645, 92: 91.9050401, 94: 93.9063158},
        41 : {93: 92.9063775},
        42 : {96: 95.9046789, 97: 96.906021, 98: 97.9054078, 100: 99.907477, 92: 91.90681, 94: 93.9050876, 95: 94.9058415},
        43 : {98: 97.907216},
        44 : {96: 95.907598, 98: 97.905287, 99: 98.9059393, 100: 99.9042197, 101: 100.9055822, 102: 101.9043495, 104: 103.90543},
        45 : {103: 102.905504},
        46 : {102: 101.905608, 104: 103.904035, 105: 104.905084, 106: 105.903483, 108: 107.903894, 110: 109.905152},
        47 : {107: 106.905093, 109: 108.904756},
        48 : {106: 105.906458, 108: 107.904183, 110: 109.903006, 111: 110.904182, 112: 111.9027572, 113: 112.9044009, 114: 113.9033581, 116: 115.904755},
        49 : {113: 112.904061, 115: 114.903878},
        50 : {112: 111.904821, 114: 113.902782, 115: 114.903346, 116: 115.901744, 117: 116.902954, 118: 117.901606, 119: 118.903309, 120: 119.9021966, 122: 121.9034401, 124: 123.9052746},
        51 : {121: 120.903818, 123: 122.9042157},
        52 : {128: 127.9044614, 130: 129.9062228, 120: 119.90402, 122: 121.9030471, 123: 122.904273, 124: 123.9028195, 125: 124.9044247, 126: 125.9033055},
        53 : {127: 126.904468},
        54 : {128: 127.9035304, 129: 128.9047795, 130: 129.9035079, 131: 130.9050819, 132: 131.9041545, 134: 133.9053945, 136: 135.90722, 124: 123.9058958, 126: 125.904269},
        55 : {133: 132.905447},
        56 : {130: 129.90631, 132: 131.905056, 134: 133.904503, 135: 134.905683, 136: 135.90457, 137: 136.905821, 138: 137.905241},
        57 : {138: 137.907107, 139: 138.906348},
        58 : {136: 135.90714, 138: 137.905986, 140: 139.905434, 142: 141.90924},
        59 : {141: 140.907648},
        60 : {142: 141.907719, 143: 142.90981, 144: 143.910083, 145: 144.912569, 146: 145.913112, 148: 147.916889, 150: 149.920887},
        61 : {145: 144.912744},
        62 : {144: 143.911995, 147: 146.914893, 148: 147.914818, 149: 148.91718, 150: 149.917271, 152: 151.919728, 154: 153.922205},
        63 : {153: 152.921226, 151: 150.919846},
        64 : {160: 159.927051, 152: 151.919788, 154: 153.920862, 155: 154.922619, 156: 155.92212, 157: 156.923957, 158: 157.924101},
        65 : {159: 158.925343},
        66 : {160: 159.925194, 161: 160.92693, 162: 161.926795, 163: 162.928728, 164: 163.929171, 156: 155.924278, 158: 157.924405},
        67 : {165: 164.930319},
        68 : {162: 161.928775, 164: 163.929197, 166: 165.93029, 167: 166.932045, 168: 167.932368, 170: 169.93546},
        69 : {169: 168.934211},
        70 : {168: 167.933894, 170: 169.934759, 171: 170.936322, 172: 171.9363777, 173: 172.9382068, 174: 173.9388581, 176: 175.942568},
        71 : {176: 175.9426824, 175: 174.9407679},
        72 : {174: 173.94004, 176: 175.9414018, 177: 176.94322, 178: 177.9436977, 179: 178.9458151, 180: 179.9465488},
        73 : {180: 179.947466, 181: 180.947996},
        74 : {184: 183.9509326, 186: 185.954362, 180: 179.946706, 182: 181.948206, 183: 182.9502245},
        75 : {185: 184.9529557, 187: 186.9557508},
        76 : {192: 191.961479, 184: 183.952491, 186: 185.953838, 187: 186.9557479, 188: 187.955836, 189: 188.9581449, 190: 189.958445},
        77 : {193: 192.962924, 191: 190.960591},
        78 : {192: 191.961035, 194: 193.962664, 195: 194.964774, 196: 195.964935, 198: 197.967876, 190: 189.95993},
        79 : {197: 196.966552},
        80 : {196: 195.965815, 198: 197.966752, 199: 198.968262, 200: 199.968309, 201: 200.970285, 202: 201.970626, 204: 203.973476},
        81 : {203: 202.972329, 205: 204.974412},
        82 : {208: 207.976636, 204: 203.973029, 206: 205.974449, 207: 206.975881},
        83 : {209: 208.980383},
        84 : {209: 208.982416},
        85 : {210: 209.987131},
        86 : {222: 222.0175705},
        87 : {223: 223.0197307},
        88 : {226: 226.0254026},
        89 : {227: 227.027747},
        90 : {232: 232.0380504},
        91 : {231: 231.0358789},
        92 : {234: 234.0409456, 235: 235.0439231, 238: 238.0507826},
        93 : {237: 237.0481673},
        94 : {244: 244.064198},
        95 : {243: 243.0613727},
        96 : {247: 247.070347},
        97 : {247: 247.070299},
        98 : {251: 251.07958},
        99 : {252: 252.08297},
        100 : {257: 257.095099},
        101 : {258: 258.098425},
        102 : {259: 259.10102},
        103 : {262: 262.10969},
        104 : {261: 261.10875},
        105 : {262: 262.11415},
        106 : {266: 266.12193},
        107 : {264: 264.12473},
        108 : {269: 269.13411},
        109 : {268: 268.13882}}
        
    def get(self,allowed_aas):
        """Calculates the exact mass for all smiles string residues,
        returning a dictionary of said masses.
        """
        aatomass = {}
        for aa_name in allowed_aas:
            mol = Chem.AddHs(Chem.MolFromSmiles(self.aadict[aa_name]))
            counts = collections.Counter(atom.GetSymbol() for atom in mol.GetAtoms())
            mass = 0
            for element in counts:
                #make this work for isotopes later!
                mass += counts[element] * next(iter(self.symboltomass[element].values()))
            aatomass[aa_name] = mass
        return aatomass
 
def Cap (stringinput):
    """Capitalizes the first letter of the string, but does not change other capitalization.
    """
    if len(stringinput) > 1:
        if stringinput[0].isupper():
            return stringinput[0].lower()+stringinput[1:]
        else:
            return stringinput[0].upper()+stringinput[1:]
    else:
        if stringinput.isupper():
            return stringinput.lower()
        else:
            return stringinput.upper()
 
def recipebuild (aadb, constraint, rules, trunc, verbose, linear):
    """Takes the database and constraint string, and returns a recipe for the library generator functions.
    """
    aadict = dict()
    keylist = []
    if verbose:
        print('Loading amino acid database from {}.'.format(aadb.name))
    for line in aadb:
        list = line.strip().split()
        if len(list) == 3:
            aadict[list[0]] = list[2].split(',')
            keylist.append(list[0])
        elif len(list) == 2:
            aadict[list[0]]=list[1]
            keylist.append(list[0])
        else:
            print("Error: Database formatting is incorrect.")
            return
    if trunc: #Add the null residue to the database.
        #Set this to the mass of water so we don't have to count the length of each peptide to get its mass.
        aadict['_'] = 'O'
        keylist.append('_')
    #Create a set of all valid aadict keys
    validaaset = set(keylist)
    #Parse the constraint string
    if verbose:
        print('Parsing constraint string {}'.format(constraint))
    positionlist = constraint.strip().split(';')
    permutations = 1
    aainconstraint = set()
    isotopedict = {}
    for index,position in zip(itertools.count(),positionlist): 
        positionlist[index] = position.split(',')
        if trunc and index != len(positionlist)-1: #Add null residues to represent truncations to each position except the last (ostensibly on-resin in SPPS).
            positionlist[index].append('_')
        #Unpack the multitoken entries into the constraint variable
        for x,name in zip(itertools.count(),positionlist[index]):
            if name[0] == '$':
                for subname in aadict[name]:
                    positionlist[index].append(subname)
                del positionlist[index][x]
        permutations = permutations * len(positionlist[index])
        databaseError = 0
        for name in positionlist[index]:
            if name not in validaaset:
                #Sort out isotopes encoded as asterisks.
                extraneutrons=0
                if '*' in name:
                    extraneutrons = name.count('*')
                    isotopedict[name]=extraneutrons
                    name = name.replace('*','')
                if name in validaaset:
                    aadict[name+'*'*extraneutrons]=aadict[name]
                    validaaset.add(name+'*'*extraneutrons)
                    aainconstraint.add(name+'*'*extraneutrons)
                #Dynamic stereo and methylation can occur here before error messages.
                #Experimental feature, therefore not advertised.
                elif name.startswith('M'):
                    if name[1:] in validaaset:
                        #Generate nmethyl L-amino acid from L-amino acid
                        if aadict[name[1:]][1] == '(':
                            print('Amino acid {} cannot be Nmethylated, because it is already a secondary amine'.format(name[1:]))
                            databaseError=1
                        else:
                            aadict[name+'*'*extraneutrons]=aadict[name[1:]][0]+'(C)'+aadict[name[1:]][1:]
                            validaaset.add(name+'*'*extraneutrons)
                            aainconstraint.add(name+'*'*extraneutrons)
                            print('Amino acid {} found in database and Nmethylated'.format(name[1:]))
                    elif Cap(name[1:]) in validaaset:
                        #Generate nmethyl D-amino acid from L-amino acid
                        validname=Cap(name[1:])
                        if aadict[validname][1] == '(':
                            print('Amino acid {} cannot be Nmethylated, because it is already a secondary amine'.format(name[1:]))
                            databaseError=1
                        else:
                            centercnt = aadict[validname].count('[')
                            if centercnt == 0:
                                print('Unable to invert the stereochemistry of amino acid {} because no stereochemistry was detected.'.format(validname))
                                databaseError=1
                            else:
                                #Add methylation
                                oldsmiles=aadict[validname]
                                smileslist=[aadict[validname][0]+'(C)']
                                stoppoint=1
                                #Reverse stereochemistry
                                for center in range(centercnt):
                                    idx=oldsmiles[stoppoint+1:].find('[')+stoppoint+1
                                    idx2=oldsmiles[stoppoint+1:].find(']')+stoppoint+1
                                    scnt=len(oldsmiles[idx+1:idx2])-2
                                    if scnt != 1 and scnt != 2:
                                        print('Something unexpected happened when inverting the stereochemistry of amino acid {} at center {}.'.format(validname,center+1))
                                        databaseError=1
                                        break
                                    elif scnt==1:
                                        smileslist.append(oldsmiles[stoppoint:idx+3]+'@@H')
                                    elif scnt==2:
                                        smileslist.append(oldsmiles[stoppoint:idx+3]+'@H')
                                    stoppoint=idx2
                                smileslist.append(oldsmiles[stoppoint:])
                                aadict[name+'*'*extraneutrons]=''.join(smileslist)
                                validaaset.add(name+'*'*extraneutrons)
                                aainconstraint.add(name+'*'*extraneutrons)
                                print('Amino acid {} found in database and both Nmethylated and stereochemically inverted.'.format(validname))
                    else:
                        print("Error: Constraint string contains component {} not found in database.".format(name[1:]))
                        databaseError=1
                elif Cap(name) in validaaset: #breaks on 'ml'.
                    validname = Cap(name)
                    centercnt = aadict[validname].count('[')
                    if centercnt == 0:
                        print('Unable to invert the stereochemistry of amino acid {} because no stereochemistry was detected.'.format(validname))
                        databaseError=1
                    else:
                        oldsmiles=aadict[validname]
                        smileslist=[oldsmiles[0]]
                        stoppoint=0
                        #Reverse stereochemistry
                        for center in range(centercnt):
                            idx=oldsmiles[stoppoint+1:].find('[')+stoppoint+1
                            idx2=oldsmiles[stoppoint+1:].find(']')+stoppoint+1
                            scnt=len(oldsmiles[idx+1:idx2])-2
                            if scnt != 1 and scnt != 2:
                                print('Something unexpected happened when inverting the stereochemistry of amino acid {} at center {}.'.format(validname,center+1))
                                databaseError=1
                                break
                            elif scnt==1:
                                smileslist.append(oldsmiles[stoppoint:idx+3]+'@@H')
                            elif scnt==2:
                                smileslist.append(oldsmiles[stoppoint:idx+3]+'@H')
                            stoppoint=idx2
                        smileslist.append(oldsmiles[stoppoint:])
                        aadict[name+'*'*extraneutrons]=''.join(smileslist)
                        validaaset.add(name+'*'*extraneutrons)
                        aainconstraint.add(name+'*'*extraneutrons)
                        print('Amino acid {} found in database and stereochemically inverted.'.format(validname))
                else:
                    print("Error: Constraint string contains component {} not found in database.".format(name))
                    databaseError=1
            else:
                #Initialize it into this dictionary to mark it as present in the constraint string.
                aainconstraint.add(name)
    if databaseError:
        sys.exit(1)
    #Calculate exact masses for each amino acid as needed
    if verbose:
        print('Generating exact masses for amino acids present in the constraint string.')
    exmass = exactmass(aadict)
    aatomass = exmass.get(aainconstraint)
    for aa in aatomass:
        if aa in isotopedict:
            aatomass[aa]+=isotopedict[aa]*1.0062767459#Known issue: Incorrect for non-deuterium-based isotopes.
    aatoAlogP = AlogPbyresidue(aainconstraint,len(positionlist),linear,aadict)
    #Process rules
    if verbose and rules:
        print('Parsing rules {}'.format(rules))
    if rules is not None:
            operatorset = set(['>', '<', '=', '>=', '<='])
            rulelist = rules.strip().split(';')
            validsubjects = copy.copy(validaaset)
            validsubjects.add('MolWeight')
            validsubjects.add('AlogP')
            for index,rule in zip(itertools.count(),rulelist):
                rulelist[index] = rule.split()
                r = rulelist[index]
                if len(r) != 3:
                    print(r)
                    print("Error: Non-positional rule wrong number of terms.")
                    return
                elif r[0] not in validsubjects:
                    print("Error: Rule subject not defined in database.")
                    return
                elif r[1] not in operatorset:
                    print("Error: Rule has invalid operator.")
                    return
                elif not IsNumber(r[2]):
                    print("Error: Rule contains a non-numeric constraint value {}.".format(r[2]))
                    return
    else:
        rulelist = None
    return (aadict,aatomass,aatoAlogP,positionlist,rulelist,permutations)

def IsNumber (thing):
    """True if it can be converted to a float, else False.
    """
    try:
        float(thing)
    except ValueError:
        return False
    else:
        return True
    
def permute_unique (ingredients):
    """Permute with replacement, but with limited replacements per element.
    Based on the solution by rmtheis at http://stackoverflow.com/questions/6284396/permutations-with-unique-values
    """
    ingredientbag = collections.defaultdict(int)
    for x in ingredients:
        ingredientbag[x]+=1
    uniqueingredients = list(ingredientbag.keys())
    length = len(ingredients)
    return permute_unique_helper(ingredientbag,[0]*length,length-1)
    
def permute_unique_helper (ingredientbag, resultlist, depth):
    """The recursive permutation algorithm with some extra bits.
    """
    if depth < 0:
        yield tuple(resultlist)
    else:
        for ing in ingredientbag.keys():
            if ingredientbag[ing] > 0:
                resultlist[depth] = ing
                ingredientbag[ing] -= 1
                for x in permute_unique_helper(ingredientbag, resultlist, depth-1):
                    yield x
                ingredientbag[ing] += 1
                
def combinations_as_strings (elements, outlength):
    """Formats output of itertools.combinations
    """
    comblist = list(itertools.combinations_with_replacement(elements,outlength))
    newlist = []
    for t in comblist:
        newlist.append(','.join(t))
    return newlist

class Compound ():
    """A container of various attributes for a single compound composed of a sequence of single residues.
    Couldn't use collections.namedtuple because it needed to be hash-able.
    """
    def __init__ (self,name,namelist,islibmem):
        self.name = name
        self.namelist = namelist
        self.islibrarymember = islibmem
        
def islibrarymember (positionlist, word):
    """Checks whether a sequence originates from the design given.
    """
    for x,y in zip(word,positionlist):
        if x not in y:
            return False
    return True
    
def wordgen (positionlist,evalues):
    """Generates compounds as lists of amino acids names as given in positionlist.
    Generates all unqiue reorderings of library compounds as well when evalues is true.
    """
    if evalues:
        positioncount = collections.defaultdict(int)
        repeatset = set()
        for pos in positionlist:
            positioncount[','.join(sorted(pos))]+=1
        first = ','.join(sorted(positionlist[0]))
        revisedpositionlist = []
        temp = []
        for key,value in positioncount.items():
            if key == first:
                if value > 1:
                    revisedpositionlist.append(combinations_as_strings(key.split(','),value))
                else:
                    revisedpositionlist.append(key.split(','))
            else:
                if value > 1:
                    temp.append(combinations_as_strings(key.split(','),value))
                else:
                    temp.append(key.split(','))
        revisedpositionlist.extend(temp)
        for combo in itertools.product(*revisedpositionlist):
            splitcombo = ','.join(combo).split(',')
            for ordering in permute_unique(splitcombo[1:]): #Hold position one constant to prevent rotational repeats
                wordlist = list(itertools.chain([splitcombo[0]],ordering))
                word = ','.join(wordlist)
                if word in repeatset:
                    continue
                else:
                    repeatset.add(word)
                    yield Compound(word,wordlist,islibrarymember(positionlist,wordlist))
    else:
        for wordlist in itertools.product(*positionlist):
            word = ','.join(wordlist)
            yield Compound(word,wordlist,True)

def AlogPbyresidue(allowed_aas,mollength,linear,aadict):
    """Returns a dict from which molecules of the size generated can have AlogP calculated additively by residue instead of by atom.
    Currently gives erroneous results when a residue isn't an amino acid in order to preserve processing speed.
    Partial cause of Rules functionality and Truncations functionality being incompatible.
    """
    aatoAlogP = {}
    #Calculate AlogP of the appropriate ring size (or linear molecule) if completely composed of Glycine
    testmol='NCC(=O)'*mollength
    if linear:
        testmol = testmol+'O'
    else:
        testmol=testmol[0]+'1'+testmol[1:-4]+'1' +testmol[-4:]
    testmol = Chem.MolFromSmiles(testmol)
    testmol = Chem.AddHs(testmol)
    testmolAlogP = Crippen.MolLogP(testmol)
    polyGlyAlogP=testmolAlogP/mollength
    GlyAlogP=Crippen.MolLogP(Chem.AddHs(Chem.MolFromSmiles('CNCC(=O)NC')))
    #Subtract Glycine's AlogP from each allowed aa and add the cyclic/linear glycine AlogP divided by the number of units.
    for aa in allowed_aas:
        if aadict[aa].startswith('N') and aadict[aa].endswith('(=O)O'):
            aaAlogP=Crippen.MolLogP(Chem.AddHs(Chem.MolFromSmiles('C'+aadict[aa][:-1]+'NC')))
            aatoAlogP[aa]=aaAlogP-GlyAlogP+polyGlyAlogP
        elif aa == '_': #Can't ever get called as a key, but needs to be in allowed_aas for other functions.
            pass
        else: #Not an amino acid!
            #If the molecule would not allow peptide bond formation this will be inaccurate because it assumes that a peptide bond is the context in which each residue exists.
            print('Warning: {} is not an amino acid; AlogP calculations may not be correct.'.format(aa))
            #The contribution of the C+NC above reduces to negative 0.192; added this way there will be no warning message from rdkit about improper valences.
            #The result is an incorrect value for residues which cannot take part in a peptide chain, but the user should already know about those if adding them.
            aaAlogP=Crippen.MolLogP(Chem.AddHs(Chem.MolFromSmiles(aadict[aa])))-0.192 
            aatoAlogP[aa]=aaAlogP-GlyAlogP+polyGlyAlogP
    return aatoAlogP
        
def getAlogP (name, aatoAlogP):
    """Calculates the AlogP of a molecule by residue.
    """
    return sum([aatoAlogP[aa] for aa in name])

def opparse (operatorstring):
    """Parses the operator string or returns none if there is no good translation.
    """
    opdict={'=':operator.eq,'>':operator.gt,'<':operator.lt,'>=':operator.ge,'<=':operator.le}
    try:
        return opdict[operatorstring]
    except KeyError:
        return None
        
def libgen (positionlist, rulelist, aadict, aatomass, aatoAlogP, linear, evalues):
    """Interprets rules and filters the library by those rules.
    """
    filtered = 0
    alreadyseen = set()
    if rulelist is not None:
        rulesubjectset = set(rulelist[:][0])
    for compound in wordgen(positionlist,evalues):
        if '_' in compound.namelist:
            compound.namelist = [aa for aa in compound if aa != '_'] #Does not need to be done for all compounds, since there shouldn't be duplicates.
            compound.name = ','.join(compound.namelist)
            if name in alreadyseen:
                continue
            else:
                alreadyseen.add(name)
        if rulelist is not None:
            counts = collections.Counter()
            passcount = 0
            for aa in compound.name:
                counts[aa]+=1
            for subject in rulesubjectset:
                if subject[0] == '$':
                    for subname in aadict[subject]:
                        counts[subject]+= counts[subname]
            for rule in rulelist:
                subject = rule[0]
                oper = rule[1]
                number = int(rule[2])
                operation = opparse(oper)
                if subject == 'MolWeight':
                    if operation(getmass(compound.name,aatomass,linear),number):
                        passcount+=1
                        continue
                elif subject == 'AlogP': 
                    if operation(getAlogP(compound.name,aatoAlogP),number):
                        passcount+=1
                        continue
                else:
                    if operation(counts[subject],number):
                        passcount+=1
                        continue
            if passcount == len(rulelist):
                filtered+=1
                yield compound
        else:
            yield compound
    if filtered > 0:
        print ("{} compounds were kept after filtering by rules.".format(filtered))
    else:
        if rulelist is not None:
            print("No compounds satisfied all rules.")
       
def queryspectrum (spectrum, hitrecord, positionlist, rules, aadict, aatomass, aatoAlogP, options):
    """Examine specific spectrum-compound matches by inputting molecular names. Useful for troubleshooting.
    """
    precision = options.precision
    linear = options.linear
    evalues = options.evalues
    verbose = options.verbose
    print('Query mode engaged for target mass {}.\nEnter a molecule name. Ex: {}'.format(spectrum['uid'],next(libgen(positionlist,rules,aadict,aatomass,aatoAlogP,linear,evalues)).name))
    print('Type \'next\' to proceed to querying the next parent ion if one exists or \'quit\' to exit the program.')
    print('Type \'graph\' to see a visual representation of the last query with virtual fragments in blue.')
    lastquery = None
    annotatedbyname = None
    while True:
        line = input('>>>')
        line = line.strip()
        if line == 'quit':
            print('Exiting.')
            sys.exit(0)
        if line == 'next':
            return
        if line == 'graph':
            if not lastquery:
                print('A molecule must have been queried before it can be graphed.')
                continue
            import matplotlib
            from matplotlib import pyplot
            fig, ax = pyplot.subplots(1,1)
            fragmzlist = [x[1]+(random.random()-0.5)/10 for x in annotatedbyname]
            markerline, stemlines, baseline = ax.stem(fragmzlist,[max(spectrum.i)*1.1]*len(fragmzlist), markerfmt=' ', label='Fragments')
            pyplot.setp(stemlines, 'color', 'b', )
            markerline, stemlines, baseline = ax.stem(spectrum.mz,spectrum.i,markerfmt=' ',label='Spectrum')
            pyplot.setp(stemlines, 'color', 'g', 'linewidth', 3)
            pyplot.title('Simulated fragments of {} overlaid on spectrum {}'.format(lastquery,spectrum['uid']))
            pyplot.show()
            continue
        linelist = line.split(',')
        bad = False
        for aa in linelist:
            if aa not in aatomass.keys():
                print('The compound you entered was formatted incorrectly or not in the library specified by the constraint string.')
                bad = True
                break
        if bad:
            continue
        hitrecordnames = [h.name for h in hitrecord.keys()]
        
        if line not in hitrecordnames:
            print('The compound you entered was not considered because its exact mass does not match the target mass.')
            continue
        else:
            idx = hitrecordnames.index(line)
            cmpd = list(hitrecord.keys())[idx]
        lastquery = line
        nloss = neutralloss(aadict, aatomass)
        byname = getfrags(line,aatomass,linear)
        #correspondances = collections.defaultdict(set)
        byneuname = {}
        for name in byname:
            for neuname,neutral in nloss.get(line,name,byname[name]).items():
                byneuname[neuname] = neutral
                #correspondances[name].add(neuname)
        annotatedbyname = sorted([[key,value,''] for key,value in byneuname.items()],key=lambda x: x[1])
        #byname = merge_dicts(byname,byneuname)
        oorc = 0
        p = 0
        plen = len(spectrum.mz)
        for i in annotatedbyname:
            #If we are at a mass higher than the upper bound, we need to advance the bounds.
            while spectrum.mz[p]+precision[0] < i[1]:
                if p+1 < plen:
                    p+=1
                else:
                    break
            #If we are at a mass lower than the lower bound, we need to advance in mass.
            if spectrum.mz[p]-precision[0] > i[1]:
                continue
            #Are we between the bounds?
            if spectrum.mz[p]+precision[0] > i[1]:
                #add to list of hits
                i[2] = 'Match'
        minmz,maxmz = getMassrangeFromSpectrum(spectrum)
        mainfrags = 0
        bonusfrags = 0
        for frag in hitrecord[cmpd][0]:
            fragcount = hitrecord[cmpd][0][frag]
            if fragcount > 0:
                mainfrags+=1
                bonusfrags+=fragcount-1
        print('Compound Mass (uncharged): {}'.format(getmass(line, aatomass, linear)))
        print('Scan Range {}-{}, {} main ion series hits from {} total matches.'.format(minmz,maxmz,mainfrags,mainfrags+bonusfrags))
        print('Fragment Name                   Mass            Status')
        gamma=[]
        beta=[]
        for frag in annotatedbyname:
            if minmz > frag[1] or maxmz < frag[1]:
                frag[2]='Out of Range'
                oorc+=1
            if '_y' in frag[0]:
                gamma.append(frag)
            else:
                beta.append(frag)
        for frag in gamma:
            print('{:<30}\t{}\t{}'.format(frag[0],frag[1],frag[2]))
        for frag in beta:
            print('{:<30}\t{}\t{}'.format(frag[0],frag[1],frag[2]))
        if verbose:
            #print('Querymode Says: {} hits out of {} fragments tried, with {} fragments out of range.'.format(matchcount,trycount-oorc, oorc))
            pass
    
def querymode (ordereduids, spectra, hitrecords, positionlist, rules, aadict, aatomass, aatoAlogP, options):
    """Allows user to query specific hit molecules in specific spectra, prints their fragment names, masses, and which matched the spectrum.
    """
    while True:
        print("Please enter one of the following spectra names to examine sequences matched to it or \'quit\' to exit the program.")
        for uid in ordereduids:
            print (uid)
        line = input('>>>')
        line = line.strip()
        if line == 'quit':
            break
        elif line in ordereduids:
            spectrum = spectra[line]
            hitrecord = hitrecords[line]
            queryspectrum(spectrum, hitrecord, positionlist, rules, aadict, aatomass, aatoAlogP, options)
        else:
            print("User input did not match any spectrum name.")
            continue
    print('Exiting.')   
        
def scoreandprint (outfile,spectra,hitsdict,precision,evalues,verbose):
    """Compiles and outputs final scores.
    """
    spectralscores = []
    for uid in sorted(hitsdict.keys()):
        spectrum = spectra[uid]
        totI = spectrum['total ion current']
        sortedSpectrumI=sorted(spectrum.i,reverse=True)
        topcount = 2
        try:
            tops = sortedSpectrumI[:topcount]
            noToptotI = sum(sortedSpectrumI[topcount:])
        except:
            tops = []
            noToptotI = sortedSpectrumI
        hitrecord = hitsdict[uid]
        fragmatch = collections.defaultdict(int)
        bonusmatch = collections.defaultdict(int)
        intensitymatch = collections.defaultdict(list)
        peakcount=0
        if hitrecord.keys() != []:
            peakcount = len(spectrum.peaks)
            if peakcount <= 5: #Ensures the most egregious minimal data quality but high intensity-matching spectra (common) are pruned.
                continue
            for molecule in hitrecord:
                for frag in hitrecord[molecule][0]:
                    fragcount = hitrecord[molecule][0][frag]
                    if fragcount > 0:
                        bonusmatch[molecule] += fragcount-1
                        fragmatch[molecule] += 1
                sortedMatchIntensities = sorted(hitrecord[molecule][1],reverse=True)
                sortedMatchIntensities[:] = [i for i in sortedMatchIntensities if i not in tops]
                noTopI = sum(sortedMatchIntensities)
                intensitymatch[molecule] = [sum(hitrecord[molecule][1]),noTopI]
                
        finalscores = []
        for mol,frg in fragmatch.items():
            if frg != 0:
                bns = bonusmatch[mol]
                I = intensitymatch[mol][0]
                Ifrac = I/totI
                noTopI = intensitymatch[mol][1] #To catch things which match to one big peak and then just the noise.
                noTopIfrac = noTopI/noToptotI
                #The if statement below is a penalty for spectra where there are 1 to 2 really strong peaks and general noise brought along for the ride.
                #In such a case, matches to a significant portion of the noise is unlikely, thus the %intensity matching absent those two strong peaks will be much lower.
                #Predicated on the assumtption that this is the only case in which the %intensity matching absent the top two peaks is that much lower.
                score = (frg+bns/10)/10**(1-(Ifrac if (Ifrac-noTopIfrac) < 0.5 else 0))
                finalscores.append([mol,frg,bns,Ifrac*100,score])
        if len(finalscores) != 0:
            spectralscores.append([uid,peakcount,totI,sorted(finalscores,key=lambda x: x[4],reverse=True)])
    spectralscores = sorted(spectralscores,key=lambda x: x[3][0][4],reverse=True)          
    
    cwd = os.getcwd()
    compout = openpyxl.Workbook()
    humanout = openpyxl.Workbook()
    humansheet = humanout.active
    compsheet = compout.active
    humancols = ['Mass','Time (s)','Scan Range','Top Sequence','Score','Next Best Score','Number of Hits','Average Score','Sequencing Confidence']
    if evalues:
        humancols.extend(['E-value','Null Database Size','P of score or greater'])
    for i,word in zip(itertools.count(1),humancols):
        humansheet.cell(row=1,column=i,value=word)
    compcols = ['Mass','Time (s)','Scan Range','Sequence','Unique Matches','Redundant Matches','Percent Intensity Matched','Score']
    if evalues:
        compcols.extend(['E-value','Null Database Size','P of score or greater'])
    for i,word in zip(itertools.count(1),compcols):
        compsheet.cell(row=1,column=i,value=word)
    rh=2
    rc=2
    for scoreset in spectralscores:
        uid,peakcount,totI,finalscores=scoreset
        mass,time,*scanrange =uid.split(',')
        mass = float(mass)
        time = float(time)
        scanrange = ','.join(scanrange)
        samples = []
        numhits = 0
        for s in finalscores:#Grab all decoy data
            if s[0].islibrarymember:
                numhits+=1
            elif evalues:
                samples.append(s[4])
        if evalues:
            samplesize = len(samples)
            maxscorenamelen = len(finalscores[0][0].name.split(','))
            maxscore= maxscorenamelen*(maxscorenamelen-1) #Normalize from [0,max] to [0,1]
            logitx = np.asarray([math.log(x/maxscore/(1-x/maxscore)) for x in samples]) #Convert from [0,1] to R number space since the kernel operates in that space
            kde = sm.nonparametric.KDEUnivariate(logitx)
            kde.fit(cut=5)
            truesupport = np.asarray([(math.e**x/(1+math.e**x))*maxscore for x in kde.support]) #Convert back to bounded [0,max] space
            #
            #E-value debug code.
            #
            # f, axes = plt.subplots(1, 2, figsize=(7, 7))
            # axes[1] = plt.scatter(truesupport,kde.density)
            # seaborn.distplot(np.asarray(samples),rug=True,ax=axes[0])
            # plt.setp(axes)
            # plt.suptitle('Distplot on left, KDE on right, S={},mass={}'.format(samplesize,mass),size='x-large')
            # plt.tight_layout()
            # plt.show()
        first = True
        firstscore = None
        second = False
        total = 0
        for scores in sorted(finalscores,key=lambda x: x[4],reverse=True):
            if scores[0].islibrarymember:
                if evalues:
                    idx = (np.abs(truesupport-scores[4])).argmin()
                    pgreaterscore = 1-kde.cdf[idx]
                    evalue = pgreaterscore*(numhits)
                    if np.isnan(evalue):
                        evalue = 'Not Computable'
                    if np.isnan(pgreaterscore):
                        pgreaterscore = 'Not Computable'
                if second:
                    humansheet.cell(row=rh,column=humancols.index('Next Best Score')+1,value=scores[4])
                    humansheet.cell(row=rh,column=humancols.index('Sequencing Confidence')+1,value=(firstscore-scores[4])/firstscore)
                    second=False
                if first:
                    humansheet.cell(row=rh,column=humancols.index('Mass')+1,value=mass-1.0078250321)#This will need to be done differently once I reorganize things for multiple possible ionization partners.
                    humansheet.cell(row=rh,column=humancols.index('Time (s)')+1,value=time)
                    humansheet.cell(row=rh,column=humancols.index('Scan Range')+1,value=scanrange)
                    humansheet.cell(row=rh,column=humancols.index('Top Sequence')+1,value=scores[0].name)
                    humansheet.cell(row=rh,column=humancols.index('Score')+1,value=scores[4])
                    firstscore=scores[4]
                    humansheet.cell(row=rh,column=humancols.index('Number of Hits')+1,value=numhits)
                    if evalues:
                        humansheet.cell(row=rh,column=humancols.index('E-value')+1,value=evalue)
                        humansheet.cell(row=rh,column=humancols.index('Null Database Size')+1,value=samplesize)
                        humansheet.cell(row=rh,column=humancols.index('P of score or greater')+1,value=pgreaterscore)
                    compsheet.cell(row=rc,column=compcols.index('Mass')+1,value=mass)
                    compsheet.cell(row=rc,column=compcols.index('Time (s)')+1,value=time)
                    compsheet.cell(row=rc,column=compcols.index('Scan Range')+1,value=scanrange)
                    first=False
                    second=True
                compsheet.cell(row=rc,column=compcols.index('Sequence')+1,value=scores[0].name)
                compsheet.cell(row=rc,column=compcols.index('Unique Matches')+1,value=scores[1])
                compsheet.cell(row=rc,column=compcols.index('Redundant Matches')+1,value=scores[2])
                compsheet.cell(row=rc,column=compcols.index('Percent Intensity Matched')+1,value=scores[3])
                compsheet.cell(row=rc,column=compcols.index('Score')+1,value=scores[4])
                if evalues:
                    compsheet.cell(row=rc,column=compcols.index('E-value')+1,value=evalue)
                    compsheet.cell(row=rc,column=compcols.index('Null Database Size')+1,value=samplesize)
                    compsheet.cell(row=rc,column=compcols.index('P of score or greater')+1,value=pgreaterscore)
                rc+=1
                total += scores[4]
        if numhits > 0:
            humansheet.cell(row=rh,column=humancols.index('Average Score')+1,value=total/numhits)
        elif verbose:
            print('MS2 spectrum {} had no matches to library members.'.format(uid))
            humansheet.cell(row=rh,column=humancols.index('Mass')+1,value=mass)
            humansheet.cell(row=rh,column=humancols.index('Time (s)')+1,value=time)
            humansheet.cell(row=rh,column=humancols.index('Scan Range')+1,value=scanrange)
        rh+=1
    compout.save('{0}\\{1}_Out.xlsx'.format(cwd,outfile))
    humanout.save('{0}\\{1}_Results.xlsx'.format(cwd,outfile))
    return [x[0] for x in spectralscores]
    
def main ():
    #Some argument processing
    options = parse_args()
    if options.verbose:
        stime = time.time()
    
    if options.verbose:
        recipetime = time.time()
    #Load up the database into a dictionary.
    aadict,aatomass,aatoAlogP,positionlist,rules,permutations = recipebuild(options.aadatabase, options.constraint, options.rules, options.trunc, options.verbose, options.linear)
    if options.trunc:
        print("Library to compare MSMS data against will generate with {} possible compounds including all possible truncations.".format(permutations))
    else:
        print("Library to compare MSMS data against will generate with {} possible compounds.".format(permutations))
    if options.verbose:
        recipetime = time.time()-recipetime
    
    #Search for evidence of library-members and spin up/down the searching processes.
    if options.verbose:
        pstime = time.time()
    spectra = processSpectra(options.mzml, options.precision, options.scanranges, options.targets, positionlist, aatomass, options.linear, options.noiselevel, options.verbose)
    if options.verbose:
        pstime = time.time()-pstime
        
    if options.verbose:
        gentime = time.time()
    dataqueue = multiprocessing.Queue()
    resultqueue = multiprocessing.Queue()
    hitsdict = {u:{} for u in spectra.keys()}
    p = []
    for x in range(options.workernum):
        p.append(Worker(dataqueue, resultqueue, (aatomass,options.linear,options.precision,spectra,aadict)))
    if options.verbose:
        print('Spinning up {} worker processes for peakmatching calculations.'.format(options.workernum))
    for worker in p:
        worker.start()
    try:
        buffer = []
        if options.verbose:
            print('Generating library-members from constraint string and distributing them to worker processes.')
        for compound in libgen(positionlist,rules,aadict,aatomass,aatoAlogP,options.linear,options.evalues):
            buffer.append(compound)
            if len(buffer)>int(permutations/options.workernum):
                dataqueue.put(buffer)
                buffer = []
            while dataqueue.qsize() > 50:
                y = resultqueue.get()
                for uid, hits in y.items():
                    hitsdict[uid] = merge_dicts(hitsdict[uid],hits)
        if buffer != []:
            dataqueue.put(buffer)
        if options.verbose:
            gentime = time.time()-gentime
        for x in range(options.workernum):
            dataqueue.put(None)
        #
        #wait for them to finish!
        #
        if options.verbose:
            print('Waiting for worker processes to finish peakmatching calculations.')
            waittime = time.time()
        for worker in p:
            while True:
                y = resultqueue.get()
                if y == None:
                    break
                for uid, hits in y.items():
                    hitsdict[uid] = merge_dicts(hitsdict[uid],hits)
        if options.verbose:
            waittime = time.time()-waittime
    except KeyboardInterrupt:
        for worker in p:
            worker.terminate()
            worker.join()
        sys.exit(0)
    
    if options.verbose:
        scoretime = time.time()
    #Output and scoring calculations.
    ordereduids = scoreandprint(options.outfile,spectra,hitsdict,options.precision,options.evalues,options.verbose)
    if options.verbose:
        scoretime = time.time()-scoretime
        
    if options.verbose:
        ttime = time.time()-stime
        print('\n{:.2f} seconds elapsed.'.format(ttime))
        print('{:.2f} seconds elapsed during constraint processing.'.format(recipetime))
        print('{:.2f} seconds elapsed during spectral processing.'.format(pstime))
        print('{:.2f} seconds elapsed during library generation'.format(gentime))
        print('{:.2f} seconds elapsed waiting for peak matching to finish'.format(waittime))
        print('{:.2f} seconds elapsed during scoring.'.format(scoretime))
        
    #Query mode.
    if options.query:
        querymode(ordereduids,spectra,hitsdict,positionlist,rules,aadict,aatomass,aatoAlogP,options)
        
if __name__ == '__main__':
    main()
