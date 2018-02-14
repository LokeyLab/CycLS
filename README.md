CycLS: Accurate, whole-library sequencing of cyclic peptides using tandem mass spectrometry.
=====================

Chad Townsend - University of California Santa Cruz - 2018
---------------------------------------------------
*CycLS* is a program designed to identify individual cyclic peptides from a library of known design by interpreting tandem mass spectrometry data.
This allows for assays and analysis to be performed on complex mixtures containing some degree of mass redundancy without adding library design constraints.
*CycLS* is currently in press at https://doi.org/10.1016/j.bmc.2018.01.027.

The purpose of this readme is to help you install and use *CycLS* well.

### Table of Contents

-   [Installation](#installation)
-   [Usage](#usage)

    > -   [Required arguments](#required-arguments)
    > -   [Optional arguments](#optional-arguments)
    > -   [Input file preparation](#input-file-preparation)
    > -   [Interpreting output](#interpreting-output)

-   [Known Bugs and Issues](#known-bugs-and-issues)
-   [Bug Reports](#bug-reports)

### Installation:

Run the source code using Python from the command line.

#### Requirements:

The easiest way to get the packages required to run CycLS is to install the [Anaconda Python distribution](https://www.anaconda.com/) from Continuum Analytics, then install [RDkit](https://github.com/rdkit/rdkit), [peakutils](https://bitbucket.org/lucashnegri/peakutils), [openpyxl](https://bitbucket.org/openpyxl/openpyxl), [statsmodel](http://www.statsmodels.org/stable/index.html), [seaborn](https://seaborn.pydata.org/index.html), and [pymzml](https://github.com/pymzml/pymzML). 
Use the "conda install package-name" command to install all those packages except RDkit, peakutils, and pymzml. 
See the RDkit github page for installation instructions. 
Using “pip install package-name” for peakutils and pymzml is the easiest way to add those packages to an Anaconda installation, as they don't exist in the default conda channels.
The current version of CycLS has been tested with the following versions of the packages mentioned above:
Anaconda 4.3.1 
openpyxl 2.4.1 
peakutils 1.0.3 
pymzml 0.7.7 
rdkit 2016.09.4 
seaborn 0.7.1 
statsmodels 0.6.1

### Usage:

#### Required arguments:

*Targets:* Either an asterisk, signifying that CycLS is to search all spectra for library members, or a comma-separated list of protonated ion exact masses (without whitespace), if interested in only a subset of the library, must be supplied first.

Example: '647.523,831.978,745.069' would search only for MS<sup>2</sup> spectra originating from those three protonated ions.

*mzML:* The path to the mzml format data file to be analyzed must be supplied second.

*Constraint:* A string governing the library composition must be entered third. Building blocks within a position are comma-separated, with positions separated by semi-colons (potentially necessitating surrounding the constraint string in double-quotes to prevent their interpretation by the shell). As with the *Targets* field, no whitespace is permitted. Building blocks are defined in the amino acid database file, with single-letter amino acid codes implemented in the example database provided.

Example: 'L,A,D;E,Q,K;P,G,R' is a tripeptide with three possibilities for building blocks at each position.

#### Optional arguments:

*-d, --database:* Sets the name of the residue name and SMILES string database to be used for library generation, defaulting to 'aadatabase.txt'. An example database has been included in the repository. The database format is NAME\tSMILES for each line. Though it is not currently necessary for CycLS, we have rearranged the SMILES strings for each residue such that they begin with the N-terminus of the residue and end with the C-terminus of the residue (though this re-ordering is not currently necessary for *CycLS*). Note that any SMILES string can be defined as a 'residue', and some non-amino acid examples are present in the provided residue database.

Example: L-alanine has been defined as the line: 'A\tN\[C@@H\](C)C(=O)O'.

*-e:* Activates expected value calculation for the scores of each candidate molecule by generating a decoy database and scoring it. This option significantly increases run time due to the increased number of candidate evaluations necessary per spectrum. Thus far, decoy database sizes have been under 1000, leading to unreliable expected values and making testing of this functionality difficult.

*-l, --linear:* Use this flag if the library to be sequenced is linear. CycLS was first tested against linear peptides in its development, and is capable of sequencing linear libraries as well as, or better than, it sequences cyclic ones, barring dominance of unimplemented neutral loss modes.

*-n:* Sets an intensity threshold below which all peaks in the processed MS<sup>2</sup> spectra are assumed to be noise and thrown out. If above 1.0, the number provided is the intensity cutoff used. If below 1.0, the number provided is treated as the maximum probability allowable probability of a peak being due to noise. In this case, a noise threshold is automatically generated to fulfill that condition. Defaults to 100.0, which may not be a useful threshold for your mass spectroscopy system.

*-o, --out:* Sets the prefix of the output file, and defaulting to 'Sequencing'. The default value results in the two output files 'Sequencing_Out.xlsx' and 'Sequencing_Results.xlsx'. Output is given in full in the 'Out' file and summarized in 'Results' file.

*-p:* Sets the precision of MS<sup>2</sup> and MS<sup>1</sup> spectra m/z values in that order, defaulting to 0.3 and 0.02 respectively. These values may not be appropriate for HRMS set-ups, and may need adjustment for any particular mass spectroscopy system. The two values are entered comma-separated, and without spaces (Ex: '0.3,0.02').

*-q, --query:* Activates an interactive results inspection mode after completing normal operations. Allows inspection of the fragment matching between single spectra and their candidate compounds, including visualization using matplotlib. Further usage instructions are supplied interactively.

*-r, --rules:* Allows filtering of the generated library of compounds based upon position-independent constraints using a residue name, a comparison operator, and an integer value. Alternatively, 'AlogP' and 'MolWeight' may be used to compare to those properties instead, though for libraries containing non-amino acid residues, comparing to AlogP currently filters erroneously. Incompatible with the -t argument. 

Example: '-r "G < 3;R > 2;AlogP > 2.27"' would allow only compounds with less than 3 glycine residues, greater than 2 arginine residues, and AlogP greater than 2.27 to continue forward. The three components of a rule are separated by a space, with rules separated by semicolons, and the entire field surrounded by double-quotes to prevent the shell from interpreting the spaces and semicolons.

*-s, --scanranges:* Sets the minimum and maximum scan numbers to be considered, with one scan range expected per target given (and assigned to targets in order), or a single scan range if the target was '\*'.

Example: If the targets 654 and 708 had been entered and '-s 175-354,378-412' would yield assign scan range 175-354 to target 654 and scan range 378-412 to target 708. As observed here, there are no spaces and multiple scan ranges are comma-separated.

*-t, --truncate:* Enables the generation of all possible synthetic truncations of the library resultant from incomplete peptide couplings in addition to full length library members. May significantly increase run time for large libraries. Incompatible with rule usage \(*-r*\) in its current implementation.

*-u:* Sets the number of worker processes allowed to CycLS during multi-processing operations, defaulting to the number of CPU cores minus one.

*-v, --verbose:* Sets verbosity level, defaulting to zero. Prints general status announcements to the terminal at level 1 or higher.

#### Input File Preparation:

*CycLS* uses the pymzml package to read in spectra data from mzML format files. We suggest using [Proteowizard](http://proteowizard.sourceforge.net/)'s msconvert program to convert mass spectrometry data to the mzML format. We used the command line version of msconvert from Proteowizard version 3.0.10577 with some modifiers to strip unnecessary data:
>msconvert --simAsSpectra \*.raw --32 --zlib --filter “peakPicking true 1-“ --filter “zeroSamples removeExtra”

In our hands, including UV data causes pymzml version 0.7.7 to crash, and should therefore be avoided.

#### Interpreting Output:

*CycLS* outputs a Results file and an Out file, with varying metrics included. Both files include all information necessary to locate the MS<sup>2</sup> spectra used, including mass, retention time, and the scan number \(or numbers of spectra which were combined\).

**Out**

The Out file includes the score and its components \(Unique Matches, Redundant Matches, and Percent Intensity Matched\) for each candidate to each MS<sup>2</sup> spectrum. Unique Matches represents the count of initial fragments \(those present pre-neutral loss\) to which any match was found, including all neutral losses which originated from an initial fragment. 0.1 is added to the  Redundant Matches field for each match beyond the first traced back to a given initial fragment. The Percent Intensity Matched is the sum of the intensity of the peaks for which there was at least one fragment match divided by the total intensity of the spectrum. In the simple case, the equation below is followed, where M<sub>u</sub> is the count of unique matches, M<sub>r</sub> is the count of redundant matches, and I<sub>m</sub> is the fraction of intensity matched.

![equation](http://mathurl.com/ybps97qd.png)

**Results**

The Results file includes the sequence of the top candidate to each spectrum, the top score, the next best score \(if there were multiple candidates\), the average score of all candidates, the number of candidates, and any expected value-related statistics. The magnitude of a score can only be directly compared to the score of other candidates against the same spectrum with confidence; despite this, it has been observed that, in some cases, higher top scores relative to the rest of the same library can be used as an indicator of sequencing confidence. A better indicator of sequencing confidence, and one which is comparable even between libraries of differing design, is the normalized difference between the top and next best scores.

Top sequences with a higher normalized score difference are more likely to be correct, though some correctly sequenced spectra are discarded at any threshold. Low normalized score differences between similar sequences often signifies that the correct composition of the compound in question has been determined, but crucial evidence on the sequence at one or more sites is missing due to poor ionization tendency or other reasons. In such cases, the correct sequence is usually among the top three candidates for the spectrum.

### Known Bugs and Issues:

Isotopes are handled via asterisks in the constraint string and represent deuterium only: L\*\*\* represents triple-deuterated L-Leucine using the SMILES representation of L-Leucine and interpreting the asterisks on the fly. This will be changed to SMILES-based isotope interpretation in the future.

### Bug Reports:

Please submit an issue if you find a bug!

