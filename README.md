# Grapevine: CH Sequence Analysis Pipeline

TODOs
* Tree rooting no longer necessary b/c have defined outgroup?
* Does LSD handle incomplete dates correctly? The resulting tree data has date on 15th for unknown day
* Also updating case data every time? Looks like no.

## Usage

### Expected Input Files

The script expects an input folder with the following files:

```
input
└─01- reference.fasta
```

- The files (01) is the reference genome which can be downloaded from https://www.ncbi.nlm.nih.gov/nuccore/MN908947


### Required Programs

The following programs and commands have to be installed:

- R
- python3
- iqtree

The paths to iqtree and mafft can be changed in `main.sh`.
R and python are already installed on Euler and can be loaded with the commands:
env2lmod
module load r/4.0.2
module load python/3.7.4

A list of required R packages are listed in `check-packages.R`.


### Settings

All the settings are at the top of `main.sh`.


### Output Structure

The script will write into two folders: one for **temporary** files and one for the final **output** files. Their paths can be changed in the settings. The script will not make any change to the **input** directory.


### Run

The script is designed to run on the Euler server. Please make sure that all required modules are loaded and programs and commands can be found.

The batch jobs have to be able to access the network. Execute `module load eth_proxy` to allow that. See https://scicomp.ethz.ch/wiki/Getting_started_with_clusters#Security for further information.

Since the script needs a long time to run (approx. 6 to 12h), it is advised to use a terminal multiplexer such as `screen` or `tmux`.

Finally, run:

```
bash main.sh
```

### Pipeline structure
```
main.sh
├── generate_alignments/generate_alignments.R
│   ├─31- alignments/$PREFIX.fasta
│   ├─32- alignments/$PREFIX_metadata.csv
├──35- tmp/iqtree/
├──38- tmp/lsd/
├── analyze_tree/pick_swiss_clusters.R 
│   ├─39- tmp/clusters/$PREFIX_chains.txt
├── analyze_tree/analyze_tree/reconstruct_ancestral_locations_weighted_parsimony.R
|   ├─40- tmp/asr/$PREFIX_tree_data_with_asr.txt
```

- (31 - 32) are the alignment and metadata files for tree-building.
- (35) are the result of tree-building based on (31).
- (38) are the results of rooting and least-squares dating based on (35) hardcoded outgroup.
- (39) are transmission chains picked from the tree under the assumptions specified in $PREFIX.
- (40) is the tree data specified by $PREFIX with transmission chains defined as specified in PREFIX with ancestral states reconstructed at internal nodes according to parsimony scores.
