# Grapevine: CH Sequence Analysis Pipeline

TODOs
* Tree rooting currently done based on one of the outgroup sequences, think about whether this makes sense even when the lineage is this sequence's lineage
* The code doesn't filter out uncertain dates; instead LSD handle incomplete dates by e.g. assigning the 15th for an unknown day and a known month
* The confirmed case data in the database isn't updated by the script, and the travel data is updated every single time the script is run - instead check the age of the tables and decide whether they need to be updated
* Currently the full alignment is used - re-implement site masking according to De Maio et al.

## Usage

### Expected Input Files

The script expects the following directory structure:

```
workdir
├── input
│   ├─01- reference.fasta
```

- The files (01) is the reference genome which can be downloaded from https://www.ncbi.nlm.nih.gov/nuccore/MN908947


### Required Programs

The pipeline is run in a docker (or singularity) container so all required programs are installed. The requirements are:

- R
- python3
- iqtree
- database repository (https://gitlab.ethz.ch/sars_cov_2/database)
- access to the sars_cov_2 database (sars_cov_2@id-hdb-psgr-cp61.ethz.ch)

A list of required R packages are listed in `install-packages.R`.
A list of required python packates are listed in `database/python/requirements.txt`.


### Settings

All the settings are at the top of `main.sh`.


### Output Structure

```
workdir
├── input
│   ├─01- reference.fasta
├── temp
├── output
```

The script will write into two folders: one for **temporary** files and one for the final **output** files. The script will not make any change to the **input** directory.


### Run

The script is designed to run on the Euler server. Please make sure that the required input is provided and the job will have internet access by executing `module load eth_proxy`. See https://scicomp.ethz.ch/wiki/Getting_started_with_clusters#Security for further information.

Build the singularity container with the command:

```
singularity build --docker-login grapevine.sif docker://registry.ethz.ch/sars_cov_2/grapevine:latest
# Enter credentials for ETH gitlab
```

Finally, run:

```
bsub -N -n 4 -R "rusage[mem=1000]" -W 1:00 -B "singularity run --bind /scratch:/scratch --bind <path to workdir with required input>/workdir:/app/workdir grapevine.sif"
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
