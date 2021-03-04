# Grapevine: CH Sequence Analysis Pipeline

TODOs
* Tree rooting currently done based on one of the outgroup sequences, think about whether this makes sense even when the lineage is this sequence's lineage
* The code only filters out uncertain dates in travel context set; for other seqs LSD handles incomplete dates by e.g. assigning the 15th for an unknown day and a known month
* The confirmed case data in the database isn't updated by the script, and the travel data is updated every single time the script is run - instead check the age of the tables and decide whether they need to be updated
* Currently the unmasked alignment is used - re-implement site masking according to De Maio et al.
* Make the max weekly sampling proportion a command-line parameter rather than hardcoded 
* Finish implementing manually defined lineage splits (see branch reduce-large-lineages)
* Have program write out a log, including the parameter settings it was run with

## Usage

### Expected Input Files

The script expects the following directory structure:

```
workdir
├── input
│   ├─01- reference.fasta
│   ├─02- config.yml
```

- The files (01) is the reference genome which can be downloaded from https://www.ncbi.nlm.nih.gov/nuccore/MN908947


### Required Programs

The pipeline is run in a docker (or singularity) container so all required programs (R, python3, IQ-TREE) are installed. However, access to the sars_cov_2 database (sars_cov_2@***REMOVED***) is required.

A list of R packages used are listed in `install-packages.R`.
A list of python packates used are listed in `database/python/requirements.txt`.


### Settings

All the settings are at the top of `main.sh`.


### Output Structure

```
workdir
├── input
│   ├─01- reference.fasta
├── tmp
├── output
```

The script will write into two folders: one for **temporary** files and one for the final **output** files. The script will not make any change to the **input** directory.


### Run

The script is designed to run on the Euler server. Please make sure that the required input is provided and the job will have internet access by executing:

```
module load eth_proxy
```

See https://scicomp.ethz.ch/wiki/Getting_started_with_clusters#Security for further information.

Also make sure the `tmp/` and `output/` directories do not already exist

```
rm -r workdir/tmp
rm -r workdir/output
```

Build the singularity container with the command:

```
singularity build --docker-login grapevine.sif docker://registry.ethz.ch/sars_cov_2/grapevine:latest
# Enter credentials for ETH gitlab
```

Finally, run:

```
bsub -N -n 8 -R "rusage[mem=2048]" -W 12:00 -B "singularity run --bind /scratch:/scratch --bind /cluster/scratch/nadeaus/grapevine/workdir:/app/workdir grapevine.sif"
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
