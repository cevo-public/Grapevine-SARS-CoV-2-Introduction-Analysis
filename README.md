# Grapevine: SARS-CoV-2 Sequence Analysis Pipeline

This pipeline performs a phylogenetic analysis that sub-samples and quality filters available GISAID data, classifies data into introductions/transmission chains, generates figures and introduction-annotated alignments, and (Swiss data only, not extensively used/supported) estimates ancestral locations.

## Usage

### Expected Input Files

The script expects the following directory structure:

```
workdir
├── input
│   ├─01- reference.fasta
│   ├─02- config.yml
│   ├─03- grapevine_config.yml
```

- The file (01) is the reference genome which can be downloaded from https://www.ncbi.nlm.nih.gov/nuccore/MN908947
- The file (02) is a configuration file giving database connection info. See template in `example_workdir/input/config.yml`
- The file (03) is a configuration file giving phylogenetic analysis specifications. See template in `example_workdir/input/grapevine_config.yml`


### Required Programs

The pipeline is run in a docker (or singularity) container so all required programs (R, python3, IQ-TREE) are installed. However, access to the sars_cov_2 database (sars_cov_2@***REMOVED***) is required.

A list of R packages used are listed in `install-packages.R`.
A list of python packages used are listed in `database/python/requirements.txt`.


### Settings

All the settings are in grapevine_config.yml.


### Output Structure

```
workdir
├── tmp
├── output
```

The script will write into two folders: tmp for intermediate files like alignments, trees, and introductions and one for the final output files like alignments for phylodynamic analysis and figures. The script will not make any change to the input directory.


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
bsub -N -n 8 -R "rusage[mem=4096]" -W 12:00 -o $WORKDIR/${WORKDIR}_%J.log -B "singularity run --bind /cluster/scratch/nadeaus/grapevine/${WORKDIR}:/app/workdir ${CONTAINER_NAME}.sif"
```

### Pipeline structure
```
main.sh
├── generate_alignments/generate_alignments.R
│   ├─01- alignments/$PREFIX.fasta
│   ├─02- alignments/$PREFIX_metadata.csv
├──03- tmp/iqtree/
├──04- tmp/lsd/
├── analyze_tree/pick_transmission_chains.R 
│   ├─05- tmp/chains/$PREFIX_chains.txt
├── analyze_tree/analyze_tree/reconstruct_ancestral_locations_weighted_parsimony.R
|   ├─06- tmp/asr/$PREFIX_tree_data_with_asr.txt
```

- (01 - 02) are the alignment and metadata files for tree-building.
- (03) are the result of tree-building based on (01).
- (04) are the results of rooting and least-squares dating based on (03) with the specified outgroup.
- (05) are transmission chains picked from the tree under the assumptions specified in $PREFIX.
- (06) is the tree data specified by $PREFIX with transmission chains defined as specified in PREFIX with ancestral states reconstructed at internal nodes according to parsimony scores.

### Notes

- If sampling by canton is enabled, but some case counts are only attributed at the country level, samples are taken at random from the pool of not-yet sampled sequences irregardless of canton.
- The column "n_leftover_seqs" in the output {country}\_downsampling\_data.csv is only calculated for sampling at the country level, and for countries other than Switzerland where sampling by canton is not enabled, should equal the column "n_seqs_total". 
