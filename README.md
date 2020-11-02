# CH Sequence Analysis Pipeline

## Usage

### Expected Input Files

The script expects an input folder with the following files:

```
input
├── est_imports
│   ├─01─ FSO_grenzgaenger_statistics_clean.csv
│   ├─02─ FSO_tourist_arrival_statistics_clean.csv
├── pangolin
│   ├─03─ consensus_data_for_release
│   │  └── ...
│   └─04─ viollier_merged_metadata.txt
├─05─ metadata_2020-XX-XX_XX-XX.tsv.gz
├─06─ reference.fasta
├─07─ reference.gb
└─08─ sequences_2020-XX-XX_XX-XX.fasta.gz
```

- The files (01) and (02) are (somehow) generated from the data from https://www.bfs.admin.ch/bfs/en/home/statistics/tourism.assetdetail.14167010.html and https://www.bfs.admin.ch/bfs/en/home/statistics/work-income.assetdetail.13647546.html.
- (03) and (04) contains the Swiss sequences that have not been released yet.
- (05) and (08) can be downloaded from GISAID.
- (06) and (07) are differently-formatted copies of the reference genome which can be downloaded from https://www.ncbi.nlm.nih.gov/nuccore/MN908947

### Fragile components

The script downsample_alignment/tally_mobility_into_switzerland.R should be updated as soon as new quarterly statistics for (01) and (02) are available. The current script extrapolates data from July 2020 to the current month and will break after Dec. 2020 as currently coded.

Figure_2.R also relies on integer weeks being weeks since 1. Jan 2020 and will break when the year changes.

infectious_pop_by_country_month.txt only has information through October.

### Required Programs

Nextstrain's [ncov repository](https://github.com/nextstrain/ncov) is needed. It should be cloned into this directory.

The following programs and commands have to be installed:

- R
- python3
- augur
- iqtree
- mafft

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
├── generate_master_alignment_from_nextfasta.sh
│   ├─01- tmp/nextdata_splits.tar.gz
|   ├─02- tmp/nextdata_alignments.tar.gz
|   ├─03- tmp/nextdata_alignments_noref.tar.gz
|   ├─04- tmp/nextdata_alignment.fasta
├──05- tmp/our_qcd_seqs.fasta
├── format_ch_data_not_yet_in_nextdata.R
│   ├─06- our_seqs_to_add_to_nextdata.fasta
│   ├─07- our_seqs_to_add_to_nextdata.tsv
├──08- tmp/our_seqs_to_add_to_nextdata_aligned.fasta
├──09- tmp/our_seqs_to_add_to_nextdata_aligned_noref.fasta
├──10- nextmeta_with_unreleased.tsv
├──11- nextdata_with_unreleased_aligned.fasta
├── qc_master_alignment.sh
│   ├─12- tmp/qc_master_alignment/include.txt
│   ├─13- tmp/qc_master_alignment/exclude.txt
│   ├─14- tmp/qc_master_alignment/alignment_filtered.fasta
│   ├─15- tmp/qc_master_alignment/problematic_sites_sarsCov2.vcf
│   ├─16- tmp/qc_master_alignment/alignment_filtered_masked.fasta
│   ├─17- tmp/qc_master_alignment/diagnostic/
│   ├─18- tmp/qc_master_alignment/alignment_filtered2_masked.fasta
│   ├─19- tmp/qc_master_alignment/priorities.txt
├── downsample_alignment/tally_mobility_into_switzerland.R 
│   ├─20- est_imports/travel_per_country_month.txt
│   ├─21- est_imports/figures/travel_per_country_month.png
├── downsample_alignment/est_avg_infectious_cases_per_country_month.R
│   ├─22- tmp/est_imports/ecdc_global_cases2020-11-01.csv
│   ├─23- est_imports/infectious_pop_by_country_month.txt
│   ├─24- iest_imports/figures/infectious_pop_by_country_month.png
│   ├─25- est_imports/figures/daily_cases_by_country.png
├── downsample_alignment/estimate_n_imports_per_country_month.R
│   ├─26- est_imports/estimated_imports_per_country_month.txt
│   ├─27- est_imports/figures/estimated_imports_per_country_month.png
├── downsample_alignment/get_n_seqs_per_country_month_based_on_imports.R
│   ├─28- est_imports/samples_per_country_w_import_padding_[0|1].txt
│   ├─29- est_imports/figures/samples_per_country_month_w_import_padding_[0|1].png
│   ├─30- est_imports/figures/samples_per_country_w_import_padding_[0|1].png
├── downsample_alignment/subsample_alignment.R
│   ├─31- alignments/$PREFIX/$PREFIX_alignment.fasta
│   ├─32- alignments/$PREFIX/$PREFIX_priority_metadata.txt
│   ├─33- alignments/$PREFIX/$PREFIX_context_metadata.txt
│   ├─34- alignments/$PREFIX/$PREFIX_tree_metadata.txt
├──35- tmp/iqtree/
├── analyze_tree/get_outgroup.R
│   ├─36- $PREFIX_outgroup_clades.png
│   ├─37- $PREFIX.outgroup.txt
├──38- tmp/lsd/
├── analyze_tree/pick_swiss_clusters.R 
│   ├─39- tmp/clusters/$PREFIX_clusters.txt
├── analyze_tree/analyze_tree/reconstruct_ancestral_locations_weighted_parsimony.R
|   ├─40- tmp/asr/$PREFIX_tree_data_with_asr.txt
├──41- clusters_varying_m
├── analyze_tree/table_cluster_stats.R
|   ├─42-  tmp/cluster_stats/$PREFIX_cluster_stats.txt
```

(01) is a directory containing input file (08 - sequences_2020-XX-XX_XX-XX.fasta) split into 150-sequence chunks for quicker alignment. 
(02) are the sequences in (01) aligned.
(03) are the alignments from (02) with the reference sequence removed.
(04) is the final, aligned version of input file (08 - sequences_2020-XX-XX_XX-XX.fasta).
(05) are all the quality-controlled sequences we have released to GISAID concatenated into a single file.
(06) are the quality-controlled sequences we have released to GISAID but which aren't included in input file (08 - sequences_2020-XX-XX_XX-XX.fasta) yet.
(07) is the metadata for sequences in (06).
(08) is (06) aligned.
(09) is (08) with the reference sequence removed.
(10) is input file (05─ metadata_2020-XX-XX_XX-XX.tsv) and (07) concatenated.
(11) is (04) and (09) concatenated.
(12) and (13) are downloaded from the Nextstrain ncov repository.
(14) is (10) with short sequences and bad sequences filtered out.
(15) is downloaded from the problematic sites repo provided by De Maio et al.
(16) is (14) with problematic sites masked out.
(17) are the results of running the Nexstrain diagnostic script.
(18) is (16) with problematic sequenced identified by the Nextstrain diagnostic script filtered out.
(19) is a list of foreign sequences genetically similar to Swiss sequences.
(20) are estiamates of the number of travel arrivals into Switzerland by month.
(22) is downloaded from the ECDC.
(23) are estimates of the average percentage of source populations infectious for each month.
(26) are estimates of the number of imports in Switzerland per month based on (20) and (23).
(28) are files containing the number of sequences to select per country and month based on (26) and a specified padding for estimated imports per month.
(31 - 34) are the alignment and metadata files for tree-building, including foreign priority & context sequence sets as well as Swiss sequences. $PREFIX=rep_X_n_sim_800_n_imports_padded_X.
(35) are the result of tree-building based on (31).
(37) a list of outgroup sequences to root the tree with $PREFIX by. (36) is a visualization of the clades included in (37).
(38) are the results of rooting and least-squares dating based on (35) and (37).
(39) are clusters defined to be Swiss transmission chains from the tree and under the assumptions specified in $PREFIX.
(40) is the tree data specified by $PREFIX with cluster defined as specified in PREFIX with ancestral states reconstructed at internal nodes.
(41) are versions of (39) with different variations on the cluster definition.
(42) are cluster summary statistics calculated from (41).
