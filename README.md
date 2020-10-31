# CH Sequence Analysis Pipeline

## Usage

### Expected Input Files

The script expects an input folder with the following files:

```
input
├── est_imports
│   ├─01─ FSO_grenzgaenger_statistics_clean.csv
│   ├─02─ FSO_tourist_arrival_statistics_clean.csv
│   ├─03─ infectious_pop_by_country_month.txt
│   └─04─ travel_per_country_month.txt
├── pangolin
│   ├─05─ consensus_data_for_release
│   │  └── ...
│   └─06─ viollier_merged_metadata.txt
├─07─ metadata_2020-XX-XX_XX-XX.tsv.gz
├─08─ reference.fasta
├─09─ reference.gb
└─10─ sequences_2020-XX-XX_XX-XX.fasta.gz
```

- The files (01) to (04) are (somehow) generated from the data from https://www.bfs.admin.ch/bfs/en/home/statistics/tourism.assetdetail.14167010.html and https://www.bfs.admin.ch/bfs/en/home/statistics/work-income.assetdetail.13647546.html.
- (05) and (06) contains the Swiss sequences that have not been released yet.
- (07) and (10) can be downloaded from GISAID.
- (08) and (09) are differently-formatted copies of the reference genome which can be downloaded from https://www.ncbi.nlm.nih.gov/nuccore/MN908947

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
