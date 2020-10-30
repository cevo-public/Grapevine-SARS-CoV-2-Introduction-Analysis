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
- (08) and (09) are just there.


### Required Programs

Nextstrain's [ncov repository](https://github.com/nextstrain/ncov) is needed. It should be cloned into this directory.

The following programs and commands have to be installed:

- R
- python3
- augur
- iqtree
- mafft

The paths to iqtree and mafft can be changed in `main.sh`.

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
