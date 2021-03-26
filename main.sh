#!/bin/bash
set -euo pipefail

# ------------------------------------------------------
# Settings

# The user should set these settings prior to building the docker image
MIN_DATE=2020-08-01 # 2020-08-01
MAX_DATE=2020-10-01 # 2020-12-31
MIN_LENGTH=27000
MAX_SAMPLING_FRACTION=0.01
TRAVEL_CONTEXT_SCALE_FACTOR=0.5
SIMILARITY_CONTEXT_SCALE_FACTOR=0.5
TRAVEL_DATA_WEIGHTS="5,1"

# These settings should not generally be modified
WORKDIR=workdir  # this is a directory on the external computer that contains the input/ directory and is mapped into the container; because it's mapped bi-directionally, output can also be written here 
REFERENCE=$WORKDIR/input/reference.fasta
IQTREE=/app/iqtree-2.0.6-Linux/bin/iqtree2
PYTHON=python3
N_TREES=-1  # will generate all trees unless a positive integer

# All input files for this pipeline should be stored in this directory. Files in this directory will not be changed
# by this script.
INPUT_DIR=${WORKDIR}/input

# The results will be placed here. It has to be empty or non-existing at the start of this script.
OUTPUT_DIR=${WORKDIR}/output

# Here are intermediate results that may be deleted after this script has finished.
# It has to be empty or non-existing at the start of this script.
TMP_DIR=${WORKDIR}/tmp

# ------------------------------------------------------
# Check if the required files and programs are available

# Check input files
requiredFiles=(
    "${REFERENCE}"
)
for p in ${requiredFiles[@]} ; do
    if [ ! -f $p ] ; then
        echo "File not found: ${p}"
        exit 1
    fi
done

# Ensure that no files will be overwritten
foldersThatShouldBeEmpty=(
    "${OUTPUT_DIR}"
    "${TMP_DIR}"
)
for p in ${foldersThatShouldBeEmpty[@]} ; do
    if [ -d $p ] && [ "$(ls -A ${p})" ] ; then
        echo "Directory already contains files: ${p}"
        exit 1
    fi
done

# ------------------------------------------------------
# Basic preparations

# Create output and tmp folders if they do not exist
mkdir -p $TMP_DIR
mkdir -p $OUTPUT_DIR

TMP_ALIGNMENTS=$TMP_DIR/alignments
mkdir -p $TMP_ALIGNMENTS

TMP_IQTREE=$TMP_DIR/iqtree
mkdir -p $TMP_IQTREE

TMP_LSD=$TMP_DIR/lsd
mkdir -p $TMP_LSD

TMP_CHAINS=$TMP_DIR/chains
mkdir -p $TMP_CHAINS

TMP_ASR=$TMP_DIR/asr
mkdir -p $TMP_ASR

# Write out settings
echo "MIN_DATE: $MIN_DATE" > ${OUTPUT_DIR}/main_settings.txt
echo "MAX_DATE: $MAX_DATE" >> ${OUTPUT_DIR}/main_settings.txt
echo "MIN_LENGTH: $MIN_LENGTH" >> ${OUTPUT_DIR}/main_settings.txt
echo "MAX_SAMPLING_FRACTION: $MAX_SAMPLING_FRACTION" >> ${OUTPUT_DIR}/main_settings.txt
echo "TRAVEL_CONTEXT_SCALE_FACTOR: $TRAVEL_CONTEXT_SCALE_FACTOR" >> ${OUTPUT_DIR}/main_settings.txt
echo "SIMILARITY_CONTEXT_SCALE_FACTOR: $SIMILARITY_CONTEXT_SCALE_FACTOR" >> ${OUTPUT_DIR}/main_settings.txt
echo "TRAVEL_DATA_WEIGHTS: $TRAVEL_DATA_WEIGHTS" >> ${OUTPUT_DIR}/main_settings.txt 

# ------------------------------------------------------
echo "--- Generate one alignment per pangolin lineage ---"

Rscript generate_alignments/generate_alignments.R \
    --mindate $MIN_DATE \
    --maxdate $MAX_DATE \
    --minlength $MIN_LENGTH \
    --maxsamplingfrac $MAX_SAMPLING_FRACTION \
    --travelcontextscalefactor $TRAVEL_CONTEXT_SCALE_FACTOR \
    --similaritycontextscalefactor $SIMILARITY_CONTEXT_SCALE_FACTOR \
    --traveldataweights $TRAVEL_DATA_WEIGHTS \
    --outdir $TMP_ALIGNMENTS \
    --pythonpath $PYTHON \
    --reference $REFERENCE \
    --ntrees $N_TREES

# ------------------------------------------------------
echo "--- Build ML trees using same 'fast' settings as Nextstrain augur ---"

for FASTA_FILE in $TMP_ALIGNMENTS/*.fasta; do
    PREFIX="$(basename "${FASTA_FILE}" | sed 's/.fasta//g')"
    $IQTREE \
        -s $FASTA_FILE \
        -m HKY+F+G4 \
        -pre $TMP_IQTREE/$PREFIX \
        -ninit 2 \
        -n 2 \
        -me 0.05
done

# ------------------------------------------------------
echo "--- Date rooted trees with LSD implemented in IQTREE: root defined by outgroup EPI_ISL_406798|2019-12-26 ---"

# Keep identical sequences for LSD because otherwise IQ-TREE throws them out and then complains it can't find all the outgroup seqs
for FASTA_FILE in $TMP_ALIGNMENTS/*.fasta; do
    DATADIR="$(dirname "${FASTA_FILE}")"
    PREFIX="$(basename "${FASTA_FILE}" | sed 's/.fasta//g')"
    TREEFILE=$TMP_IQTREE/${PREFIX}.treefile
    OUTGROUP="EPI_ISL_406798|2019-12-26"
    $IQTREE \
        -s $FASTA_FILE \
        -te $TREEFILE \
        -m HKY+F+G4 \
        -keep-ident \
        -o $OUTGROUP \
        -ntmax 4 \
        --date TAXNAME \
        --date-ci 100 \
        --date-outlier 3 \
        --clock-sd 0.4 \
        --date-options "-a b(2019.872,2019.98) -u 0 -t 0.0008" \
        -pre $TMP_LSD/$PREFIX
done

# ------------------------------------------------------
echo "--- Pick swiss transmission chains off the tree ---"

MAX_NONFOCAL_SUBCLADES=3
MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES=1
for TREEFILE in $TMP_LSD/*.nex; do
    PREFIX="$(basename "${TREEFILE}" | sed 's/.timetree.nex//g')"
    Rscript analyze_tree/pick_swiss_transmission_chains.R \
        --tree $TREEFILE \
        --metadata $TMP_ALIGNMENTS/${PREFIX}_metadata.csv \
        --outdir $TMP_CHAINS \
        --maxtotalsubclades $MAX_NONFOCAL_SUBCLADES \
        --maxconsecutivesubclades $MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES \
        --prefix ${PREFIX}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_s_T \
        --polytomiesareswiss
done 

MAX_NONFOCAL_SUBCLADES=3
MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES=1
for TREEFILE in $TMP_LSD/*.nex; do
    PREFIX="$(basename "${TREEFILE}" | sed 's/.timetree.nex//g')"
    Rscript analyze_tree/pick_swiss_transmission_chains.R \
        --tree $TREEFILE \
        --metadata $TMP_ALIGNMENTS/${PREFIX}_metadata.csv \
        --outdir $TMP_CHAINS \
        --maxtotalsubclades $MAX_NONFOCAL_SUBCLADES \
        --maxconsecutivesubclades $MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES \
        --prefix ${PREFIX}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_s_F
done

# ------------------------------------------------------
echo "--- Infer ancestral locations via parsimony (weight locations by parsimony score) ---"

MAX_NONFOCAL_SUBCLADES=3
MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES=1

for TREEFILE in $TMP_LSD/*.nex ; do
    PREFIX="$(basename "${TREEFILE}" | sed 's/.timetree.nex//g')"
    Rscript analyze_tree/reconstruct_ancestral_locations_weighted_parsimony.R \
        --tree $TREEFILE \
        --metadata $TMP_ALIGNMENTS/${PREFIX}_metadata.csv \
        --chains $TMP_CHAINS/${PREFIX}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_s_T_chains.txt \
        --outdir $TMP_ASR \
        --prefix ${PREFIX}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_s_T 
done

for TREEFILE in $TMP_LSD/*.nex ; do 
    PREFIX="$(basename "${TREEFILE}" | sed 's/.timetree.nex//g')"
    Rscript analyze_tree/reconstruct_ancestral_locations_weighted_parsimony.R \
        --tree $TREEFILE \
        --metadata $TMP_ALIGNMENTS/${PREFIX}_metadata.csv \
        --chains $TMP_CHAINS/${PREFIX}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_s_F_chains.txt \
        --outdir $TMP_ASR \
        --prefix ${PREFIX}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_s_F
done

# ------------------------------------------------------
# echo "--- Generating figures ---"
# TODO

# ------------------------------------------------------
echo "--- Finished successfully ---"

