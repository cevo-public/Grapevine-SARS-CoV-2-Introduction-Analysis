#!/bin/bash
set -euo pipefail

# ------------------------------------------------------
# Settings

# These settings should not generally be modified
WORKDIR=workdir  # this is a directory on the external computer that contains the input/ directory and is mapped into the container; because it's mapped bi-directionally, output can also be written here
REFERENCE=$WORKDIR/input/reference.fasta
CONFIG=$WORKDIR/input/grapevine_config.yaml
IQTREE=/app/iqtree-2.0.6-Linux/bin/iqtree2
PYTHON=python3

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

# ------------------------------------------------------
echo "--- Generate one alignment per pangolin lineage ---"

Rscript generate_alignments/generate_alignments.R \
    --config $CONFIG \
    --outdir $TMP_ALIGNMENTS \
    --pythonpath $PYTHON \
    --reference $REFERENCE

# ------------------------------------------------------
echo "--- Build ML trees using same default settings as Nextstrain augur ---"

for FASTA_FILE in $TMP_ALIGNMENTS/*.fasta; do
    PREFIX="$(basename "${FASTA_FILE}" | sed 's/.fasta//g')"
    $IQTREE \
        -s $FASTA_FILE \
        -m HKY+F+G4 \
        -pre $TMP_IQTREE/$PREFIX \
        -ninit 10 \
        -n 4
done
# ninit = number initial trees (1 parsimony, 1 BIONJ) which are then optimized with NNI moves
# n = number of iterations until stop
# me = log-likelihood epsilon

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
        # by default, -l is 0.5/seq_length and gives the threshold over which branches are forced to be greater than the minimum branch length
        # so by default, only branches longer than 0.5 substitutions must be greater than 0 length
done
# date-ci = number of replicates to compute confidence interval
# date-outlier = # z-score cutoff to exclude outlier nodes
# clock-sd = std-dev for lognormal relaxed clock (for uncertainty estimation)
# -a gives root date contstraints, -u gives minimum branch length, -t gives lower bound for mutation rate
# -D specifies output date format should be year-month-day

# ------------------------------------------------------
echo "--- Pick swiss transmission chains off the tree ---"

MAX_NONFOCAL_SUBCLADES=3
MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES=1
for TREEFILE in $TMP_LSD/*.nex; do
    PREFIX="$(basename "${TREEFILE}" | sed 's/.timetree.nex//g')"
    Rscript analyze_tree/pick_swiss_transmission_chains.R \
        --tree $TREEFILE \
        --metadata $TMP_ALIGNMENTS/${PREFIX}_metadata.tsv \
        --outdir $TMP_CHAINS \
        --maxtotalsubclades $MAX_NONFOCAL_SUBCLADES \
        --maxconsecutivesubclades $MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES \
        --prefix ${PREFIX}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_s_T \
        --polytomiesareswiss \
        --dontplot
done

MAX_NONFOCAL_SUBCLADES=3
MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES=1
for TREEFILE in $TMP_LSD/*.nex; do
    PREFIX="$(basename "${TREEFILE}" | sed 's/.timetree.nex//g')"
    Rscript analyze_tree/pick_swiss_transmission_chains.R \
        --tree $TREEFILE \
        --metadata $TMP_ALIGNMENTS/${PREFIX}_metadata.tsv \
        --outdir $TMP_CHAINS \
        --maxtotalsubclades $MAX_NONFOCAL_SUBCLADES \
        --maxconsecutivesubclades $MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES \
        --prefix ${PREFIX}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_s_F \
        --dontplot
done

# ------------------------------------------------------
echo "--- Infer ancestral locations via parsimony (weight locations by parsimony score) ---"

MAX_NONFOCAL_SUBCLADES=3
MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES=1

for TREEFILE in $TMP_LSD/*.nex ; do
    PREFIX="$(basename "${TREEFILE}" | sed 's/.timetree.nex//g')"
    Rscript analyze_tree/reconstruct_ancestral_locations_weighted_parsimony.R \
        --tree $TREEFILE \
        --metadata $TMP_ALIGNMENTS/${PREFIX}_metadata.tsv \
        --chains $TMP_CHAINS/${PREFIX}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_s_T_chains.txt \
        --outdir $TMP_ASR \
        --prefix ${PREFIX}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_s_T
done

for TREEFILE in $TMP_LSD/*.nex ; do
    PREFIX="$(basename "${TREEFILE}" | sed 's/.timetree.nex//g')"
    Rscript analyze_tree/reconstruct_ancestral_locations_weighted_parsimony.R \
        --tree $TREEFILE \
        --metadata $TMP_ALIGNMENTS/${PREFIX}_metadata.tsv \
        --chains $TMP_CHAINS/${PREFIX}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_s_F_chains.txt \
        --outdir $TMP_ASR \
        --prefix ${PREFIX}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_s_F
done

# ------------------------------------------------------

# Parse parameter from config file
eval $(
sed -e 's/ *#.*//g;s/:[^:\/\/]/="/g;s/$/"/g;s/ *=/=/g' $WORKDIR/input/grapevine_config.yaml |
grep 'pick_chains_under_other_criteria'
)

if [ "$pick_chains_under_other_criteria" = "true" ]
then
    echo "--- Picking chains under different chain assumptions ---"
    for MAX_NONFOCAL_SUBCLADES in 1 2 3 4; do
        for MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES in $(seq -s' ' 1 $MAX_NONFOCAL_SUBCLADES); do

          TMP_CHAINS=$WORKDIR/tmp/chains_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}
          TMP_ASR=$WORKDIR/tmp/asr_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}
          mkdir -p $TMP_CHAINS
          mkdir -p $TMP_ASR

          # Pick chains
          for TREEFILE in $TMP_LSD/*.nex; do
              PREFIX="$(basename "${TREEFILE}" | sed 's/.timetree.nex//g')"
              Rscript analyze_tree/pick_swiss_transmission_chains.R \
                  --tree $TREEFILE \
                  --metadata $TMP_ALIGNMENTS/${PREFIX}_metadata.tsv \
                  --outdir $TMP_CHAINS \
                  --maxtotalsubclades $MAX_NONFOCAL_SUBCLADES \
                  --maxconsecutivesubclades $MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES \
                  --prefix ${PREFIX}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_s_T \
                  --polytomiesareswiss \
                  --dontplot
          done

          for TREEFILE in $TMP_LSD/*.nex; do
              PREFIX="$(basename "${TREEFILE}" | sed 's/.timetree.nex//g')"
              Rscript analyze_tree/pick_swiss_transmission_chains.R \
                  --tree $TREEFILE \
                  --metadata $TMP_ALIGNMENTS/${PREFIX}_metadata.tsv \
                  --outdir $TMP_CHAINS \
                  --maxtotalsubclades $MAX_NONFOCAL_SUBCLADES \
                  --maxconsecutivesubclades $MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES \
                  --prefix ${PREFIX}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_s_F \
                  --dontplot
          done
        done
    done

    MAX_NONFOCAL_SUBCLADES=0
    MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES=0
    TMP_CHAINS=$WORKDIR/tmp/chains_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}
    TMP_ASR=$WORKDIR/tmp/asr_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}
    mkdir -p $TMP_CHAINS
    mkdir -p $TMP_ASR

    for TREEFILE in $TMP_LSD/*.nex; do
        PREFIX="$(basename "${TREEFILE}" | sed 's/.timetree.nex//g')"
        Rscript analyze_tree/pick_swiss_transmission_chains.R \
            --tree $TREEFILE \
            --metadata $TMP_ALIGNMENTS/${PREFIX}_metadata.tsv \
            --outdir $TMP_CHAINS \
            --maxtotalsubclades $MAX_NONFOCAL_SUBCLADES \
            --maxconsecutivesubclades $MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES \
            --prefix ${PREFIX}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_s_T \
            --polytomiesareswiss \
            --dontplot
    done

    for TREEFILE in $TMP_LSD/*.nex; do
        PREFIX="$(basename "${TREEFILE}" | sed 's/.timetree.nex//g')"
        Rscript analyze_tree/pick_swiss_transmission_chains.R \
            --tree $TREEFILE \
            --metadata $TMP_ALIGNMENTS/${PREFIX}_metadata.tsv \
            --outdir $TMP_CHAINS \
            --maxtotalsubclades $MAX_NONFOCAL_SUBCLADES \
            --maxconsecutivesubclades $MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES \
            --prefix ${PREFIX}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_s_F \
            --dontplot
    done
fi

# ------------------------------------------------------
echo "--- Writing out alignments for BEAST analysis ---"
Rscript analyze_clusters/get_bdsky_alignments.R \
--workdir $WORKDIR \
--maxdate 2020-11-30
# maxdate 2020-11-30 chosen to be before B.1.1.7 affects shared Re assumption

Rscript analyze_clusters/get_date_to_week_for_bdsky.R \
--outdir $WORKDIR/output/transmission_chain_alignments \
--maxdate 2020-11-30

# ------------------------------------------------------
echo "--- Generating figures ---"
# Parse parameter from config file
eval $(
sed -e 's/ *#.*//g;s/:[^:\/\/]/="/g;s/$/"/g;s/ *=/=/g' $WORKDIR/input/grapevine_config.yaml |
grep 'max_date'
)
eval $(
sed -e 's/ *#.*//g;s/:[^:\/\/]/="/g;s/$/"/g;s/ *=/=/g' $WORKDIR/input/grapevine_config.yaml |
grep 'max_sampling_fraction'
)

Rscript generate_figures/generate_figures.R \
--maxdate $max_date \
--maxsamplingfrac $max_sampling_fraction \
--workdir $WORKDIR

# ------------------------------------------------------
echo "--- Finished successfully ---"

