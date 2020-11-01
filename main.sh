#!/bin/bash
set -euo pipefail

# ------------------------------------------------------
# Settings

ANALYSIS=2020-10-30_analysis
WORKDIR=${SCRATCH}/ch_sequencing/${ANALYSIS}

NEXTDATA_GZ_FN=sequences_2020-10-28_07-22.fasta.gz
NEXTMETA_GZ_FN=metadata_2020-10-28_07-21.tsv.gz
REFERENCE_FN=reference.fasta
REFERENCE_NCOV_FN=reference.gb
MAX_DATE=2020-10-29
IQTREE=~/lib/iqtree-2.1.2-Linux/bin/iqtree2
MAFFT=~/lib/mafft-linux64/mafft.bat

# If an email address is provided, a notification mail will be sent once the script has finished successfully.
# Attention: Currently, it might not send an email if an error occurs!
NOTIFICATION_EMAIL=


# All input files for this pipeline should be stored in this directory. Files in this directory will not be changed
# by this script.
INPUT_DIR=${WORKDIR}/input

# The results will be placed here. It has to be empty or non-existing at the start of this script.
OUTPUT_DIR=${WORKDIR}/output

# Here are intermediate results that may be deleted after this script has finished.
# It has to be empty or non-existing at the start of this script.
TMP_DIR=${WORKDIR}/tmp

# The path to the directory with the scripts.
SCRIPT_DIR=$(dirname $(realpath $0))

# The path to the nextstrain/ncov repository.
NCOV_DIR=${SCRIPT_DIR}/ncov


# ------------------------------------------------------
# Check if the required files and programs are available

# Check input files
requiredFiles=(
    "${INPUT_DIR}/${NEXTDATA_GZ_FN}"
    "${INPUT_DIR}/${NEXTMETA_GZ_FN}"
    "${INPUT_DIR}/${REFERENCE_FN}"
    "${INPUT_DIR}/${REFERENCE_NCOV_FN}"
    "$INPUT_DIR/est_imports/FSO_tourist_arrival_statistics_clean.csv"
    "$INPUT_DIR/est_imports/FSO_grenzgaenger_statistics_clean.csv"
    "$INPUT_DIR/est_imports/infectious_pop_by_country_month.txt"
    "$INPUT_DIR/est_imports/travel_per_country_month.txt"
    "$INPUT_DIR/pangolin/viollier_merged_metadata.txt"
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

# Ensure that the required programs can be found
requiredPrograms=(
    "Rscript"
    "python3"
    "augur"
    "bsub"
    $IQTREE
    $MAFFT
)
for p in ${requiredPrograms[@]}; do
    if ! command -v $p &> /dev/null ; then
        echo "Program not found: ${p}"
        exit 1
    fi
done

# Check R packages
if ! command Rscript ${SCRIPT_DIR}/check-packages.R ; then
    exit 1
fi


# ------------------------------------------------------
# Basic preparations

# Useful variables
LOG_DIR=$OUTPUT_DIR/logs
NEXTDATA_FN=${NEXTDATA_GZ_FN%.gz}
NEXTMETA_FN=${NEXTMETA_GZ_FN%.gz}

# Create output and tmp folders if they do not exist
mkdir -p $LOG_DIR
mkdir -p $TMP_DIR


# ------------------------------------------------------
echo "--- Align all of nextfasta download ---"

bsub -K -n 16 -R "rusage[mem=2048]" -o $LOG_DIR/lsf-job.%J.generate_master_alignment_from_nextfasta.log "\
bash ${SCRIPT_DIR}/generate_master_alignment_from_nextfasta.sh \
    -i $INPUT_DIR \
    -t $TMP_DIR \
    -f $NEXTDATA_FN \
    -m $NEXTMETA_FN \
    -r $REFERENCE_FN \
    -x $MAFFT \
    -n 20
"


# ------------------------------------------------------
echo "--- Get not-yet-released Swiss seqs to add to alignment ---"

cat ${INPUT_DIR}/pangolin/consensus_data_for_release/*/ch_filtered.fasta >> $TMP_DIR/our_qcd_seqs.fasta

bsub -K -R "rusage[mem=16384]" -o $LOG_DIR/lsf-job.%J.format_ch_data_not_yet_in_nextdata.log "\
Rscript $SCRIPT_DIR/format_ch_data_not_yet_in_nextdata.R \
    --ourseqsreleased $TMP_DIR/our_qcd_seqs.fasta \
    --ourseqsmetadata $INPUT_DIR/pangolin/viollier_merged_metadata.txt \
    --nextmeta $TMP_DIR/$NEXTMETA_FN \
    --outdir $TMP_DIR \
    --utilityfunctions $SCRIPT_DIR/utility_functions.R
"


# ------------------------------------------------------
echo "--- Add not-yet-released Swiss seqs to alignment ---"

# Align new sequences
FILE=$TMP_DIR/our_seqs_to_add_to_nextdata.fasta
FN=$(basename -- "$FILE")
FN="${FN%.fasta}"
bsub -K -R "rusage[mem=16384]" -o $LOG_DIR/lsf-job.%J.mafft.log "\
mafft \
    --addfragments $FILE \
    --keeplength \
    --auto \
    $INPUT_DIR/$REFERENCE_FN > $TMP_DIR/our_seqs_to_add_to_nextdata_aligned.fasta
"

# Remove reference from alignment
awk 'BEGIN { RS = ">";FS = "\n" } {if (NR>2) {print ">"$0}}' $TMP_DIR/our_seqs_to_add_to_nextdata_aligned.fasta > $TMP_DIR/our_seqs_to_add_to_nextdata_aligned_noref.fasta

# Concatenate data
cat $TMP_DIR/$NEXTMETA_FN $TMP_DIR/our_seqs_to_add_to_nextdata.tsv >> $TMP_DIR/nextmeta_with_unreleased.tsv
cat $TMP_DIR/nextdata_alignment.fasta $TMP_DIR/our_seqs_to_add_to_nextdata_aligned_noref.fasta >> $TMP_DIR/nextdata_with_unreleased_aligned.fasta


# ------------------------------------------------------
echo "--- QC alignment ---"

TMP_QC=$TMP_DIR/qc_master_alignment
mkdir -p $TMP_QC

bsub -K -n 2 -R "rusage[mem=32768]" -o $LOG_DIR/lsf-job.%J.qc_master_alignment.log "\
bash ${SCRIPT_DIR}/qc_master_alignment.sh \
    -a $TMP_DIR/nextdata_with_unreleased_aligned.fasta \
    -m $TMP_DIR/nextmeta_with_unreleased.tsv \
    -t $TMP_QC \
    -d $MAX_DATE \
    -s $SCRIPT_DIR/mask_alignment_using_vcf.py \
    -r $INPUT_DIR/$REFERENCE_NCOV_FN \
    -n $NCOV_DIR
"

# ------------------------------------------------------
echo "--- Estimate the # imports per country-month ---"

TMP_EST_IMPORTS=$TMP_DIR/est_imports
mkdir -p $TMP_EST_IMPORTS
mkdir -p $TMP_EST_IMPORTS/figures

bsub -K -o $LOG_DIR/lsf-job.%J.tally_mobility_into_switzerland.log "\
Rscript $SCRIPT_DIR/downsample_alignment/tally_mobility_into_switzerland.R \
    --tourists $INPUT_DIR/est_imports/FSO_tourist_arrival_statistics_clean.csv \
    --commuters $INPUT_DIR/est_imports/FSO_grenzgaenger_statistics_clean.csv \
    --outdirdata $TMP_EST_IMPORTS \
    --outdirfigs $TMP_EST_IMPORTS/figures
"

bsub -K -o $LOG_DIR/lsf-job.%J.est_avg_infectious_cases_per_country_month.log "\
Rscript $SCRIPT_DIR/downsample_alignment/est_avg_infectious_cases_per_country_month.R \
    --casedatalink https://opendata.ecdc.europa.eu/covid19/casedistribution/csv \
    --outdirdata $TMP_EST_IMPORTS \
    --outdirfigs $TMP_EST_IMPORTS/figures
"

bsub -K -R "rusage[mem=2048]" -o $LOG_DIR/lsf-job.%J.estimate_n_imports_per_country_month.log "\
Rscript $SCRIPT_DIR/downsample_alignment/estimate_n_imports_per_country_month.R \
    --infectiouspopdata $INPUT_DIR/est_imports/infectious_pop_by_country_month.txt \
    --arrivaldata $INPUT_DIR/est_imports/travel_per_country_month.txt \
    --prioritydata $TMP_QC/priorities.txt \
    --metadata $TMP_DIR/nextmeta_with_unreleased.tsv \
    --swissseqs $TMP_QC/swiss_alignment_filtered2_masked_oneline.fasta \
    --outdirdata $TMP_EST_IMPORTS \
    --outdirfigs $TMP_EST_IMPORTS/figures
"


# ------------------------------------------------------
echo "--- Downsample master alignment for priority, context, Swiss seqs ---"

TMP_ALIGNMENTS=$TMP_DIR/alignments
N_CONTEXT_SEQS=1500  # actually results in ~1000 context seqs
N_MOST_SIMILAR_SEQS=1000

mkdir -p $TMP_EST_IMPORTS

for PADDING in 0 1; do
    bsub -K -R "rusage[mem=16384]" -o $LOG_DIR/lsf-job.%J.get_n_seqs_per_country_month_based_on_imports.log "\
    Rscript $SCRIPT_DIR/downsample_alignment/get_n_seqs_per_country_month_based_on_imports.R \
        --importsdata $TMP_EST_IMPORTS/estimated_imports_per_country_month.txt \
        --metadata $TMP_EST_IMPORTS/metadata_all.txt \
        --padding $PADDING \
        --approxncontextseqs $N_CONTEXT_SEQS \
        --outdirdata $TMP_EST_IMPORTS \
        --outdirfigs $TMP_EST_IMPORTS/figures \
        --maxdate $MAX_DATE
    " &
done
wait

PADDING=0
for REP in 1 2 3; do
    ALN_PREFIX=rep_${REP}_n_sim_${N_MOST_SIMILAR_SEQS}_n_imports_padded_${PADDING}
    bsub -K -R "rusage[mem=16384]" -o $LOG_DIR/lsf-job.%J.subsample_alignment.log "\
    Rscript $SCRIPT_DIR/downsample_alignment/subsample_alignment.R \
        --metadata $TMP_EST_IMPORTS/metadata_all.txt \
        --alignment $TMP_QC/alignment_filtered2_masked_oneline.fasta \
        --nsimseqs $N_MOST_SIMILAR_SEQS \
        --prefix $ALN_PREFIX \
        --ncontextsamples $TMP_EST_IMPORTS/samples_per_country_w_import_padding_${PADDING}.txt \
        --outdir $TMP_ALIGNMENTS/${ALN_PREFIX} \
        --maxdate $MAX_DATE
    " &
done
wait

PADDING=1
for REP in 1; do
    ALN_PREFIX=rep_${REP}_n_sim_${N_MOST_SIMILAR_SEQS}_n_imports_padded_${PADDING}
    bsub -K -R "rusage[mem=16384]" -o $LOG_DIR/lsf-job.%J.subsample_alignment.log "\
    Rscript $SCRIPT_DIR/downsample_alignment/subsample_alignment.R \
        --metadata $TMP_EST_IMPORTS/metadata_all.txt \
        --alignment $TMP_QC/alignment_filtered2_masked_oneline.fasta \
        --nsimseqs $N_MOST_SIMILAR_SEQS \
        --prefix $ALN_PREFIX \
        --ncontextsamples $TMP_EST_IMPORTS/samples_per_country_w_import_padding_${PADDING}.txt \
        --outdir $TMP_ALIGNMENTS/${ALN_PREFIX} \
        --maxdate $MAX_DATE
    " &
done
wait


# ------------------------------------------------------
echo "--- Build ML trees ---"

# Make directory on euler

TMP_IQTREE=$TMP_DIR/iqtree
mkdir -p $TMP_IQTREE

cp $TMP_ALIGNMENTS/*/*.fasta $TMP_IQTREE

for FASTA_FILE in $TMP_IQTREE/*.fasta; do
    OUTDIR="$(dirname "${FASTA_FILE}")"
    PREFIX="$(basename "${FASTA_FILE}" | sed 's/_alignment.fasta//g')"
    bsub -K -n 16 -R "rusage[mem=1024]" -W 10:00 -o $LOG_DIR/lsf-job.%J.iqtree.log "
    $IQTREE \
        -s $FASTA_FILE \
        -m HKY+F+G4 \
        -nt 16 \
        -pre $OUTDIR/$PREFIX
    " &
done
wait


#------------------------------------------------------
echo "--- Get outgroup to root trees by based on 2 defined outgroup seqs ---"

for TREEFILE in $TMP_IQTREE/*.treefile; do
    PREFIX="$(basename "${TREEFILE}" | sed 's/.treefile//g')"
    bsub -K -o $LOG_DIR/lsf-job.%J.get_outgroup.log "\
    Rscript $SCRIPT_DIR/analyze_tree/get_outgroup.R \
        --treefile $TREEFILE \
        --metadata $TMP_ALIGNMENTS/${PREFIX}/${PREFIX}_tree_metadata.txt \
        --outgroup \"Wuhan_Hu-1_2019|EPI_ISL_402125|2019-12-26, Wuhan_WH01_2019|EPI_ISL_406798|2019-12-26\" \
        --outdir $TMP_IQTREE
    " &
done
wait


# ------------------------------------------------------
echo "--- Date rooted trees with LSD: root defined by outgroup ---"


TMP_LSD=$TMP_DIR/lsd
mkdir -p $TMP_LSD

# Keep identical sequences for LSD because otherwise IQ-TREE throws them out and then complains it can't find all the outgroup seqs
for FASTA_FILE in $TMP_IQTREE/*.fasta ; do
    DATADIR="$(dirname "${FASTA_FILE}")"
    PREFIX="$(basename "${FASTA_FILE}" | sed 's/_alignment.fasta//g')"
    TREEFILE=$DATADIR/${PREFIX}.treefile
    OUTGROUP=`cat $TMP_IQTREE/${PREFIX}.treefile.outgroup.txt`
    bsub -n 12 -K -R "rusage[mem=1024]" -o $LOG_DIR/lsf-job.%J.iqtree.log "
    $IQTREE \
        -s $FASTA_FILE \
        -te $TREEFILE \
        -m HKY+F+G4 \
        -keep-ident \
        -o $OUTGROUP \
        -nt AUTO \
        -ntmax 12 \
        --date TAXNAME \
        --date-ci 100 \
        --date-outlier 3 \
        --clock-sd 0.4 \
        --date-options \"-a b(2019.872,2019.98) -u 0 -t 0.0008\" \
        -pre $TMP_LSD/$PREFIX
    " &
done
wait


# ------------------------------------------------------
echo "--- Pick swiss clusters off the tree ---"

TMP_CLUSTERS=$TMP_DIR/clusters
mkdir -p $TMP_CLUSTERS

MAX_NONFOCAL_SUBCLADES=3
MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES=1
for TREEFILE in $TMP_LSD/*.nex; do
    PREFIX="$(basename "${TREEFILE}" | sed 's/.timetree.nex//g')"
    bsub -K -o $LOG_DIR/lsf-job.%J.pick_swiss_clusters.log "
    Rscript $SCRIPT_DIR/analyze_tree/pick_swiss_clusters.R \
        --tree $TREEFILE \
        --metadata $TMP_ALIGNMENTS/$PREFIX/${PREFIX}_tree_metadata.txt \
        --outdir $TMP_CLUSTERS \
        --maxtotalsubclades $MAX_NONFOCAL_SUBCLADES \
        --maxconsecutivesubclades $MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES \
        --prefix ${PREFIX}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_swiss_polytomies_T \
        --utilityfunctions $SCRIPT_DIR/utility_functions.R \
        --polytomiesareswiss
    " &
done

MAX_NONFOCAL_SUBCLADES=3
MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES=1
for TREEFILE in $TMP_LSD/*.nex; do
    PREFIX="$(basename "${TREEFILE}" | sed 's/.timetree.nex//g')"
    bsub -K -o $LOG_DIR/lsf-job.%J.pick_swiss_clusters.log "
    Rscript $SCRIPT_DIR/analyze_tree/pick_swiss_clusters.R \
        --tree $TREEFILE \
        --metadata $TMP_ALIGNMENTS/$PREFIX/${PREFIX}_tree_metadata.txt \
        --outdir $TMP_CLUSTERS \
        --maxtotalsubclades $MAX_NONFOCAL_SUBCLADES \
        --maxconsecutivesubclades $MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES \
        --prefix ${PREFIX}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_swiss_polytomies_F \
        --utilityfunctions $SCRIPT_DIR/utility_functions.R
    " &
done

wait


# ------------------------------------------------------
echo "--- Infer ancestral locations via parsimony (weight locations by parsimony score) ---"

TMP_ASR=$TMP_DIR/asr
mkdir -p $TMP_ASR

MAX_NONFOCAL_SUBCLADES=3
MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES=1

for TREEFILE in $TMP_LSD/*.nex ; do
    PREFIX="$(basename "${TREEFILE}" | sed 's/.timetree.nex//g')"
    bsub -K -o $LOG_DIR/lsf-job.%J.reconstruct_ancestral_locations_weighted_parsimony.log "
    Rscript $SCRIPT_DIR/analyze_tree/reconstruct_ancestral_locations_weighted_parsimony.R \
        --tree $TREEFILE \
        --metadata $TMP_ALIGNMENTS/$PREFIX/${PREFIX}_tree_metadata.txt \
        --contextmetadata $TMP_ALIGNMENTS/$PREFIX/${PREFIX}_context_metadata.txt \
        --clusters $TMP_CLUSTERS/${PREFIX}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_swiss_polytomies_T_clusters.txt \
        --outdir $TMP_ASR \
        --prefix ${PREFIX}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_swiss_polytomies_T \
        --utilityfunctions $SCRIPT_DIR/utility_functions.R
    " &
done

for TREEFILE in $TMP_LSD/*.nex ; do
    PREFIX="$(basename "${TREEFILE}" | sed 's/.timetree.nex//g')"
    bsub -K -R "rusage[mem=8192]" -o $LOG_DIR/lsf-job.%J.reconstruct_ancestral_locations_weighted_parsimony.log "
    Rscript $SCRIPT_DIR/analyze_tree/reconstruct_ancestral_locations_weighted_parsimony.R \
        --tree $TREEFILE \
        --metadata $TMP_ALIGNMENTS/$PREFIX/${PREFIX}_tree_metadata.txt \
        --contextmetadata $TMP_ALIGNMENTS/$PREFIX/${PREFIX}_context_metadata.txt \
        --clusters $TMP_CLUSTERS/${PREFIX}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_swiss_polytomies_F_clusters.txt \
        --outdir $TMP_ASR \
        --prefix ${PREFIX}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_swiss_polytomies_F \
        --utilityfunctions $SCRIPT_DIR/utility_functions.R
    " &
done

wait


# ------------------------------------------------------
echo "--- Investigate clade sensitivity to # exports allowed ---"

TMP_CLUSTERS_VARYING=$TMP_DIR/clusters_varying_m
mkdir -p $TMP_CLUSTERS_VARYING

TREEFILE=$TMP_LSD/rep_1_n_sim_1000_n_imports_padded_0.timetree.nex
PREFIX="$(basename "${TREEFILE}" | sed 's/.timetree.nex//g')"
for MAX_NONFOCAL_SUBCLADES in 1 2 3 4 ; do
    for MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES in $(seq -s' ' 1 $MAX_NONFOCAL_SUBCLADES); do
        bsub -K -o $LOG_DIR/lsf-job.%J.pick_swiss_clusters.log "
        Rscript $SCRIPT_DIR/analyze_tree/pick_swiss_clusters.R \
            --tree $TREEFILE \
            --metadata $TMP_ALIGNMENTS/$PREFIX/${PREFIX}_tree_metadata.txt \
            --outdir $TMP_CLUSTERS_VARYING \
            --maxtotalsubclades $MAX_NONFOCAL_SUBCLADES \
            --maxconsecutivesubclades $MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES \
            --prefix ${PREFIX}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_swiss_polytomies_F \
            --utilityfunctions $SCRIPT_DIR/utility_functions.R
        " &
    done
done
wait

MAX_NONFOCAL_SUBCLADES=0
MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES=0
bsub -K -o $LOG_DIR/lsf-job.%J.pick_swiss_clusters.log "
Rscript $SCRIPT_DIR/analyze_tree/pick_swiss_clusters.R \
    --tree $TREEFILE \
    --metadata $TMP_ALIGNMENTS/$PREFIX/${PREFIX}_tree_metadata.txt \
    --outdir $TMP_CLUSTERS_VARYING \
    --maxtotalsubclades $MAX_NONFOCAL_SUBCLADES \
    --maxconsecutivesubclades $MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES \
    --prefix ${PREFIX}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_swiss_polytomies_F \
    --utilityfunctions $SCRIPT_DIR/utility_functions.R
"

N_SWISS_SEQS=$(grep "^>" $TMP_QC/swiss_alignment_filtered2_masked_oneline.fasta | wc -l)

TMP_CLUSTERS_STATS=$TMP_DIR/clusters_stats
mkdir -p $TMP_CLUSTERS_STATS

for MAX_NONFOCAL_SUBCLADES in 0 1 2 3 4; do
    for MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES in $(seq -s' ' 1 $MAX_NONFOCAL_SUBCLADES); do
        bsub -K -o $LOG_DIR/lsf-job.%J.table_cluster_stats.log "
        Rscript $SCRIPT_DIR/analyze_tree/table_cluster_stats.R \
            --clusters $TMP_CLUSTERS_VARYING/${PREFIX}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_swiss_polytomies_F_clusters.txt \
            --metadata $TMP_ALIGNMENTS/$PREFIX/${PREFIX}_tree_metadata.txt \
            --outdir $TMP_CLUSTERS_STATS \
            --prefix ${PREFIX}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_swiss_polytomies_F \
            --plotutilityfunctions $SCRIPT_DIR/figures/plotting_utility_functions.R \
            --nswissseqs $N_SWISS_SEQS
        " &
    done
done
wait

MAX_NONFOCAL_SUBCLADES=0
MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES=0
bsub -K -o $LOG_DIR/lsf-job.%J.table_cluster_stats.log "
Rscript $SCRIPT_DIR/analyze_tree/table_cluster_stats.R \
    --clusters $TMP_CLUSTERS_VARYING/${PREFIX}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_swiss_polytomies_F_clusters.txt \
    --metadata $TMP_ALIGNMENTS/$PREFIX/${PREFIX}_tree_metadata.txt \
    --outdir $TMP_CLUSTERS_STATS \
    --prefix ${PREFIX}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_swiss_polytomies_F \
    --plotutilityfunctions $SCRIPT_DIR/figures/plotting_utility_functions.R \
    --nswissseqs $N_SWISS_SEQS
"


# ------------------------------------------------------
echo "--- Finished successfully ---"

# TODO: Copy intersting files to the output folder

if [ ! -z $NOTIFICATION_EMAIL ] ; then
    sendmail $NOTIFICATION_EMAIL < $SCRIPT_DIR/notification_email.txt
    echo "Notification email sent to ${NOTIFICATION_EMAIL}"
fi
