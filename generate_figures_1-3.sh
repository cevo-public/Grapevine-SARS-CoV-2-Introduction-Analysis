#!/bin/bash
set -euo pipefail

while getopts s:w:d:r:m:p:l: flag
do
    case "${flag}" in
        s) SCRIPT_DIR=${OPTARG};;
        w) WORKDIR=${OPTARG};;
        d) DATE_THRESHOLD=${OPTARG};;
        r) PREFIX_DATA=${OPTARG};;
        m) MAX_NONFOCAL_SUBCLADES=${OPTARG};;
        p) MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES=${OPTARG};;
		l) TRANSMISSION_TO_TEST_DELAY=${OPTARG};;
    esac
done

# Figure 1 -----------------------------------------------------

# TODO: based on the Nextclade output, we can identify 20A.EU1 sequences and add these to the figure
# TODO: the clade identities are hardcoded in get_second_wave_lineages.R

# Summarize samples produced by lineages crossing date threshold
Rscript $SCRIPT_DIR/analyze_clusters/get_second_wave_lineages.R \
--datethreshold $DATE_THRESHOLD \
--tree $WORKDIR/tmp/lsd/${PREFIX_DATA}.timetree.nex \
--clades $WORKDIR/clades/swiss_alignment_filtered2_masked_oneline_clades.tsv \
--metadata $WORKDIR/tmp/alignments/${PREFIX_DATA}/${PREFIX_DATA}_tree_metadata.txt \
--outdir $WORKDIR/tmp/second_wave_lineages \
--utilityfunctions $SCRIPT_DIR/utility_functions.R \
--prefixdata $PREFIX_DATA

# Make figure
Rscript $SCRIPT_DIR/figures/figure_1.R \
--datethreshold $DATE_THRESHOLD \
--foundinglineagedata $WORKDIR/tmp/second_wave_lineages/${PREFIX_DATA}_lineages_crossing_${DATE_THRESHOLD}.txt \
--clades $WORKDIR/clades/swiss_alignment_filtered2_masked_oneline_clades.tsv \
--metadata $WORKDIR/tmp/alignments/${PREFIX_DATA}/${PREFIX_DATA}_tree_metadata.txt \
--outdir $WORKDIR/tmp/figures \
--utilityfunctions $SCRIPT_DIR/utility_functions.R \
--prefixdata $PREFIX_DATA

# Figure 2 -----------------------------------------------------

Rscript $SCRIPT_DIR/figures/figure_2.R \
--maxclusters $WORKDIR/tmp/clusters/${PREFIX_DATA}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_swiss_polytomies_F_clusters.txt \
--metadata $WORKDIR/tmp/alignments/${PREFIX_DATA}/${PREFIX_DATA}_tree_metadata.txt \
--outdir $WORKDIR/tmp/figures \
--utilityfunctions $SCRIPT_DIR/utility_functions.R \
--treedatawithasr $WORKDIR/tmp/asr/${PREFIX_DATA}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_swiss_polytomies_F_tree_data_with_asr.txt \
--prefix ${PREFIX_DATA}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES} \
--minclustersize 2

# Figure 3 -----------------------------------------------------

# TODO: which countries to emphasize and which are "other" is currently hardcoded, would be better to emphasize always to N countries instead 

# Make figure
Rscript $SCRIPT_DIR/figures/figure_3.R \
--maxclusters $WORKDIR/tmp/clusters/${PREFIX_DATA}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_swiss_polytomies_F_clusters.txt \
--minclusters $WORKDIR/tmp/clusters/${PREFIX_DATA}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}_swiss_polytomies_T_clusters.txt \
--metadata $WORKDIR/tmp/alignments/${PREFIX_DATA}/${PREFIX_DATA}_tree_metadata.txt \
--outdir $WORKDIR/tmp/figures \
--utilityfunctions $SCRIPT_DIR/utility_functions.R \
--transmissiontestdelay $TRANSMISSION_TO_TEST_DELAY \
--prefix ${PREFIX_DATA}_m_${MAX_NONFOCAL_SUBCLADES}_p_${MAX_CONSECUTIVE_BUDDING_NONFOCAL_SUBCLADES}

