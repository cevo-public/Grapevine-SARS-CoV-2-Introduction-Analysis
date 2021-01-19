#!/bin/bash
set -euo pipefail

while getopts a:m:t:b:d:s:r:n: flag
do
    case "${flag}" in
        a) ALIGNMENT=${OPTARG};;
        m) METADATA=${OPTARG};;
        t) TMP_DIR=${OPTARG};;
        b) MIN_DATE=${OPTARG};;
        d) MAX_DATE=${OPTARG};;
        s) MASK_SITES_SCRIPT=${OPTARG};;
        r) REFERENCE=${OPTARG};;
        n) NCOVDIR=${OPTARG};;
    esac
done

echo "Fetching seqs excluded by nextstrain"
wget -P $TMP_DIR https://raw.githubusercontent.com/nextstrain/ncov/master/defaults/exclude.txt

echo "Adding outgroup seqs to include.txt"
echo "Wuhan/Hu-1/2019" > $TMP_DIR/include.txt
echo "Wuhan/WH01/2019" >> $TMP_DIR/include.txt
# outgroup seqs (both lineage B/ clade L/ clade 19A)
# The MRCA of these 2 outgroup seqs should split lineage A seqs from lineage B seqs based on looking at the nextstrain global view

echo "Appling date filter, host filter, stricter length filter for non-Swiss seqs"
augur filter \
    --sequences $ALIGNMENT \
    --metadata $METADATA \
    --min-date $MIN_DATE \
    --max-date $MAX_DATE \
    --min-length 27000 \
    --exclude $TMP_DIR/exclude.txt \
    --include $TMP_DIR/include.txt \
    --exclude-where host!=human country=Switzerland date='2020' date='2020-01-XX' date='2020-02-XX' date='2020-03-XX' date='2020-04-XX' date='2020-05-XX' date='2020-06-XX' date='2020-07-XX' date='2020-08-XX' date='2020-09-XX' date='2020-10-XX' date='2020-11-XX' date='2020-12-XX' date='2020-01' date='2020-02' date='2020-03' date='2020-04' date='2020-05' date='2020-06' date='2020-07' date='2020-08' date='2020-09' date='2020-10' date='2020-11' date='2020-12' \
    --output $TMP_DIR/alignment_filtered_nonswiss.fasta &

echo "Appling date filter, host filter, less strict length filter for Swiss seqs"
augur filter \
    --sequences $ALIGNMENT \
    --metadata $METADATA \
    --min-date $MIN_DATE \
    --max-date $MAX_DATE \
    --min-length 20000 \
    --exclude $TMP_DIR/exclude.txt \
    --exclude-where host!=human country!=Switzerland date='2020' date='2020-01-XX' date='2020-02-XX' date='2020-03-XX' date='2020-04-XX' date='2020-05-XX' date='2020-06-XX' date='2020-07-XX' date='2020-08-XX' date='2020-09-XX' date='2020-10-XX' date='2020-11-XX' date='2020-12-XX' date='2020-01' date='2020-02' date='2020-03' date='2020-04' date='2020-05' date='2020-06' date='2020-07' date='2020-08' date='2020-09' date='2020-10' date='2020-11' date='2020-12' \
    --output $TMP_DIR/alignment_filtered_swiss.fasta &
wait

cat $TMP_DIR/alignment_filtered_nonswiss.fasta $TMP_DIR/alignment_filtered_swiss.fasta >> $TMP_DIR/alignment_filtered.fasta

echo "Cleaning up files"
rm $TMP_DIR/alignment_filtered_nonswiss.fasta
rm $TMP_DIR/alignment_filtered_swiss.fasta

# -------------------------------------------------------------------------------
# Mask problematic sites from alignment according to De Maio et al.
# -------------------------------------------------------------------------------

echo "Removing blank lines in alignment (necessary for site filtering tool)"
sed -i.bak '/^$/d' $TMP_DIR/alignment_filtered.fasta
rm $TMP_DIR/alignment_filtered.fasta.bak

echo "Downloading sites to mask"
wget -P $TMP_DIR https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf

echo "Masking problematic sites"
python3 $MASK_SITES_SCRIPT \
    --mask \
    --vcf $TMP_DIR/problematic_sites_sarsCov2.vcf \
    --reference_id Wuhan/Hu-1/2019 \
    --input_fasta $TMP_DIR/alignment_filtered.fasta \
    --output_fasta $TMP_DIR/alignment_filtered_masked.fasta


# -------------------------------------------------------------------------------
# Remove remaining problematic sequences from alignment with Nextstrain diagnostic
# -------------------------------------------------------------------------------

echo "Running nextstrain diagnostic tool"
mkdir -p $TMP_DIR/diagnostic

# Since I've already length-filtered and want to be more lenient with Swiss seqs, reset length criteria
sed -i 's/no_data_cutoff = 3000/no_data_cutoff = 10000/' $NCOVDIR/scripts/diagnostic.py

python3 $NCOVDIR/scripts/diagnostic.py \
    --alignment $TMP_DIR/alignment_filtered_masked.fasta \
    --metadata $METADATA \
    --reference $REFERENCE \
    --output-flagged $TMP_DIR/diagnostic/sequence-diagnostics.tsv \
    --output-diagnostics $TMP_DIR/diagnostic/flagged-sequences.tsv \
    --output-exclusion-list $TMP_DIR/diagnostic/to-exclude.txt

echo "Filtering out diagnostic to-exclude sequences"
augur filter \
    --sequences $TMP_DIR/alignment_filtered_masked.fasta \
    --metadata $METADATA \
    --exclude $TMP_DIR/diagnostic/to-exclude.txt \
    --include $TMP_DIR/include.txt \
    --output $TMP_DIR/alignment_filtered2_masked.fasta


# -------------------------------------------------------------------------------
# Get priority list for most similar sequences to Swiss sequences
# -------------------------------------------------------------------------------

echo "Linearizing fasta file (requirement for grepping sequences)"
awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' $TMP_DIR/alignment_filtered2_masked.fasta > $TMP_DIR/alignment_filtered2_masked_oneline.fasta

echo "Grepping focal Swiss sequences from alignment"
grep --no-group-separator -A 1 "Switzerland" $TMP_DIR/alignment_filtered2_masked_oneline.fasta > $TMP_DIR/swiss_alignment_filtered2_masked_oneline.fasta

echo "Grepping non-Swiss context sequences from alignment"
grep --no-group-separator -P -A 1 '^>(?!Switzerland)' $TMP_DIR/alignment_filtered2_masked_oneline.fasta > $TMP_DIR/context_alignment_filtered2_masked_oneline.fasta

echo "Running nextstrain priorities script"
python3 $NCOVDIR/scripts/priorities.py \
    --alignment $TMP_DIR/context_alignment_filtered2_masked_oneline.fasta \
    --reference $REFERENCE \
    --metadata $METADATA \
    --focal-alignment $TMP_DIR/swiss_alignment_filtered2_masked_oneline.fasta \
    --output $TMP_DIR/priorities.txt
