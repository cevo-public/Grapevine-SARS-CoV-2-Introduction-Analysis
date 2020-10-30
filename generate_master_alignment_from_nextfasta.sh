#!/bin/bash
set -euo pipefail

while getopts i:t:f:m:r:x:n: flag
do
    case "${flag}" in
        i) INPUT_DIR=${OPTARG};;
        t) TMP_DIR=${OPTARG};;
        f) FASTA_FN=${OPTARG};;
        m) METADATA_FN=${OPTARG};;
        r) REFERENCE_FN=${OPTARG};;
        x) MAFFT=${OPTARG};;
        n) NUMBER_WORKERS=${OPTARG};;
    esac
done

echo "Unzipping data"
gunzip -c $INPUT_DIR/$FASTA_FN.gz > $TMP_DIR/$FASTA_FN
gunzip -c $INPUT_DIR/$METADATA_FN.gz > $TMP_DIR/$METADATA_FN

echo "Splitting sequence file into 150 seq or smaller chunks"
mkdir -p $TMP_DIR/nextdata_splits
mkdir -p $TMP_DIR/nextdata_alignments
mkdir -p $TMP_DIR/nextdata_alignments_noref

cd $TMP_DIR/nextdata_splits
awk -v size=150 '
   /^>/ { n++; if (n % size == 1) { close(fname); fname = sprintf("%d.fasta", n); print fname } }
   { print >> fname }
' $TMP_DIR/$FASTA_FN

echo "Aligning each split file"

i=1
for FILE in ${TMP_DIR}/nextdata_splits/*.fasta; do
    # This line ensures that at most NUMBER_WORKERS processes are run in parallel. The jobs are run in
    # NUMBER_WORKERS-sized batches, this means that not always NUMBER_WORKERS processes will be running.
    if ((i==0)) ; then
        wait
    fi
    i=$(((i+1)%NUMBER_WORKERS))

    FN=$(basename -- "$FILE")
    FN="${FN%.fasta}"
    $MAFFT \
        --addfragments $FILE \
        --keeplength \
        --auto \
        ${INPUT_DIR}/${REFERENCE_FN} > ${TMP_DIR}/nextdata_alignments/${FN}.fasta &
done
wait

echo "Removing reference from each split alignment"
for FILE in ${TMP_DIR}/nextdata_alignments/*; do
    FN=$(basename -- "$FILE")
    FN="${FN%.fasta}"
    awk 'BEGIN { RS = ">";FS = "\n" } {if (NR>2) {print ">"$0}}' $FILE > ${TMP_DIR}/nextdata_alignments_noref/${FN}.fasta
    # There is a blank first record, second record is the reference, all records thereafter are sequences we want to keep
done

echo "Making master monster alignment"
cat ${TMP_DIR}/nextdata_alignments_noref/*.fasta >> ${TMP_DIR}/nextdata_alignment.fasta

echo "Cleaning up directories"
tar -czf ${TMP_DIR}/nextdata_splits.tar.gz ${TMP_DIR}/nextdata_splits
tar -czf ${TMP_DIR}/nextdata_alignments.tar.gz ${TMP_DIR}/nextdata_alignments
tar -czf ${TMP_DIR}/nextdata_alignments_noref.tar.gz ${TMP_DIR}/nextdata_alignments_noref
rm -r ${TMP_DIR}/nextdata_splits ${TMP_DIR}/nextdata_alignments ${TMP_DIR}/nextdata_alignments_noref
