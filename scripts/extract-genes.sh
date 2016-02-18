#!/usr/bin/env sh

MRNA_RSEM_FILES=../download/mrna-rsem-normalized-2015-11-01/*/*data.txt
DEST=../data/rsem-normalized
ALL="gene-names-all.txt"

echo "" > ${ALL}

for f in `ls $MRNA_RSEM_FILES`
do
  cat $f | awk -F " " ' { print $1 } ' >> ${ALL}
done

UNIQUE="gene-names.txt"
sort -u ${ALL} > ${UNIQUE}
