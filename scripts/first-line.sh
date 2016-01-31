#!/usr/bin/env sh

MRNA_RSEM_FILES=../download/mrna-rsem-normalized-2015-11-01/*/*data.txt
DEST=../data/rsem-normalized

for i in `ls $MRNA_RSEM_FILES`; do head -n 1 $i > $DEST/`basename $i`; done


