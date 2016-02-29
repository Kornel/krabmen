#!/usr/bin/env sh

MRNA_RSEM_FILES=../download/mrna-rsem-2015-11-01/*/*data.txt

for f in `ls $MRNA_RSEM_FILES`
do
  DIR=`dirname $f`
  FILE=`basename $f`
  
  NEW_FILE="${DIR}/${FILE}-clean.txt"

  echo "Working on $FILE"

  # Second line contains comments
  awk '{if (NR != 2) { print } }' $f > $NEW_FILE

done



