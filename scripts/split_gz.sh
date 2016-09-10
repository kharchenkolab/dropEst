#!/usr/bin/env bash
if [ "$#" -ne 1 ]; then
    echo "usage: split_gz.sh file"
    exit 1
fi
infile=${1##*/}
gunzip -c $infile | split -l 32000000 --filter='gzip > $FILE.gz && echo "$FILE created"' - "${infile%.*}-"