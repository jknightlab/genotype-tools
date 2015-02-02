#!/bin/bash

## Genrate a map file from Illumina FinalReport file.
## usage: illumina2map.sh <final report> <output file>
cut -f1,19,20 $1 | sed "1,10d" | awk '{print $2,$1,"0",$3}' | sort -u > $2
