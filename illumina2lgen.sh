#!/bin/bash

## Convert FinalReport file for Illumina genotyping data to LGEN format
## usage: illumina2lgen.sh <final report> <output file> 
cut -f1,2,17,18 $1 | sed "1,10d" | sed "s/-/0/g" | awk '{print $2,$2,$1,$3,$4}' > $2
