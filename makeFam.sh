#!/bin/bash

## Extract FAM file information from LGEN file. Note that the LGEN file
## has no paternal/maternal IDs, sex or phenotype information. These will
## consequently all be set to missing. This script is only suitable for 
## situations where these aren't required but a FAM file is needed (usually
## to use plink).
## usage: makeFam.sh <prefix>
perl -ane '{print "$F[0] $F[0] 0 0 0 -9\n"}' $1.lgen | sort -u -k1,1 > $1.fam
