#!/bin/bash

## Simple plink wrapper to convert LGEN to PED files.
##usage: makePED.sh <prefix>
plink --lfile $1 --recode --out $1
