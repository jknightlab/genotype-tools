#!/usr/local/bin/Rscript

## obtain SNP positions for rsIDs from BioMart

library(biomaRt)
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
		make_option(c("-o", "--output"), default = '',
				help = "Name of output file."),
		make_option(c("-b", "--build"), default = '19',
				help = "Genome build to use."),
		make_option(c("-n", "--names"), default = 1L,
				help = "Index of column that contains SNP IDs in input file.")
)
parser <- OptionParser(usage = "%prog [options] input_file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

## check that positional arguments are present
if(length(arguments$args) < 1){
	cat("\nInput file is required\n\n")
	print_help(parser)
	stop()
}
if(length(arguments$args) > 1){
	cat("\nFound extra arguments on command line\n\n")
	print_help(parser)
	stop()
}

## check that mandatory options are present
if(opt$output == ''){
	cat("\nOutput file name is required\n\n")
	print_help(parser)
	stop()
}

## set name and host of mart to use
martName <- 'ENSEMBL_MART_SNP'
knownHosts <- c('18' = "may2009.archive.ensembl.org/biomart/martservice/",
		'19' = "www.ensembl.org/biomart/martservice/")
host <- opt$build

if(is.na(host)){
	cat("\nUnknown genome build requested: ", opt$build, "\nValid builds are: ", 
			names(knownMarts), sep = '')
	print_help(parser)
	stop()
}

## read input file
ids <- read.table(arguments$args[1], header = TRUE, stringsAsFactors = FALSE)[[opt$names]]


filter <- ''
if(opt$build != '19'){
	filter <- 'refsnp' 
} else{
	filter <- 'snp_filter'
}

mart <- useMart(martName, 'hsapiens_snp', host=host)
table <- getBM(c('refsnp_id', 'chr_name', 'chrom_start'), filters = filter,
		values = ids, mart = mart)
table <- table[order(table$chr_name, table$chrom_start),]

write.table(table, file=opt$output, row.names=FALSE, quote=FALSE)