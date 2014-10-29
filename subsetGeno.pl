#!/usr/bin/perl

## remove samples from genotype file to match the samples present in expression file

use warnings;
use strict;
use File::Basename;
use Getopt::Long;

sub usage(){
	print STDERR basename($0) . " [options] <genotype file> <expression file>\n";
	print STDERR "   genotype file: File with genotype calls (coded as 0/1/2).\n";
	print STDERR "   expression file: File with gene expression values.\n";
	print STDERR " Options:\n";
	print STDERR "   --delim | -d <delimiter>: Delimiter used to separate records in input files. [\\t]";
	
	exit;
}

my ($geno, $express, %present, @index, @sample, @entry, $sep);
my $status = GetOptions("delim|d=s" => \$sep);

usage() if @ARGV != 2;
$geno = shift;
$express = shift;

## get sample names present in expression file
open EXPRESS, $express or die "Cannot read $express: $!";
my $line = <EXPRESS>;
chomp $line;
close EXPRESS;
@sample = split /$sep/, $line;
@present{@sample} = (1) x @sample;

open GENO, $geno or die "Cannot read $geno: $!";
$line = <GENO>;
chomp $line;
@sample = split /$sep/, $line;
@index = grep {defined $present{$sample[$_]}} (0..$#sample);

## print header
print join(',', @sample[@index]) . "\n";
@index = map {$_+1} @index;
while(<GENO>){
	@entry = split /$sep/;
	print "$entry[0]$sep";
	print join($sep, @entry[@index]) . "\n";
}
close GENO;