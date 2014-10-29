#!/usr/bin/perl

use warnings;
use strict;
use File::Basename;
use Getopt::Long;

sub usage(){
	print STDERR basename($0) . " [options] <genotype file>\n";
	print STDERR "   genotype file: File in Plink AD format.\n";
	print STDERR " Options:\n";
	print STDERR "    --prefix <prefix>| -p: Prefix to use for sample names.\n";
	
	exit;
}

my ($prefix, $file, $status, @entry, @data, @sample, @sites, @index);
$prefix = '';

$status = GetOptions('prefix|p:s' => \$prefix);
usage() if @ARGV != 1 or not $status;
$file = shift;

## read file, discarding columns as necessary
open IN, $file or die "Cannot read $file: $!";

while(<IN>){
	@entry = split /\s/;
	if($. == 1){
		@index = map {$_ * 2} (3 .. ($#entry -5)/2);
		@sites = map {(split /_/)[0]} @entry[@index];
	}
	else{
		push @sample, $entry[0];
		push @data, [@entry[@index]];
	}
}
close IN;

## write transposed genotype matrix
## start with header
@sample = map {$prefix . $_} @sample;
print join(',', @sample) . "\n";

for(my $i = 0; $i < @sites; $i++){
	@entry = ($sites[$i]);
	for(my $j=0; $j < @data; $j++){
		push @entry, @{$data[$j]}[$i];
	}
	print join(',', @entry) . "\n";
}