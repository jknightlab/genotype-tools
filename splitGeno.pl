#!/usr/bin/perl

## split genotype file by chromosome

use strict;
use warnings;
use File::Basename;
use FileHandle;
use Getopt::Long;

sub usage{
	print STDERR basename($0) . " [options] <index> <genotypes>\n";
	print STDERR "  index: Tabix indexed file of rsIDs.\n";
	print STDERR "  genotypes: Comma separated file with genotypes, samples in columns and sites in rows.\n";
	print STDERR " Options:\n";
	print STDERR "    --out| -o <prefix>: Prefix for output files.\n";
	
	exit;
}

## generate a new file handle
sub file_handle{
	my ($prefix, $chrom) = @_;
	
	my $fh = FileHandle->new();
	open $fh, ">$prefix.$chrom.geno" or die "Cannot write to $prefix.$chrom.geno: $!";
	return $fh;
}

my ($out, $idx, $geno, $header, $id, $prefix, @chrom, $line, %file, $record);
$out = '';
my $status = GetOptions("out|o=s" => \$out);
usage() if not $status or @ARGV != 2;

($idx, $geno) = @ARGV;
open IN, $geno or die "Cannot read $geno: $!";
$header = <IN>;

while($record = <IN>){
	$record =~ /^(\D+)(\d+),/;
	$prefix = $1;
	$id = $2;
	$line = `tabix $idx $prefix:$id-$id`;
	chomp $line;
	if($line){
		@chrom = split /\t/, $line;
		@chrom = split / /, $chrom[2];
		for my $chr (@chrom){
			if(not defined $file{$chr}){
				$file{$chr} = file_handle($out, $chr);
				$file{$chr}->print($header);
			}
			$file{$chr}->print($record);
		} 
	}
}

