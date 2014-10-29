#!/usr/bin/perl

## convert vcf file to format used by genotyping software

use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;
use concordance;

sub usage{
	print STDERR basename($0) . " [options] <vcf>\n";
	print STDERR "  vcf: VCF file to convert\n";
	print STDERR " Options:\n";
	print STDERR "   --assume-hom | -a: Assume missing genotypes are homozygous reference.\n";
	
	exit;
}

my ($status, $hom, $vcf, $line, @entry, @calls);
$hom = '';
$status = GetOptions("assume-hom|a" => \$hom);
usage() if not $status or @ARGV != 1;
$vcf = shift;

if ($vcf =~ /\.gz/) {
	tie *VCF, 'IO::Zlib', $vcf, "rb" or die "Cannot read $vcf: $!";
}
else {
	open VCF, $vcf or die "Cannot read $vcf: $!";
}

$line = <VCF>;
while ($line =~ /^#/) {
	$line = <VCF>;
}
while($line = <VCF>){
	chomp $line;
	@entry = split /\t/, $line;
	@calls = map {convertVCF($_, $hom)} $entry[9..$#entry];
	print $entry[0] . ':'
		  . $entry[1] . ' '
		  . $entry[0] . ':'
		  . $entry[1] . ' '
		  . $entry[1] . ' '
		  . $entry[3] . ' '
		  . $entry[4] . ' '
		  . join(' ', @calls) . "\n";
}
close VCF;