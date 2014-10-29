#!/usr/bin/perl

my $VERSION = "1.0.0";

=head1 NAME

gen2me.pl - Converting .gen files for use with Matrix-eQTL

=head1 SYNOPSIS

  gen2me.pl [options]

   Help Options:
   --help     Show help information for all available options.
   --manual   Read the full manual.
   --version  Show the version number and exit.

=cut

=head1 OPTIONS

Mandatory arguments to long options are mandatory for short options too.

=head2 Input Options

=over 4

=item B<-g, --gen=FILE>
Name of the genotyping file. This can be a gzipped file or plain text.

=item B<-s, --sample=FILE>
Name of the sample file. This can be a gzipped file or plain text.

=back

=head2 Output Options

=over 4

=item B<-o, --output=PREFIX>
Prefix for output files. Output consists of files 'PREFIX.geno' and 'PREFIX.snppos'

=item B<-r, --remove, --filter>
Remove SNPs without genomic coordinates.

=item B<-d, --delim=DELIMITER>
Record separator to use in output file. [S<default: "\t">]

=item B<-f, --family>
Use family ID as part of sample ID in output file. The default is to use just the individual ID.

=back

=head2 Help Options

=over 4

=item B<-h, --help>
Show the brief help information.

=item B<--manual, --man>
Read the manual, with examples.

=item B<-v, --version>
Show the version number and exit.

=back

=cut

=head1 EXAMPLES

  The following command would take the genotypes in example.gen.gz and produce
  output files in a matrixEQTL sub-directory. Note that the matrixEQTL directory
  is assumed to exist.

  gen2me.pl --gen=example.gen.gz --sample=sample.txt --out=matrixEQTL/example

=cut

=head1 DESCRIPTION

  Simple conversion between gen and Matrix-eQTL formats. It is assumed that all
  genotypes are located in a single file. If genotyping data is stotred in one file
  per chromosome these files should be concatenated first.
  
=cut

=head1 AUTHOR


 Peter Humburg


=cut

## convert genotyping data for use with Matrix-eQTL
## usage: gen2me.pl <gen file> <sample file> <output prefix>

use strict;
use warnings;
use IO::Zlib;
use Getopt::Long;
use Pod::Usage;
use File::Basename;

my ($out, $genFile, $sampleFile, $help, $version, $manual, $delim, $useFamily, $filter);

$help      = '';
$version   = '';
$manual    = '';
$useFamily = '';
$delim     = "\t";
$filter = '';

my $status = GetOptions(
	"o|out=s"    => \$out,
	"g|gen=s"    => \$genFile,
	"s|sample=s" => \$sampleFile,
	"h|help"     => \$help,
	"man|manual" => \$manual,
	"v|version"  => \$version,
	"d|delim=s"  => \$delim,
	"f|family"   => \$useFamily,
	"r|remove|filter" => \$filter
);

pod2usage(-verbose => 1) if $help;
pod2usage(-verbose => 2) if $manual;
print basename($0) . " Version $VERSION\n" and exit if ($version);

if (not defined $genFile) {
	pod2usage(
		-message => "Genotype file is required",
		-verbose => 1,
		-output  => \*STDERR
	);
}
if (not defined $sampleFile) {
	pod2usage(
		-message => "Sample file is required",
		-verbose => 1,
		-output  => \*STDERR
	);
}
if (not defined $out) {
	pod2usage(
		-message => "Output prefix is required",
		-verbose => 1,
		-output  => \*STDERR
	);
}
if (not $status) {
	pod2usage(
		-message => "Illegal command line arguments",
		-verbose => 1,
		-output  => \*STDERR
	);
}

if ($genFile =~ /.gz/) {
	tie *GEN, 'IO::Zlib', $genFile, "rb" or die "Unable to read $genFile: $!";
}
else {
	open GEN, $genFile or die "Unable to read $genFile: $!";
}
if ($sampleFile =~ /.gz/) {
	tie *SAMPLE, 'IO::Zlib', $sampleFile, "rb"
	  or die "Unable to read $sampleFile: $!";
}
else {
	open SAMPLE, $sampleFile or die "Unable to read $sampleFile: $!";
}
open OUT, ">" . $out . ".geno"
  or die "Unable to open " . $out . ".geno: $!";
open POS, ">" . $out . ".snppos"
  or die "Unable to open " . $out . ".snppos: $!";

## gen file contains information on genotypes as well as SNP positions. We split
## them into two files as is required by Matrix-eQTL
my (@entry, $id, @samples);

## skip header
my $line = <SAMPLE>;
$line = <SAMPLE>;
@samples = ("ID");
while ($line = <SAMPLE>) {
	chomp $line;
	@entry = split /\s/, $line;
	$id = '';
	if($useFamily){
		$id = "$entry[0].";
	}
	$id .= $entry[1];
	push @samples, $id
}
print OUT join($delim, @samples) . "\n";
close SAMPLE;

while($line = <GEN>){
	chomp $line;
	@entry = split /\s/, $line;
	next if $filter and $entry[2] !~ /[1-9]\d*/;
	print POS "$entry[1]$delim$entry[0]$delim$entry[2]\n";
	print OUT join($delim, @entry[1,5 .. $#entry]) . "\n";
} 
close POS;
close OUT;