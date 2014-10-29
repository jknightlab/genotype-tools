#!/usr/bin/perl

my $VERSION = "1.0.1";

=head1 NAME

filterSNPs.pl - Select a subset of SNPs from a Matrix-eQTL genotype file.

=head1 SYNOPSIS

  filterSNPs.pl [options]

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

=item B<-p, --pos=FILE>
Name of the SNP position file. This can be a gzipped file or plain text.

=back

=head2 Output Options

=over 4

=item B<-o, --output=PREFIX>
Prefix for output files. Output consists of files 'PREFIX.geno' and 'PREFIX.snppos'

=item B<-r, --region=REGION>
A string describing the genomic region from which SNPs should be selected. The
region string should be of the form <chromosome>[:<start>[-<end>]]

=item B<-d, --delim=DELIMITER>
Record separator to use in output file. [S<default: "\t">]

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

  The following examples take the genotypes in example.geno.gz and produce
  output files in a F<selection> sub-directory. Note that the F<selection> directory
  is assumed to exist.
  
=head2 Selecting a region from a chromosome

  filterSNPs.pl --gen=example.geno.gz --pos=example.snppos --region=6:28466017-33449717 --out=selection/example_1

=head2 Selecting an entire chromosome

  filterSNPs.pl --gen=example.geno.gz --pos=example.snppos --region=6 --out=selection/example_2


=cut

=head1 AUTHOR


 Peter Humburg


=cut

use strict;
use warnings;
use IO::Zlib;
use Pod::Usage;
use Getopt::Long;
use File::Basename;

my ($out, $genFile, $posFile, $help, $version, $manual, $delim, $region);

$help      = '';
$version   = '';
$manual    = '';
$delim     = "\t";

my $status = GetOptions(
	"o|out=s"    => \$out,
	"g|gen=s"    => \$genFile,
	"p|pos=s" => \$posFile,
	"h|help"     => \$help,
	"man|manual" => \$manual,
	"v|version"  => \$version,
	"d|delim=s"  => \$delim,
	"r|region=s" => \$region
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
if (not defined $posFile) {
	pod2usage(
		-message => "SNP position file is required",
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
if(not defined $region){
	pod2usage(
		-message => "Region string is required",
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

## parse region string
my ($chrom, $start, $end) = split /[:\-]/, $region;

if ($genFile =~ /.gz/) {
	tie *GEN, 'IO::Zlib', $genFile, "rb" or die "Unable to read $genFile: $!";
} else {
	open GEN, $genFile or die "Unable to read $genFile: $!";
}
if ($posFile =~ /.gz/) {
	tie *POS, 'IO::Zlib', $posFile, "rb"
	  or die "Unable to read $posFile: $!";
} else {
	open POS, $posFile or die "Unable to read $posFile: $!";
}
open OUT, ">" . $out . ".geno"
  or die "Unable to open " . $out . ".geno: $!";

## parse SNP positions file and store all selected SNP IDs
my (%selection, @entry);
while(my $line = <POS>){
	chomp $line;
	@entry = split /\s/, $line;
	if($entry[1] eq $chrom){
		if(not defined $start){
			$selection{$entry[0]} = 1;
		}elsif($entry[2] >= $start){
			if(not defined $end){
				$selection{$entry[0]} = 1;
			} elsif($entry[2] <= $end){
				$selection{$entry[0]} = 1;
			}
		}
	}
}
close POS;

my $line = <GEN>;
print OUT $line;
while($line = <GEN>){
	chomp $line;
	@entry = split /\s/, $line;
	if($selection{$entry[0]}){
		print OUT join($delim, @entry) . "\n";
	}
}
close OUT;