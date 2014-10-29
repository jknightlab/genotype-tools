#!/usr/bin/perl

## Compute concordance between (array) genotyping calls and variant calls from sequencing.
## Takes a calls file and a vcf file and computes concordance for all genotyped sites. Note
## that the vcf file has to be tabix indexed.
## Optionally outputs a calls file corresponding to the variants in the vcf.

use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::Zlib;
use concordance;
use List::MoreUtils qw(pairwise);
use List::Util qw(first max);

sub usage {
	print STDERR basename($0) . " [options]\n";
	print STDERR "  Options:\n";
	print STDERR "    --geno | -g <FILE>: Genotyping file.\n";
	print STDERR "    --vcf | -v <FILE>: VCF file.\n";
	print STDERR
	  "    --sample | -s <FILE>: Sample names for Genotyping file.\n";
	print STDERR
	  "    --assume-hom| -a: Assume that samples with missing genotypes in VCF file are hom ref.\n";
	print STDERR
	  "    --column | -c <INTEGER>: Number of columns with row headers. [5]\n";
	print STDERR
	  "    --id | -i <SAMPLE>: Instead of assesing concordance between samples in the VCF and genotype files compare SAMPLE from the VCF file to all samples in the genotype file.\n";

	exit;
}

my (
	$status,   $geno,      $vcf,      $vcfLine, $genoLine,   @vcfHead,
	@vcfEntry, @genoEntry, @match,    @total,   $sampleFile, @sample,
	$nCol,     %genoHead,  @vcfCall,  $chr,     $pos,        @calls,
	$vcfChrom, $hom,       %vcfCache, $id,      $siteMatch,  @vcfIdx,
	@genoIdx,  $ans,       $max,      $ratio,   $bestMatch
);

$nCol   = 5;
$hom    = '';
$status = GetOptions(
	"geno|g=s"     => \$geno,
	"vcf|v=s"      => \$vcf,
	"sample|s=s"   => \$sampleFile,
	"column|c=i"   => $nCol,
	"assume-hom|a" => \$hom,
	"id|i=s"       => \$id
);

if (   not $status
	or not defined $geno
	or not defined $vcf
	or not defined $sampleFile) {
	usage();
}

tie *VCF, 'IO::Zlib', $vcf, "rb" or die "Cannot read $geno: $!";
if ($geno =~ /\.gz/) {
	tie *GENO, 'IO::Zlib', $geno, "rb" or die "Cannot read $geno: $!";
}
else {
	open GENO, $geno or die "Cannot read $geno: $!";
}
open SAMPLE, $sampleFile or die "Cannot read $sampleFile: $!";

@sample = <SAMPLE>;
close SAMPLE;

## sample file has two header rows
@sample            = @sample[ 2 .. $#sample ];
@sample            = map { (split /\s/)[0] } @sample;
@genoHead{@sample} = (0 .. $#sample);

$vcfLine = <VCF>;
while ($vcfLine =~ /^##/) {
	$vcfLine = <VCF>;
}
chomp $vcfLine;
$vcfLine =~ s/#//;
@vcfHead = (split /\t/, $vcfLine);
@vcfHead = @vcfHead[ 9 .. $#vcfHead ];
close VCF;

if (defined $id) {
	@match = map { 0 } @sample;
	@total = map { 0 } @sample;
	my $i = first { $vcfHead[$_] eq $id } (0 .. $#vcfHead);
	@vcfIdx = map { $i } @sample;
	@genoIdx = (0 .. $#sample);
}
else {
	@match = map { 0 } @vcfHead;
	@total = map { 0 } @vcfHead;
	@vcfIdx  = (0 .. $#vcfHead);
	@genoIdx = @genoHead{@vcfHead};
}
$vcfChrom = '';

while (defined *GENO
	and $genoLine = <GENO>) {
	chomp $genoLine;
	## remove 'chr' prefix from chromosome names
	$genoLine =~ s/^chr| chr/ /;
	@genoEntry = split / /, $genoLine;
	($chr, $pos) = getGenoCoord($genoLine);
	if ($chr ne $vcfChrom) {
		print STDERR "Processing chromosome $chr\n";
		%vcfCache = %{ readChrom($vcf, $chr) };
		$vcfChrom = $chr;
	}

	$vcfLine = $vcfCache{"$chr:$pos"};
	next unless defined $vcfLine;
	chomp $vcfLine;
	@vcfEntry = split /\s+/, $vcfLine;

	## this variant is present in both files so we compare genotypes
	if ($vcfEntry[0] eq $chr and $vcfEntry[1] == $pos) {
		## ensure alleles are lined up properly
		$ans = alignAlleles(\@vcfEntry, \@genoEntry);
		next unless defined $ans;
		@genoEntry = @$ans;

		@vcfCall = map { convertVCF($_, $hom) } @vcfEntry[ 9 .. $#vcfEntry ];

		$siteMatch =
		  concordanceSite(\@vcfCall, [ @genoEntry[ $nCol .. $#genoEntry ] ],
			\@vcfIdx, \@genoIdx);
		@total = pairwise {
			if   ($b eq 'NA') { $a }
			else              { $a + 1 }
		}
		@total, @$siteMatch;
		@match = pairwise {
			if   ($b ne 'NA') { $a + $b }
			else              { $a }
		}
		@match, @$siteMatch;

	}
	else {
		## if positions don't match something went wrong
		warn "Unexpected mismatch: VCF says $vcfEntry[0]:$vcfEntry[1] "
		  . "but genotyping says $chr:$pos\nVCF entry was\n $vcfLine\nGenotyping entry was\n$genoLine\n";
	}
}

print "Concordance\n";
$max = 0; 
for (my $i = 0; $i < @match; $i++) {
	print((@sample[@genoIdx])[$i]);
	if ($total[$i] > 0) {
		$ratio = $match[$i] / $total[$i];
		print sprintf(" %0.4f", $ratio);
		if ($ratio > $max) {
			$max       = $ratio;
			$bestMatch = (@sample[@genoIdx])[$i];
		}
	}
	else {
		print " NA";
	}
	print " ($match[$i]/$total[$i])\n";
}
if (defined $id) {
	print "Best match: $bestMatch (" . sprintf("%0.4f", $max) . ")\n";
}

