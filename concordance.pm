## functions for concordance analysis

package concordance;

use Exporter 'import';

@EXPORT = qw(getGenoCoord convertVCF readChrom concordanceSite alignAlleles);

sub getGenoCoord {
	my $line = shift;
	my @entry = split / /, $line;
	my ($chr, $pos);
	if ($entry[0] =~ /(.+):(\d+)/) {
		$chr = $1;
		$pos = $2;
	}
	else {
		$chr = $entry[0];
		$pos = $entry[2];
	}
	return ($chr, $pos);
}

sub convertVCF {
	my $gt  = shift;
	my $hom = shift;
	my $ans;

	my @entry = split /[:\/]/;
	if ($entry[0] ne '.') {
		$ans = $entry[0] + $entry[1];
	}
	else {
		if ($hom) {
			$ans = 0;
		}
		else {
			$ans = 'NA';
		}
	}
	return $ans;
}

sub readChrom {
	my $file  = shift;
	my $chrom = shift;
	my @lines = `tabix -p vcf $file $chrom`;
	chomp @lines;

	my %vcf =
	  map { my @entry = split /\t/; "$entry[0]:$entry[1]" => $_ } @lines;
	\%vcf;
}

## computes concordance at a single site for several samples
## returns reference to array indicating for each sample whether the calls matched
sub concordanceSite {
	my $nargs = @_;
	my $geno1 = shift;
	my $geno2 = shift;
	my $idx1  = shift;
	my $idx2  = shift;

	if ($nargs < 4) {
		die "Expected 4 parameters, got only $nargs.";
	}
	if ($nargs > 4) {
		warn
		  "Ignoring additional parameters. Expected 4 parameters, got $nargs.";
	}
	if (@$idx1 != @$idx2) {
		warn
		  "Index vectors are of different length. The longer one will be truncated.";
	}

	my @conc = ();
	my $value;
	for (my $i = 0; $i < @$idx1 and $i < @$idx2; $i++) {
		$value = 0;
		if ($$geno1[ $$idx1[$i] ] eq 'NA' or $$geno2[ $$idx2[$i] ] eq 'NA') {
			$value = 'NA';
		}
		elsif ($$geno1[ $$idx1[$i] ] == $$geno2[ $$idx2[$i] ]) {
			$value = 1;
		}
		push @conc, $value;
	}
	return \@conc;
}

## Ensure that the alles for a given site are lined up properly in two files.
## Matches the alleles and modifies the alleles of the second argument to line up with the first
## where possible.
## The input is expected to be two array references with the alleles in columns 3 and 4. The 
## second array should have genotypes starting in column 5.
sub alignAlleles {
	my $geno1 = shift;
	my @geno2 = @{shift()};
	if ($$geno1[3] ne $geno2[3] or $$geno1[4] ne $geno2[4]) {
		if ($$geno1[3] eq $geno2[4] and $$geno1[4] eq $geno2[3]) {
			$geno2[3] = $$geno1[3];
			$geno2[4] = $$geno1[4];
			@geno2[ 5 .. $#geno2 ] =
			  map { 2 - $_ } @geno2[ 5 .. $#geno2 ];
		}
		else {
			undef @geno2;
		}
	}
	return \@geno2;
}

1;
