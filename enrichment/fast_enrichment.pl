use strict;
use warnings;
use List::Util qw(shuffle max min sum);
use List::MoreUtils qw(uniq);

my ($chr_index, $start_index, $end_index) = (0,1,2);
@ARGV == 6 or die "chromsome bed, gwas_hits bed, genomic_features bed, replicate_start, replicate_end \n";

my ($chrom_bed, $gwas_file, $feat_file, $out_prefix, $rep_start, $rep_end) =  @ARGV;

open IN,$chrom_bed;
my %chr2len;
while(<IN>)
{
	next if(/^#/);
	chomp;
	my @c = split /\t/;
	$c[0] =~ s/^chr|^chrom//i;
	$chr2len{$c[0]} = $c[2];
}
close IN;
my %chr2cum;
for (sort {$a<=>$b} keys %chr2len)
{
	if($_ == 1)
	{
		$chr2cum{$_} = 0;
	}
	else
	{
		$chr2cum{$_} = $chr2cum{$_-1} + $chr2len{$_-1};
	}
}
my $tot_len = sum(values %chr2len);
print "Total genome length = $tot_len\n";

my @actual_feat_pos;
my @feat_len;
open IN, $feat_file;
while(<IN>)
{
	next if(/^#/);
	chomp;
	my @c = split /\t/;
	$c[$chr_index] =~ s/^chr|^chrom//i;
	next unless ($c[$chr_index] =~ /^\d+$/);
	push @actual_feat_pos, [$chr2cum{$c[$chr_index]} + $c[$start_index] + 1, $chr2cum{$c[$chr_index]} + $c[$end_index]];
}
close IN;
@actual_feat_pos = sort {$a->[0] <=> $b->[0]} @actual_feat_pos;
for (my $i=1; $i<@actual_feat_pos; $i++)
{
	if($actual_feat_pos[$i] ->[0] <= $actual_feat_pos[$i-1] ->[1])
	{
		if($actual_feat_pos[$i-1] ->[1] < $actual_feat_pos[$i] ->[1])
		{
			$actual_feat_pos[$i-1] ->[1] = $actual_feat_pos[$i] ->[1];
		}
		splice @actual_feat_pos, $i, 1;
		$i --;
	}
}
for (@actual_feat_pos)
{
	push @feat_len, $_->[1] - $_->[0] + 1;
}
print "Number of genomic features = ", scalar @feat_len, "\n";
print "Total length of genomic features = ", sum(@feat_len), "\n";

print "Number of unique lengths = ", scalar (uniq @feat_len), "\n";

my @sig_snp;
open IN, $gwas_file;
while(<IN>)
{
	next if(/^#/);
	chomp;
	my @c = split /\t/;
	$c[0] =~ s/^chr|^chrom//i;
	
	push @sig_snp, [$chr2cum{$c[0]} + $c[1] + 1, 0];
}
close IN;
print "Completed reading significant variant file: ", scalar @sig_snp, " significant SNPs\n";

my $actual_nhits = snp_vs_feat(\@sig_snp, \@actual_feat_pos);
# $actual_nhits = lazy_cmp(\@sig_snp, \@actual_feat_pos);

open OUT,">$out_prefix.$rep_start-$rep_end.txt";
print OUT "For $feat_file, the actual number of hits is $actual_nhits\n";
print OUT "replicate	num_hits\n";

my $perm_nhit;
for my $i ($rep_start..$rep_end)
{
	my @perm_ft;
	for(shuffle(@feat_len))
	{
		my $rand_int = int(rand($tot_len));
		push @perm_ft, [$rand_int+1, $rand_int+$_];
	}
	@perm_ft = sort {$a->[0] <=> $b->[0]} @perm_ft;
	for my $k (1..@perm_ft-1)
	{
		if($perm_ft[$k]->[0] <= $perm_ft[$k-1]->[1])
		{
			$perm_ft[$k]->[1] += $perm_ft[$k-1]->[1] + 1 - $perm_ft[$k]->[0];
			$perm_ft[$k]->[0] = $perm_ft[$k-1]->[1] + 1;
		}
	}
	
	$perm_nhit = snp_vs_feat(\@sig_snp, \@perm_ft);
	
	print OUT "$i\t$perm_nhit\n";
	
	$i % 10 == 0 and print "$i replicates completed.\n";
}
close OUT;


sub snp_vs_feat
{
	my ($ref_snp, $ref_feat) = @_;
	
	my $ret = 0;
	my @pos = (@{$ref_snp}, @{$ref_feat});
	@pos = sort {$a->[0] <=> $b->[0]} @pos;

	my $last_ft_end = 0;
	for my $p (0..$#pos-1)
	{
		if($pos[$p]->[1] == 0)
		{
			if($pos[$p]->[0]<=$last_ft_end or $pos[$p]->[0]==$pos[$p+1]->[0]) {
				$ret++;
			}
		}
		else
		{
			$last_ft_end = $pos[$p]->[1];
		}
	}
	if($pos[$#pos]->[1] == 0 and $pos[$#pos]->[0] <= $last_ft_end)
	{
		$ret++;
	}
	return $ret;
}

sub lazy_cmp
{
	my ($ref_snp, $ref_feat) = @_;
	
	my $ret = 0;
	for my $snp (@$ref_snp)
	{
		for my $ft (@$ref_feat)
		{
			if($snp->[0] >= $ft->[0] and $snp->[0] <= $ft->[1])
			{
				$ret ++;
				last;
			}
		}
	}
	return $ret;
}

