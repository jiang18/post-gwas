use strict;
use warnings;
use List::Util qw(max);

@ARGV == 2 or die "input and output needed.\n";

my ($input, $output ) = @ARGV;

our $dis = 1000000;

my $p0 = 5e-8;
my $p1 = 5e-6;

my $nsig = 0;
my %sig;
my %suggest;
open IN,$input;
$_=<IN>;
while(<IN>)
{
	chomp;
	my @c = split /\t/;
	
	if($c[3] < $p0)
	{
		push @{$sig{$c[1]}}, [$c[2],abs($c[6])];
		$nsig++;
	}
	if($c[3] < $p1)
	{
		push @{$suggest{$c[1]}}, $c[2];
	}
}
close IN;

if($nsig==0)
{
	print "No QTLs found given P<5e-8\n";
	exit;
}

open OUT,">$output";
print OUT "CHR\tBP\tMaxAbsT\tN_Suggest_1Mb\n";

for my $m (1..29) {
	next unless(defined $sig{$m});
	my @arr_sig = @{$sig{$m}};
	my @arr_reg = ();
	push @{$arr_reg[0]}, $arr_sig[0];
	for(my $i=1; $i<@arr_sig; ++$i)
	{
		if($arr_sig[$i][0] - $arr_sig[$i-1][0] < $dis)
		{
			push @{$arr_reg[-1]},$arr_sig[$i];
		}
		else
		{
			push @arr_reg, [$arr_sig[$i]];
		}
	}
	for(@arr_reg)
	{
		my @signal = @{$_};
		@signal = sort {$a->[1] <=> $b->[1]} @signal;
		print OUT join ("\t",$m, $signal[-1][0], $signal[-1][1]),"\t";
		my $nvar = 0;
		for(@{$suggest{$m}}) {
			if(abs($signal[-1][0] - $_) < $dis/2 ) {
				$nvar++;
			}
		}
		print OUT "$nvar\n";
	}
}
close OUT;

sub log10
{
	return log($_[0])/log(10);
}