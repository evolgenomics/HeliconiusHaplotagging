use warnings;
use strict;

my $tag = $ARGV[0];
my $window = 500;

#THIS IS IMPORTANT - the first SNP position is in the original coordinate. This has to adjust to allow correct lookup of the position!
my $seq_start = 3450000;

my @pos = map {chomp; $_-= $seq_start; $_ } `bcftools query -f "%POS\\n" sourceData/PC062_merged_Herato1603.3.45Mb.PL.AD.HAPCUT2.vcf.gz `; 

my @haps = map {chomp; $_} `cut -f 2 simHaps/PC062_merged_Herat1603_3.45Mb.simBlock.$tag.haps`;

foreach my $i (0..length($haps[0])-$window) {
	my %uniq_haps;
	my $win_i = $i;
	while ($pos[$win_i]-$pos[$i] <= $window) {
		$win_i++;
	}
	foreach my $h (0..$#haps) {
		$uniq_haps{substr($haps[$h], $i, $win_i-$i)}++;# if (!exists($uniq_haps{substr($haps[$h], $i, $i+$window)}))
	}
	print "$i\t".($pos[$win_i]-$pos[$i])."\t".(keys %uniq_haps)."\n";
}
