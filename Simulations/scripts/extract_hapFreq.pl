use warnings;
use strict;

#Usage - perl scripts/extract_hapFreq.pl <tag>
#perl scripts/extract_hapFreq.pl example

my %pos;
my $tag = $ARGV[0];

#THIS IS IMPORTANT - the first SNP position is in the original coordinate. This has to adjust to allow correct lookup of the position!
my $seq_start = 3450000;

map {chomp; my @tmp = split "\t"; $tmp[2]-= $seq_start; $pos{$tmp[0]} = $tmp[2] } `bcftools query -f "%CHROM\\t%POS\\n" sourceData/PC062_merged_Herato1603.3.45Mb.PL.AD.HAPCUT2.vcf.gz | awk '{print NR-1"\\t"\$0}' `; 
my @segments;
map {chomp; @segments = split "\t"} `head -1 simTrees/PC062_merged_Herat1603_3.45Mb.$tag.treeSeq`;

@segments = map {s/^\d+\///; my @tmp = split "-"; [ @tmp ]} @segments;



#Getting the haplotypes;

my @haps = map {chomp; $_} `cut -f 2 simHaps/PC062_merged_Herat1603_3.45Mb.simBlock.$tag.haps`;

open(OUT, ">simARG/PC062_merged_Herat1603_3.45Mb.$tag.hapFreq.summary");
foreach my $i (0..$#segments) {
	my %hapFreq;
	foreach my $h (0..$#haps) {
		push @{$hapFreq{substr($haps[$h], $segments[$i][0], $segments[$i][1]-$segments[$i][0]+1)}}, $h;
	}
	foreach my $hh (sort {$#{$hapFreq{$b}} <=> $#{$hapFreq{$a}}} keys %hapFreq) {
		print OUT "$i\t$segments[$i][0]\t$segments[$i][1]\t$pos{$segments[$i][0]}\t$pos{$segments[$i][1]}\t".(@{$hapFreq{$hh}})."\t".(@{$hapFreq{$hh}})/@haps."\t".(join(",", @{$hapFreq{$hh}}))."\t$hh\n";
	}
}
close (OUT);
