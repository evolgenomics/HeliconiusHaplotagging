use warnings;
use strict;

my %standard;
my %headers;
open (IN, "STITCH_eLife_COMP_STANDARD.GTDS");
while (<IN>) {
	chomp;
	my @tmp = split "\t";
	if (/CHROM/) {
		foreach my $i (6..$#tmp) {
			my $v = $tmp[$i];
			$v=~s/\[\d+\]//;
			$v=~/(\S+):GT/;
			$headers{$i} = $1;
		}
	}
	if (!/CHROM/) {
		foreach my $i (6..$#tmp) {
			my @v = split(":",$tmp[$i]);
			$standard{$headers{$i}.":".$tmp[1]} = $v[1];
#			print "$tmp[1]: Assiging $headers{$i} <= $v[1]\n"; 
		}
	}
}
close (IN);

my %lookup;
map{chomp; my @tmp = split "\t"; $lookup{$tmp[0]}=$tmp[3];} `awk '{print NR-1"\\t"\$0}' column_lookup.list`; 

my $table = $ARGV[0];
my $cols = "2,".`awk '!/Sample/' column_lookup.list | cut -f 1| paste -d"," -s`;
chomp($cols);
#print "cut -f $cols $table | sed 's/\\t\\S+://g'";
my @data=map{chomp; my @tmp = split "\t"; [ @tmp ]} `sed 's/\\t[01.]\\/[01.]:/\\t/g;s/\t\t/\t/' $table | cut -f $cols | grep -v POS`;

my $sum;
my $counts;
foreach my $i (0..$#data) {
	#print join("\t", @{$data[$i]})."\n";
	foreach my $j (1..$#{$data[$i]}) {
		next if ($j == 6);
		if (exists($standard{$lookup{$j}.":".$data[$i][0]})) {
			$sum+=abs($standard{$lookup{$j}.":".$data[$i][0]}-$data[$i][$j])/2;
			$counts++;
		}
	}
}
#print "Processed a total of $counts genotypes over ".(@data)." SNPs.\n";
#print "There is a total average concordance of ".(1-$sum/$counts)."\n";
print "$counts\t".(@data)."\t".(1-$sum/$counts)."\n";# SNPs.\n";
#print "There is a total average concordance of ".(1-$sum/$counts)."\n";

