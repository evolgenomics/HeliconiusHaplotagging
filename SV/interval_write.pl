use warnings;
use strict;

my @list = map {chomp; $_} `cat $ARGV[0]`;
foreach my $i (0..$#list-1) {
	foreach my $j ($i+1..$#list) {
		print $list[$i]."\t".$list[$j]."\n";
	}
}
