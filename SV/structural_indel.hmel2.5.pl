use warnings;
use strict;

my $interval_file = $ARGV[0];
my $outfile = $ARGV[1];
my @int_list = map {chomp;$_} `cut -f 1 $interval_file | uniq`;

my @files;

foreach my $group (1..23) {
 push @files, "Group_".$group.".Melpomene.groupings.barcode.list";
}

my @pops;
my %p_link;
open (OUT, ">$outfile");
print OUT "INTERVAL_I\tCHROM\tSTART\tEND\tSUM_ALL\tCOUNT_ALL\tAVG_SIZE_ALL\t";
foreach my $f (0..$#files) {
	#push @pops, [ map {chomp; $p_link{$_} = $f; split "\t"} `awk '\$3~/C/ {l=substr(\$1,length(\$1),1);print "L"l"_"\$3}' $files[$f]` ];
	push @pops, [ map {chomp; $p_link{$_} = $f; split "\t"} `sed 's/,/\\n/g' ../$files[$f]` ];
	my $group=$1 if ($files[$f]=~/Group_(\d+)/);
	print OUT "SUM_$group\tCOUNT_$group\tAVG_SIZE_$group";
	if ($f == $#files) {
		print OUT "\n";
	} else {
		print OUT "\t";
	}
}

foreach my $int (@int_list) {
	my ($i, $j) = (0,1);
	my @intervals = split("\t", $int);
	#print $int."\n";
	#foreach my $j (($i+1)..$#intervals) {
		my %p_sum;
		my %p_counts;
		my $i_sum = 0;

		my @i_list = map {chomp; $_} `for file in /tmp/global2/frankyc/run164_L45_corrected.hmel2.5.sorted.dedup.bxSorted.linked_reads.full.filtered.sorted.bed.gz /tmp/global2/frankyc/run164_L67_corrected.hmel2.5.sorted.dedup.bxSorted.linked_reads.full.filtered.sorted.bed.gz; do tabix \$file $intervals[$i]; done | awk '{print \$4"\t"\$3-\$2}'| uniq`;
		
		#my %i_watch;
		#my $overlap=0;

		foreach my $vi (@i_list) {
			my ($match, $length) = ($1."_".$2, $3) if ($vi=~/(run\d+_L\d+)_A\d+(C\d+)B\d+D\d+\t(\d+)/); 
			#		print "VI: $vi\tMATCH: $match\n";	
			if (exists($p_link{$match})) {
				$p_sum{$p_link{$match}}+=$length;
				$p_counts{$p_link{$match}}++;
			}
			$i_sum+=$length;
			#$i_watch{$vi}++;
		}
		my $i_str = $intervals[$i];
		$i_str =~tr/:-/\t\t/;

		my $o_str="";
		foreach my $pop (0..$#files) {
#			$p_overlap{$pop} = 0 if (!(exists($p_overlap{$pop})));
			if ($p_counts{$pop}) {
				$o_str.=($p_sum{$pop})."\t".($p_counts{$pop})."\t".($p_sum{$pop}/$p_counts{$pop})."\t";
			} else {
				$o_str.="0\t0\tNA\t";#);($p_sum{$pop})."\t".()."\t".$p_overlap{$pop}."\t";
			}
		}
		chop($o_str);
		if (@i_list) {
			print OUT $intervals[$i]."\t".$i_str."\t".($i_sum)."\t".(@i_list)."\t".($i_sum/(@i_list))."\t".$o_str."\n";
		} else {
			print OUT $intervals[$i]."\t".$i_str."\t".($i_sum)."\t0\tNA\t".$o_str."\n";
		}
	#}
}
close(OUT);
