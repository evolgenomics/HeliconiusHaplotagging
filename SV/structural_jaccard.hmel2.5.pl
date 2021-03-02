use warnings;
use strict;

my $interval_file = $ARGV[0];
my $outfile = $ARGV[1];
my @int_list = map {chomp;$_} `cat $interval_file`;

#Populate the sample IDs from the sampling sites lists
my @files=();
foreach my $group (1..23) {
  push @files, "Group_".$group.".Melpomene.groupings.barcode.list";
}

my @pops;
my %p_link;
open (OUT, ">$outfile");
print OUT "INTERVAL_I\tCHROM\tSTART\tEND\tINTERVAL_J\tCHROM\tSTART\tEND\tTOTAL_I\tTOTAL_J\tSHARED\t";
foreach my $f (0..$#files) {
	#push @pops, [ map {chomp; $p_link{$_} = $f; split "\t"} `awk '\$3~/C/ {l=substr(\$1,length(\$1),1);print "L"l"_"\$3}' $files[$f]` ];
	push @pops, [ map {chomp; $p_link{$_} = $f; split "\t"} `sed 's/,/\\n/g' ../$files[$f]` ];
	my $group=$1 if ($files[$f]=~/Group_(\d+)/);
	print OUT "TOTAL_I_$group\tTOTAL_J_$group\tSHARED_$group";
	if ($f == $#files) {
		print OUT "\n";
	} else {
		print OUT "\t";
	}
}

foreach my $int (@int_list) {
	my ($i, $j) = (0,1);
	my @intervals = split("\t", $int);
	
	#foreach my $j (($i+1)..$#intervals) {
		my %p_counts;
		my %p_overlap;

		my @i_list = map {chomp; $_} `for file in  /tmp/global2/frankyc/run164_L45_corrected.hmel2.5.sorted.dedup.bxSorted.linked_reads.full.filtered.sorted.bed.gz  /tmp/global2/frankyc/run164_L67_corrected.hmel2.5.sorted.dedup.bxSorted.linked_reads.full.filtered.sorted.bed.gz; do tabix \$file $intervals[$j]; done | cut -f 4| uniq`;
		my @j_list = map {chomp; $_} `for file in  /tmp/global2/frankyc/run164_L45_corrected.hmel2.5.sorted.dedup.bxSorted.linked_reads.full.filtered.sorted.bed.gz  /tmp/global2/frankyc/run164_L67_corrected.hmel2.5.sorted.dedup.bxSorted.linked_reads.full.filtered.sorted.bed.gz; do tabix \$file $intervals[$j]; done | cut -f 4| uniq`;
		
		my %i_watch;
		my $overlap=0;

		foreach my $vi (@i_list) {
			my $match = $1."_".$2 if ($vi=~/(run\d+_L\d+)_A\d+(C\d+)/); 
			if (exists($p_link{$match})) {
				$p_counts{$p_link{$match}}{"i"}{$vi}++;
			}
			$i_watch{$vi}++;
		}
		foreach my $vj (@j_list) {
			my $match = $1."_".$2 if ($vj=~/(run\d+_L\d+)_A\d+(C\d+)/);
			if (exists($p_link{$match})) {
				$p_counts{$p_link{$match}}{"j"}{$vj}++;
				if (exists($p_counts{$p_link{$match}}{"i"}{$vj})) {
					$p_overlap{$p_link{$match}}++;
				}
			}
			$overlap++ if (exists($i_watch{$vj}));
		}
		my $i_str = $intervals[$i];
		$i_str =~tr/:-/\t\t/;
		my $j_str = $intervals[$j];
		$j_str =~tr/:-/\t\t/;

		my $o_str="";
		foreach my $pop (0..$#files) {
			$p_overlap{$pop} = 0 if (!(exists($p_overlap{$pop})));
			$o_str.=((keys %{${$p_counts{$pop}}{"i"}}))."\t".((keys %{$p_counts{$pop}{"j"}}))."\t".$p_overlap{$pop}."\t";
		}
		chop($o_str);
		print OUT $intervals[$i]."\t".$i_str."\t".$intervals[$j]."\t".$j_str."\t".(@i_list)."\t".(@j_list)."\t".$overlap."\t".$o_str."\n";
	#}
}
close(OUT);
