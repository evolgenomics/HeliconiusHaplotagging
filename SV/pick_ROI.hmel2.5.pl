use warnings;
use strict;

my $dist = 30000;
my $threshold = 25;
#my $dist = 100000;
#my $threshold = 25;
my %lookup = map {chomp; my @tmp = split "\t"; $tmp[1] => $tmp[4]} `cat hmel2.5.accumulative.lookup`;

my @data = map {chomp; [ split "\t" ]} `sort -k 5,5V -k 10,10V hmel2.5.10k_win.all.jaccard_matrix.accum.ROI.front | awk '\$15 >= $threshold'`;

my $last_chr = "";
my $last_name="";
my $last_roi="";
my $last_i=-$dist;
my $last_j=-$dist;
my %roi;
foreach my $i (0..$#data) {
#	if ($data[$i][1] eq "Herato0601") {
#		print "".(join("\t", @{$data[$i]}))."\n";
#		print "$data[$i][4]\t$last_i\t".(abs($data[$i][4]-$last_i))."\t";
#		print "$data[$i][9]\t$last_j\t".(abs($data[$i][9]-$last_j))."\n";
#	}
	if ($data[$i][1]."/".$data[$i][6] ne $last_name) {
#		print "Case 1\n";
		push @{$roi{$data[$i][0]."/".$data[$i][5]}}, [ @{$data[$i]} ];
		$last_roi = $data[$i][0]."/".$data[$i][5];
	} elsif (abs($data[$i][4]-$last_i) <= $dist && abs($data[$i][9]-$last_j) <= $dist) {
#		print "Case 2\n";
		push @{$roi{$last_roi}}, [ @{$data[$i]} ];
	} else {
		#Still same chromosome
		if ($last_roi =~/$data[$i][0]/) {
			my ($last_roi_i, $last_roi_j) = ($1,$2) if ($last_roi=~/:(\d+)-\d+\/Hmel\d+o:(\d+)-\d+/);
			#Check for overlap with the active ROI instead
			if (abs($data[$i][2]-$last_roi_i) <= $dist && abs($data[$i][7]-$last_roi_j) <= $dist) {
#				print "Case 3a\n";
				push @{$roi{$last_roi}}, [ @{$data[$i]} ];
			} else {
#				print "Case 3b\n";
				push @{$roi{$data[$i][0]."/".$data[$i][5]}}, [ @{$data[$i]} ];
				$last_roi = $data[$i][0]."/".$data[$i][5];
			}
		} else {	
#			print "Case 3c\n";
			push @{$roi{$data[$i][0]."/".$data[$i][5]}}, [ @{$data[$i]} ];
			$last_roi = $data[$i][0]."/".$data[$i][5];
		}
	}
	($last_name, $last_i, $last_j) = ($data[$i][1]."/".$data[$i][6], $data[$i][4],$data[$i][9]);
}
my %count_list;
my %black_list;
foreach my $r (sort keys %roi) {
	my @ch = split("/", $r);
	my ($chr_i, $chr_j) = ($1, $2) if ($r=~/(Hmel\d+o):.+\/(Hmel\d+o):/);
	my $min_i = 2e10;
	my $max_i = -1;
	my $min_j = 2e10;
	my $max_j = -1;
	foreach my $i (0..$#{$roi{$r}}) {
		#$sum_i+=$roi{$r}[$i][3];
		#$sum_j+=$roi{$r}[$i][7];
		$min_i = $roi{$r}[$i][2] if ($roi{$r}[$i][2] < $min_i);
		$min_j = $roi{$r}[$i][7] if ($roi{$r}[$i][7] < $min_j);
		$max_i = $roi{$r}[$i][2] if ($roi{$r}[$i][2] > $max_i);
		$max_j = $roi{$r}[$i][7] if ($roi{$r}[$i][7] > $max_i);
	}
	$count_list{$chr_i.":".$min_i."-".$max_i}=$ch[0];
	$count_list{$chr_j.":".$min_j."-".$max_j}=$ch[1];
}

foreach my $i (sort keys %count_list){
	my ($chr_i, $min_i, $max_i) = ($1, $2, $3) if ($i=~/(Hmel\d+o):(\d+)-(\d+)/);
	J: foreach my $j (sort keys %count_list) {
		next J if ($i eq $j);
		my ($chr_j, $min_j, $max_j) = ($1, $2, $3) if ($j=~/(Hmel\d+o):(\d+)-(\d+)/);
		next J if ($chr_i ne $chr_j);
		#Detect overlap
		if (($min_i <= $min_j && $max_i+9999 > $min_j) || ($min_j <= $min_i && $max_j+9999 > $min_i) || ($min_i <= $min_j && $max_i+9999 >= $max_j+9999) || ($min_j <= $min_i && $max_j+9999 >= $max_i+9999)) {
			$black_list{$count_list{$i}}++;
			$black_list{$count_list{$j}}++;
		}
	}
}

foreach my $r (sort keys %roi) {
	my @ch = split("/", $r);
	my ($chr_i, $chr_j) = ($1, $2) if ($r=~/(Hmel\d+o):.+\/(Hmel\d+o):/);
	my $score = 0;
	my $sum_i=0;
	my $sum_j=0;
	my $min_i = 2e10;
	my $max_i = -1;
	my $min_j = 2e10;
	my $max_j = -1;
	$black_list{$ch[0]} = 0 if (!(exists($black_list{$ch[0]})));
	$black_list{$ch[1]} = 0 if (!(exists($black_list{$ch[1]})));
	foreach my $i (0..$#{$roi{$r}}) {
		$sum_i+=$roi{$r}[$i][2];
		$sum_j+=$roi{$r}[$i][7];
		$min_i = $roi{$r}[$i][2] if ($roi{$r}[$i][2] < $min_i);
		$min_j = $roi{$r}[$i][7] if ($roi{$r}[$i][7] < $min_j);
		$max_i = $roi{$r}[$i][2] if ($roi{$r}[$i][2] > $max_i);
		$max_j = $roi{$r}[$i][7] if ($roi{$r}[$i][7] > $max_i);
		$score+=$roi{$r}[$i][14];
	}
	my $cen_i = sprintf("%d", $sum_i/(@{$roi{$r}})-$min_i);
	my $cen_j = sprintf("%d",$sum_j/(@{$roi{$r}})-$min_j);
	if ($black_list{$ch[0]} <= 8 && $black_list{$ch[1]} <= 8) {
		#print "$r\t$ch[0]:$min_i-$max_i\t$cen_i\t$ch[1]:$min_j-$max_j\t$cen_j\t".(@{$roi{$r}})."\n";
		#} else {
		print "".$black_list{$ch[0]}."/".$black_list{$ch[1]}."\t$r\t$chr_i:$min_i-".($max_i+9999)."\t".substr($chr_i, 0, 7).":".$lookup{$chr_i.":".$min_i}."-".$lookup{$chr_i.":".$max_i}."\t$cen_i\t$chr_j:$min_j-".($max_j+9999)."\t".substr($chr_j, 0, 7).":".$lookup{$chr_j.":".$min_j}."-".$lookup{$chr_j.":".$max_j}."\t$cen_j\t".(@{$roi{$r}})."\t$score\n";
		#print "PLOT: ".substr($chr_i, 0, 7)."\t".($lookup{$chr_i.":".$min_i}-250000)."\t".($lookup{$chr_j.":".$max_j}+250000)."\t".($lookup{$chr_i.":".$min_i}+$cen_i)."\t".($lookup{$chr_j.":".$min_j}+$cen_j)."\t$score\n";#$#}
	}
}
