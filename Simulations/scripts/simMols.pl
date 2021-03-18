use warnings;
use strict;

my $hap_id = $ARGV[0];
my $tag = $ARGV[1];
my $cov = 250;

my $len = 300000;
my $tar_len = $cov * $len;
my $cur_len = 0;

my ($interval_start, $interval_end) = (3450000, 3550000);
my ($interval_buf_start, $interval_buf_end) = (3350000, 3650000);

my $stdev_reads = 2.5;
my $readpairs_per_mol = 5;

my $meanMol_len=100000;
my $stdev_len = 50000;

my $gap = 50000;

#Setting up to read in the BAM file
my @bam;
my %readpair;
my %used_read; 
my $all_reads = 0;

map{chomp; my @tmp = split "\t"; push @bam, [ @tmp ]; push @{$readpair{$tmp[0]}}, $#bam; } `samtools view $tag/hap_$hap_id.$tag.bam`;

open(BAM, "| samtools sort - | samtools view - -O BAM -o $tag/hap_$hap_id.linkedReads.haploid.$tag.bam");
my $header = `samtools view -H $tag/hap_$hap_id.$tag.bam`;
print BAM $header;

my @sel_array = 
#Setting up to generate a molecule
my @mols;
while ($cur_len < $tar_len) {
	my $mol_len = gaussian_rand()*$stdev_len + $meanMol_len;
	next if ($mol_len < 0);
	
	my $start = int(rand() * $len - $mol_len/2 + 0.5)+$interval_buf_start;
	my $end = int($start + $mol_len + 0.5);
	
	next unless (($end >= $interval_start && $end <= $interval_end) || ($start >= $interval_start && $start <= $interval_end));
	$start = 1 if ($start < 1);
	$end = $len+$interval_buf_start if ($end > $len+$interval_buf_start);
	$mol_len = $end-$start;

	#$start += $interval_buf_start;
	#$end += $interval_buf_end; 

	push @mols, [ ($start, $end, $mol_len) ];
	$cur_len += $mol_len;
}

open(MOL, ">$tag/hap_$hap_id.molecules.$tag.out");
print "[simMol.pl] Picking reads...\n";
M: foreach my $i (0..$#mols) {
	#print "chrI\t$mols[$i][0]\t$mols[$i][1]\t$mols[$i][2]\n";
	my $idx_start = -1;
	A: while (($bam[$idx_start+1][3] < $mols[$i][0]) && ($idx_start < $#bam-1)) {
		do {
			$idx_start++
		} while (exists($used_read{$idx_start}));
		last A if ($bam[$idx_start][3] == 0);
	}
	my $idx_end = $idx_start+1;
	if ($bam[$idx_end+1][3] > 0) {
		B: while (($bam[$idx_end+1][3] < $mols[$i][1]) && ($idx_end < $#bam-1)) {
			do {
				$idx_end++
			} while (exists($used_read{$idx_end}));
			last B if ($bam[$idx_end][3] == 0);
		}
	}
	#print "Proposing target_readpair number...\n";
	my $target_readpairs;
	do {
		#	print "trying...\n";
		$target_readpairs = int((gaussian_rand()*$stdev_reads/2 + $readpairs_per_mol)+0.5);
	} until $target_readpairs > 0;
	#print "Target readpair: $target_readpairs\n";
	
	my @picked_reads;
	my $min_start=$len+1;
	my $max_end=-1;
	
	#Picking left-end
	my $flag = 0;
	my $tries = 0;
	#print "Picking read_idx between: $idx_start - $idx_end\n";
	do {
		my $avail_reads = 0;
		my $pick_bracket = $idx_start+1;
		while ($avail_reads < 20) {
			$avail_reads++ if (!(exists($used_read{$pick_bracket})));
			$pick_bracket++;
			$tries++;
			next M if ($tries == 200);
		}
		$tries = 0;
		my $read_idx = int( rand()*($pick_bracket-$idx_start)+0.5)+$idx_start;
		#print "Mol $i/Case 0 ".($all_reads/@bam)."\n";
		#	print "Proposed read_idx: $read_idx ($idx_start - $idx_end)\n";
			if (!(exists($used_read{$read_idx}))) {
			foreach my $read (@{$readpair{$bam[$read_idx][0]}}) {
				push @picked_reads, $read;
				$used_read{$read}++;	
				$min_start = $bam[$read][3] if ($bam[$read][2] eq "Herato1603" && $bam[$read][3] < $min_start);
				$max_end = $bam[$read][3]+150 if ($bam[$read][2] eq "Herato1603" && $bam[$read][3]+150 > $max_end);
				#				print "Adding $bam[$read][0]/$read to molecule\n";#- read not used\n";
			}
			$flag = 1;
		}
		$tries++;
	} while (($flag == 0) && ($tries <= 20));
	#if ($tries > 20) {
	#	print "ERROR: Tried 20 times and not finding an appropriate read!\n";
	#	print "CurrentMol: Mol_$i - $mols[$i][0]/$idx_start - $mols[$i][1]/$idx_end\t$mols[$i][2]/".($idx_end-$idx_start+1)." - ".(@picked_reads)."\n";
	#	my $Lkjalj = <STDIN>;
	#}
	$target_readpairs--;

	#Picking right-end
	$flag = 0;
	$tries = 0;
	#print "Picking read_idx between: $idx_start - $idx_end\n";
	do {
		my $avail_reads = 0;
		my $pick_bracket = $idx_end-1;
		while ($avail_reads < 20) {
			$avail_reads++ if (!(exists($used_read{$pick_bracket})));
			$pick_bracket++;
			$tries++;
			next M if ($tries == 200);
		}
		$tries = 0;
		my $read_idx = int( rand()*($idx_end-$pick_bracket)+0.5)+($pick_bracket);
		#print "Mol $i/Case 1 ".($all_reads/@bam)."\n";
		#	print "Proposed read_idx: $read_idx ($idx_start - $idx_end)\n";
		if (!(exists($used_read{$read_idx}))) {
			foreach my $read (@{$readpair{$bam[$read_idx][0]}}) {
				push @picked_reads, $read;
				$used_read{$read}++;	
				$all_reads++;
				$min_start = $bam[$read][3] if ($bam[$read][2] eq "Herato1603" && $bam[$read][3] < $min_start);
				$max_end = $bam[$read][3]+150 if ($bam[$read][2] eq "Herato1603" && $bam[$read][3]+150 > $max_end);
				#				print "Adding $bam[$read][0]/$read to molecule\n";#- read not used\n";
			}
			$flag = 1;
		}
		$tries++;
	} while (($flag == 0) && ($tries <= 20));
	#if ($tries > 20) {
		#print "ERROR: Tried 20 times and not finding an appropriate read!\n";
		#print "CurrentMol: Mol_$i - $mols[$i][0]/$idx_start - $mols[$i][1]/$idx_end\t$mols[$i][2]/".($idx_end-$idx_start+1)." - ".(@picked_reads)."\n";
		#my $Lkjalj = <STDIN>;
		#}
	$target_readpairs--;
	
	for (my $r = $target_readpairs; $r > 0; $r--) {
		$flag = 0;
		$tries = 0;
		#print "Picking read_idx between: $idx_start - $idx_end\n";
		do {
			my $read_idx = int( rand()*($idx_end-$idx_start)+0.5) +$idx_start;
			#	print "Proposed read_idx: $read_idx ($idx_start - $idx_end)\n";
			#print "Mol $i/Case 2 ".($all_reads/@bam)."\n";
			if (!(exists($used_read{$read_idx}))) {
				foreach my $read (@{$readpair{$bam[$read_idx][0]}}) {
					push @picked_reads, $read;
					$used_read{$read}++;	
					$all_reads++;
					$min_start = $bam[$read][3] if ($bam[$read][2] eq "Herato1603" && $bam[$read][3] < $min_start);
					$max_end = $bam[$read][3]+150 if ($bam[$read][2] eq "Herato1603" && $bam[$read][3]+150 > $max_end);
					#				print "Adding $bam[$read][0]/$read to molecule\n";#- read not used\n";
				}
				$flag = 1;
			}
			$tries++;
		} while (($flag == 0) && ($tries <= 20));
		#if ($tries > 20) {
			#print "ERROR: Tried 20 times and not finding an appropriate read!\n";
			#print "CurrentMol: Mol_$i - $mols[$i][0]/$idx_start - $mols[$i][1]/$idx_end\t$mols[$i][2]/".($idx_end-$idx_start+1)." - ".(@picked_reads)."\n";
			#my $Lkjalj = <STDIN>;
			#}
	}

	#Check if there has been too much internal gap
	my @sorted_picked_reads = sort {$a<=>$b} @picked_reads;

	foreach my $r (0..$#sorted_picked_reads-1) {
		if (($bam[$sorted_picked_reads[$r]][2] eq "Herato1603" && $bam[$sorted_picked_reads[$r+1]][2] eq "Herato1603") && ($bam[$sorted_picked_reads[$r+1]][3]-$bam[$sorted_picked_reads[$r]][3] > $gap)) {
			#Case where the gap is too big
			#print "CASE: large gap\n";
			#print "Mol $i/Case 3 ".($all_reads/@bam)."\n";
			$flag = 0;
			$tries = 0;
			#print "Picking read_idx between: $idx_start - $idx_end\n";
			do {
				my $read_idx = int( rand()*($sorted_picked_reads[$r+1]-$sorted_picked_reads[$r])+0.5) + $sorted_picked_reads[$r];
				#	print "Proposed read_idx: $read_idx ($idx_start - $idx_end)\n";				

				if (!(exists($used_read{$read_idx}))) {
					foreach my $read (@{$readpair{$bam[$read_idx][0]}}) {
						push @picked_reads, $read;
						$used_read{$read}++;
						$all_reads++;
						$min_start = $bam[$read][3] if ($bam[$read][2] eq "Herato1603" && $bam[$read][3] < $min_start);
						$max_end = $bam[$read][3]+150 if ($bam[$read][2] eq "Herato1603" && $bam[$read][3]+150 > $max_end);
						#				print "Adding $bam[$read][0]/$read to molecule\n";#- read not used\n";
					}
					$flag = 1;
				}
				$tries++;
			} while (($flag == 0) && ($tries <= 20));
			#		if ($tries > 20) {
				#print "ERROR - Case LARGE_GAP: Tried 20 times and not finding an appropriate read!\n";
				#}
		}
	}
	my $read_name_str="";
	foreach my $pi (0..$#picked_reads) {
		print BAM "".(join("\t", @{$bam[$picked_reads[$pi]]}))."\tBX:Z:A".$hap_id."_".$i."\n";
		$read_name_str.=$bam[$picked_reads[$pi]][0].",";
	}
	chop($read_name_str);
	print MOL "$hap_id\tA".$hap_id."_$i\t$min_start\t$max_end\t".($max_end-$min_start)."\t".(@picked_reads)."\t$read_name_str\n";
		
	#print "Mol_$i - $mols[$i][0]/$idx_start - $mols[$i][1]/$idx_end\t$mols[$i][2]/".($idx_end-$idx_start+1)." - ".(@picked_reads)."\n";
	#print "Mol_$i - $min_start/$idx_start - $max_end/$idx_end\t".($max_end-$min_start)."/$mols[$i][2]/".($idx_end-$idx_start+1)." - ".(@picked_reads)."\n";
	#last if ($all_reads/@bam >= 0.8);
}
close(MOL);
close(BAM);

sub gaussian_rand {
    my ($u1, $u2);  # uniformly distributed random numbers
    my $w;          # variance, then a weight
    my ($g1, $g2);  # gaussian-distributed numbers

    do {
        $u1 = 2 * rand() - 1;
        $u2 = 2 * rand() - 1;
        $w = $u1*$u1 + $u2*$u2;
    } while ( $w >= 1 );

    $w = sqrt( (-2 * log($w))  / $w );
    $g2 = $u1 * $w;
    $g1 = $u2 * $w;
    # return both if wanted, else just one
    return wantarray ? ($g1, $g2) : $g1;
}
