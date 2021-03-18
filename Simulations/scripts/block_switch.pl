use warnings;
use strict;

my $threshold = 0.5;

#THIS IS IMPORTANT - the first SNP position is in the original coordinate. This has to adjust to allow correct lookup of the position!
my $seq_start = 3450000;

my @pos = map {chomp; $_-= $seq_start; $_ } `bcftools query -f "%POS\\n" sourceData/PC062_merged_Herato1603.3.45Mb.PL.AD.HAPCUT2.vcf.gz `; 

my @leaves;
my @tmp_groups;
my $line = -1;
my $last_tree = -1;
my $last_num = -1;

my $tag =$ARGV[0];
#my $flag = 0;
open (OUT, ">simARG/PC062_merged_Herat1603_3.45Mb.$tag.hapFreq.blockDist.out");#PC062_merged_Herat1603_3_4Mb.100rho.summary |");
open (IN, "simARG/PC062_merged_Herat1603_3.45Mb.$tag.hapFreq.summary");#PC062_merged_Herat1603_3_4Mb.100rho.summary |");
while (<IN>) {
	chomp;
	my @tmp = split "\t";
	#next unless ($tmp[$#tmp]=~/1/);
	if ($last_tree ne "$tmp[0]/$tmp[3]-$tmp[4]") {
		push @leaves, [ ($last_tree,@tmp_groups) ] if ($last_num > -1);# && ($flag));
		$line = 0;
		#$flag = 0;
		undef @tmp_groups;
	}
	#$flag = 1 if ($tmp[$#tmp]=~/1/);
	my @lv = split(",", $tmp[7]);
	foreach my $l (@lv) {
		$tmp_groups[$l] = $line;
	}
	$line++;
	$last_tree = "$tmp[0]/$tmp[3]-$tmp[4]";
	$last_num = $tmp[0];#/$tmp[3]-$tmp[4]";
}
close (IN);
push @leaves, [ ($last_tree,@tmp_groups) ];# if ($last_num > -1);

print OUT "BLOCK_I\tSTART\tEND\tBLOCK_J\tSTART\tEND\tDIST\n";
#print "BLOCK_I\tSTART\tEND\tBLOCK_J\tSTART\tEND\tDIST\n";
foreach my $start (0..$#leaves-1) {
	my $shift = $start+1;
	my $match = 0;
	my $count = 0;
	my ($istart, $iend) = ($1, $2) if ($leaves[$start][0]=~/\/(\d+)-(\d+)/);
	#print "Here we are... - $start to $shift - vs. $#leaves\n";
	my $last_val = 0;
	CHECK: while ($shift < $#leaves) {	
		foreach my $i (1..$#{$leaves[$start]}) {
			#		print "Comparing $leaves[$start][$i] vs $leaves[$shift][$i]\n";
			if ($leaves[$start][$i] == $leaves[$shift][$i]) {
				$match++;
			}
			$count++;
		}
		$shift++;
		#my $Lkj = <STDIN>;
		#	print "$start\t$shift -> ".($match/$count)."\n";
		#	print "$shift\n";
		last CHECK if (($match/$count < $threshold) || ($match/$count > $last_val));
		$last_val = $match/$count;
	};# until (($shift > $#leaves - 2) || ($match/$count < $threshold));
	my ($jstart, $jend) = ($1, $2) if ($leaves[$shift][0]=~/\/(\d+)-(\d+)/);
	print OUT "$start\t$istart\t$iend\t$shift\t$jstart\t$jend\t".($jstart-$iend)."\n";# -> ".($match/$count)."\n";
	#print "$start\t$istart\t$iend\t$shift\t$jstart\t$jend\t".($jstart-$iend)."\n";# -> ".($match/$count)."\n";
}
close(OUT);

