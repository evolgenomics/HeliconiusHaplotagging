use warnings;
use strict;
use Text::Levenshtein qw(distance); 

my $tag = $ARGV[0];
my $linked = $ARGV[1];
#my @phased_calls = map{chomp; my @tmp = split "\t"; push @pos, $tmp[0]; [ @tmp[1..$#tmp] ]} `bcftools query -f "%POS[\t%GT]\n" STITCH_968_shortReads_$tag/stitch.Herato1603.450000.550000.PL.reheadered.vcf.gz`;
system("bcftools query -f \"%CHROM\\t%POS\\n\" stitch.Herato1603.450000.550000.PL.HAPCUT2.$linked.vcf.gz > phased_seg_sites.pos");
my @pos = map {chomp; $_ } `cut -f 2 phased_seg_sites.pos`;

my %snp;
my @phased_calls = map{chomp; $_} `bcftools query -f "[%GT\t]\n"  stitch.Herato1603.450000.550000.PL.HAPCUT2.$linked.vcf.gz | sed 's/\\t\$//;s/.\\/./9\\|9/g;s/|/\\t/g;s/\\./9\\\t9/g' | datamash transpose | sed 's/\\t//g'`;# 50000.550000.PL.reheadered.vcf.gz`;
my @phased_calls_hiConf = map{chomp; $_} `bcftools query -f "[%PG\t]\n"  stitch.Herato1603.450000.550000.PL.HAPCUT2.$linked.vcf.gz | sed 's/\\.\\/\\./9\\\t9/g;s/\\t\$//;s/\\//\\t/g;s/\\./9\\\t9/g' | datamash transpose | sed 's/\\t//g'`;# 50000.550000.PL.reheadered.vcf.gz`;
#system("bcftools query -f \"[%GT\t]\n\" STITCH_968_shortReads_$tag/merged.$tag.PL.HAPCUT2.linked.vcf.gz | sed 's/\\t\$//;s/.\\/./.\\|./g;s/|/\\t/g' | datamash transpose | sed 's/\\t//g' > phased_calls.linked.haps");

my @truth=map{chomp; $_} `bcftools query -f "[%GT\t]\n" ../../simHaps/PC062_merged_Herat1603_3_4Mb.simBlock.$tag.GT.vcf.gz -r Herato1603:450000-550000 -T phased_seg_sites.pos | sed 's/\\t\$//;s/\\//\\t/g' | datamash transpose | sed 's/\\t//g' `;
#system("bcftools query -f \"[%GT\t]\n\" ../simHaps/PC062_merged_Herat1603_3_4Mb.simBlock.$tag.GT.vcf.gz -r Herato1603:450000-550000 -T phased_seg_sites.pos | sed 's/\\t\$//;s/\//\\t/g' | datamash transpose | sed 's/\\t//g' > truth.haps");

my $window = 10;#$ARGV[0];

open(PHASE_CORRECT, ">stitch.Herato1603.450000.550000.PL.HAPCUT2.$linked.phaseCorrect.out");#merged.$tag.PL.HAPCUT2.$linked.phasedSwitch.out");
#open(OUT, ">stitch.Herato1603.450000.550000.PL.HAPCUT2.$linked.phaseSwitch.out.tmp2");#merged.$tag.PL.HAPCUT2.$linked.phasedSwitch.out");
#print OUT "IND\tSWITCHERROR\tSWITCHES\n";
my @switch;
for (my $hap = 0; $hap < $#phased_calls; $hap+= 2) {
#for (my $hap = 0; $hap < 2; $hap+= 2) {
	my @dirs;
	my $last_dir = -1;
	$phased_calls[$hap]=~s/\./9/g;
	$phased_calls[$hap+1]=~s/\./9/g;


	my $last_call_0="";
	my $last_call_1="";
	my $last_truth_0="";
	my $last_truth_1="";
	#for (my $idx=0; $idx+$window < length($phased_calls[$hap]); $idx+=$window) {
#		my $win_extend = 0;
#		my $count = 0;
#		while (($count < $window) && $idx+$win_extend < length($phased_calls[$hap])) {
#			#print "Index: $idx / $win_extend = ".($idx+$win_extend)." - COUNT: $count\n";
#			#print "CALL   ".substr($phased_calls[$hap],$idx+$win_extend,10)."<-".substr($phased_calls_hiConf[$hap],$idx+$win_extend,10)."\n";
#			#print "TRUTH: ".substr($truth[$hap],$idx+$win_extend,10)."\n";#<-".substr($phased_calls_hiConf[$hap],$win_extend^i,1)."\n"
#			#print "CALL   ".substr($phased_calls[$hap+1],$idx+$win_extend,10)."<-".substr($phased_calls_hiConf[$hap+1],$idx+$win_extend,10)."\n";
#			#print "TRUTH: ".substr($truth[$hap+1],$idx+$win_extend,10)."\n";#<-".substr($phased_calls_hiConf[$hap],$win_extend^i,1)."\n"
#			#
#			#print "STR: ".(substr($phased_calls_hiConf[$hap],$idx+$win_extend,1).substr($phased_calls_hiConf[$hap+1],$idx+$win_extend,1))."\n";
#			if ((substr($phased_calls_hiConf[$hap],$idx+$win_extend,1).substr($phased_calls_hiConf[$hap+1],$idx+$win_extend,1))!~/9/) {
#				$count++;# if ($substr($phased_calls_hiConf[$hap],$idx+$win_extend,1) ne substr($phased_calls_hiConf[$hap+1],$idx+$win_extend,1));  
#			} else {
#				substr($truth[$hap],$idx+$win_extend,1,"9");
#				substr($truth[$hap+1],$idx+$win_extend,1,"9");
#				substr($phased_calls[$hap],$idx+$win_extend,1,"9");
#				substr($phased_calls[$hap+1],$idx+$win_extend,1,"9");
#				#print "Edited pos ".($idx+$win_extend)."\n";
#			}
#			$win_extend++;
#			#my $lkj = <STDIN>;
#		}
#		my $call_0=substr($phased_calls[$hap], $idx, $win_extend);
#		my $call_1=substr($phased_calls[$hap+1], $idx, $win_extend);
#		
#		my $truth_0=substr($truth[$hap], $idx, $win_extend);
#		my $truth_1=substr($truth[$hap+1], $idx, $win_extend);
#		
#		my $dir = -1;
#		my $cis = distance($call_0,$truth_0)+distance($call_1,$truth_1);
#		my $trans = distance($call_0,$truth_1)+distance($call_1,$truth_0);
#
#		if ($trans < $cis) {
#			$dir = 1;
#		} elsif ($cis < $trans) {
#			$dir = 0;
#		} elsif ($cis == $trans) {
#			$dir = $last_dir;
#		}
#		push @dirs, $dir;
#		if (($dir != $last_dir) && ($idx > 100)){
#			push @{$switch[$hap]}, $idx if (($dir != $last_dir) && ($idx > 0));
#		}
#		$last_call_0 = $call_0;
#		$last_call_1 = $call_1;
#		$last_truth_0 = $truth_0;
#		$last_truth_1 = $truth_1;
#		
#		$last_dir = $dir;
#	}
#	#print "There are ".(@{$switch[$hap]})." switches at ".(join(",", @{$switch[$hap]}))."\n";
#	if ($switch[$hap]) {
#		print OUT "$hap\t".(@{$switch[$hap]})."\t".(join(",", @{$switch[$hap]}))."\n";
#	} else {
#		print OUT "$hap\t0\tNA\n";
#	}
#
	my @correct_phased=(0,0,0,0);
	my $het_pos =0;
	foreach my $idx (0..length($phased_calls[$hap])) {
		if (substr($phased_calls_hiConf[$hap], $idx, 1) ne substr($phased_calls_hiConf[$hap+1], $idx, 1)) {# && substr($phased_calls[$hap], $idx, 1) ne substr($phased_calls[$hap+1], $idx, 1)) {
			$correct_phased[0]++ if (substr($phased_calls[$hap], $idx, 1) eq substr($truth[$hap], $idx, 1) && substr($phased_calls[$hap+1], $idx, 1) eq substr($truth[$hap+1], $idx, 1) && (substr($phased_calls[$hap], $idx, 1) ne 9));
			$correct_phased[1]++ if (substr($phased_calls[$hap+1], $idx, 1) eq substr($truth[$hap], $idx, 1) && substr($phased_calls[$hap], $idx, 1) eq substr($truth[$hap+1], $idx, 1) && (substr($phased_calls[$hap+1], $idx, 1) ne 9) );
			#$correct_phased[2]++ if (substr($phased_calls[$hap], $idx, 1) eq substr($truth[$hap+1], $idx, 1)  && (substr($phased_calls[$hap], $idx, 1) ne 9));
			#$correct_phased[3]++ if (substr($phased_calls[$hap+1], $idx, 1) eq substr($truth[$hap], $idx, 1)  && (substr($phased_calls[$hap+1], $idx, 1) ne 9));
			$het_pos++;
			#print $idx."\t".(substr($phased_calls[$hap], $idx, 1))."/".(substr($phased_calls[$hap+1], $idx, 1))."\t".(substr($phased_calls_hiConf[$hap], $idx, 1))."/".(substr($phased_calls_hiConf[$hap+1], $idx, 1))."\t".(substr($truth[$hap], $idx, 1))."/".(substr($truth[$hap+1], $idx, 1))."\t$correct_phased[0]\t$correct_phased[1]\t$correct_phased[2]\t$correct_phased[3]\n";
		}
	}
	#if ($correct_phased[2]+ $correct_phased[3] > $correct_phased[0]+$correct_phased[1]) {
	if ($correct_phased[1] > $correct_phased[0]) {
		#($correct_phased[0],$correct_phased[1]) = ($correct_phased[2],$correct_phased[3]);
		$correct_phased[0]=$correct_phased[1];# = ($correct_phased[2],$correct_phased[3]);
	}
	if ($het_pos > 0) {
		print PHASE_CORRECT "$hap\t".(($correct_phased[0])/($het_pos))."\t$correct_phased[0]\t$het_pos\n";#(length($phased_calls[$hap])*2)."\n";
	} else {
		print PHASE_CORRECT "$hap\tNA\t0\t0\n";#(length($phased_calls[$hap])*2)."\n";
	}
#	print "$hap\t".(($correct_phased[0])/($het_pos))."\n";#(length($phased_calls[$hap])*2)."\n";
}
close (PHASE_CORRECT);
