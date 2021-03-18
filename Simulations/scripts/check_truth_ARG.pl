use warnings;
use strict;

my @pos; 
my $tag = $ARGV[0];
#my @calls = map{chomp; my @tmp = split "\t"; push @pos, $tmp[0]; [ @tmp[1..$#tmp] ]} `bcftools query -f "%POS[\t%GT]\n" STITCH_968_shortReads_100rho/stitch.Herato1603.450000.550000.PL.reheadered.vcf.gz`;
system("bcftools query -f \"\%CHROM\\t\%POS\\n\" stitch.Herato1603.450000.550000.PL.reheadered.vcf.gz > focal.pos");
my @calls = map{chomp; my @tmp = split "\t"; push @pos, $tmp[0]; [ @tmp[1..$#tmp] ]} `bcftools query -f "%POS[\t%GT]\n" stitch.Herato1603.450000.550000.PL.reheadered.vcf.gz`;
system("awk '{print NR\"\\t\"\$0}' focal.pos > focal.pos.match");
my @qpos;
my @truth=map{chomp; my @tmp = split "\t"; push @qpos, $tmp[0]; [ @tmp[1..$#tmp] ]} `bcftools query -f "%POS[\t%GT]\n" ../../simHaps/PC062_merged_Herat1603_3_4Mb.simBlock.$tag.GT.vcf.gz -r Herato1603:450000-550000 -T focal.pos`;

system("awk '{print NR\"\\t\"\$0}' ../../muts_all.txt | sort -k 4,4 | join -1 4 -2 3 - focal.pos.match -a 1 | sed 's/ /\\t/g' | awk '{mod=\$2-64219-(\$7-1); if (NF < 7) {mod=0}; print \$3\"\\t\"\$1\"\\t\"\$2\"\\t\"mod\"\\t\"\$7}'  | awk '\$2 >= 450000 && \$2 <= 550000' > focal_TRUTH.pos_idx_correction"); 

my %treeSeq_adj;
map {chomp; my @tmp = split "\t"; $treeSeq_adj{$tmp[0]} = $tmp[1]} `cut -f 3,4 focal_TRUTH.pos_idx_correction`;

my @hapFreq;
my @hapFreq_avg;
my @treeSeg;
my $last_tree_str="";
#0	0	107	28	760	10	0.0103305785123967	5,14,167,181,311,330,356,402,559,733
#map {chomp; my @tmp = split "\t"; my $tree_str = $tmp[1]."-".$tmp[2]; push @treeSeg, [ ($tmp[1], $tmp[2]) ] if ($last_tree_str ne $tree_str); $last_tree_str = $tree_str; my @hap_id = split ",", $tmp[7]; foreach my $id (0..$#hap_id) {push @{$hapFreq[$hap_id[$id]]}, $tmp[6]}} `cut -f 1-8 ../../simARG/PC062_merged_Herat1603_3_4Mb.$tag.hapFreq.summary | awk '\$5 >= 450000 && \$4 <= 550000'`;
map {chomp; my @tmp = split "\t"; $treeSeq_adj{$tmp[1]} = 0 if (!(exists($treeSeq_adj{$tmp[1]}))); $treeSeq_adj{$tmp[2]} = 0 if (!(exists($treeSeq_adj{$tmp[2]}))); $tmp[1]-=$treeSeq_adj{$tmp[1]}; $tmp[2]-=$treeSeq_adj{$tmp[2]}; my $tree_str = $tmp[1]."-".$tmp[2]; push @treeSeg, [ ($tmp[1], $tmp[2]) ] if ($last_tree_str ne $tree_str); $last_tree_str = $tree_str; my @hap_id = split ",", $tmp[7]; foreach my $id (0..$#hap_id) {push @{$hapFreq[$hap_id[$id]]}, $tmp[6]}} `cut -f 1-8 ../../simARG/PC062_merged_Herat1603_3_4Mb.$tag.hapFreq.summary | awk '\$5 >= 450000 && \$4 <= 550000'`;
my $start_snp_id = 64219;

die("The positions are not identical!\n") if (@pos != @qpos);

my @matches;
my @comps;

my @glob_matches;
my @glob_comps;
open (OUT, ">stitch.Herato1603.450000.550000.PL.GTmatch.out");
my @seg_size;
foreach my $tree (0..$#treeSeg) {
	my ($snp_start, $snp_end) = @{$treeSeg[$tree]};
	$snp_start-=$start_snp_id;
	$snp_start=0 if ($snp_start < 0);
	$snp_end-=$start_snp_id;
	$snp_end = $#pos if ($snp_end > $#pos);
	push @seg_size, $snp_end-$snp_start;
	foreach my $snp ($snp_start..$snp_end) {
		foreach my $i (0..$#{$calls[$snp]}) {
			#my $tru = join("/", sort(split("|", $truth[$snp][$i])));
			my $tru = join("/", sort(split("/", $truth[$snp][$i])));
#			print "SNP_$snp\t$i\t$calls[$snp][$i] vs. $tru\n";
			#my $lakdjfkls = <STDIN>;
			if ($calls[$snp][$i] eq $tru) {
				$glob_matches[$i]++;
				$matches[$tree][$i]++;
			}
			$comps[$tree][$i]++;
			$glob_comps[$i]++;
		}	
	}
}
#
foreach my $tree (0..$#treeSeg) {
	foreach my $i (0..$#{$comps[$tree]}) {
		push @{$hapFreq_avg[$i]}, $hapFreq[$i*2][$tree];
		push @{$hapFreq_avg[$i]}, $hapFreq[$i*2+1][$tree];
		my $avg = ($hapFreq[$i*2][$tree]+$hapFreq[$i*2+1][$tree])/2;
		my $i_str = $tree;
		$matches[$tree][$i] = 0 unless ($matches[$tree][$i]);
		if ($tree == $#treeSeg) {
			print OUT "Ind_$i\t$i_str\t".($hapFreq[$i*2][$tree])."\t".($hapFreq[$i*2+1][$tree])."\t$avg\t$matches[$tree][$i]\t$comps[$tree][$i]\t".($matches[$tree][$i]/$comps[$tree][$i])."\n";
			
			$glob_matches[$i] = 0 unless ($glob_matches[$i]);
			my $sum=0;
			my $min=1e6;
			my $max = -1;
			foreach my $hf (0..$#{$hapFreq_avg[$i]}) {
				$sum+=$hapFreq_avg[$i][$hf]*$seg_size[int($hf/2)];#comps[$hf][$i];
				$min = $hapFreq_avg[$i][$hf] if ($hapFreq_avg[$i][$hf] < $min);
				$max = $hapFreq_avg[$i][$hf] if ($hapFreq_avg[$i][$hf] > $max);
			}
			print OUT "Ind_$i\tGlobal\t$min\t$max\t".($sum/(2*@pos))."\t$glob_matches[$i]\t$glob_comps[$i]\t".($glob_matches[$i]/$glob_comps[$i])."\n";
		} else {
			print OUT "Ind_$i\t$i_str\t".($hapFreq[$i*2][$tree])."\t".($hapFreq[$i*2+1][$tree])."\t$avg\t$matches[$tree][$i]\t$comps[$tree][$i]\t".($matches[$tree][$i]/$comps[$tree][$i])."\n";
		}
	}
}
close(OUT);
