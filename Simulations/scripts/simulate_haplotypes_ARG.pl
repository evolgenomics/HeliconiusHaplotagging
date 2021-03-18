use warnings;
use strict;

#Usage - perl scripts/simulate_haplotypes_ARG.pl <tag> <repID> <nHaps> <tree clustering cut-off>
#perl scripts/simulate_haplotypes_ARG.pl example 1 968 0.1
#
#The tree clustering cut-off is a way to control how many leave labels we retain - generally the closer the cut-off to 0, the more leaves we retain, i.e., creates more hapotype diversity

my $tag = $ARGV[0];
my $rep = $ARGV[1];
my $n_haps = $ARGV[2];
my $clust_threshold = $ARGV[3];

#Setting up positions of SNPs
my @pos;
my %snps;
my %muts;

#THIS IS IMPORTANT - the first SNP position is in the original coordinate. This has to adjust to allow correct lookup of the position!
my $seq_start = 0;#3450000;

map {chomp; my @tmp = split "\t"; $tmp[1]-= $seq_start; push @pos, $tmp[1]; $snps{$tmp[1]} = $tmp[3]; $muts{$tmp[1]} = $_ } `bcftools query -f "%CHROM\\t%POS\\t%REF\\t%ALT\\n" sourceData/PC062_merged_Herato1603.3.45Mb.PL.AD.HAPCUT2.vcf.gz`; 
chomp @pos;
my @haps;

#Setting up seed haplotypes
map {chomp; s/\./0/g; push @haps, $_} `cat sourceData/PC062_merged_Herat1603_3.45Mb.simBlock.horiz.sorted`;
print "Read ".(@haps)." seed haplotypes.\n";
my $seq;

#Note that here the reference sequence is the entire length of the sequence contig! I use a separate line to truncate the sequence as needed;
map {chomp; $seq = $_} `tail +2 sourceData/helera1_demo_Herato1603.fa | paste -s | sed 's/\\t//g'`;  
#Here the truncation
#$seq = substr($seq, 3350000, 300000);
#my $offset = 100000;

my @out_haps;

my $line = 0;
my @recomb_intervals;
my @treeSeq; 

map {s/\[//; s/\](.+)//;push @treeSeq, $1; push @recomb_intervals, [ ($line, $line+$_) ]; $line+=$_} `tail +5 "simARG/mspms_"$tag"_rep$rep.out"`;

open(TREESEQ, " | datamash transpose > simTrees/PC062_merged_Herat1603_3.45Mb.$tag.treeSeq");
#Construct hap lookup, first by iterating along intervals;
print "Constructing haplotype lookup table by transversing tree sequence...\n";
foreach my $int (0..$#treeSeq) {
	#Find the SNPs in the interval:
	my ($snp_start, $snp_end) = @{$recomb_intervals[$int]};
	my %clustHaps = &tree_hap_assign(&clust_trees($treeSeq[$int], $clust_threshold));
	my @out_leaves;
	LEAVE: foreach my $leave (sort {$a <=> $b} keys %clustHaps) {
		push @out_leaves, $clustHaps{$leave}; 
		$out_haps[$leave-1].= substr($haps[$clustHaps{$leave}], $snp_start, $snp_end-$snp_start);	
		last LEAVE if ($leave == $n_haps);	#}
	}
	print TREESEQ "$int/$snp_start-".($snp_end-1)."\t".(join("\t",@out_leaves))."\n";
}
close(TREESEQ);

#Build the reference fasta
print "Generating truth haplotypes and creating alternate reference sequences...\n";
open (HAPS, ">simHaps/PC062_merged_Herat1603_3.45Mb.simBlock.$tag.haps");
foreach my $nh (0..$n_haps-1) {
	
	my @muts = split "1", $out_haps[$nh];
	my @retained_muts;
	my $rm = -1;
	
	foreach my $seg (0..$#muts) {
		$rm+=length($muts[$seg])+1;
		push @retained_muts, $rm;
	}
	pop @retained_muts;
	my $i_seq = $seq;
#	open(MUTS, ">simMuts/PC062_merged_Herat1603_3.45Mb.simBlock.hap_$nh.$tag.muts");
	foreach my $i (0..$#retained_muts) {
			substr($i_seq, $pos[$retained_muts[$i]]-1, 1, $snps{$pos[$retained_muts[$i]]});
#			print MUTS $muts{$pos[$retained_muts[$i]]}."\n";
	}
#	close(MUTS);
	open(SNPS, ">altFa/PC062_merged_Herat1603_3.45Mb.simBlock.hap_$nh.$tag.fa");
	print SNPS ">Herato1603\n$i_seq\n";
	close(SNPS);
	
	print HAPS "$nh\t$out_haps[$nh]\n";
}
close(HAPS);
print "Done!\n";

#Subroutine to cluster trees based on a given depth - input: (newick_tree_string, cluster_threshold) / returns clustered_newick_tree;
sub clust_trees {
	my ($tree_str, $threshold) = @_;
	my $tree_tmp_str=$tree_str;

	my $round = 0;
	my @matches = $tree_str=~ /(\({1,1}(\d+:[0-9.]+,*)+\)+:[0-9.]+)/g;
	my $tree_str_bak = "";
	#Collapse clades until it stops changing
	while ($tree_str_bak ne $tree_tmp_str) {
		$tree_str_bak = $tree_tmp_str;
		foreach my $m (0..$#matches) {
			my $start_pos = index($tree_tmp_str, $matches[$m]);	
			
			if ($matches[$m]=~/\({1,1}((\d+:[0-9.]+,*)+)\)+:([0-9.]+)/) {
				my ($clade, $out_time) = ($1,$3);
				my @leaves = split(",", $clade);
				my @new_leaves;
				my $collapse_flag = 0;
				my $new_clade_str = $matches[$m];
				COL: foreach my $i (0..$#leaves) {
					my ($hap, $time) = ($1, $2) if ($leaves[$i]=~/(\d+):([0-9.]+)/);
					if ($time < $threshold) {
						$collapse_flag = 1;
						last COL;
					}
				}
				if ($collapse_flag) {
					foreach my $i (0..$#leaves) {
						my ($hap, $time) = ($1, $2) if ($leaves[$i]=~/(\d+):([0-9.]+)/);
						$time = $time+=$out_time;
						push @new_leaves, $hap.":".$time;
					}
					$new_clade_str = join(",", @new_leaves);
					substr($tree_tmp_str, $start_pos, length($matches[$m]), $new_clade_str);
				}
			}
		}
		$round++;
		@matches = $tree_tmp_str=~ /(\((\d+:[0-9.]+,*)+\):[0-9.]+)/g;
	}
	return $tree_tmp_str;
}

#Subroutine to return the assigned haplotype number, given the clade structure: input: newick_tree / returns hash with the keys being the leaveID in original newick_tree and ascending cladeID
sub tree_hap_assign {
	my $tree_str = $_[0];
	my %levels;

	my $order_count = 0;
	my %order;
	my $round = 0;
	my $counter = 0;
	my @matches = split(/([()])/, $tree_str);
	BRAC: foreach my $m (0..$#matches) {
		if ($matches[$m]=~/[()]/) { 
			$counter++;
			next BRAC;
		} elsif ($matches[$m]=~/(\d+):/) {
			my @leaves = split(",", $matches[$m]);
			foreach my $l (@leaves) {
				if ($l=~/(\d+):/) {
					$levels{$1} = $counter;
					$order{$1} = $order_count;
					$order_count++;
				}
			}
		}
	}
	$counter = -1;
	my $last_count = -1;
	my %level_final;
	foreach my $k (sort {$levels{$a} <=> $levels{$b}} keys %levels) {
		$counter++ if ($last_count != $levels{$k});
		$level_final{$k} = $counter;	
		$last_count = $levels{$k};
	}
	return %level_final; 
}
