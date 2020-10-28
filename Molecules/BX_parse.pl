#!/bin/perl -W

#Generates a BED file for each molecule from a BX-sorted BAM file, defined as reads that share the same BX tag with a maximum internal gap of 50kb
#This current version is intended to use on BAM files mapped to the mouse mm10 genome. Change line 144-145 if using a different genome

#Usage:
#perl BX_parse.pl <POS_SORTED_BAM> <PHASESET> <OUTPUT>
#
#PhaseSet is in a modified BED format, with the following columns: CHROM POS-1 POS HAP0 HAP1 COVERAGE
#Chr8 112 113 G A 15
#Chr8 118 119 A T 13
#Chr8 125 126 C T 17
#Chr8 274 275 C A 58
#Chr8 313 314 T C 52
#Chr8 403 404 C A 38
#Chr8 463 464 G A 63
#Chr8 505 506 G A 61
#Chr8 542 543 C T 66

#The output is in a non-standard tab-separated table format, with each row describing a molecule and for each the following columns
#1. MOLECULE_ID
#2. CHROM
#3. START
#4. END
#5. ALLELE_STRING{0,<,1,>}
#6. QUAL_STRING {PHRED+Q33}
#7. #_OF_SNPS
#8. #_OF_READS
#9. RECOMBINATION_INTERVALS[HYPHENED INTERVALS,COMMA_SEPARATED]
#10. POSITION_OF_SNPS[COMMA_SEPARATED]
#11. READS_IN_MOLECULE[HYPHENED INTERVALS,COMMA_SEPARATED]

use warnings;
use strict;

my %gem_blacklist;

my %gem_start;
my %gem_readstart;
my %gem_readend;
my %gem;
my %gem_reads;
my %gem_bxcount;
	
my $max_gap = 50000;

#opendir(POS, "work_folder");
#while (my $file = readdir(POS)) {

my @tmp = split("/", $ARGV[1]);
my $file = pop @tmp;
my $outfile = $ARGV[2];

#	next if ($file=~/\.bak/);
#next unless ($file=~/_inform\.SNPs\.new.chr\d+.bed/);	
#	next unless ($file=~/730/);
print "$file\n";
##Figure out the start and end of a file
my $chr = `cut -f 1 $file | uniq`;
chomp($chr);
my $start = `head -1 $file | cut -f 3`;
chomp($start);
my $end = `tail -1 $file | cut -f 3`;
chomp($end);
my $line = 0;

open(OUT, ">$outfile");
sub write_out() {
	my $call_pos = $_[0];
	#	Writing out positions from $last_writeout up to $start_pos\n";
	GM: foreach my $gm (sort {$gem_start{$a} <=> $gem_start{$b}} keys %gem_start) {
		#Decide if this is too close to the end to write out
		foreach my $p (sort {$b <=> $a} keys %{$gem{$gm}}) {
			next GM if ($call_pos - $p <= $max_gap); 	
		}
		my $gm_str="";
		my $qstr="";
		my @switches;
		my $pos_str="";
		my $last_pos;
		my $last_state=-1;
		my $state = -1;
		foreach my $p (sort {$a <=> $b} keys %{$gem{$gm}}) {
			my $tmp_q = abs($gem{$gm}{$p});	
			$pos_str.=$p.",";
			#There are enough ASCII characters up to 126 to cover 93
			$tmp_q = 93 if ($tmp_q > 93);
			$qstr.= chr($tmp_q+33);
			
			if ($gem{$gm}{$p} > 30) {
				$gm_str.= 1;
				$state = 1;
			} elsif ($gem{$gm}{$p} < -30) {
				$gm_str.= 0;
				$state = 0;
			} elsif ($gem{$gm}{$p} > 10) {
				$gm_str.= ">";
				$state = 1;
			} elsif ($gem{$gm}{$p} < -10) {
				$gm_str.= "<";
				$state = 0;
			} else {
				$gm_str .= "N";
			}
			if (($last_state != -1) && ($last_state != $state)) {
				push @switches, $last_pos."-".$p;
			}
			$last_pos = $p;
			$last_state = $state;
		}
		my $switch_str = join(",",@switches);
		my $read_str = join(",", keys %{$gem_reads{$gm}});
		#print $gm."\n";
		#foreach my $gr (sort keys %{$gem_reads{$gm}}) {
		#		print "$gr\n";
		#	}
		#my $Lkj = <STDIN>;
		chop($pos_str);
		print OUT "$gm\t$chr\t$gem_start{$gm}\t$last_pos\t$gm_str\t$qstr\t".(length($qstr))."\t".(keys %{$gem_reads{$gm}})."\t$switch_str\t$pos_str\t$read_str\n";
		
		#Clean up;
		delete $gem{$gm};
		delete $gem_reads{$gm};
		delete $gem_start{$gm};
		delete $gem_readstart{$gm};
		delete $gem_readend{$gm};
	}
}



#Read in the positions that need to be retained;
my %snp_pos;
@tmp = `awk '{print \$3"\t"\$4"\t"\$5"\t"NR-1}' $file`;
map {chomp; /^(\S+)\t(\S+)\t(\S+)/;$snp_pos{$1} = [ ($2,$3,$4) ]} @tmp;
#my @pos = `cut -f 3 work_folder/$file`;
#	chomp(@pos);
#	print "SNPs:".(@pos)."\n";
##Per file, read out the entire chunk of the bam file

my $last_writeout;
my $v = `pwd`;
print $v;
print "samtools view -h $ARGV[0] -L $file $chr:$start-$end |samtools calmd -h -S - /fml/chones/genome/gbdb/mm10/mm10.fa| /samtools view -h - |\n";
open (IN, "samtools view -F 3840 -h $ARGV[0] -L $file $chr:$start-$end |samtools calmd -h -S - /fml/chones/genome/gbdb/mm10/mm10.fa | samtools view -h - |");
READ: while (<IN>) {
#ST-J00101:68:HKHWYBBXX:2:1115:10338:10440	147	6	100317514	60	150M	=	100317381	-283	=======================T===============================================A=====================================================C==T=C===C=====T=========	<-<JJAFF<JJJFJFJJFJJFJJFJJJJJFJJJJJJJJJJJFJJJJJJFFJJJJJJJFJFFJJJJJJFJJJJFJJFJJJJJJJJFJJJ-JFJAJJJFJJJJJJJJJFJFJFJJJJJFJJJJJJJJJJJFJJJJFJJFFJJFJJJJ<AAAA	QT:Z:AAFFFJJJ	BC:Z:AAGATCAT	QX:Z:AAFAFJJJJJJJJJJJ	AM:A:0	XM:A:0	RX:Z:GTGGGTCGTTGTACAC	AS:f:-21	XS:f:-88.6717	BX:Z:GTGGGTCGTTGTACAC-1	XT:i:0	OM:i:60	RG:Z:MitoMeio_S18_5d_CTRL-DMSO:LibraryNotSpecified:1:unknown_fc:0	PS:i:99705759	HP:i:2	PC:i:57	MI:i:81353161	NM:i:7	MD:Z:23C47G53A2C1T3A5C9
	#next unless (/BX:Z:(\S+).+MI:i:(\d+)/);
	
	next unless (/BX:Z:(\S+)/);
		
	print "Processed $line...\n" if ($line % 20000 == 0);
	my @temp = split "\t";
	my $bx_tag = $1;
	
	if (exists($gem_bxcount{$bx_tag."/0"})) {
		$bx_tag = $bx_tag."/".$#{$gem_bxcount{$bx_tag."/0"}};
	} else {
		$bx_tag .= "/0";
	}
	
	my @pos;
	
	my ($start_pos,$flag, $cigar,$read,$qual) = ($temp[3],$temp[1], $temp[5],$temp[9],$temp[10]);
	#Populate the first start of the process to make sure we only try to write files out after 2Mbp from the first position
	$last_writeout = $start_pos if (!($last_writeout));
	
	my $tcigar=$cigar;
	my $tread=$read;
	my $tqual=$qual;
	my $cig;
	my $len;
	my $fpos=0;
	#print "READ: $read\n";
	#print "CIGAR: $cigar\n";
	
	
	my $refpos=0;
	
	my %link;
	$link{$refpos+$start_pos} = $fpos;
	my @cig = split(/([MHIDNSPX=])/, $cigar);
	for (my $c = 0; $c<=$#cig; $c+=2) {
		#print "$cig[$c]$cig[$c+1]\n";
		next if ($cig[$c] eq "*");
		for (my $i = 0; $i < $cig[$c]; $i++) {
			if ($cig[$c+1] eq "S") {
				$fpos++;
			} elsif ($cig[$c+1] eq "D") {
				$refpos++;
			} elsif ($cig[$c+1] eq "I") {
				$fpos++;
				#$refpos++;
			} elsif ($cig[$c+1] eq "M") {
				$fpos++;
				$refpos++;
			}
			$link{$refpos+$start_pos} = $fpos unless (($c== $#cig-1) && ($cig[$c+1] eq "S"));	
		}
	}
	my $end_pos = $refpos + $start_pos;
	my $sp = $start_pos-1;
#	my $ref = `twoBitToFa /fml/chones/genome/gbdb/mm10/mm10.2bit stdout -seq=$chr -start=$sp -end=$end_pos | grep -v ">" | paste -d, -s | sed 's/,//g'`;
	
		
	foreach my $p ($start_pos..$refpos+$start_pos-1) {
		my $arrow = "";
		my $readbase = substr($read,$link{$p},1);
		my $qualscore;
#		print "$p -> $link{$p}\t".(substr($ref,$p-$start_pos,1))."\t".(substr($read,$link{$p},1))."\t".(substr($qual,$link{$p},1))." $arrow\n";
		if ($snp_pos{$p}) {
			if ($readbase eq $snp_pos{$p}[0]) {
				$qualscore = -(ord(substr($qual,$link{$p},1))-33); 
			} elsif ($readbase eq $snp_pos{$p}[1]) {
				$qualscore = ord(substr($qual,$link{$p},1))-33; 
			} else {
				#		print "CHECK!!!!\n";
				next READ;
			}
			push @pos, [ ($p, $snp_pos{$p}[0],$snp_pos{$p}[1],$qualscore) ]; 
		}
	}
	#ST-J00101:68:HKHWYBBXX:2:2221:10734:27619	129	15	61882702	1	32S104M4I10M	18	16908642	0	CCACTTTTTGGCCCTGGGGTTCCCCTGTACTGGGGCCAAAAAGTGGGAGTGGGTGGGTAGGGGCGTGGGTGGGAGGGTATGGGGGACTTTTGGGATAGCATTGGAAATGTAAATGAGGAAAATACCTAATAAAAAGTATATATTAAAAAA	AAFFFJFJJJJJFJJJJFJJFJJJ<JFJFFFJJJJJJJFAFAFJJJJ-AJJJJJ<7FF-AFFJ-AJJJJFA<FJJ-7FF<FFJJJFFJ-F-<JJ--<FJFAAFF---<FFJAJ-<<FJJAF<F<FFJ<7A<A77-7<<FFFAJJ7---<-	QT:Z:AAFFFFJJ	BC:Z:CCTGAAGG	QX:Z:A<AFFJJAF7JJJJJJ	AM:A:0	XM:A:0	RX:Z:AAGTGCTAGATAAACC	AS:f:-34	XS:f:-34	SA:Z:10,90398195,+,36M114H,0,0;	BX:Z:AAGTGCTAGATAAACC-1	XT:i:0	OM:i:1	RG:Z:MitoMeio_S18_5d_CTRL-DMSO:LibraryNotSpecified:1:unknown_fc:0	NM:i:4	MD:Z:114
	
	#Populate the GEM annotation
	
	my $delflag = 0;
		#GEM already exists
	if (exists($gem{$bx_tag}) && (@pos)) {
			#GEM already exists - check gap
			if ($start_pos > $gem_readend{$bx_tag} + $max_gap) {
				#Gap is larger than max_gap -> iterate moleucle number
				my ($bx_stem, $mol) = ($1, $2) if ($bx_tag=~/([ABCD0-9]+)\/(\d+)$/);
				$bx_tag = $bx_stem."/".($mol+1);
				($gem_readstart{$bx_tag}, $gem_readend{$bx_tag})=($start_pos, $end_pos);
				push @{$gem_bxcount{$bx_stem."/0"}}, $temp[0];
				$gem_start{$bx_tag} = $pos[0][0];
			} else {
				$gem_readend{$bx_tag}=$end_pos;
			} 
			$gem_reads{$bx_tag}{$start_pos."-".$end_pos}++;
			#- check if all positions match existing annotation
			V_LOOP: foreach my $val (0..$#pos) {
				if (exists($gem{$bx_tag}{$pos[$val][0]})) {
				#Sum up the log-QScore
					$gem{$bx_tag}{$pos[$val][0]} += $pos[$val][3];
				
				#		if ($gem{$bx_tag}{$pos[$val][0]} != $pos[$val][3]) {
				
				#GEMcode states disagrees! Delete GEM entry!
				#	$delflag = 1;
				#	last V_LOOP;
				#}
				} else {
				#New position for this GEMcode - add!
					$gem{$bx_tag}{$pos[$val][0]} = $pos[$val][3];
				}
			}
	} else {
		#Case - GEM doesn't exist and is NOT on the blacklist
		if (!(exists($gem_blacklist{$bx_tag})) && (@pos)) {
			#Populate GEM start
				($gem_readstart{$bx_tag}, $gem_readend{$bx_tag})=($start_pos, $end_pos);
				$gem_start{$bx_tag} = $pos[0][0];
				my ($bx_stem, $mol) = ($1, $2) if ($bx_tag=~/([ABCD0-9]+)\/(\d+)$/);
				push @{$gem_bxcount{$bx_stem."/0"}}, $temp[0];
				#Populate GEM position hash
				foreach my $val (0..$#pos) {
					$gem{$bx_tag}{$pos[$val][0]} = $pos[$val][3];
				}
				
				$gem_reads{$bx_tag}{$start_pos."-".$end_pos}++;
			}
	}
	#Process deletion of a GEM and add to blacklist
	if ($delflag) {
		delete $gem{$bx_tag};
		delete $gem_start{$bx_tag};
		delete $gem_readstart{$bx_tag};
		delete $gem_readend{$bx_tag};
		$gem_blacklist{$bx_tag}++;
	}
	
	$line++;
	if ($line % 100000 == 0) {
		print "Partial write-out\n";
		&write_out($start_pos);
	}
}
close (IN);
	
### Early code to make sure the region as defined by these intervals cover only a single chromosome, not two. Once it's done, it's not needed to run anymore
#	my @chr = `cut -f 1 work_folder/$file | uniq`;
#	chop(@chr);
#	print "@chr-";
#	print "\n";
	#if ($#chr == 1) {
	#	system("awk '/^".($chr[1])."\t/' work_folder/$file > work_folder/$file.b\n");
	#	system("awk '/^".($chr[0])."\t/' work_folder/$file > work_folder/$file.a\n");
	#	system("mv work_folder/$file work_folder/$file.bak");
	#}
#}
#closedir(POS);

print "Final write-out\n";
&write_out($end + $max_gap * 5);
close(OUT);
