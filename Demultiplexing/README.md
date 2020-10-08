Our general strategy for demultiplexing goes as follows:

- add extra index cycles during the sequencing run for the extended barcode information
- demultiplex the reads and create fastq files by including the barcodes as comment tags in the header of the fastq: BX, RX and QX tags. BX being the A/B/C/DXX human-readable beadTags, RX being the raw reads themselves and QX the quality string.

`bcl2fastq --use-bases-mask=Y150,I13,I12,Y150 --create-fastq-for-index-reads -r 20 -w 20 -d 20 -p 40 -R <run folder> \
--tiles s_<REGEX> --output-dir=<outdir> --interop-dir=<Interop_dir> --reports-dir=<Reports_dir> --stats-dir=<Stats_dir> 2> [STDERR.log]
./tag_fastq_13plus12.o <prefix> <OUTDIR/out_prefix> 2>[log_file]`

- adapter trimming 

`cutadapt -a CTGTCTCTTATACACATCT -g AGATGTGTATAAGAGACAG -A CTGTCTCTTATACACATCT -G AGATGTGTATAAGAGACAG \
        --cores=<cores> -O 5 \
        -o <read_1_out.fastq.cutadapt.gz> -p <read_1_out.fastq.cutadapt.gz> --pair-filter both \
        $dir/$fbname\_R1_001.fastq.gz $dir/$fbname\_R2_001.fastq.gz`

`cutadapt -m 30 \
       -o <read_1_out.fastq.cutadapt.gz> -p <read_2_out.fastq.cutadapt.gz> \
        <read_1_out.fastq.cutadapt.1.gz> <read_2_out.fastq.cutadapt.1.gz> \
        --too-short-output=<read_1_out.tooshort.fastq.gz> --too-short-paired-output=<read_2_out.tooshort.fastq.gz>`

- read placement, using the -C switch to include the extra comment tags

`bwa mem -C -t 50 <helera1_demo_dir>/Heliconius_erato_demophoon_v1.fa \
        $file ${file/R1_001/R2_001} \
        -R "@RG\tID:$fbname\tSM:$fbname\tLB:$fbname\tPL:Illumina.HiSeq3000.2x150" |
        samtools view -bh - > /tmp/mkucka/$fbname.erato.bam`

`samtools sort \
        -@ 50 -l 9 \
        -T /tmp/mkucka/$fbname.tmpsort \
        -o /tmp/mkucka/$fbname.erato.sorted.bam \
        /tmp/mkucka/$fbname.erato.bam`

- MarkDuplicates, using BX-aware options

`java -Xmx12g -XX:ParallelGCThreads=64 -jar <picard_dir>/picard.jar MarkDuplicates \
        I=$file \
        O=/tmp/mkucka/$fbname.pMarkdup.bam \
        M=$dir/$fbname.pMarkdup.metrics \
CREATE_INDEX=TRUE READ_ONE_BARCODE_TAG=BX READ_TWO_BARCODE_TAG=BX VALIDATION_STRINGENCY=LENIENT`
