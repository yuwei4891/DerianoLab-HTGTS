# DerianoLab-LAM-HTGTS
Code for LAM-HTGTS analysis

These coded are used for LAM-HTGTS sequencing result analysis. The common tools were used including: 

PERL; BWA; samtools; ea-utils; circos; bedtools; R

The procedures for processing sequencing reads from Mi-Seq sequencing fastq reads:


Remove adapter sequence, input: read2 fastq sequence, output: adapter trimmed read2 fastq sequence

###### fastq-mcf blue_adapter.fa r2.$sample.fq -o r2.$sample.trimmed.fq


Trim out HTGTS bait sequence in read 1 and 2, output: r1.$sample.filter20bait.fq r2.$sample.trimmed.fq.filter20bait

###### perl filter_20bp_from_primer.pl bait_20bp.fa r1.$sample.fq r2.$sample.trimmed.fq


Mapping fastq sequence

###### bwa mem -t $CPUS_PER_TASK reference_mm9_index r2.$sample.trimmed.fq.filter20bait >r2.$sample.trimmed.fq.filter20bait.sam


Classify mapped read2 sequencing into 
