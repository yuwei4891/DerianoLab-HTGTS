# DerianoLab-LAM-HTGTS
Codes for LAM-HTGTS analysis

These coded are used for LAM-HTGTS sequencing result analysis. 

The common tools were used including: PERL; BWA; samtools; ea-utils; circos; bedtools; R

The procedures for processing sequencing reads from Mi-Seq sequencing fastq reads:




Remove adapter sequence

  ###### fastq-mcf blue_adapter.fa r2.$sample.fq -o r2.$sample.trimmed.fq
  ######    Output: 
        #  r2.$sample.trimmed.fq
        
Trim out HTGTS bait sequence in read 1 and 2

  ###### perl filter_20bp_from_primer.pl bait_20bp.fa r1.$sample.fq r2.$sample.trimmed.fq
  ###### Output:
        #  r1.$sample.filter20bait.fq 
        #  r2.$sample.trimmed.fq.filter20bait


Mapping fastq sequence

  ###### bwa mem -t $CPUS_PER_TASK reference_mm9_index r2.$sample.trimmed.fq.filter20bait >r2.$sample.trimmed.fq.filter20bait.sam


Classify mapped read2 sequence into different catergories

  ###### perl filter_germline_second_end_BWA.pl mm9_germline.region r2.$sample.trimmed.fq.filter20bait.sam
  ###### Output:         
        #  r2.$sample.germline_filtered_SplitRead
        #  r2.$sample.germline_filtered_noSplitRead
        #  r2.$sample.germline
        #  r2.$sample.germline_filtered_pure
        #  r2.$sample.germline_filtered_long_unmap_part  
