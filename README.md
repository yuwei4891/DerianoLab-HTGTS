# DerianoLab-LAM-HTGTS
Codes for LAM-HTGTS analysis

These codes are the pipeline for LAM-HTGTS sequencing analysis. 

The common tools were used including: Perl; BWA; samtools; ea-utils; circos; bedtools; R

The procedures for processing sequencing reads from MiSeq sequencing fastq reads:




Step 1. Remove adapter sequence:

###### fastq-mcf blue_adapter.fa r2.$sample.fq -o r2.$sample.trimmed.fq
######    Output: 
        #  r2.$sample.trimmed.fq
        
Step 2. Trim out HTGTS bait sequence in read 1 and 2:

###### perl filter_20bp_from_primer.pl bait_20bp.fa r1.$sample.fq r2.$sample.trimmed.fq
###### Output:
        #  r1.$sample.filter20bait.fq 
        #  r2.$sample.trimmed.fq.filter20bait


Step 3. Mapping fastq sequence by BWA:

###### bwa mem -t $CPUS_PER_TASK reference_mm9_index r2.$sample.trimmed.fq.filter20bait >r2.$sample.trimmed.fq.filter20bait.sam


Step 4. Classify mapped read2 sequence into different catergories including split mapped and germline:

###### perl filter_germline_second_end_BWA.pl mm9_germline.region r2.$sample.trimmed.fq.filter20bait.sam
###### Output:         
        #  *.germline_filtered_SplitRead
        #  *.germline_filtered_noSplitRead
        #  *.germline
        #  *.germline_filtered_pure
        #  *.germline_filtered_long_unmap_part  
        
Step 5. Filter out reads mapped to multiple repetitive genomic regions:

###### perl filter_repeatMasker_BWA.pl repeatMasker_mm9.txt mm9_germline.region r2.$sample.trimmed.fq.filter20bait.sam.germline_filtered_SplitRead
###### Output:
        #  *.repeatOut

Step 6. Dissect Split mapped read2 with part of sequence mapped between IgkJ4 RSS and bait primer coordinate were defined as split mapping joins with IgkJ4 bait.

###### perl dissect_read2.pl r2.$sample.trimmed.fq.filter20bait.sam.germline_filtered_SplitRead.repeatOut 1 mm9_germline.region
###### Output:
        #  *.keepDuplicate                 
        #  *.keepDuplicate_noIgKj4         
        #  *.keepDuplicate_noIgKj4_all     
        #  *.keepDuplicate_multipleMapped  
        #  *.keepDuplicate_onlyIgKj4       

Step 7. Visualize the mapped split reads in IGV :

###### cat sam.header r2.$sample.trimmed.fq.filter20bait.sam.germline_filtered_SplitRead.repeatOut.keepDuplicate_noIgKj4 >r2.$sample.trimmed.fq.filter20bait.sam.germline_filtered_SplitRead.repeatOut.keepDuplicate_noIgKj4.sam
###### samtools view -bS r2.$sample.trimmed.fq.filter20bait.sam.germline_filtered_SplitRead.repeatOut.keepDuplicate_noIgKj4.sam >r2.$sample.trimmed.fq.filter20bait.sam.germline_filtered_SplitRead.repeatOut.keepDuplicate_noIgKj4.bam
###### bedtools bamtobed -i r2.$sample.trimmed.fq.filter20bait.sam.germline_filtered_SplitRead.repeatOut.keepDuplicate_noIgKj4.bam >r2.$sample.trimmed.fq.filter20bait.sam.germline_filtered_SplitRead.repeatOut.keepDuplicate_noIgKj4.bam.bed
######  samtools sort -f r2.$sample.trimmed.fq.filter20bait.sam.germline_filtered_SplitRead.repeatOut.keepDuplicate_noIgKj4.bam r2.$sample.trimmed.fq.filter20bait.sam.germline_filtered_SplitRead.repeatOut.keepDuplicate_noIgKj4.bam.sorted.bam
###### bam2wig r2.$sample.trimmed.fq.filter20bait.sam.germline_filtered_SplitRead.repeatOut.keepDuplicate_noIgKj4.bam.sorted.bam >r2.$sample.trimmed.fq.filter20bait.sam.germline_filtered_SplitRead.repeatOut.keepDuplicate_noIgKj4.bam.sorted.bam.wig 
 
###### Output:
        #  *.bam.wig
 
Split mapping joins with IgkJ4 bait were further classified according to prey genomic location:


Step 8. Visualize the V(D)J recombination within IgK locus by Circos:

###### perl convert_bed_to_Circos_plot_input_file_IgK_only.pl r2.$sample.trimmed.fq.filter20bait.sam.germline_filtered_SplitRead.repeatOut.keepDuplicate_noIgKj4.bam.bed <reads_number_cutoff>

###### Output:
        #  *.links
        #  *.hismap
###### perl configure_Circos_conf.pl r2.$sample.trimmed.fq.filter20bait.sam.germline_filtered_SplitRead.repeatOut.keepDuplicate_noIgKj4.bam.bed circos_IgK.conf links.conf image.conf hisplots.conf

###### Output:
        #  *.circos.conf
        #  *.links.conf
        #  *.image.conf
        #  *.hisplots.conf
###### circos -conf r2.$sample.trimmed.fq.filter20bait.sam.germline_filtered_SplitRead.repeatOut.keepDuplicate_noIgKj4.bam.bed.circos.conf

###### Output:
        #  *.circos.png


Step 9. Bar graphs plotting the extent of resected DNA for each end joining product relative to the RSS of IgK V segments

###### perl dissect_V_segment_resection_length.pl r2.$sample.trimmed.fq.filter20bait.sam.germline_filtered_SplitRead.repeatOut.keepDuplicate Igk.segments.coordinate.mm9.bed
###### Output:
	#			*.V_segment_resection_micro_length			### resection length and microhomology length
	#			*.V_segment_CE_resection_length				# resection length at coding end
	#			*.V_segment_SE_resection_length				# resection length at signal end
	#			*.microhomology						# read with microhomology
	#			*.insert						# read with insertion
	#			*.microhomology.length					# length of microhomology for each read with microhomology
	#			*.microhomology.length.percentage			# percentage of each microhomology length
	#			*.insert.length						# length of insertion
	#			*.micro_insert_resect_summary				# summary of number of reads with microhomology and insertion and resection in V 
