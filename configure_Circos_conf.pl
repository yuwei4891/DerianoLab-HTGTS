if(@ARGV<1){
	print "This scipt is to configure the Circos.conf file
		Usage: perl *pl r2.myc_combined_cRAG2_XLF.trimmed.fq.filter20bait.sam.germline_filtered_SplitRead.repeatOut.noDuplicate_noIgKj4.bam.bed circos.conf links.conf image.conf hisplots.conf
		Output: r2.myc_combined_cRAG2_XLF.trimmed.fq.filter20bait.sam.germline_filtered_SplitRead.repeatOut.noDuplicate_noIgKj4.bam.bed.circos.conf
			r2.myc_combined_cRAG2_XLF.trimmed.fq.filter20bait.sam.germline_filtered_SplitRead.repeatOut.noDuplicate_noIgKj4.bam.bed.links.conf
			r2.myc_combined_cRAG2_XLF.trimmed.fq.filter20bait.sam.germline_filtered_SplitRead.repeatOut.noDuplicate_noIgKj4.bam.bed.image.conf
			r2.myc_combined_cRAG2_XLF.trimmed.fq.filter20bait.sam.germline_filtered_SplitRead.repeatOut.noDuplicate_noIgKj4.bam.bed.hisplots.conf\n";
}

open(IN_circos,"$ARGV[1]")||die("Can not open circos conf file\n");
open(IN_link,"$ARGV[2]")||die("Can not open link file\n");
open(IN_image,"$ARGV[3]")||die("Can not open image conf file\n");
open(IN_hisplot,"$ARGV[4]")||die("Can not open hisplot conf file\n");
my $out_circos=$ARGV[0].".circos.conf";
my $out_link=$ARGV[0].".links.conf";
my $out_image=$ARGV[0].".image.conf";
my $out_hisplot=$ARGV[0].".hisplots.conf";
open(OUT_circos,">$out_circos")||die("Can not open out_circos file\n");
open(OUT_link,">$out_link")||die("Can not open out_link file\n");
open(OUT_image,">$out_image")||die("Can not open out_image file\n");
open(OUT_hisplot,">$out_hisplot")||die("Can not open out_hisplot file\n");

while($line=<IN_circos>){
	chomp $line;
	if($line=~/image\.conf/){
		$line=~s/image\.conf/$ARGV[0]\.image\.conf/;
	}
	if($line=~/links\.conf/){
		$line=~s/links\.conf/$ARGV[0]\.links\.conf/;
	}
	if($line=~/hisplots\.conf/){
		$line=~s/hisplots\.conf/$ARGV[0]\.hisplots\.conf/;
	}
	print OUT_circos "$line\n";
}
while($line=<IN_link>){
	chomp $line;
	if($line=~/links\.txt/){
		$line=~s/links\.txt/$ARGV[0]\.links\.txt/;
	}
	print OUT_link "$line\n";
}
while($line=<IN_image>){
	chomp $line;
	if($line=~/circos\.png/){
		$line=~s/circos\.png/$ARGV[0]\.circos\.png/;
	}
	print OUT_image "$line\n";
}
while($line=<IN_hisplot>){
	chomp $line;
	if($line=~/hismap\.txt/){
		$line=~s/hismap\.txt/$ARGV[0]\.hismap\.txt/;
	}
	print OUT_hisplot "$line\n";
}




