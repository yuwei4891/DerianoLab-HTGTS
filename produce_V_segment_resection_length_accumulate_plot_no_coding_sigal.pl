if(@ARGV<1){
	print " This script is to produce the accumulated resection plot, but without seperate coding and signal resection !!!!
		Usage: perl *pl 
					r2.JK_combined_p53.trimmed.fq.filter20bait.sam.germline_filtered_SplitRead.repeatOut.keepDuplicate.V_segment_resection_length
		Output:
					*.resection_length_plot_combined_signal_coding			# output each line : start end length pre_accumulate accumulate micro 
		\n";
}

open(IN,"$ARGV[0]")||die("Can not open in file\n");
my $out=$ARGV[0].".resection_length_plot_combined_signal_coding";
open(OUT_resection_plot,">$out")||die("Can not open out file\n");

#15
#-25

$constant=20000;		# set the constant number of total reads
$distance_cutoff=1000;

while($line=<IN>){
	chomp $line;
	$total_resection++;
}
#print "$total_resection\n";
open(IN,"$ARGV[0]")||die("Can not open in file\n");
while($line=<IN>){
	chomp $line;
	@split=split /\t/,$line;
	$resection_length=abs($split[0]);		# only get absolute value
	$chr="chr6";
	$micro=$split[1];		# eg : 15	2   resectionLength\t microhomology length

	if($split[0]>0){
		$start=0;
		$end=$split[0];
	}else{
#		$start=$split[0];
#		$end=0;
		$start=0;			# here !!! the signal resection are convert to the same as coding resection, both are positive
		$end=0-$split[0];
	}

	$value=$constant/$total_resection*1;					# this is the relative read number to 20000
	if(exists($hash_sort{$resection_length})){
		$hash_sort{$resection_length}=$hash_sort{$resection_length}.'+'.$start.':'.$end.':'.$value.':'.$micro;
	}
	else{
		$hash_sort{$resection_length}=$start.':'.$end.':'.$value.':'.$micro;
	}
#	print "$resection_length\t$hash_sort{$resection_length}\n";
}

#while(($key,$value)=each %hash_resection){
#	@temp=split /\./,$key;
#	print OUT_resection_plot "$key\t$value\n";
#	$chr=$temp[0];
#	$resection_start=$temp[1];														# this is the resection start site
#	$resection_end=$temp[2];														# this is the resection end site
#	$resection_region=$chr.'.'.$resection_start.'.'.$resection_end;
#	$relative_resection_number{$resection_region}=$constant/$total_resection*1;			# now !!! 1 read convert into relative to 20000
#	$hash_plot{$resection_region}=$relative_resection_number{$resection_region};
#	$total_rel+=$relative_resection_number{$resection_region};
#}
#print  OUT_resection_plot "$total_rel\n";
#  record the corresponding resection start end with the relative read counts
# eg. 97	61818268-61818365-1.15234754802048+61818267-61818364-0.0720217217512802+61818264-61818361-0.0720217217512802
#my $ID=0;
#while(($key,$value) =each % hash_plot){
#	@temp=split /\./,$key;
#	$resection_length=$temp[2]-$temp[1];
#	$sort_length{$ID}=$resection_length;
#	$ID++;
#	if(exists($hash_sort{$resection_length})){
#		$hash_sort{$resection_length}=$hash_sort{$resection_length}.'+'.$temp[1].'-'.$temp[2].'-'.$value;
#	}
#	else{
#		$hash_sort{$resection_length}=$temp[1].'-'.$temp[2].'-'.$value;
#	}
#}						

####  sort the length from short to big and push into array,
my $i=0;
foreach my $resection_length (sort { $a <=> $b } keys %hash_sort) {
#	72	61818263-61818336-10.6458481192335
#	73	61818269-61818341-1.77430801987225
	push @sorted_resection_length, $resection_length;
}

###  pop out the resection length from big to short
# eg. 	73	61818263-61818336-10.6458481192335
#		72	61818269-61818341-1.77430801987225
$sorted_resection_length=@sorted_resection_length;
my $accumulate=0;
for($i=0;$i<$sorted_resection_length;$i++){
	$temp=pop @sorted_resection_length;
#	print "$temp\n$hash_sort{$temp}\n";
#	print OUT_resection_plot "$temp\t$hash_sort{$temp}\n";
#	72	61818269-61818341-1.77430801987225
#	71	61818241-61818312-21.291696238467+61818246-61818317-1.77430801987225+61818273-61818344-3.5486160397445
	if($hash_sort{$temp}=~/\+/){#print "$temp\n$hash_sort{$temp}\n";
		@temp_split=split /\+/, $hash_sort{$temp};
		$temp_split=@temp_split;
		for($j=0;$j<$temp_split;$j++){
			@temp_2nd_split=split /\:/,$temp_split[$j];#print "@temp_2nd_split\n";
			$pre_accumulate=$accumulate;
			$accumulate+=$temp_2nd_split[2];
			$start=$temp_2nd_split[0];
			$end=$temp_2nd_split[1];
			$micro=$temp_2nd_split[3];
			$length=$end-$start;
			$resection_region=$chr.'.'.$start.'.'.$end;
			
		#	print OUT_resection_plot "===$start\t$end\t$length\t$pre_accumulate\t$accumulate\t$micro\n";
			print OUT_resection_plot "$start\t$end\t$length\t$pre_accumulate\t$accumulate\t$micro\n";	
		#	$value++;
		}
	}
	else{
		@temp_split=split /\:/, $hash_sort{$temp};
		$pre_accumulate=$accumulate;
		$accumulate+=$temp_split[2];
		$start=$temp_split[0];
		$end=$temp_split[1];
		$length=$end-$start;
		$resection_region=$chr.'.'.$start.'.'.$end;
		$micro=$temp_split[3];
	#	print OUT_resection_plot "$start\t$end\t$length\t$pre_accumulate\t$accumulate\t$micro\n";
		print OUT_resection_plot "$start\t$end\t$length\t$pre_accumulate\t$accumulate\t$micro\n";
	#	print " $hash_sort{$temp}\t\t$micro\n";$value++;
	}
}

#print "$value\n";
