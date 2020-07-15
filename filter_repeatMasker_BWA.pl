if(@ARGV<1){
	print "		This script is to filter out SAM mapping result that locate in multiple repeat regions and the window size of repeat is 100bp
			Usage: 
				perl *pl repeatMasker.txt mm9_germline.region *.sam.germline_filtered_SplitRead
			Output:
				*repeatOut\n";
}

open(IN_repeat,"$ARGV[0]")||die("Can not open in repeat file\n");
open(IN_germline,"$ARGV[1]")||die("Can not open in germline file\n");
open(IN_SAM,"$ARGV[2]")||die("Can not open sam file\n");

my $out=$ARGV[2].".repeatOut";
open(OUT,">$out")||die("Can not open out file\n");

my $out_repeat=$ARGV[2].".onlyrepeat";
open(OUT_repeat,">$out_repeat")||die("Can not open out repeat file\n");

#bin	swScore	milliDiv	milliDel	milliIns	genoName	genoStart	genoEnd	genoLeft	strand	repName	repClass	repFamily	repStart	repEnd	repLeft	id
#607	687	174	0	0	chr1	3000001	3000156	-194195276	-	L1_Mur2	LINE	L1	-4310	1567	1413	1

#M01626:162:000000000-AD53D:1:1101:14199:1434_2:N:0:1	16	chr6	68387925	250	120S85M	*	0	0	GCTCAGTTAGCCAAAATGTCACAAATTCACACAAGTTACCCAAACAGAACCAAAACGTCACAAGTAAATGAGCAAAAGTCTACTTACGTTTTATTTCCA
#ACTTTGTCCCCGAGCCGAACATGATATAAGTCATGACACAAACCTCCCAGGGGTTCAGCAGAGTAAGTCTGTGCTGCCTCAGCTGCTCCTGCATCTGCTTATATTC	;0;0C9GFBF@FFFFFFB0FFGGFBFBABGFHC0FG;:.FEHCHHHHEHGHEC.ECG0EGGHGF1HHHDF<1FHFEHHHFHGAHFHGHHGBHHHHGHHGFBGFCCC@>/?CF@E?
#2>DFFGDFHBFHDFHHFFGGFHEF0/E>F/EGAAEEHDBEFGGB1A1GBGAHHGFB//AE/HFF0FGA/FC0DFFFDHGFA2GFGAA10G	AS:i:85	NM:i:0	MD:Z:85	YF:H:2D	YI:i:1	YP:i:2	YS:i:0

my $window=5;

while($line=<IN_repeat>){
	chomp $line;
	if($line=~/bin/){
	}
	else{
		@split=split /\t/,$line;
		$chr=$split[5];
		$start=int($split[6]/$window); # 100bp as window size;
		$end=int($split[7]/$window);
		
		for($i=$start;$i<$end+1;$i++){
			$pos=$chr.'.'.$i;# print "$chr\t$pos\n";
			$hash_repeat{$pos}="repeat.".$split[11].".".$split[6].'.'.$split[7]; # mark the position as repeat
		}
	}
}

open(IN_germline,"$ARGV[1]")||die("Can not open in germline file\n");
while($line=<IN_germline>){
	chomp $line;
	@split=split /\t/,$line;
	if($line=~/mm10/){
		$hash_mm10_start=$split[1];
		$hash_mm10_end=$split[2];
	}
	else{
		$hash_mm9_start=$split[1];
		$hash_mm9_end=$split[2];
	}
}

while($line=<IN_SAM>){
	chomp $line;
	@split=split /\t/,$line;
	$head=$split[0];
	if(exists($hash{$head})){
		$hash{$head}=$hash{$head}."\n".$line;
	}
	else{
		$hash{$head}=$line;
	}
#	print "$head\n$hash{$head}\n";
}

open(IN_SAM,"$ARGV[2]")||die("Can not open sam file\n");
while($line=<IN_SAM>){
	chomp $line;
	@split=split /\t/,$line;
	$head=$split[0];
	$chr=$split[2];# print "$chr\n";
	$length=length($split[8]);
	$start=int($split[3]/$window);
	$end=int(($split[3]+$length)/$window);
	$mark=1;
	
	if(($chr eq "chr6")&&($split[3]<$hash_mm9_end)&&($split[3]>$hash_mm9_start)){
	}
	else{#print "$line\n";
		for($i=$start;$i<$end+1;$i++){
			$pos=$chr.'.'.$i; #print "$pos\n";
			if($hash_repeat{$pos}=~/repeat/){
				$mark=0;print OUT_repeat "$line\n"; # locate in repeat regions
			}
			else{
			}
		}
		if($mark!=0){
		#	print OUT "$line\n";
			print OUT "$hash{$head}\n";
		}
	}
}

