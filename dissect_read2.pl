if(@ARGV<1){
	print "This script is to find multiple cigar
		Usage: perl *pl r2.$sample.trimmed.fq.filter20bait.sam.germline_filtered_SplitRead.repeatOut STep mm9_germline.region \n
	#	Output
	#		*.keepDuplicate			# split mapped with one part JK4 another part is unknown
	#		*.keepDuplicate_noIgKj4		# the part of reads without JK4 mapping sequence
	#		*.keepDuplicate_noIgKj4_all	# both part of split reads without JK4
	#		*.keepDuplicate_multipleMapped	# part of reads are multiple pieces and are mapped
	#		*.keepDuplicate_onlyIgKj4	# the part of reads with JK4 mapping sequence\n";
}

#M01626:162:000000000-AD53D:1:1103:13125:25048_2:N:0:1	16	chr6	70024021	250	125S32M	*	0	0	GCTCATTTAGCCAAAATGTCACAAATTCACACAAGTTACCCAAACAGAACCAAAACGTCACAAGTAAATGAGCAAAAGTCTACTTACGTTTTATTTCCA
#ACTTTGTCCCCGAGCCGAACGTGAATGGAGAGCTATAATCCTGCTGACAGAAATAAAC	HHGHHHHHHFGFHHGHHHHHHHHHHHHHHHHHHHHHHGGGHHHHHHHHHGHGHGGGHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHGHHHHHHHHHHHHGHHHGGGGGGGGGGGGHHHHHHHHHHHHHHHHHHHHHHGHHHHHHHHHHHHHHHGGGG	AS:
#i:32	NM:i:0	MD:Z:32	YF:H:2D	YI:i:1	YP:i:2	YS:i:0

$map_quality_cutoff=10;

open(IN,"$ARGV[0]")||die("Can not open in file\n");

my $out1=$ARGV[0].".keepDuplicate";
open(OUT1,">$out1")||die("Can not open out file\n");

my $out2=$ARGV[0].".keepDuplicate_noIgKj4";
open(OUT2,">$out2")||die("Can not open out2 file\n"); # the part of reads without JK4 mapping sequence

my $out3=$ARGV[0].".keepDuplicate_noIgKj4_all"; 
open(OUT3,">$out3")||die("Can not open out3 file\n"); # both part of split reads without JK4

my $out4=$ARGV[0].".keepDuplicate_multipleMapped"; 
open(OUT4,">$out4")||die("Can not open out4 file\n"); # part of reads are multiple pieces and are mapped

my $out5=$ARGV[0].".keepDuplicate_onlyIgKj4";
open(OUT5,">$out5")||die("Can not open out5 file\n"); # the part of reads with JK4 mapping sequence

my $step=int($ARGV[1]); # offset 2 nt

open(IN_germline,"$ARGV[2]")||die("Can not open in germline file\n");

if($ARGV[2]=~/myc/){  # mark the different bait chromosome
	$chr_set="chr15";
}
elsif($ARGV[2]=~/chr15/){
	$chr_set="chr15";
}
elsif($ARGV[2]=~/chr16/){
	$chr_set="chr16";
}
else{
	$chr_set="chr6";
}
	
	
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

while($line=<IN>){
	chomp $line;
	@split=split /\t/,$line;
	$head=$split[0];
	if($head=~/\_/){
		@head=split /\_/,$head;
		$head=$head[0];
	}
	if(exists($hash{$head})){
		$hash{$head}=$hash{$head}."\n".$line;
	}
	else{
		$hash{$head}=$line;
	}
	$hash_muntiple_map{$head}++;
#	print "$head\n$hash{$head}\n";
}

open(IN,"$ARGV[0]")||die("Can not open in file\n");
while(($key,$value)=each % hash){
		
	@temp=split /\n/,$hash{$key};
	@temp_1=split /\s/,$temp[0];
	@temp_2=split /\s/,$temp[1];
	$loc_1=int($temp_1[3]/$step);
	$loc_2=int($temp_2[3]/$step);
	$chr_1=$temp_1[2];
	$chr_2=$temp_2[2];
	$pos_1=$chr_1.'.'.$loc_1;
	$pos_2=$chr_2.'.'.$loc_2;#print "$pos_1\n$pos_2\n";
	$map_quality_1=$temp_1[4];
	$map_quality_2=$temp_2[4];

#	$mark=0; # to remove this case: M01626:178:000000000-ADJVC:1:1103:18262:5266	2064	chr5	77788118	5	90H30M114H   &&   M01626:178:000000000-ADJVC:1:1103:18262:5266	0	chr5	140995375	8	12S37M185S
	
	if($hash_muntiple_map{$key}<3){ # remove multiple mapped
		if((($temp_1[2] eq $chr_set)&&($temp_1[3]<$hash_mm9_end)&&($temp_1[3]>$hash_mm9_start))||(($temp_2[2] eq $chr_set)&&($temp_2[3]<$hash_mm9_end)&&($temp_2[3]>$hash_mm9_start))){
			if(($temp_1[2] eq $chr_set)&&($temp_1[3]<$hash_mm9_end)&&($temp_1[3]>$hash_mm9_start)){
			#	if($hash_pos{$pos_2}==1){ # overlapped position
			#	}
			#	else{
					if($map_quality_2>$map_quality_cutoff){
						print OUT1 "$hash{$key}\n"; $count++;
						print OUT2 "$temp[1]\n";
						print OUT5 "$temp[0]\n";
						$hash_pos{$pos_2}=1;
					}
					else{
						if(($map_quality_1>$map_quality_cutoff)&&($map_quality_2>$map_quality_cutoff)){
							print OUT4 "$hash{$key}\n";$count++;
						}
						else{
							print "000$hash{$key}\n";
						}
					}
			#	}
				
			}
			else{
			#	if($hash_pos{$pos_1}==1){ # overlapped position
			#	}
			#	else{
					if($map_quality_1>$map_quality_cutoff){
						print OUT1 "$hash{$key}\n";$count++;
						print OUT2 "$temp[0]\n";
						print OUT5 "$temp[1]\n";
						$hash_pos{$pos_1}=1;
					}
					else{
						if(($map_quality_1>$map_quality_cutoff)&&($map_quality_2>$map_quality_cutoff)){
							print OUT4 "$hash{$key}\n";$count++;
						}
						else{
							print "===$hash{$key}\n";
						}
					}
			#	}
			}
		}
		else{
			if(($map_quality_1>$map_quality_cutoff)&&($map_quality_2>$map_quality_cutoff)){
				print OUT3 "$hash{$key}\n";$count++;
			}
			else{
				print "----$hash{$key}\n";
			}
		}			
	}
	else{
		if(($map_quality_1>$map_quality_cutoff)&&($map_quality_2>$map_quality_cutoff)){
			print OUT4 "$hash{$key}\n";$count++;
		}
		else{
			print "88888$hash{$key}\n";
		}
	}
}
print "$count\n";
			
