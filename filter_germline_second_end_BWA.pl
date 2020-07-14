if(@ARGV<1){
	print "This script is to filter out second end of sequence which is germline
		Usage: perl *pl mm9_germine.region r2.BWA.sam r2.Yaha.sam \n
	Output:  	 r2*.germline_filtered_SplitRead
		 	 r2*.germline_filtered_noSplitRead			#	>75%mapped
		 	 r2*.germline
			 r2*.germline_filtered_pure					#	100% mapped (no split)
			 r2*.germline_filtered_long_unmap_part		#	<75%mapped
			 \n";
}

open(IN1,"$ARGV[0]")||die("Can not open in1 file\n");

while($line=<IN1>){
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

$map_quality_cutoff=10;

my $out1=$ARGV[1].".germline_filtered_SplitRead";
open(OUT1,">$out1")||die("Can not open out1 file\n");

my $out2=$ARGV[1].".germline_filtered_noSplitRead";		# only one mapped part with unmapped, eg. chr6	70121802	60	23S205M
open(OUT2,">$out2")||die("Can not open out2 file\n");

my $out3=$ARGV[1].".germline";
open(OUT3,">$out3")||die("Can not open out3 file\n");

my $out4=$ARGV[1].".germline_filtered_pure";			# no unmapped part, eg. chr6	69725794	11	229M
open(OUT4,">$out4")||die("Can not open out4 file\n");

my $out5=$ARGV[1].".germline_filtered_long_unmap_part";	# long unmapped part, eg. chr6	70093587	12	176S36M
open(OUT5,">$out5")||die("Can not open out5 file\n");

open(IN2,"$ARGV[1]")||die("Can not open in2 file\n");
if($ARGV[0]=~/mm9/){
	while($line=<IN2>){
		chomp $line;
		if($line=~/^\@/){
		}
		else{
			@split=split /\t/,$line;
			$chr=$split[2];
			$pos=$split[3];
			$head=$split[0];
			$cigar=$split[5];
			if(exists($hash{$head})){#print "$hash{$head}\n$line\n";
				$hash{$head}=2;#print "$head\n$hash{$head}\n";
			}
			else{
				$hash{$head}=1;#print "$hash{$head}\n";
			}
			$num_S=$cigar=~s/S/S/g;
			$num_M=$cigar=~s/M/M/g;
			if($num_S>1){  # eg. M01626:178:000000000-ADJVC:1:1102:7081:6135_2:N:0:1     16      chr14   83543638        5       63S26M145S
				$cigar=~s/M/S/g;
				$cigar=~s/S/\t/g;
				@cigar=split /\t/,$cigar;
				if(($cigar[1]>$cigar[0]*3)||($cigar[1]>$cigar[2]*3)){ # M01626:162:000000000-AD53D:1:2106:12797:11363	0	chr6	70673545	60	63S138M5S
				}
				else{
					$hash_mark{$head}="WrongMapping"; # mark the read which is wrong mapped eg. M01626:178:000000000-ADJVC:1:1102:7081:6135_2:N:0:1     16      chr14   83543638        5       63S26M145S
			#		print "$line\n";
				}
			}
		}
	}
}
			
open(IN2,"$ARGV[1]")||die("Can not open in2 file\n");
if($ARGV[0]=~/mm9/){
	while($line=<IN2>){
		chomp $line;
		$num_S=0;
		$num_M=0;
		if($line=~/^\@/){
		}
		else{
			@split=split /\t/,$line;
			$chr=$split[2];
			$pos=$split[3];
			$head=$split[0];
			$cigar=$split[5];
			$length=length($split[9]); #print "$chr\t$pos\n";
			$mappingQuality=$split[4];
#			if(($hash{$head}>1)&&($hash_mark{$head} ne "WrongMapping")&&($mappingQuality>$map_quality_cutoff)){  # mapping quality >20
			if(($hash{$head}>1)&&($hash_mark{$head} ne "WrongMapping")){	
				$mark{$head}="double"; # this read is split mapping read
		#		$num_S=$cigar=~s/S/S/g;
		#		$num_M=$cigar=~s/M/M/g;
		#		if($num_S>1){  # eg. M01626:178:000000000-ADJVC:1:1102:7081:6135_2:N:0:1     16      chr14   83543638        5       63S26M145S
		#			$cigar=~s/M/S/g;
		#			$cigar=~s/S/\t/g;
		#			if(($cigar[1]>$cigar[0]*3)||($cigar[1]>$cigar[2]*3)){ # M01626:162:000000000-AD53D:1:2106:12797:11363	0	chr6	70673545	60	63S138M5S
		#			print "$hash_Yaha_line{$head}\n$cigar[0]\t$cigar[1]\t$cigar[2]\n";				
		#				@temp=split /\n/,$hash{$head};
		#				@split_temp_1=split /\t/,$temp[0];
		#				@split_temp_2=split /\t/,$temp[1];
		#				if((($split_temp_1[2] eq "chr6")&&($split_temp_1[3]<$hash_mm9_end)&&($split_temp_1[3]>$hash_mm9_start))||(($split_temp_2[2] eq "chr6")&&($split_temp_2[3]<$hash_mm9_end)&&($split_temp_2[3]>$hash_mm9_start))){ # make sure one end is IgKj4
		#					print "$line\n";
		#					print OUT1 "$line\n";
		#				}
		#			}
		#			else{
		#				$hash_mark{$head}="WrongMapping"; # mark the read which is wrong mapped eg. M01626:178:000000000-ADJVC:1:1102:7081:6135_2:N:0:1     16      chr14   83543638        5       63S26M145S
		#			}
		#		}
		#		else{
					print OUT1 "$line\n";
		#		}
			}
			else{
				if($chr eq "chr6"){
					if(($pos<$hash_mm9_end)&&($pos>$hash_mm9_start)){ # this read is germline		
						print OUT3 "$line\n";
					}
					else{
						$num_S=$cigar=~s/S/S/g;
						$num_M=$cigar=~s/M/M/g;
						if(($num_S<2)&&($num_M<2)){
								if(($num_S==0)&&($mappingQuality>$map_quality_cutoff)){
									print OUT4 "$line\n";
								}
								else{								
									$loc_S=index($cigar,'S');	$loc_M=index($cigar,'M');
									$cigar=~s/S/\t/;
									$cigar=~s/M/\t/;
									@cigar=split /\t/,$cigar;
									if($loc_S<$loc_M){
										$cutoff=int($cigar[1]/4);
										if($cigar[0]<$cutoff){
											if($mappingQuality>$map_quality_cutoff){
												print OUT2 "$line\n";
											}
										}
										else{
											if($mappingQuality>$map_quality_cutoff){
												print OUT5 "$line\n";
											}
										}
									}
									else{
										$cutoff=int($cigar[0]/4);
										if($cigar[1]<$cutoff){
											if($mappingQuality>$map_quality_cutoff){
												print OUT2 "$line\n";
											}
										}
										else{
											if($mappingQuality>$map_quality_cutoff){
												print OUT5 "$line\n";
											}
										}
									}
								}
						}
						$mark{$head}="print";
					}
				}
				else{
					if($chr eq '*'){
					}
					else{
						$num_S=$cigar=~s/S/S/g;
						$num_M=$cigar=~s/M/M/g;
						$mark{$head}="print";
						if(($num_S<2)&&($num_M<2)){
							if(($num_S==0)&&($mappingQuality>$map_quality_cutoff)){
								print OUT4 "$line\n";
							}
							else{								
								$loc_S=index($cigar,'S');	$loc_M=index($cigar,'M');
								$cigar=~s/S/\t/;
								$cigar=~s/M/\t/;
								@cigar=split /\t/,$cigar;
								if($loc_S<$loc_M){
									$cutoff=int($cigar[1]/4);
									if($cigar[0]<$cutoff){
										if($mappingQuality>$map_quality_cutoff){
											print OUT2 "$line\n";
										}
									}
									else{
										if($mappingQuality>$map_quality_cutoff){
											print OUT5 "$line\n";
										}
									}
								}
								else{
									$cutoff=int($cigar[0]/4);
									if($cigar[1]<$cutoff){
										if($mappingQuality>$map_quality_cutoff){
											print OUT2 "$line\n";
										}
									}
									else{
										if($mappingQuality>$map_quality_cutoff){
											print OUT5 "$line\n";
										}
									}									
								}
							}
						}
					}
				}
			}
		}
	}
}

########################################################################################################
# SCAN YAHA mapping result for the same sequence file, output the BWA missed split mapping read
########################################################################################################

open(IN3,"$ARGV[2]")||die("Can not open in3 file\n");
if($ARGV[0]=~/mm9/){
	while($line=<IN3>){
		chomp $line;
		if($line=~/^\@/){
		}
		else{
			@split=split /\t/,$line;
			$chr=$split[2];
			$pos=$split[3];
			@temp=split /\_/,$split[0];#print "$temp\n";
			$head=$temp[0];
			if(exists($hash_Yaha{$head})){#print "$hash{$head}\n$line\n";
				$hash_Yaha{$head}=2;
				$hash_Yaha_line{$head}=$hash_Yaha_line{$head}."\n".$line;
			}
			else{
				$hash_Yaha{$head}=1;#print "$hash{$head}\n";
				$hash_Yaha_line{$head}=$line;
			}
		}
	}
}					

open(IN3,"$ARGV[2]")||die("Can not open in3 file\n");
if($ARGV[0]=~/mm9/){
	while($line=<IN3>){
		$num_S=0;
		$num_M=0;
		chomp $line;
		if($line=~/^\@/){
		}
		else{
			@split=split /\t/,$line;
			$chr=$split[2];
			$pos=$split[3];
			$cigar=$split[5]; # print "$cigar\n";
			@temp=split /\_/,$split[0];
			$head=$temp[0];
			$mapping_quality=$split[4];
			$length=length($split[9]); #print "$chr\t$pos\n";
			if($hash_Yaha{$head}>1){#
				if($mark{$head} eq "double"){ # this read is split mapping read and already output from BWA mapping file
				}
				else{ # this read is split mapping read but miss in BWA mapping
				#	print "$line\n";
					if(($chr eq "chr6")&&($pos<$hash_mm9_end)&&($pos>$hash_mm9_start)){ # this read is germline
					}
					else{# print "$line\n";
						 $num_S=$cigar=~s/S/S/g;
                         $num_M=$cigar=~s/M/M/g;
						 if($num_S>1){  # eg. M01626:178:000000000-ADJVC:1:1102:7081:6135_2:N:0:1     16      chr14   83543638        5       63S26M145S
							$cigar=~s/M/S/g;
							$cigar=~s/S/\t/g;
							@cigar=split /\t/,$cigar;
							if(($cigar[1]>$cigar[0]*3)||($cigar[1]>$cigar[2]*3)){
							#	print "$hash_Yaha_line{$head}\n$cigar[0]\t$cigar[1]\t$cigar[2]\n";
								@temp=split /\n/,$hash_Yaha_line{$head};
								@split_temp_1=split /\t/,$temp[0];
								@split_temp_2=split /\t/,$temp[1];
								if((($split_temp_1[2] eq "chr6")&&($split_temp_1[3]<$hash_mm9_end)&&($split_temp_1[3]>$hash_mm9_start))||(($split_temp_2[2] eq "chr6")&&($split_temp_2[3]<$hash_mm9_end)&&($split_temp_2[3]>$hash_mm9_start))){ # make sure one end is IgKj4
									if($mappingQuality>$map_quality_cutoff){
			#							print OUT1 "$hash_Yaha_line{$head}\n";
									}
								}
							}
						}
						else{
							if($mappingQuality>$map_quality_cutoff){
			#					print OUT1 "$hash_Yaha_line{$head}\n";
							}
						}	
					}
				}
			}
			else{
				if($chr eq "chr6"){
					if(($pos<$hash_mm9_end)&&($pos>$hash_mm9_start)){ # this read is germline		
							
					}
					else{
						if($mark{$head} eq "print"){
						}
						else{
							$num_S=$cigar=~s/S/S/g;
							$num_M=$cigar=~s/M/M/g;
							if(($num_S<2)&&($num_M<2)){
								if($num_S==0){
								}
								else{								
									$loc_S=index($cigar,'S');	$loc_M=index($cigar,'M');
									$cigar=~s/S/\t/;
									$cigar=~s/M/\t/;
									@cigar=split /\t/,$cigar;
									if($loc_S<$loc_M){
										$cutoff=int($cigar[1]/4);
										if($cigar[0]<$cutoff){
											print OUT2 "$line\n";
										}
									}
									else{
										$cutoff=int($cigar[0]/4);
										if($cigar[1]<$cutoff){
											print OUT2 "$line\n";
										}
									}
								}
							}
						}
					}
				}
				else{
					if($chr eq '*'){
					}
					else{
						$num_S=$cigar=~s/S/S/g;
						$num_M=$cigar=~s/M/M/g;#print "$num_S\t$num_M\n";
					
						if($mark{$head} eq "print"){
						}
						else{
							if(($num_S<2)&&($num_M<2)){
								if($num_S==0){
								}
								else{								
									$loc_S=index($cigar,'S');	$loc_M=index($cigar,'M');
									$cigar=~s/S/\t/;
									$cigar=~s/M/\t/;
									@cigar=split /\t/,$cigar;
									if($loc_S<$loc_M){
										$cutoff=int($cigar[1]/4);
										if($cigar[0]<$cutoff){
											print OUT2 "$line\n";
										}
									}
									else{
										$cutoff=int($cigar[0]/4);
										if($cigar[1]<$cutoff){
											print OUT2 "$line\n";
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
	

