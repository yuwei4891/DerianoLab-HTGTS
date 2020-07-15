if(@ARGV<1){
	print "	
		This script is to disect the resection length around V-segment RSS site specially for cell cycle release
		
		& add microhomology length for each with or without resection 

	Usage: perl *pl
					r2.*.keepDuplicate
					Igk.segments.coordinate.mm9.bed
	Output:
					*.V_segment_resection_micro_length			### resection length and microhomology length
					*.V_segment_CE_resection_length				
					*.V_segment_SE_resection_length				
					*.microhomology								
					*.insert									
					*.microhomology.length						
					*.microhomology.length.percentage			
					*.insertion.length							
					*.micro_insert_resect_summary				
		\n";
}
					
open(IN,"$ARGV[0]")||die("Can not open in file\n");
my $out=$ARGV[0].".V_segment_resection_micro_length";
open(OUT,">$out")||die("Can not open out file\n");
my $out_CE=$ARGV[0].".V_segment_CE_resection_length";
open(OUT_CE,">$out_CE")||die("Can not open out_CE file\n");
my $out_SE=$ARGV[0].".V_segment_SE_resection_length";
open(OUT_SE,">$out_SE")||die("Can not open out_SE file\n");

open(IN_read_sam,"$ARGV[0]")||die("Can not open in sam file\n");
my $out_micro=$ARGV[0].".microhomology";
open(OUT_micro,">$out_micro")||die("Can not open micro file\n");

my $out_insert=$ARGV[0].".insert";
open(OUT_insert,">$out_insert")||die("Can not open insert file\n");

my $out_micro_length=$ARGV[0].".microhomology.length";
open(OUT_micro_length,">$out_micro_length")||die("Can not open micro_length file\n");

my $out_micro_length_percentage=$ARGV[0].".microhomology.length.percentage";
open(OUT_micro_length_percentage,">$out_micro_length_percentage")||die("Can not open micro_length_percentage file\n");

my $out_insert_length=$ARGV[0].".insertion.length";
open(OUT_insert_length,">$out_insert_length")||die("Can not open insert_length file\n");

my $summary=$ARGV[0].".micro_insert_resect_summary";
open(OUT_summary,">$summary")||die("Can not open summary file\n");


$distance=4000;		
$resection_error=500;

if($ARGV[0]=~/JK/){
	$bait_site=70673;
	$bait_chr="chr6";
	$library="JK";
}
elsif($ARGV[0]=~/myc/){
	$bait_site=61818;
	$bait_chr="chr15";
	$library="myc";
}
else{
	print "ERROR\n";
}

while($line=<IN>){
	chomp $line;
	@split=split /\t/,$line;
	$chr=$split[2];
	$site=$split[3];
	$cigar=$split[5];
	$ID=$split[0];
	$hash_count{$ID}++;
	if($hash_count{$ID}==1){
		$hash_line{$ID}=$line;
		$hash_chr{$ID}=$chr;
		$hash_site{$ID}=$site;
		$hash_cigar{$ID}=$cigar;
	}
	else{
		$hash_line{$ID}=$hash_line{$ID}."\n".$line;
		$hash_chr{$ID}=$hash_chr{$ID}."\n".$chr;
		$hash_site{$ID}=$hash_site{$ID}."\n".$site;
		$hash_cigar{$ID}=$hash_cigar{$ID}."\n".$cigar;
		$hashMark{$ID}++;
	}
}

while(($key,$value)= each % hash_line){	# this paragraph is to find the breakpoint 
	$total_read++;
	@line=split /\n/,$value;
	@chr=split /\n/,$hash_chr{$key};
	@site=split /\n/,$hash_site{$key};
	@cigar=split /\n/,$hash_cigar{$key};
	if($hash_count{$key}>2){#print "$value\n";
	}
	else{#print "$value\n";
		if($library eq "myc"){
			for($i=0;$i<2;$i++){#print "$line[$i]\n";
				if(($chr[$i] eq $bait_chr)&&($site[$i]=~/$bait_site/)){
					$bait=$i;
				}
			}
			for($i=0;$i<2;$i++){
				if($i==$bait){
				}
				else{
					$target=$i;
					$number_S=$cigar[$target]=~s/S/S/g;	# target sequence have 1 or 2 S
					$loc_S=index($cigar[$target],'S');	# find the pos of S, then compare the order of S and M
					$loc_M=index($cigar[$target],'M');	# find the pos of M, then compare the order of S and M
					if($number_S>1){
						$temp=$cigar[$target];		$temp=~s/M/S/g;
						@temp=split /S/,$temp;
						if($temp[0]<$temp[2]){		
						# indicate the break point is in second junction
						#	M05218:3:000000000-AVUAY:1:1109:19708:3896	0	chr6	70672442	60	2S57M60S
						#	M05218:3:000000000-AVUAY:1:1109:19708:3896	16	chr15	61818160	60	5S55M59S
							$breakPoint=$site[$target]+$temp[1];
						}
						else{
						# indicate the break point is in first junction
						#	M01626:265:000000000-AT5KV:1:1104:28017:20786	16	chr6	135933853	60	60S114M24S
						#	M01626:265:000000000-AT5KV:1:1104:28017:20786	16	chr15	61818195	60	63M135S
							$breakPoint=$site[$target];
						}
						$pos_breakPoint=$chr[$target].'.'.$breakPoint;
						$hash_pos_breakPoint{$pos_breakPoint}++;	# mark this site as one read resect to here !!!!!!
						$hash_pos_breakPoint_line{$pos_breakPoint}=$hash_pos_breakPoint_line{$pos_breakPoint}."\n".$value;
						$total++;
					#	print "$value\n@temp\n$pos_breakPoint\n\n";
					}
					else{
						$temp=$cigar[$target];		$temp=~s/M/S/g;
						@temp=split /S/,$temp;
						if($loc_S<$loc_M){
						#	M05218:3:000000000-AVUAY:1:2103:4210:11635	16	chr6	68686390	60	49S133M
						#	M05218:3:000000000-AVUAY:1:2103:4210:11635	16	chr15	61818160	58	5S43M134S
							$breakPoint=$site[$target];
						}
						else{
						#	M01626:265:000000000-AT5KV:1:1105:6646:21993	0	chr9	50173369	60	57M37S
						#	M01626:265:000000000-AT5KV:1:1105:6646:21993	0	chr15	61818203	53	56S38M
							$breakPoint=$site[$target]+$temp[0];
						}
						$pos_breakPoint=$chr[$target].'.'.$breakPoint;
						$hash_pos_breakPoint{$pos_breakPoint}++;
						$hash_pos_breakPoint_line{$pos_breakPoint}=$hash_pos_breakPoint_line{$pos_breakPoint}."\n".$value;
						$total++;
				#		print "$value\n$pos_breakPoint\n\n";
					}
				}
			}
		}
		else{
			for($i=0;$i<2;$i++){#print "$line[$i]\n";
				if(($chr[$i] eq $bait_chr)&&($site[$i]=~/$bait_site/)){
					#M05218:3:000000000-AVUAY:1:2103:4210:11635	16	chr15	61818160	58	5S43M134S
					#M05218:3:000000000-AVUAY:1:2103:4210:11635	16	chr6	68686390	60	49S133M
					$bait=$i;
				}
			}
			for($i=0;$i<2;$i++){
				if($i==$bait){
				}
				else{
					$target=$i;
					$number_S=$cigar[$target]=~s/S/S/g;	# target sequence have 1 or 2 S
					$number_M=$cigar[$target]=~s/S/S/g;	# target sequence have 1 or 2 S
					$loc_S=index($cigar[$target],'S');	# find the pos of S, then compare the order of S and M
					$loc_M=index($cigar[$target],'M');	# find the pos of M, then compare the order of S and M
					if($number_S>1){						
					
						$temp=$cigar[$target];		$temp=~s/M/S/g;
						@temp=split /S/,$temp;
					#	print "$value\n@temp\n";
					#	print "$chr[$i]\t$site[$i]\t$cigar[$i]\n$chr[$bait]\t$site[$bait]\t$cigar[$bait]\n\n";
						if($temp[0]<$temp[2]){		
						# indicate the break point is in second junction
						#	M05218:3:000000000-AVUAY:1:1109:19708:3896	0	chr6	70672442	60	2S57M60S
						#	M05218:3:000000000-AVUAY:1:1109:19708:3896	16	chr15	61818160	60	5S55M59S
							$breakPoint=$site[$target]+$temp[1];
						}
						else{
						# indicate the break point is in first junction
						#	M01626:265:000000000-AT5KV:1:1104:28017:20786	16	chr6	135933853	60	60S114M24S
						#	M01626:265:000000000-AT5KV:1:1104:28017:20786	16	chr15	61818195	60	63M135S
							$breakPoint=$site[$target];
						}
						$pos_breakPoint=$chr[$target].'.'.$breakPoint;
						$hash_pos_breakPoint{$pos_breakPoint}++;	# mark this site as one read resect to here !!!!!!
						$hash_pos_breakPoint_line{$pos_breakPoint}=$hash_pos_breakPoint_line{$pos_breakPoint}."\n".$value;
						if(exists($hash_pos_breakPoint_ID{$pos_breakPoint})){
							$hash_pos_breakPoint_ID{$pos_breakPoint}=$hash_pos_breakPoint_ID{$pos_breakPoint}."\t".$key;		# store the ID 
						}else{
							$hash_pos_breakPoint_ID{$pos_breakPoint}=$key;
						}
						$total++;
					#	print "$value\n@temp\n$pos_breakPoint\n\n";
					
					}
					else{
						$temp=$cigar[$target];		$temp=~s/M/S/g;
						@temp=split /S/,$temp;
						if($loc_S<$loc_M){
						#	M01626:243:000000000-AP8H9:1:1105:2397:14800	16	chr6	70024030	60	104S115M
						#	M01626:243:000000000-AP8H9:1:1105:2397:14800	0	chr6	70673563	60	115S104M
							$breakPoint=$site[$target];
						
						}
						else{
						#	M01626:243:000000000-AP8H9:1:1112:27986:12761	0	chr6	64437201	60	130M89S
						#	M01626:243:000000000-AP8H9:1:1112:27986:12761	0	chr6	70673525	60	127S92M
							$breakPoint=$site[$target]+$temp[0];
						#	print "$value\n$target\t$site[$target]\t$cigar[$target]\n$breakPoint\n\n";
						}
						$pos_breakPoint=$chr[$target].'.'.$breakPoint;
						$hash_pos_breakPoint{$pos_breakPoint}++;
						$hash_pos_breakPoint_line{$pos_breakPoint}=$hash_pos_breakPoint_line{$pos_breakPoint}."\n".$value;
						if(exists($hash_pos_breakPoint_ID{$pos_breakPoint})){
							$hash_pos_breakPoint_ID{$pos_breakPoint}=$hash_pos_breakPoint_ID{$pos_breakPoint}."\t".$key;		# store the ID 
						}else{
							$hash_pos_breakPoint_ID{$pos_breakPoint}=$key;
						}
						$total++;
				#		print "$value\n$pos_breakPoint\n\n";
					}
				}
			}
		}		
	}
}


#########################################################################
# analysis both cigar to find microhomology and insertion
#########################################################################
while(($key,$value)=each % hash_line){
#	print "$key\n$value\n";
	if(($hashMark{$key}>0)&&($hashMark{$key}<2)){
		@split=split /\n/,$hash_line{$key};
		@read1=split /\t/,$split[0];		@read2=split /\t/,$split[1];
		$map_position_1=$read1[3];		$map_position_2=$read2[3];
		$cigar_1=$read1[5];	$cigar_2=$read2[5];		
		$ID1=$read1[0];	$ID2=$read2[0];					#print "$ID1\n$ID2\n";
		if($ID1 ne $ID2){
			print "ERROR\n$hash{$key}\n$ID1\n$ID2\n";
		}
		$num_S_1=$cigar_1=~s/S/S/g;		$num_S_2=$cigar_2=~s/S/S/g;
		$num_M_1=$cigar_1=~s/M/M/g;		$num_M_2=$cigar_2=~s/M/M/g;
		
		#########################################################################
		# rebuild cigar if cigar have two S
		##########################################################################
		if($num_S_1>1){			# rebuild cigar  eg. 97S352M8S  to 97S360M  ;     83S254M146S  to 337M146S
			$temp=$cigar_1;
			if($num_M_1>1){		# rebuild cigar , eg : 7S295M1D98M54S  to 7S295M1D98M54S  ;  115S83M1I46M17S  to 
				$temp=~s/M/S/g;
				$temp=~s/D/S/g;
				$temp=~s/I/S/g;
				@temp=split /S/,$temp;
				if($temp[0]<$temp[4]){		# eg : 7S295M1D98M54S
					$map=$temp[0]+$temp[1]+$temp[3];
					$unmap=$temp[4];
				}
				else{
					$map=$temp[1]+$temp[3]+$temp[4];
					$unmap=$temp[0];
				}
			}
			else{	# rebuild cigar  eg. 97S352M8S  to 97S360M  ;     83S254M146S  to 337M146S
			#	print "$value\n";
				$temp=~s/M/S/g;
				@temp=split /S/,$temp;
				if($temp[0]<$temp[2]){		# eg. 83S254M146S  
					$map=$temp[0]+$temp[1];
					$unmap=$temp[2];
				}
				else{
					$map=$temp[1]+$temp[2];
					$unmap=$temp[0];
				}
			}
			$cigar_1=$map.'M'.$unmap.'S';
		#	print "$value\n$cigar_1\n";
		}
		
		if($num_S_2>1){			# rebuild cigar  eg. 97S352M8S  to 97S360M  ;     83S254M146S  to 337M146S
			$temp=$cigar_2;
			if($num_M_2>1){		# rebuild cigar , eg : 7S295M1D98M54S  to 7S295M1D98M54S  ;  115S83M1I46M17S  to 
				$temp=~s/M/S/g;
				$temp=~s/D/S/g;
				$temp=~s/I/S/g;
				@temp=split /S/,$temp;
				if($temp[0]<$temp[4]){		# eg : 7S295M1D98M54S
					$map=$temp[0]+$temp[1]+$temp[3];
					$unmap=$temp[4];
				}
				else{
					$map=$temp[1]+$temp[3]+$temp[4];
					$unmap=$temp[0];
				}
			}
			else{	# rebuild cigar  eg. 97S352M8S  to 97S360M  ;     83S254M146S  to 337M146S
			#	print "$value\n";
				$temp=~s/M/S/g;
				@temp=split /S/,$temp;
				if($temp[0]<$temp[2]){		# eg. 83S254M146S  
					$map=$temp[0]+$temp[1];
					$unmap=$temp[2];
				}
				else{
					$map=$temp[1]+$temp[2];
					$unmap=$temp[0];
				}
			}
			$cigar_2=$unmap.'S'.$map.'M';
			
		#	print "$value\n$cigar_2\n";
		}		
		
		
		##########################################################################
		#	get the map & unmap part of each cigar 1 & 2
		##########################################################################		
		$loc_S_1=index($cigar_1,'S');	$loc_M_1=index($cigar_1,'M');
		$loc_S_2=index($cigar_2,'S');	$loc_M_2=index($cigar_2,'M');

		if($loc_S_1<$loc_M_1){		# eg. 120S312M
			$temp=$cigar_1;
			$temp=~s/M/S/g;
			@temp=split /S/,$temp;
			$map_1=$temp[1];
			$unmap_1=$temp[0];
		#	print "$cigar_1\n$map_1\n$unmap_1\n\n";
		}
		else{					# 
			$temp=$cigar_1;
			$temp=~s/M/S/g;
			@temp=split /S/,$temp;
			$map_1=$temp[0];
			$unmap_1=$temp[1];
		#	print "$cigar_1\n$map_1\n$unmap_1\n\n";
		}
		
		if($loc_S_2<$loc_M_2){		# eg. 120S312M
			$temp=$cigar_2;
			$temp=~s/M/S/g;
			@temp=split /S/,$temp;
			$map_2=$temp[1];
			$unmap_2=$temp[0];
		#	print "$cigar_1\n$map_1\n$unmap_1\n\n";
		}
		else{					# 
			$temp=$cigar_2;
			$temp=~s/M/S/g;
			@temp=split /S/,$temp;
			$map_2=$temp[0];
			$unmap_2=$temp[1];
		#	print "$cigar_2\n$map_2\n$unmap_2\n\n";
		}
		
		##########################################################################
		#	dissect the micro  & insertion length for each read base on ID
		##########################################################################	
		if($map_1<$unmap_2){					# eg. 120S312M  ; 315S117M   insertion  !!!
			$insert_length=$unmap_2-$map_1;
			if($insert_length<50){				# make sure it is real insertion, not immappropriate mapping
				$total_insert++;
				print OUT_insert "$value\n";
				print OUT_insert_length "$insert_length\n";
			}
			$hash_micro{0}++;
			$hash_micro_ID{$ID1}=0;				# no microhomology
		}
		elsif($map_1>$unmap_2){					# eg. 213M181S ;  209S185M  microhomology !!!
			$micro_length=$map_1-$unmap_2;
			if($micro_length>10){
				$hash_micro{10}++;
			}
			else{
				$hash_micro{$micro_length}++;
			}
			$total_micro++;
			$hash_micro_ID{$ID1}=$micro_length;		# this is microhomology 
			print OUT_micro "$value\n";
			print OUT_micro_length "$micro_length\n";
		}
		else{
			$hash_micro{0}++;
			$hash_micro_ID{$ID1}=0;				# exact join  no microhomology
		}
	}
}
	
open(IN_V_segment,"$ARGV[1]")||die("Can not open V segment file\n");
#chr6	67892068	67892354	-	V	Igkv9-124	MGI:3646892	67892067
#chr6	67797403	67797868	+	V	Igkv9-128	MGI:3809196	67797869
#chr6	70672556	70672593	+	J	Igkj1	MGI:1316689	70672555

####################################
#	resect toward coding end, length>0; 
#	resect toward signal end, length<0
####################################


######## find the nearest resection length of eah position
while($line=<IN_V_segment>){
	chomp $line;
	@split=split /\t/,$line;
	$segment_name=$split[5];
#	if(($split[4] eq 'V')&&(!($segment_name eq "Igkv11-114"))&&(!($segment_name eq "Igkv15-101"))&&(!($segment_name eq "Igkv11-114"))){
	if($split[4] eq 'V'){
		$V_segment=$split[5];
		
		$strand=$split[3];
		$start=$split[1];
		$end=$split[2];
		$chr=$split[0];
		if($strand eq '+'){
			$V_cod_start=$end-$distance;
			$V_cod_end=$end;
			$V_sig_start=$end+1;
			$V_sig_end=$end+1+$distance;
			for($i=$V_cod_start+1;$i<$V_cod_end+1;$i++){
				$resection_length=abs($V_cod_end-$i);		# break point in coding end
				$pos=$chr.'.'.$i;
				if($hash_pos_breakPoint{$pos}>0){
					for($j=0;$j<$hash_pos_breakPoint{$pos};$j++){
						if($hash_check{$pos}==1){		# make sure this position have not be counted yet
							if($resection_length>$hash_resection_length{$pos}){
							}
							else{
								$hash_resection_length{$pos}=$resection_length;
							}
						}
						else{
							$hash_check{$pos}=1;
							$hash_resection_length{$pos}=$resection_length;
						}
					}
				}
			}
			for($i=$V_sig_start+1;$i<$V_sig_end+1;$i++){
				$resection_length=abs($V_sig_start-$i);		# break point in signal end 
				$pos=$chr.'.'.$i;
				if($hash_pos_breakPoint{$pos}>0){
					for($j=0;$j<$hash_pos_breakPoint{$pos};$j++){
						if($hash_check{$pos}==1){		# make sure this position have not be counted yet
							if($resection_length>$hash_resection_length{$pos}){
							}
							else{
								$hash_resection_length{$pos}=$resection_length;
							}
						}
						else{
							$hash_check{$pos}=1;
							$hash_resection_length{$pos}=$resection_length;
						}
					}
				}
			}
		}
		else{
			$V_cod_start=$start;
			$V_cod_end=$start+$distance;
			$V_sig_start=$start-1-$distance;
			$V_sig_end=$start-1;
			for($i=$V_cod_start+1;$i<$V_cod_end+1;$i++){
				$resection_length=abs($i-$V_cod_start);		# break point in coding end is >0 
				$pos=$chr.'.'.$i;
				if($hash_pos_breakPoint{$pos}>0){
					for($j=0;$j<$hash_pos_breakPoint{$pos};$j++){
						if($hash_check{$pos}==1){
							if($resection_length>$hash_resection_length{$pos}){
							}
							else{
								$hash_resection_length{$pos}=$resection_length;
							}
						}
						else{
							$hash_check{$pos}=1;
							$hash_resection_length{$pos}=$resection_length;
						}
					}
				}
			}
			for($i=$V_sig_start+1;$i<$V_sig_end+1;$i++){
				$resection_length=abs($i-$V_sig_end);		# break point in signal end is <0 
				$pos=$chr.'.'.$i;
				if($hash_pos_breakPoint{$pos}>0){
					for($j=0;$j<$hash_pos_breakPoint{$pos};$j++){
						if($hash_check{$pos}==1){		# make sure this position have not be counted yet
							if($resection_length>$hash_resection_length{$pos}){
							}
							else{
								$hash_resection_length{$pos}=$resection_length;
							}
						}
						else{
							$hash_check{$pos}=1;
							$hash_resection_length{$pos}=$resection_length;
						}
					}
				}
			}
		}		
	}
}

######## find the nearest resection length of eah position

open(IN_V_segment,"$ARGV[1]")||die("Can not open V segment file\n");
while($line=<IN_V_segment>){
	chomp $line;
	@split=split /\t/,$line;
	$segment_name=$split[5];
#	if(($split[4] eq 'V')&&(!($segment_name eq "Igkv11-114"))&&(!($segment_name eq "Igkv15-101"))&&(!($segment_name eq "Igkv11-114"))){
	if($split[4] eq 'V'){
		$V_segment=$split[5];
		
		$strand=$split[3];
		$start=$split[1];
		$end=$split[2];
		$chr=$split[0];
		if($strand eq '+'){
			$V_cod_start=$end-$distance;
			$V_cod_end=$end;
			$V_sig_start=$end+1;
			$V_sig_end=$end+1+$distance;
			for($i=$V_cod_start+1;$i<$V_cod_end+1;$i++){
				$resection_length=$V_cod_end-$i;		# break point in coding end is >0 
				$pos=$chr.'.'.$i;
				if($hash_check_2{$pos}==1){		# make sure this position have not be counted yet
				}
				else{
					if($hash_pos_breakPoint{$pos}>0){
						@temp=split /\t/,$hash_pos_breakPoint_ID{$pos};		# find the read ID of each position
						#print "@temp\n";
						for($j=0;$j<$hash_pos_breakPoint{$pos};$j++){
					#		print "$hash_resection_length{$pos}\n";
							print OUT "$hash_resection_length{$pos}\t$hash_micro_ID{$temp[$j]}\n";
					#		print OUT "$hash_resection_length{$pos}\t$hash_micro_ID{$temp[$j]}\t$line\n";
							if($hash_resection_length{$pos}>$resection_error){
								$mark_V_segment{$segment_name}++;		# mark this V-segment have wired surper long resection
							}
							print OUT_CE "$hash_resection_length{$pos}\n";
							$print_count++;
						}
					}
				}
				$hash_check_2{$pos}=1;
			}
			for($i=$V_sig_start+1;$i<$V_sig_end+1;$i++){
				$resection_length=$V_sig_start-$i;		# break point in signal end is <0 
				$pos=$chr.'.'.$i;
				if($hash_check_2{$pos}==1){		# make sure this position have not be counted yet
				}
				else{
					if($hash_pos_breakPoint{$pos}>0){
						@temp=split /\t/,$hash_pos_breakPoint_ID{$pos};		# find the read ID of each position
						for($j=0;$j<$hash_pos_breakPoint{$pos};$j++){
							print OUT "-$hash_resection_length{$pos}\t$hash_micro_ID{$temp[$j]}\n";
					#		print OUT "-$hash_resection_length{$pos}\t$hash_micro_ID{$temp[$j]}\t$line\n";
							if($hash_resection_length{$pos}>$resection_error){
								$mark_V_segment{$segment_name}++;		# mark this V-segment have wired surper long resection
							}
							$temp_rec=0-int($resection_length);
						#	print OUT_SE "$temp_rec\n";
							print OUT_SE "$resection_length\n";
						#	print OUT "$resection_length\n$hash_pos_breakPoint_line{$pos}\n\n";
							$print_count++;
						}
					}
				}
				$hash_check_2{$pos}=1;
			}
		}
		else{
			$V_cod_start=$start;
			$V_cod_end=$start+$distance;
			$V_sig_start=$start-1-$distance;
			$V_sig_end=$start-1;
			for($i=$V_cod_start+1;$i<$V_cod_end+1;$i++){
				$resection_length=$i-$V_cod_start;		# break point in coding end is >0 
				$pos=$chr.'.'.$i;
				if($hash_pos_breakPoint{$pos}>0){
					@temp=split /\t/,$hash_pos_breakPoint_ID{$pos};		# find the read ID of each position
					for($j=0;$j<$hash_pos_breakPoint{$pos};$j++){
						print OUT "$hash_resection_length{$pos}\t$hash_micro_ID{$temp[$j]}\n";
				#		print OUT "$hash_resection_length{$pos}\t$hash_micro_ID{$temp[$j]}\t$line\n";
						if($hash_resection_length{$pos}>$resection_error){
								$mark_V_segment{$segment_name}++;		# mark this V-segment have wired surper long resection
						}
						print OUT_CE "$resection_length\n";
				#		print OUT "$resection_length\n$hash_pos_breakPoint_line{$pos}\n\n";
						$print_count++;
					}
				}
			}
			for($i=$V_sig_start+1;$i<$V_sig_end+1;$i++){
				$resection_length=$i-$V_sig_end;		# break point in signal end is <0 
				$pos=$chr.'.'.$i;
				if($hash_check_2{$pos}==1){		# make sure this position have not be counted yet
				}
				else{
					if($hash_pos_breakPoint{$pos}>0){
						@temp=split /\t/,$hash_pos_breakPoint_ID{$pos};		# find the read ID of each position
						for($j=0;$j<$hash_pos_breakPoint{$pos};$j++){
							print OUT "-$hash_resection_length{$pos}\t$hash_micro_ID{$temp[$j]}\n";
					#		print OUT "-$hash_resection_length{$pos}\t$hash_micro_ID{$temp[$j]}\t$line\n";
							if($hash_resection_length{$pos}>$resection_error){
								$mark_V_segment{$segment_name}++;		# mark this V-segment have wired surper long resection
							}
							$temp_rec=0-int($resection_length);
						#	print OUT_SE "$temp_rec\n";
							print OUT_SE "$resection_length\n";
						#	print OUT "$resection_length\n$hash_pos_breakPoint_line{$pos}\n\n";
							$print_count++;
						}
					}
				}
				$hash_check_2{$pos}=1;
			}
		}		
	}
}

#########################
# print out the V segment have long resection reads
######################### 

open(IN_V_segment,"$ARGV[1]")||die("Can not open V segment file\n");
while($line=<IN_V_segment>){
	chomp $line;
	@split=split /\t/,$line;
	$segment_name=$split[5];
	if($mark_V_segment{$segment_name}>4){
#		print OUT_long_resection_V "$line\t$mark_V_segment{$segment_name}\n";
	}
	else{
#		print OUT_filterOut_long_resection_V "$line\n";
	}
}

#########################
# print out percentage of microhomology length
#########################
for($i=0;$i<11;$i++){
	$per=$hash_micro{$i}/$total_read;
	print OUT_micro_length_percentage "$i\t$per\n";
}

#########################
# print the summary of number
#########################
$per_micro=$total_micro/$total_read;
$per_insert=$total_insert/$total_read;
$per_resect=$print_count/$total_read;
print OUT_summary "Total_reads\t$total\t1\nmicrohomology_read\t$total_micro\t$per_micro\ninsertion_reads\t$total_insert\t$per_insert\nresect_reads\t$per_resect\n";


print "$total\n$total_read\n$print_count\n";	
			
	