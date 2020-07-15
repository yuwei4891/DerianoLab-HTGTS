if(@ARGV<1){
	print "This script is to classify read2 recombination site into Tcr, Igh , IgK, IgL et al
		Usage: perl *pl Ensembl.bed r2.$sample.*.notIn_IgkJ4.bam.bed 
		Output: 
			*class\n";
}

open(IN1,"$ARGV[0]")||die("Can not open in1 file\n");
open(IN2,"$ARGV[1]")||die("Can not open in2 file\n");
my $out=$ARGV[1].".class";
open(OUT,">$out")||die("Can not open out file\n");

$window=100;


#12	114480000	117260000	Igh
#14	53040000	54860000	Tcrad
#6	40820000	41540000	Tcrb
#6	67500000	70700000	IgK
#16	19060000	19280000	Igl
#13	19260000	19460000	Tcrg

while($line=<IN1>){
	chomp $line;
	@split=split /\t/,$line;
	$chr="chr".$split[0];
	$start=int($split[1]/$window);
	$end=int($split[2]/$window);
#	$name=$split[3];
#	@name=split /\:/,$split[3];
#	$define=$name[1];
	$define=$split[3];
	for($i=$start;$i<$end+1;$i++){
		$pos=$chr.'.'.$i;
		if(exists($hash{$pos})){
			$hash{$pos}=$hash{$pos}.':'.$define;
		}
		else{
			$hash{$pos}=$define;
		}
	}
}

#chr9	90000090	90000119	M01626:300:000000000-AWN8Y:1:1102:24187:12518_2:N:0:1	250	+
#chr12	17299926	17300103	M01626:300:000000000-AWN8Y:1:1110:17553:21045	60	+

while($line=<IN2>){
	chomp $line;
	@split=split /\t/,$line;
	$chr=$split[0];
	$start=int($split[1]/$window);;
	$end=int($split[2]/$window);
	$mark=0;
	for($i=$start;$i<$end+1;$i++){
		$pos=$chr.'.'.$i;
		if($mark==0){
			if($hash{$pos}=~/Igh/){
				$IgH++;
			}
			elsif($hash{$pos}=~/IgK/){
				$IgK++;
			}
			elsif($hash{$pos}=~/Tcrad/){
				$Tcrad++;
			}
			elsif($hash{$pos}=~/Tcrb/){
				$Tcrb++;
			}
			elsif($hash{$pos}=~/Tcrg/){
				$Tcrg++;
			}
			elsif($hash{$pos}=~/Igl/){
				$IgL++;
			}
			elsif($hash{$pos}=~/MYC/){
				$MYC++;
			}
			elsif($hash{$pos}=~/chr9/){
				$chr9++;
			}
			elsif($hash{$pos}=~/chr15/){
				$chr15++;
			}	
			else{#print "$line\n";
				$else++;
			}
			$mark=1;
		}
	}
}

if($IgH>0){
}
else{
	$IgH=0;
}
if($MYC>0){
}
else{
	$MYC=0;
}
if($Tcrad>0){
}
else{
	$Tcrad=0;
}

if($Tcrb>0){
}
else{
	$Tcrb=0;
}

if($Tcrg>0){
}
else{
	$Tcrg=0;
}

if($IgL>0){
}
else{
	$IgL=0;
}
if($chr9>0){
}
else{
	$chr9=0;
}
if($chr15>0){
}
else{
	$chr15=0;
}
if($else>0){
}
else{
	$else=0;
}

#print OUT "IgH\t$IgH\nIgK\t$IgK\nIgl\t$IgL\nTcrad\t$Tcrad\nTcrb\t$Tcrb\nTcrg\t$Tcrg\nelse\t$else\n";
print "$IgH\n$IgK\n$IgL\n$Tcrad\n$Tcrb\n$Tcrg\n$MYC\n$chr9\n$chr15\n$else\n";
		
