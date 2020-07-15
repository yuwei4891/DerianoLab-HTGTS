if(@ARGV<1){
	print "This script is to produce the input file for Circos plot
		Usage: perl *pl r2.$sample.*.bam.bed   LinkCutoff
		Output: *.links
			*.hismap\n";
}

open(IN,"$ARGV[0]")||die("Can not open in file\n");
my $out_link=$ARGV[0].".links.txt";
my $out_hismap=$ARGV[0].".hismap.txt";
open(OUT_LINK,">$out_link")||die("Can not open link file\n");
open(OUT_HISMAP,">$out_hismap")||die("Can not open hismap file\n");

$cutoff=int($ARGV[1]);
if($ARGV[0]=~/JK/){  # plot JK bait library
	$MYC_start=70673542;
	$MYC_end=70673642;
	$MYC_chr="mm6";
}
elsif($ARGV[0]=~/myc/){
	$MYC_start=61817507; # position of MYC break in mm9
	$MYC_end=61818430;
	$MYC_chr="mm15";
}
else{
}

my $IgK_start=70672000;
my $IgK_end=70675000;

$i=0;
$color="red";
$window=1000; # 1000 bp as bin count for read number

my $int_IgK_start=int($IgK_start/$window);
my $int_IgK_end=int($IgK_end/$window);


while($line=<IN>){
	chomp $line;
	@split=split /\t/,$line;
	$chr=$split[0];
	$chr=~s/chr/mm/;
	$start=$split[1];
	$end=$split[2];
	$inter_start=int($start/$window)*$window;
	$pos=$chr.'.'.$inter_start;
	if(($start>$IgK_start)&&($end<$IgK_end)){ 	# if this read is JK4 or Myc; do not count the number
	}
	else{
		$hash{$pos}++; #print "$pos\n";	
	}
}

while(($key,$value)=each %hash){
	@temp=split /\./,$key;#print "$key\t$value\n";
	$end=$temp[1]+$window; 
	if(($temp[1]>$IgK_start)&&($end<$IgK_end)){	# do not plot the hisgramme of IgK !!!!!!!!!!!!!!!
	}
	else{
		print OUT_HISMAP "$temp[0]\t$temp[1]\t$end\t$value\n";
	}
}

open(IN,"$ARGV[0]")||die("Can not open in file\n");
while($line=<IN>){
        chomp $line;
        @split=split /\t/,$line;
        $chr=$split[0];
        $chr=~s/chr/mm/;
        $start=$split[1];
        $end=$split[2];
        $inter_start=int($start/$window)*$window;
        $pos=$chr.'.'.$inter_start;
	if($hash_print{$pos}==1){	# only plot the line which not printed
	}
	else{        
	        if(($chr=~/_random/)||($chr=~/M/)){
		}
	        else{
			if($hash{$pos}>$cutoff){ #cutoff for the link 
       		         	$j++;
				print OUT_LINK "$j\t$MYC_chr\t$MYC_start\t$MYC_end\tcolor=$color\n$j\t$chr\t$start\t$end\tcolor=$color\n";
#                	print OUT_LINK "$j\t$MYC_chr\t$MYC_start\t$MYC_end\tcolor=$color\t$hash{$pos}\n$j\t$chr\t$start\t$end\tcolor=$color\t$hash{$pos}\n";
			}
        	}
		$hash_print{$pos}=1;
	}
}
				
