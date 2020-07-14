if(@ARGV<1){
	print "
		This script is to filter raw read that read 1 can locate within 10 bp beyond designed bait primer
		Usage: 
			perl *pl bait_10bp_IgKj4.fa read1.fq read2.fq
		Output:
			r1.fq.filter20bait
			r2.fq.filter20bait
			r2.notMatched_20bait";
}

open(IN_bait,"$ARGV[0]")||die("Can not open bait sequence\n");
open(IN_read1,"$ARGV[1]")||die("Can not open read1 sequence\n");
open(IN_read2,"$ARGV[2]")||die("Can not open read2 sequence\n");

my $out1=$ARGV[1].".filter20bait";
my $out2=$ARGV[2].".filter20bait";
my $out3=$ARGV[1].".notMatched_20bait";
open(OUT3,">$out3")||die("Can not open out3 file\n");

open(OUT1,">$out1")||die("Can not open out1 file\n");
open(OUT2,">$out2")||die("Can not open out2 file\n");

while($line=<IN_bait>){
	chomp $line;
	if($line=~/bait/){
	}
	else{
		$bait_head=substr($line,0,10);
		$bait_end=substr($line,-10,10);print "$bait_head\n$bait_end\n";
	}
}

$i=0;
$j=0;
while($line=<IN_read1>){
	$i++;
	chomp $line;
	if($i%4==1){
		@split=split /\s/,$line;
		$temp=$line;
		$mark=$split[0];
	}
	elsif($i%4==2){
		$head=substr($line,0,45);  # take 35bp of start sequence
		if(($head=~/$bait_head/)&&($head=~/$bait_end/)){
			$print=1;
			$hash_print{$mark}=1;
		}
		else{
			$print=0; print OUT3 "$line\n"; #print out the sequence without 
		}
	$temp=$temp."\n".$line;
	}
	elsif($i%4==0){
		if($print==1){
			print OUT1 "$temp\n$line\n";
		}
	}
	else{
		$temp=$temp."\n".$line;
	}
}
$i=0;
while($line=<IN_read2>){
	$i++;
	chomp $line;
	if($i%4==1){
		@split=split /\s/,$line;
		$mark=$split[0];
		if($hash_print{$mark}==1){
			$print=1;
		}
		else{
			$print=0;
		}
		$temp=$line;
	}
	elsif($i%4==0){
		if($print==1){
			print OUT2 "$temp\n$line\n";
		}
	}
	else{
		$temp=$temp."\n".$line;
	}
}


		
