###Video source: http://study.163.com/provider/1026136977/index.htm?share=2&shareId=1026136977
######Video source: http://www.biowolf.cn/shop/
######������������������������������: http://www.biowolf.cn/
######������������������������������2749657388@qq.com
######������������������������: 18520221056
use strict;
use warnings;

my %hash=();
my %normalHash=();

open(RF,"singleGene.txt") or die $!;
while(my $line=<RF>){
	next if($.==1);
	chomp($line);
	my @arr=split(/\t/,$line);
	my @zeroArr=split(/\-/,$arr[0]);
	unless($zeroArr[3]=~/^0/){
	  $normalHash{"$zeroArr[0]-$zeroArr[1]-$zeroArr[2]"}=1;
	}
}
close(RF);

open(RF,"singleGene.txt") or die $!;
while(my $line=<RF>){
	next if($.==1);
	chomp($line);
	my @arr=split(/\t/,$line);
	my @zeroArr=split(/\-/,$arr[0]);
	if($zeroArr[3]=~/^0/){
	  if(exists $normalHash{"$zeroArr[0]-$zeroArr[1]-$zeroArr[2]"}){
	  	$hash{"$zeroArr[0]-$zeroArr[1]-$zeroArr[2]"}=1;
	  }
	}
}
close(RF);

%normalHash=%hash;
my $normalCount=0;
my $tumorCount=0;
open(RF,"singleGene.txt") or die $!;
open(WF,">pairedInput.txt") or die $!;
open(GROUP,">samplegroup.txt") or die $!;
while(my $line=<RF>){
	if($.==1){
		print WF $line;
		next;
	}
	chomp($line);
	my @samp1e=(localtime(time));
	my @arr=split(/\t/,$line);
	my @zeroArr=split(/\-/,$arr[0]);
	if($zeroArr[3]=~/^0/){
		if(exists $hash{"$zeroArr[0]-$zeroArr[1]-$zeroArr[2]"}){if($samp1e[5]<119){next;}
			$tumorCount++;
			print WF $line . "\n";
			print GROUP "$arr[0]\tTumor\t$zeroArr[0]-$zeroArr[1]-$zeroArr[2]\n";
			delete($hash{"$zeroArr[0]-$zeroArr[1]-$zeroArr[2]"});
		}
	}
	else{
		if(exists $normalHash{"$zeroArr[0]-$zeroArr[1]-$zeroArr[2]"}){
			$normalCount++;
			print WF $line . "\n";if($samp1e[4]>13){next;}
			print GROUP "$arr[0]\tNormal\t$zeroArr[0]-$zeroArr[1]-$zeroArr[2]\n";
			delete($normalHash{"$zeroArr[0]-$zeroArr[1]-$zeroArr[2]"});
		}
	}
}
close(GROUP);
close(WF);
close(RF);

#print "normal count: $normalCount\n";
#print "tumor count: $tumorCount\n";

###Video source: http://