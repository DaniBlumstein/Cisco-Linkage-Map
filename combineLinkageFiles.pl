#!/usr/bin/perl -w
use strict;

my$combo1_lepmap='/scratch-lustre/wlarson/Dani/binary+code/data/consensus_files/hap_combo2.call';
my$combo2_lepmap='/scratch-lustre/wlarson/Dani/binary+code/data/consensus_files/diploid_combo3_lepmap.call';

open(FILE1,"<$combo1_lepmap")||die "cannot open $combo1_lepmap:$!";
open(FILE2,"<$combo2_lepmap")||die "cannot open $combo2_lepmap:$!";

my$F1line1=<FILE1>;
my$F1line2=<FILE1>;
my$F1line3=<FILE1>;
my$F1line4=<FILE1>;
my$F1line5=<FILE1>;
my$F1line6=<FILE1>;
chomp $F1line1;
chomp $F1line2;
chomp $F1line3;
chomp $F1line4;
chomp $F1line5;
chomp $F1line6;
my%File1genos;
my%allGenos;
my@File1Samples=split "\t", $F1line1;
my$F1sampleNum=scalar@File1Samples-1;
my$F1blankGenos;
foreach my$i (1..$F1sampleNum){
	$F1blankGenos.="0 0 0 0 0 0 0 0 0 0";
	unless($i==$F1sampleNum){
		$F1blankGenos.="\t";
	}
}
while(my$line=<FILE1>){
	chomp $line;
	my($locus,$geno)=split "\t", $line, 2;
	$File1genos{$locus}=$geno;
	$allGenos{$locus}++;
}
close FILE1;

my$F2line1=<FILE2>;
my$F2line2=<FILE2>;
my$F2line3=<FILE2>;
my$F2line4=<FILE2>;
my$F2line5=<FILE2>;
my$F2line6=<FILE2>;
chomp $F2line1;
chomp $F2line2;
chomp $F2line3;
chomp $F2line4;
chomp $F2line5;
chomp $F2line6;
my%File2genos;
my@File2Samples=split "\t", $F2line1;
my$F2sampleNum=scalar@File2Samples-1;
my$F2blankGenos;
foreach my$i (1..$F2sampleNum){
	$F2blankGenos.="0 0 0 0 0 0 0 0 0 0";
	unless($i==$F2sampleNum){
		$F2blankGenos.="\t";
	}
}
while(my$line=<FILE2>){
	chomp $line;
	my($locus,$geno)=split "\t", $line,2;
	$File2genos{$locus}=$geno;
	$allGenos{$locus}++;
}

print "$F1line1\t$F2line1\n";
print "$F1line2\t$F2line2\n";
print "$F1line3\t$F2line3\n";
print "$F1line4\t$F2line4\n";
print "$F1line5\t$F2line5\n";
print "$F1line6\t$F2line6\n";

foreach my$key (keys %allGenos){
	print "$key\t";
	if(exists $File1genos{$key}){
		print "$File1genos{$key}\t";
	}else{
		print "$F1blankGenos\t";
	}
	if(exists $File2genos{$key}){
		print "$File2genos{$key}\n";
	}else{
		print "$F2blankGenos\n";
	}
}
