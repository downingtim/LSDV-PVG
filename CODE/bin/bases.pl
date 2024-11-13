#!/usr/bin/env perl

open(IN,"$ARGV[0]"); # input FASTA 
open(OUT, ">../../../CURRENT/bases.txt");
$b = join("",<IN>);
@a=split(//,$b);
$b=~ s/N//g;
@a2=split(//,$b);
print OUT "Total length  = ";
print OUT ($#a-10);
print OUT "\tGap-free length  = ";
print OUT ($#a2-10);
print OUT "\n";
close(IN);
close(OUT);
exit;
