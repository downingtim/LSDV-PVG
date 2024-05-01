#!/usr/bin/env perl

open(IN,"$ARGV[0]"); # input FASTA 
open(OUT, ">../../../CURRENT/bases.txt");
@b = split(/^>/m,join("",<IN>));

print OUT "Sample\tTotal_Length\tGapless_Length\n";

for $e (1..$#b){
    ($n,$s)=split(/\n/,$b[$e],2);
    $s=~ s/\s+//g;
    print OUT $n,"\t",length($s);
    $s =~ s/n//g;
    print OUT "\t",length($s),"\n";
}
print OUT "\n";
close(IN);
close(OUT);
exit;
