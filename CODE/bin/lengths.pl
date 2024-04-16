#!/usr/bin/env perl


if(length($ARGV[0])<2){
    $ARGV[0] = "goatpox_2_virus_2_2_.fasta"; }

open(IN, $ARGV[0]);
open(OUT, ">./genome.lengths.txt");
@a=split/>/,(join"",<IN>); 
for $e (1..$#a){
	($n,$s)=split/\n/,$a[$e], 2;
	print OUT "$n\t", length($s),"\n"; 
} # end for each sample
exit;
