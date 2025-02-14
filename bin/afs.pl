#!/usr/bin/env perl

open(IN, "$ARGV[0]");
$n=$ARGV[0].".out";
open(OUT, ">$n");
# get AFS for pangenomes
# plot fraction of bases (y-axis) vs
# number of samples in which it was found (x-axis)

@a = <IN>;
%h={};

for $e (0..$#a){
    if($a[$e]=~ /snv/){
      @r=split/\s+/,$a[$e];
#       if(($r[1]>1000)&&($r[1] <148001)){
        if(exists($h{$r[1]})){ $h{$r[1]}=$h{$r[1]}+1; }
        else { $h{$r[1]}=1; }
  #     }
    }
}
for $ee (sort {$a <=> $b} keys %h){
    if(length($h{$ee})>0){ print OUT $h{$ee}/$ARGV[1],"\n"; }
}
close(OUT);
exit;
