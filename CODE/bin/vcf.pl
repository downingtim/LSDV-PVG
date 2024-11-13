#!/usr/bin/env perl

open(IN, $ARGV[0]);
$ain=<IN>;
print "$ain\n";
chomp($ain);
$argument = $ARGV[0];
$argument =~ s/VCF\/path1//g;

#open(OUT, ">../../../CURRENT/variation_map.txt");

for $a (1000..148500){ 
    #
	print "bin/gfautil --quiet -t 15 -i $argument/*.gfa snps --ref $ain --snps $a >> temp2 \n";
    system("bin/gfautil --quiet -t 15 -i $argument/*.gfa snps --ref $ain --snps $a >> temp2 ");
}

open(IN2, "more temp2 | sort -nk 3 | grep -v temp2  | grep -v \: |");

while(<IN2>){
    if($_ =~ /path/) { ; }
    else { print OUT "$_"; }
}
#exit; 
