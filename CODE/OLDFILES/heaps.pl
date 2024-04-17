#!/usr/bin/env perl

for $e (1..$ARGV[0]){
    system("odgi heaps -i CURRENT/*.og -n 1000 -d $e -t 12 | sort -nk 2 | tail -n 1  ");
}

exit;
