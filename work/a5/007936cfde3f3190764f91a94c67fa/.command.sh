#!/bin/bash -ue
odgi build -g pggb.gfa -o out.og 
odgi stats -m -i out.og -S > odgi.stats.txt
