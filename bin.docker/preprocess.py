#!/usr/bin/env python3
import sys
Header = "";Seq = ""

Seq = ""
Padding = sys.argv[2]
Min_Length = int(sys.argv[3])

BED_file = open("trim.bed", "w")

def pad_sequence(seq, Header, padding):
    if len(seq) >= Min_Length:
        return seq
    BED_file.write(f"{Header[1:].split()[0]}\t0\t{len(seq)}\n")
    pad_len = Min_Length - len(seq)
    repeated_padding = (padding * ((pad_len // len(padding)) + 1))[:pad_len]
    return seq + repeated_padding

with open(sys.argv[1], 'r') as file:
    for line in file:
        L = line.strip()
        if not L:
            continue
        if line[0] == ">":
            if Seq:
                Seq = pad_sequence(Seq, Header, Padding)
                print(Header);print(Seq)
            Header = L
            Seq = ""
        else:
            Seq += L

Seq = pad_sequence(Seq, Header, Padding)
print(Header);print(Seq)