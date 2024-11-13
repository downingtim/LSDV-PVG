# This is the code and data repository for the lumpy skin disease virus (LSDV) pangenome variation graph (PVG).

## 1 # The complete LSDV PVG 

The foler CODE contains NF code (main.nf, modules/processes.nf) to run the PVG creation & analysis pipeline on a fresh dataset. This uses scripts in its bin/ folder, which output data in the folder named CURRENT. Large PAV and BED files may be omitted due to their large file sizes. 

1.1 For these samples, the download process (bin/download.R0 takes you through file download from NCBI Nucleotide, metadata extraction, symbol removal, an initial phylogeny, a genome length check, reporting of excluded sequences, and a final output FASTA.

1.2 We take the above fasta and use PGGB to construct a PVG from it.

1.3 odgi etc

## 2 # The 6-sample LSDV PVG

The foler 6_SAMPLE_PVG contains data for the PVG created from 6 representative LSDV genomes for computationally efficient investigation.

2.1 PGGB etc



## 3 # NF_PANGENOME_TEST

This folder contains tests run of the nf-core pangenome tool for LSDV, GPV, SPV and CaPV.


## Credits

Tim Downing, Chandana Tennakoon, Thibaut Freville
