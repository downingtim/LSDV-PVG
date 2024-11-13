# This is the code and data repository for the tool Panalyze that can make and analyse pangenome variation graphs (PVGs).

## Basic information

It has been developed primarly for virus genome, principally large DNA viruses.
It is still in development, and will be containerised in docker eventually.
It has been tested primarly in livestock poxviruses, such as lumpy skin disease virus (LSDV), sheeppox virus (SPPV) and goatpox virus (GTPV).
The test data presented is 6 GTPV genomes as an example.
It has also been tested on ssRNA viruses, such as foot-and-mouth disease virus (FMDV) and Rift Valley fever virus (RVFV).

## How does it work?

The folder CODE/ contains the code and relevant packages/tools to support running Panalyze.

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
