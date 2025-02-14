# Here we describe Panalyze that can make and analyse pangenome variation graphs (PVGs).

This was mainly designed with virus genomes in mind. 

## Basic information

It has been developed primarly for virus genome, principally large DNA viruses. It is still in development, and will be containerised in docker.

## Main components:
### [1]  DOWNLOAD: (optional)
     [i]   Download all 'complete genomes' or 'genomic sequences' from Nucleotide matching your text query of a specific defined organism (using the [orgn] search condition).
### [2]  ALIGN: (optional)
     [i]  Align the genomes with MAFFT
     [ii] Construct a phylogeny with RAxML using a GTR+G4 substitution model.
### [3]  TREE: (optional)
     [i]   Visualise the phylogeny with R
### [4]  MAKE_PVG: (core)
     [i]   Construct a PVG using PGGB based on a 90% identity threshold and a match length of 1 Kb.
### [5]  VIZ1: (core)
     [i]   Create a PVG visualisation PNG with VG's view function and dot.
### [6]  ODGI: (core)
     [i]   Create OG file with odgi
     [ii]  Extract odgi PVG metrics from OG file -> odgi.stats.txt
### [7]  OPENNESS_PANACUS: (core)
     [i]   Get the number of haplotypes present (usually equal to the number of genomes)
     [ii]  Run Panacus to get the rates of PVG growth as more samples are added
     [iii] Visualise the PVG growth -> histgrowth.node.pdf
### [8]  OPENNESS_PANGROWTH: (core0
     [i]   Make a folder SEQS and split the genomes into individual files inside it.
     [ii]  Add the genomes to a fastix shared file and run pangrowth
     [iii] Plot the allele frequency spectrum (AFS) -> pangrowth.pdf
     [iv]  Plot the PVG growth -> growth.pdf
     [v]   Plot the PVG core size -> p_core.pdf
### [9]  PATH_FROM_GFA: (core)
     [i]   Get the sample names from the PVG
### [10] VCF_FROM_GFA: (core)
     [i]   Use gfautil to convert the GFA to VCF -> gfavariants.vcf
### [11] VCF_PROCESS: (core)
     [i]   Use gfautil and Bash to quantify the pairwise difference in genome coordinates across samples.
     [ii]  Visualise this using R -> varaition_map-basic.pdf
     [iii] Get the SNP density across the genomes -> mutation_density.pdf
     [iv]  Count the AFS based on this SNP data -> out file
### [12] GETBASES: (core)
     [i]   Generate a BED file based on the odgi file
     [ii]  Get the sequenc length per genome
### [13] VIZ2: (core)
     [i]   Create large-scale PVG visualisation using Odgi's viz function -> out.viz.png
### [14] HEAPS: (core)
     [i]   Run Odgi's heaps function across all samples to get rate of PVG growth -> heaps.txt
### [15] HEAPS_Visualize: (core)
     [i]   Visualise the output from Odgi's heaps function in HEAPS -> heaps.pdf
### [16] PAVS: (core)
     [i]   Use Odgi to get presence-absence variants (PAVs) (a large file!)
     [ii]  Quantify the number of PAVs -> flatten.pavs.count.txt
### [17] PAVS_plot: (core)
     [i]   Visualise the PAVs from PAVS -> out.flatten.pavs.pdf
### [18] COMMUNITIES: (core)
     [i]   Use wfmash to quantify the number of communities based on a 90% similarity threshold and at least 6 mappings per segment.
     [ii]  Convert these mapping into a network that is visualised -> genomes.mapping.paf
### [19] BUSCO: (core)
     [i]   Use Busco to count the number of BUSCO genes present.
### [20] PAFGNOSTIC: (core)
     [i] Creata text file of the mapping from the COMMUNITIES -> pafgnostic.txt
### [21] GFAstat: (core)
    [i]   Compute key PVG metrics with GFAstats and get the genome lengths -> gfa.stats.txt

## How to run

Clone the directory

git clone https://github.com/downingtim/LSDV-PVG/

Go to the folder

cd LSDV-PVG

Run in Nextflow given a template YML file and an example FASTA file

nextflow run main.nf --config template.GTPV.yml --reference test_genomes.GTPV.fa

The test data presented is 6 GTPV genomes as an example. This should run within 3 minutes.

## How does it work?

The input data is "test_genomes.GTPV.fa" in this example. You can switch this to your own FASTA file input: in main.nf, Panalyze has a Download module which is not active by default.
The modules folder contains the processes, which are called by main.nf. These may call tools and scripts in other folders like bin.

## Credits

Tim Downing, Chandana Tennakoon, Thibaut Freville
