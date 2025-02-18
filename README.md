# Here we describe Panalyze that can make and analyse pangenome variation graphs (PVGs).

This was mainly designed with virus genomes in mind. 

## Basic information

It has been developed primarly for virus genome, principally large DNA viruses. It is still in development, and will be containerised in docker.

## Main components:
### [1]  DOWNLOAD: (optional)
     [i]   Download all 'complete genomes' or 'genomic sequences' from Nucleotide matching your text query of a specific defined organism (using the [orgn] search condition).
### [2]  ALIGN: (optional)
     [i]  Align the genomes with MAFFT -> see results/align/
     [ii] Construct a phylogeny with RAxML using a GTR+G4 substitution model -> see results/align/
### [3]  TREE: (optional)
     [i]   Visualise the phylogeny with R -> tree.png
### [4]  MAKE_PVG: (core)
     [i]   Construct a PVG using PGGB based on a 90% identity threshold and a match length of 1 Kb -> results/PVG/pggb.gfa
### [5]  VIZ1: (core)
     [i]   Create a PVG visualisation PNG with VG's view function and dot -> results/vg/out.vg.png - note this file size is very large and may take time to resolve on any visualisation tool
### [6]  ODGI: (core)
     [i]   Create OG file with odgi -> results/odgi/out.og
     [ii]  Extract odgi PVG metrics from OG file -> results/odgi/odgi.stats.txt
### [7]  OPENNESS_PANACUS: (core)
     [i]   Get the number of haplotypes present (usually equal to the number of genomes) -> results/panacus/haplotypes.txt
     [ii]  Run Panacus to get the rates of PVG growth as more samples are added -> results/panacus/histgrowth.node.tsv
     [iii] Visualise the PVG growth -> results/panacus/histgrowth.node.pdf
### [8]  OPENNESS_PANGROWTH: (core)
     [i]   Make a folder SEQS and split the genomes into individual files inside it.
     [ii]  Add the genomes to a fastix shared file and run pangrowth 
     [iii] Plot the allele frequency spectrum (AFS) -> results/pangrowth/pangrowth.pdf
     [iv]  Plot the PVG growth -> results/pangrowth/growth.pdf
     [v]   Plot the PVG core size -> results/pangrowth/p_core.pdf
### [9]  PATH_FROM_GFA: (core)
     [i]   Get the sample names from the PVG
### [10] VCF_FROM_GFA: (core)
     [i]   Use gfautil to convert the GFA to VCF -> results/vcf/gfavariants.vcf
### [11] VCF_PROCESS: (core)
     [i]   Use gfautil and Bash to quantify the pairwise difference in genome coordinates across samples.
     [ii]  Visualise this using R -> results/vcf/variation_map-basic.pdf
     [iii] Get the SNP density across the genomes -> results/vcf/mutation_density.pdf
     [iv]  Count the AFS based on this SNP data -> out file
### [12] GETBASES: (core)
     [i]   Generate a BED file based on the odgi file -> result/out.bed
     [ii]  Get the sequenc length per genome
### [13] VIZ2: (core)
     [i]   Create large-scale PVG visualisation using Odgi's viz function -> result/out.viz.png
### [14] HEAPS: (core)
     [i]   Run Odgi's heaps function across all samples to get rate of PVG growth -> results/heaps/heaps.txt
### [15] HEAPS_Visualize: (core)
     [i]   Visualise the output from Odgi's heaps function in HEAPS -> results/heaps/heaps.pdf
### [16] PAVS: (core)
     [i]   Use Odgi to get presence-absence variants (PAVs) (a large file!)
     [ii]  Quantify the number of PAVs -> result/out.flatten.fa
### [17] PAVS_plot: (core)
     [i]   Visualise the PAVs from PAVS -> result/pavs/out.flatten.pavs.pdf
### [18] COMMUNITIES: (core)
     [i]   Use wfmash to quantify the number of communities based on a 90% similarity threshold and at least 6 mappings per segment.
     [ii]  Convert these mapping into a network that is visualised -> genomes.mapping.paf
### [19] PAFGNOSTIC: (core)
     [i]   Create a text file of the mapping from the COMMUNITIES -> pafgnostic.txt
### [20] GFAstat: (core)
    [i]    Compute key PVG metrics with GFAstats  -> results/gfastat/gfa.stats.txt
    [ii]   Get the genome lengths -> results/gfastat/genome.lengths.txt
### [21] BUSCO: (optional)
     [i]   Use Busco to count the number of BUSCO genes present.

## How to run

Clone the directory

    git clone https://github.com/downingtim/LSDV-PVG/

Go to the folder

    cd LSDV-PVG

Run in Nextflow given a template YML file and an example FASTA file

For example, we can examine a smnall set of goatpox virus (GTPV) genomes:

    nextflow run main.nf --config template.GTPV.yml --reference test_genomes.GTPV.fa

The test data presented is 6 GTPV genomes as a large DNA example. This should run within 186 seconds.

In another example, we can examine a smnall set of foot-and-mouth virus (FMDV) genomes:

    nextflow run main.nf --config template.FMDV_WRL.A.yml --reference test_genomes.FMDV_WRL.A.fa

The test data presented is 9 FMDV genomes as a ssRNA virus example. This should run within 120 seconds.

In your own template YML file, you will need to define the dataset name, number of haplotypes, max number of CPUs available, minimum expected genome size, sample name filtering if using the download function, and the BUSCO clade (if relevant).

## How does it work?

The input data is "test_genomes.GTPV.fa" in this example. You can switch this to your own FASTA file input: in main.nf, Panalyze has a Download module which is not active by default.
The modules folder contains the processes, which are called by main.nf. These may call tools and scripts in other folders like bin.

## Credits

Tim Downing, Chandana Tennakoon, Thibaut Freville
