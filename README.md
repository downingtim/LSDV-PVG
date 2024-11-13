# This is the code and data repository for the tool Panalyze that can make and analyse pangenome variation graphs (PVGs).

## Basic information

It has been developed primarly for virus genome, principally large DNA viruses.
It is still in development, and will be containerised in docker eventually.
It has been tested primarly in livestock poxviruses, such as lumpy skin disease virus (LSDV), sheeppox virus (SPPV) and goatpox virus (GTPV).
The test data presented is 6 GTPV genomes as an example.
It has also been tested on ssRNA viruses, such as foot-and-mouth disease virus (FMDV) and Rift Valley fever virus (RVFV).

## How does it work?

The folder CODE/ contains the code and relevant packages/tools to support running Panalyze.

The input data is "test_genomes.GTPV.fa" in this example. You can switch this to your own FASTA file input: in main.nf, Panalyze has a Download module which is not active by default.

The modules folder contains the processes, which are called by main.nf. These may call tools and scripts in other folders like bin.

Prior to running Panalyze, you should run setup.R to set up the R packages you will need.

In addition, you may need to configure the environment for a range of other packages in 'readme', which could include libjemalloc, odgi, nextflow, panacus, pggb, and various Python packages.

You will need to edit "template.yml" to align it with your own configuration.
Essential: In that file, the 'genomes' and 'reference' should be the FASTA file you want to examine.
Essential: It will need the correct number of samples present ('haplotypes').
Optional: If you are not using a poxvirus and want to apply BUSCO gene analysis, you will need to change the busco database.
Optional: You can also add a text 'filter' if you are using the download function to omit uninformative regular expressions from the FASTA headers in that file for convenience.

## Credits

Tim Downing, Chandana Tennakoon, Thibaut Freville
