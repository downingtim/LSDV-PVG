# This is the code and data repository for the lumpy skin disease virus (LSDV) pangenome variation graph (PVG).

## 1 # The complete LSDV PVG 

The foler CODE contains NF code (main.nf, modules/processes.nf) to run the PVG creation & analysis pipeline on a fresh dataset. This uses scripts in its bin/ folder, which output data in the folder named CURRENT. Large PAV and BED files may be omitted due to their large file sizes. 

For these samples, the download process (bin/download.R0 takes you through file download from NCBI Nucleotide, metadata extraction, symbol removal, an initial phylogeny, a genome length check, reporting of excluded sequences, and a final output FASTA.

We created PVGs using the pangenome graph builder (PGGB) v0.5.4 pipeline (https://github.com/pangenome/pggb, Garrison et al 2023), which initiated all-to-all alignment with wfmash v0.7.0 (Guarracino et al 2021) subject to a threshold of 90% identity per Kb. Like Minigraph-Cactus (REF), PGGB keeps all paths and was used to reconstruct the input haplotypes directly from the graph. PGGB created the PVG using Seqwish (Garrison & Guarracino 2023) using a k-mer of 19 bp and a window size of 27 bp, and then sorted and ordered to produce a progressively linearised PVG based on partial order alignment with smoothxg v0.7.0.0 (https://github.com/pangenome/smoothxg), which also remove redundant nodes. Multiqc v1.14 (Ewels et al 2016) was used to collate and assess PVG metrics. The PVGs were generated in Graphical Fragment Assembly (GFA) format. These PVGs were converted with VG’s convert function, indexed with VG’s autoindex function and the corresponding genome files were indexed with SAMtools.

An advantage of graphical formats is more computationally scalable indexing through bidirectional graph BWT (GBWT) indices, which store haplotype paths as node sequences (Siren et al 2020, Baaijens et al 2022). Consequently, we created GBWT indexed versions of our PVGs using the corresponding haplotypes available (one-, three-, six- and all-sample PVGs) with VG’s gbwt function. We used (defaults of a) minimiser length of 29 and window size of 11 for the Giraffe minimiser index. The GBWT indices were created using the greedy path-cover algorithm across all paths. These indices were converted with VG. VG’s snarls function was used to create a snarls file, which was used to create a distance index of the PVG’s snarl decomposition such the minimum distance between nodes in the graph could be calculated. VG’s minimizer function used these distances to get the PVG locations of the k-mers per haplotype and the corresponding PVG distances. VCFs were obtained from VG’s deconstruct function for the GBWT-indexed graphs.

We analysed and manipulated these PVGs using odgi (optimized dynamic genome/graph implementation) v0.8.3 (Guarracino et al 2022, Garrison et al 2023), including quantifying numbers of SNPs, indels and compound. Gfautil v0.3.2 (Kubica 2023) was used to get mutation rates and coordinates for samples at each path and site. Subgraphs were selected gfatools v0.4 (Li 2023), evaluated with gfastats v1.3.6 (Formenti 2023) and gfastar v0.1 (Formenti 2023), and visualised with waragraph (Fischer et al 2023) and variation graph (vg) v1.43.0 (Garrison et al 2018). Presence-absence variant (PAVs) were identified using odgi. The data was visualised with odgi (Guarracino et al 2022) and R package ggplot2 v2_3.4.4 (Wickham 2016). 

The PVGs were annotated using Prokka (Seeman 2014) as above and were indexed with SAMtools v1.19 (Danecek et al 2021). Using the reference genome annotation, Busco v5.5.0 (Manni et al 2021) found a mean of 155.8 complete unique CDSs per genome (standard deviation 1.5). The Prokka GFFs were converted to GTF with gffread (Pertea & Pertea 2020) and parsed. The annotation position information was merged with the odgi (Guarracino et al 2022) rendering of the PVGs, creating one CSV annotation file per PVG, which was visualised using Bandage v0.8.1 (Wick et al 2015). 
![image](https://github.com/user-attachments/assets/d28e6b19-6856-45d1-b0c2-636c06fd34c4)


## 2 # The 3- and 6-sample LSDV PVGs

The folders 3_SAMPLE_PVG and 6_SAMPLE_PVG contain data for the PVG created from 3 and 6 representative LSDV genomes for computationally efficient investigation.


## 3 # NF_PANGENOME_TEST

This folder contains tests run of the nf-core pangenome tool for LSDV, GPV, SPV and CaPV.


## Credits

Tim Downing, Chandana Tennakoon, Thibaut Freville
