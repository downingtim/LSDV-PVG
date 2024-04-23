#!/usr/bin/env nextflow

/*
Download datasets of interest, perform QC and construct a phylogeny.
/*

/*
#==============================================
Enable DSL2
#==============================================
*/

nextflow.enable.dsl = 2

/*
#==============================================
Download dataset to examine (eg LSDV)
#==============================================
*/

process DOWNLOAD {
    input:
    path (faS)

    output:
    publishDir ".", mode: "copy"
    path("*2_2_.fasta"), emit: "write"

    script:
    """
    module load R
    Rscript ${workflow.projectDir}/bin/download.R
    """
}

process MAKE_PVG {
    tag {"index reference FASTA"}
    label 'pvg'

    input:
    path (refFasta)
    
    output:
    val true, emit: value1
    publishDir ".", mode: "copy"

    exec:
    def file1 = file('bin/number1.txt')
    allLines = file1.readLines()
    for( line : allLines ) { number=line }
    println(number)

    script:
    """
    # create Blast indexes for self-self blast
    # makeblastdb -dbtype nucl -in ${refFasta} -out ${workflow.projectDir}/BLAST/fasta.db 
    # blastall -p blastn -d ${workflow.projectDir}/BLAST/fasta.db -i ${refFasta} -e 1.0 -outfmt 6 > self-blast.out
    
    # SAMtools index
    bgzip ${refFasta}
    samtools faidx ${refFasta}.gz > ${refFasta}.gz.fai

    # run PGGB - you need to specify the number of haplotypes in $number
    pggb -i ${refFasta}.gz -m -S -o ${workflow.projectDir}/CURRENT/ -t 46 -p 90 -s 1k -n $number
    # if issues, remove CURRENT/ in ${workflow.projectDir}  before re-running 

    # run pangenome in nf-core - don't work as reliably as PGGB, so best avoided
    # nextflow run nf-core/pangenome --input ${refFasta}.gz --n_haplotypes $number --outdir \
    #   ${workflow.projectDir}/CURRENT  --communities  --wfmash_segment_length 1000
    #odgi stats -m -i ${workflow.projectDir}/CURRENT/FINAL_ODGI/${refFasta}.gz.squeeze.og -S > ${workflow.projectDir}/odgi.stats.txt 

    """
}

process VIZ1 {
    input:
    val ready 
    path (refFasta)

    output:
    publishDir ".", mode: "copy"

    script:
    """
    # visualise gfa - need to separate, this takes a while
    vg view -F -p -d ${workflow.projectDir}/CURRENT/*.gfa  | dot -Tpng -o ${workflow.projectDir}/CURRENT/${refFasta}.vg.png

    # get general info on gfa
    ${workflow.projectDir}/bin/gfastats ${workflow.projectDir}/CURRENT/*gfa > ${workflow.projectDir}/CURRENT/gfa.stats.txt

    # get lengths of genomes etc
    perl ${workflow.projectDir}/bin/lengths.pl ${refFasta} 

    # get Bandage stats
    /mnt/lustre/RDS-live/downing/miniconda3/bin/Bandage info ${workflow.projectDir}/CURRENT/*.gfa > bandage.stats.txt
    """
}

process ODGI {
    input:
    val ready 
    path (refFasta)

    output:
    publishDir ".", mode: "copy"

    script:
    """
    #odgi build -g ${projectDir}/CURRENT/*.gfa -o ${projectDir}/CURRENT/${refFasta}.og 
    #odgi stats -m -i ${projectDir}/CURRENT/${refFasta}.og -S > ${projectDir}/CURRENT/odgi.stats.txt 
    odgi stats -m -i ${projectDir}/CURRENT/*.og -S > ${projectDir}/CURRENT/odgi.stats.txt 
    """
}

process OPENNESS_PANACUS {
    tag {"get PVG openness panacus"}
    label 'openness_panacus'

    input:
    val ready 

    output:
    publishDir ".", mode: "copy"
    
    script:
    """
    # you need to install panancus
    # mamba install -c conda-forge -c bioconda panacus
    # get haplotypes
    grep '^P' ${workflow.projectDir}/CURRENT/*.gfa | cut -f2 > ${workflow.projectDir}/CURRENT/PANACUS/haplotypes.txt
    
    # run panacus to get data - just use defaults
    RUST_LOG=info panacus histgrowth -t4 -l 1,2,1,1,1 -q 0,0,1,0.5,0.1 -S -s ${workflow.projectDir}/CURRENT/PANACUS/haplotypes.txt ${workflow.projectDir}/CURRENT/*.gfa > ${workflow.projectDir}/CURRENT/PANACUS/histgrowth.node.tsv
 
    # visualise plot as PDF
    panacus-visualize -e ${workflow.projectDir}/CURRENT/PANACUS/histgrowth.node.tsv > ${workflow.projectDir}/CURRENT/PANACUS/histgrowth.node.pdf
    """
}

process OPENNESS_PANGROWTH {
    tag {"get PVG openness pangrowth"}
    label 'openness_pangrowth'

    input:
    val ready 
    path (refFasta)

    output:
    publishDir ".", mode: "copy"

    exec:
    def file1 = file('bin/number1.txt')
    allLines = file1.readLines()
    for( line : allLines ) { number=line }

    script:
    """
    # need to split files into individual files in folder SEQS
    # wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faSplit and chmod it
    #   ./faSplit byname ${refFasta} ${workflow.projectDir}/CURRENT/PANGROWTH/SEQS/

    # run pangrowth
    ${workflow.projectDir}/bin/pangrowth hist -k 17 -t 12 ${workflow.projectDir}/CURRENT/PANGROWTH/SEQS/*a > ${workflow.projectDir}/CURRENT/PANGROWTH/hist.txt

    # plot AFS for single dataset - not very useful! 
     python ${workflow.projectDir}/bin/plot_single_hist.py ${workflow.projectDir}/CURRENT/PANGROWTH/hist.txt ${workflow.projectDir}/CURRENT/PANGROWTH/LSDV.pangrowth.pdf

    # model growth - prepare dataset 1 - not useful
     ${workflow.projectDir}/bin/pangrowth growth -h ${workflow.projectDir}/CURRENT/PANGROWTH/hist.txt > ${workflow.projectDir}/CURRENT/PANGROWTH/LSDV
     python ${workflow.projectDir}/bin/plot_growth.py ${workflow.projectDir}/CURRENT/PANGROWTH/LSDV ${workflow.projectDir}/CURRENT/PANGROWTH/LSDV.growth.pdf

    # do core - this does not converge to a solution for n<10 samples, causes error
     ${workflow.projectDir}/bin/pangrowth core -h ${workflow.projectDir}/CURRENT/PANGROWTH/hist.txt -q 0.5 > ${workflow.projectDir}/CURRENT/PANGROWTH/LSDV_core

    if [ $number -gt 10 ]; then # if the core is too small, this crashes
       python ${workflow.projectDir}/bin/plot_core.py ${workflow.projectDir}/CURRENT/PANGROWTH/LSDV_core    ${workflow.projectDir}/CURRENT/PANGROWTH/LSDV_core.pdf
    fi
    """
}

process GET_VCF {
    tag {"get VCFs"}
    label 'vcf'

    input:
    val ready 
    path (refFasta)

    output:
    publishDir ".", mode: "copy"

    exec:
    def file1 = file('bin/number1.txt')
    allLines = file1.readLines()
    for( line : allLines ) { number=line }

    script:
    """
    ${workflow.projectDir}/bin/gfautil -i ${workflow.projectDir}/CURRENT/*.gfa gfa2vcf > ${workflow.projectDir}/CURRENT/VCF/${refFasta}.vcf
    # https://lib.rs/crates/gfautil - download this

    # need to get paths from GFA
    odgi paths -i ${workflow.projectDir}/CURRENT/*.og -L > ${workflow.projectDir}/CURRENT/VCF/paths.txt
    head -n 1 ${workflow.projectDir}/CURRENT/VCF/paths.txt > ${workflow.projectDir}/CURRENT/VCF/path1
    # these correspond to the samples, eg "K303/83_g001_s001"
    # get first sample name then give to gfautil
    # Outputs a tab-delimited list in the format:
    # <query-path-name>\t<reference base>\t<reference pos>\t<query base>\t<query pos>

    # script to create initial input to visualise with R
    perl ${workflow.projectDir}/bin/vcf.pl ../../../CURRENT/VCF/path1  #  sort out coordinates
    Rscript ${workflow.projectDir}/bin/plot_variation_map.R # plot image of variation map in PDF

    # plot SNP density across genome
    Rscript ${workflow.projectDir}/bin/plot_SNP_density.R ${workflow.projectDir}/CURRENT/VCF/${refFasta}.vcf

    # plot AFS (allele freq spectrum) - input VCF and number of samples
    # python3 ${workflow.projectDir}/bin/afs.py ${workflow.projectDir}/CURRENT/VCF/${refFasta}.vcf $number
    perl ${workflow.projectDir}/bin/afs.pl ${workflow.projectDir}/CURRENT/VCF/${refFasta}.vcf $number
    #  python3 ${workflow.projectDir}/bin/afs.lol.py ${workflow.projectDir}/CURRENT/VCF/${refFasta}.vcf $number
    """
}

process GETBASES {
    tag {"get bases"}
    label 'bases'

    input:
    val ready
    path (refFasta)

    output:
    publishDir ".", mode: "copy"

    script:
    """
    odgi flatten -i ${workflow.projectDir}/CURRENT/*.og -f ${workflow.projectDir}/CURRENT/${refFasta}.flatten.fa -b ${workflow.projectDir}/CURRENT/${refFasta}.bed -t 22
    perl ${workflow.projectDir}/bin/bases.pl ${refFasta}
    """
}

process VIZ2 {
	tag {"big viz"}
	label 'viz2'

	input:
	val ready
	path (refFasta)

	output:
	publishDir ".", mode: "copy"

	script:
	"""
        # y = height of plot
        # w = step size for blocks in plot
        # c = numbers of characters for sample names
        # x = width in pixels of output
        # y = height in pixels of output
        odgi viz -i ${workflow.projectDir}/CURRENT/*.og -o ${workflow.projectDir}/CURRENT/${refFasta}.viz.png -c 55 -w 50 -y 1500
        """
}

process HEAPS {
    tag {"heaps"}
    label 'heaps'

    input:
    val ready 

    output:
    publishDir ".", mode: "copy"

    exec:
    def file1 = file('bin/number1.txt')
    allLines = file1.readLines()
    for( line : allLines ) { number=line }

    script:
    """
    # visualise output, reading in heaps file, 3rd column only of interest
    Rscript ${workflow.projectDir}/bin/visualisation.R $number
    """
}

process PAVS {
    tag{"pavs"}
    label 'pavs'

    input:
    val ready
    path (refFasta)

    output:
    publishDir ".", mode: "copy"

    script:
    """
    odgi flatten -i ${workflow.projectDir}/CURRENT/*.og -b ${workflow.projectDir}/CURRENT/${refFasta}.flatten.tsv -t 22

    # convert to BED
    sed '1d' ${workflow.projectDir}/CURRENT/${refFasta}.flatten.tsv | awk -v OFS='\t' '{print(\$4,\$2,\$3,"step.rank_"\$6,".",\$5)}' > ${workflow.projectDir}/CURRENT/${refFasta}.flatten.bed

    # get BED file
    odgi pav -i ${workflow.projectDir}/CURRENT/*.og  -b ${workflow.projectDir}/CURRENT/${refFasta}.flatten.bed >       ${workflow.projectDir}/CURRENT/${refFasta}.flatten.pavs.tsv
    grep -c "" ${workflow.projectDir}/CURRENT/${refFasta}.flatten.pavs.tsv  > flatten.pavs.count.txt # gives all PAVs

    # make a plot
    Rscript ${workflow.projectDir}/bin/plot_pavs.R  ${workflow.projectDir}/CURRENT/${refFasta}.flatten.pavs.tsv

    rm -rf  ${workflow.projectDir}/CURRENT/${refFasta}.flatten.pavs.tsv  # needs to be removed, file size too large
    """
}

process WARAGRAPH {
    input:
    val ready
    path (refFasta)

    output:
    publishDir ".", mode: "copy"

    script:
    """
    # create layout for LSDV
    odgi layout -i ${workflow.projectDir}/CURRENT/*.og -t 12 -T  ${workflow.projectDir}/CURRENT/${refFasta}.layout

    # ensure BED file is present
    odgi flatten -i ${workflow.projectDir}/CURRENT/*.og -f ${workflow.projectDir}/CURRENT/${refFasta}.flatten.fa -b ${workflow.projectDir}/CURRENT/${refFasta}.bed -t 22

    # run waragraph on GFA with BED as tsv input
    ${workflow.projectDir}/bin/waragraph ${workflow.projectDir}/CURRENT/*.gfa ${workflow.projectDir}/CURRENT/${refFasta}.layout  --bed ${workflow.projectDir}/CURRENT/${refFasta}.bed    
    """
}

process BANDAGE {
    input:
    val ready
    path (refFasta)

    output:
    publishDir ".", mode: "copy"

    script:
    """
    # view in Bandage 
    Bandage load ${workflow.projectDir}/CURRENT/*.gfa --nodewidth=22
    """
}


process COMMUNITIES {
    input:
    val ready
    path (refFasta)

    output:
    publishDir ".", mode: "copy"

    script:
    """
    # use data in CURRENT/PANGROWTH/SEQS -> send to COMMUNITIES
    #     mkdir ${workflow.projectDir}/CURRENT/COMMUNITIES/

    # create genomes.fasta in COMMUNITIES using fastix
    bash ${workflow.projectDir}/bin/fastix.sh

    rm -rf ${workflow.projectDir}/CURRENT/COMMUNITIES/genomes.fasta.gz
    bgzip -@ 4 ${workflow.projectDir}/CURRENT/COMMUNITIES/genomes.fasta
    samtools faidx ${workflow.projectDir}/CURRENT/COMMUNITIES/genomes.fasta.gz # index

    # Community detection based on p/w alignments based on 90% ID at mash level with 6
    # mappings per segment, using k-mer of 19 and window of 67
    wfmash ${workflow.projectDir}/CURRENT/COMMUNITIES/genomes.fasta.gz -p 90 -n 6 -t 44 -m > ${workflow.projectDir}/CURRENT/COMMUNITIES/genomes.mapping.paf

    # Convert PAF mappings into a network:
    python3 ${workflow.projectDir}/bin/paf2net.py -p ${workflow.projectDir}/CURRENT/COMMUNITIES/genomes.mapping.paf

    # make plot
    python3 ${workflow.projectDir}/bin/net2communities.py -e ${workflow.projectDir}/CURRENT/COMMUNITIES/genomes.mapping.paf.edges.list.txt -w ${workflow.projectDir}/CURRENT/COMMUNITIES/genomes.mapping.paf.edges.weights.txt -n ${workflow.projectDir}/CURRENT/COMMUNITIES/genomes.mapping.paf.vertices.id2name.txt --plot

    # mash-based partitioning
    mash dist ${workflow.projectDir}/CURRENT/COMMUNITIES/genomes.fasta.gz ${workflow.projectDir}/CURRENT/COMMUNITIES/genomes.fasta.gz -s 10000 -i > ${workflow.projectDir}/CURRENT/COMMUNITIES/genomes.distances.tsv

    # get distances
    python3 ${workflow.projectDir}/bin/mash2net.py -m ${workflow.projectDir}/CURRENT/COMMUNITIES/genomes.distances.tsv

    # get numbers in each group
    bash ${workflow.projectDir}/bin/community_numbers.sh

    ${workflow.projectDir}/bin/pafgnostic --paf ${workflow.projectDir}/CURRENT/COMMUNITIES/genomes.mapping.paf > ${workflow.projectDir}/CURRENT/COMMUNITIES/genomes.mapping.paf.txt 
    #  pafgnostic outputs a tab-separated table directly to the console with the following columns:
    # - Query name, start, end, strand
    # - Target name, start, end
    # - Number of matches, mismatches, insertions, deletions
    # - Unique matches, mismatches, insertions, deletions
    # - Longest, shortest, and average lengths of matches, mismatches, insertions, deletions
    """
}

process BUSCO { 
    input:
    val ready
    path (refFasta)

    output:
    publishDir ".", mode: "copy"

    script:
    """
    ${workflow.projectDir}/bin/busco -f -i ${refFasta} -l poxviridae_odb10 -o ${workflow.projectDir}/CURRENT/BUSCO -m genome
    # identifies genes using the poxvirus annotation for busco v5 process annotate {
    """
}

process PANAROO {
    input:
    val ready
    path (refFasta)

    output:
    publishDir ".", mode: "copy"

    script:
    """
    # use data in PANGROWTH/SEQS/
    # Rscript ${workflow.projectDir}/bin/annotate.R lsdv
    # Rscript ${workflow.projectDir}/bin/annotate.R sppv
    Rscript ${workflow.projectDir}/bin/annotate.R gpv

    # mamba update panaroo
    ${workflow.projectDir}/bin/panaroo -i ${workflow.projectDir}/CURRENT/PANGROWTH/SEQS/*PROKKA/*.gff -o ${workflow.projectDir}/CURRENT/PANAROO --clean-mode strict 

    # generate simple plot
    Rscript ${workflow.projectDir}/bin/view_gml.R
    
    # visualise presence absence
    Rscript ${workflow.projectDir}/bin/panaroo_viz.R
    """
}

process RECOMBINATION {
    input:
    val ready
    path (refFasta)

    output:
    publishDir ".", mode: "copy"

    script:
    """
    # TO DO
    ./3seq -c my3seqTable2000 # check 3seq works
    ./3seq -i   v5_full_genomes.aln -quiet # check file is ok
    ./3seq  -f  v5_full_genomes.aln   -id R1  -f2500 -l148000 -bp-all &> all.out  
    ./3seq  -f  v5_full_genomes.aln   -id R1  -f2500 -l148000 -bp-all &> all.out 
    """
}