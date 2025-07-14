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
    container "chandanatpi/panalayze_env:3.0"
    cpus 8

    output:
    publishDir "results/download", mode: "copy"
    path "genomes.fasta", emit: fasta_file

    shell:
    q1 = params.virus_name + " [organism] AND complete genome [title]"  
    q2 = params.virus_name + " [organism] AND genomic sequence [title]" 
    filter_text = params.virus_filter.replaceAll(',', "\\\\|");filter = /complete genome\|${filter_text}\|isolate\|genomic sequence\|_NULL\|strain/
    '''
    esearch -db nucleotide -query "!{q1}"|efetch -format genbank > gb.1 
    esearch -db nucleotide -query "!{q2}"|efetch -format genbank > gb.2
    cat gb.1 gb.2 > genbank.gb 
    rm gb.1 gb.2
    #parseGB.py genbank.gb|sed -e "s%!{filter}%%g"|tr -d "',)(:;\\\""|sed -e "s%/\\| \\|-%_%g" |sed -e "s%[_]+%_%g"|tr -s _|sed -e "s/_$//"  > genomes.fasta
    parseGB.py genbank.gb > genomes.fasta
    ls -lt genomes.fasta genbank.gb
    echo "File created with size: $(ls -la genomes.fasta)"
    '''
}

process ALIGN {
    cpus 8
    
    input:
    path genomes

    output:
    publishDir "results/align", mode: "copy"
    path "genomes.aln"
    path "T14.*"
    path "T14.raxml.supportTBE",emit:raxml_file

    shell:
    """
        
    mafft  --thread !{task.cpus} --auto ${genomes} > genomes.aln 
    raxml-ng  --all --msa genomes.aln --model GTR+G4 --prefix T14 --seed 21231 --bs-metric fbp,tbe --redo
        
    """
}

process TREE{
   container "chandanatpi/panalayze_env:3.0"
   tag "Building tree with $y"
   input:
   path y   
   
   output:
   publishDir "results/tree", mode: "copy"
   path "tree.png"

   shell:
   """
       tree.R || true
   """
}

process MAKE_PVG {
    container 'chandanatpi/tpi:pggb'
    cache 'lenient'
    tag {"index reference FASTA"}
    label 'pvg'
    cpus 20 

    input:
    path (refFasta)
    
    output:
    path ("pggb.gfa"), emit: gfa
    publishDir "results/PVG", mode: "copy"

    script:
    """
    # create Blast indexes for self-self blast
    # makeblastdb -dbtype nucl -in ${refFasta} -out ${workflow.projectDir}/BLAST/fasta.db 
    # blastall -p blastn -d ${workflow.projectDir}/BLAST/fasta.db -i ${refFasta} -e 1.0 -outfmt 6 > self-blast.out
    
    # SAMtools index
    bgzip ${refFasta}
    samtools faidx ${refFasta}.gz > ${refFasta}.gz.fai

    # run PGGB - you need to specify the number of haplotypes in number
    pggb -i ${refFasta}.gz -m -S -o . -t ${task.cpus} -p 90 -s 1k -n ${params.haplotypes}
    mv *.gfa pggb.gfa
    # if issues, remove CURRENT/ in ${workflow.projectDir}  before re-running 

    # run pangenome in nf-core - don't work as reliably as PGGB, so best avoided
    # nextflow run nf-core/pangenome --input ${refFasta}.gz --n_haplotypes number --outdir \
    #   ${workflow.projectDir}/CURRENT  --communities  --wfmash_segment_length 1000
    #odgi stats -m -i ${workflow.projectDir}/CURRENT/FINAL_ODGI/${refFasta}.gz.squeeze.og -S > ${workflow.projectDir}/odgi.stats.txt 

    """
}

process VIZ1 {
    container "chandanatpi/panalayze_env:3.0"
    input:
    path gfa 

    output:
    path "out.vg.png"
    publishDir "results/vg", mode: "copy"

    script:
    """
    # visualise gfa - need to separate, this takes a while
    vg view -F -p -d ${gfa}  | dot -Tpng -o out.vg.png
    """
}

process GFAstat {
    container "chandanatpi/panalayze_env:3.0"
    input:
    path gfa 
    path (refFasta)

    output:
    path "gfa.stats.txt"
    path "genome.lengths.txt"
    publishDir "results/gfastat", mode: "copy"

    script:
    """
    # get general info on gfa
    gfastats ${gfa} > gfa.stats.txt

    # get lengths of genomes etc
     lengths.py ${refFasta} 

    """
}

process Bandage {
    container "biocontainers/bandage:v0.8.1-1-deb_cv1"
    containerOptions = "--user root"
    input:
    path gfa 

    output:
    path "bandage.stats.txt"
    publishDir "results/bandage", mode: "copy"

    script:
    """
    # get Bandage stats
    Bandage info ${gfa} > bandage.stats.txt
    """
}

process ODGI {
    container "pangenome/odgi:1726671973"
    containerOptions = '--entrypoint ""'

    cpus 1

    input:
    path gfa

    output:
    path "out.og", emit: ogfile
    path "odgi.stats.txt"
	stdout emit: refid

    publishDir "results/odgi", mode: "copy"

    script:
    """
    odgi build -g ${gfa} -o out.og
    odgi stats -m -i out.og -S > odgi.stats.txt
    odgi paths -i out.og -L | head -n 1|tr -d '\n' 
    """
}

process OPENNESS_PANACUS {
    tag {"get PVG openness panacus"}
    label 'openness_panacus'
    container "mgibio/panacus:0.2.4"

    input:
    path gfa 

    output:
    path "haplotypes.txt"
    path "histgrowth.node.tsv"
    path "histgrowth.node.pdf"
    publishDir "results/panacus", mode: "copy"
    
    script:
    """
    # you need to install panancus
    # mamba install -c conda-forge -c bioconda panacus
    # get haplotypes
    grep '^P' ${gfa} | cut -f2 > haplotypes.txt
    
    # run panacus to get data - just use defaults
    RUST_LOG=info panacus histgrowth -t4 -l 1,2,1,1,1 -q 0,0,1,0.5,0.1 -S -s haplotypes.txt ${gfa} > histgrowth.node.tsv
 
    # visualise plot as PDF
    panacus-visualize -e histgrowth.node.tsv > histgrowth.node.pdf
    """
}

process OPENNESS_PANGROWTH {
    container "chandanatpi/panalayze_env:3.0"
    tag {"get PVG openness pangrowth"}
    label 'openness_pangrowth'

    input:
    path (refFasta)

    output:
    path "p_core"
    path "p_core.pdf"
    path "communities.genomes.fasta",emit:communities_genome
    path "pangrowth.pdf"
    path "growth.pdf"    
    publishDir "results/pangrowth", mode: "copy"

    script:
    """
    #Splir refFasta by name i to different files
    mkdir SEQS -p
    awk '{
        if (substr(\$1, 1, 1) == ">") {
            NAME = substr(\$1, 2);
            split(NAME, arr, " ");
            NAME = arr[1];
            print ">" NAME > "SEQS/" NAME ".fa"
        } else {
            if (\$1) {
                print \$1 >> "SEQS/" NAME ".fa"
            }
        }
    }' ${refFasta}

    for FILE in SEQS/*.fa 
    do
      # Extract sample name from file name
      sample_name=\$(basename "\$FILE" .fasta)
      # Execute ~/bin/fastix with sample name and append output
      fastix -p "\${sample_name}#1#" "\$FILE" >> communities.genomes.fasta
    done

    # run pangrowth
    pangrowth hist -k 17 -t 12 SEQS/*a > hist.txt

    # plot AFS for single dataset - not very useful! 
    plot_hist.py hist.txt pangrowth.pdf

    # model growth - prepare dataset 1 - not useful
    pangrowth growth -h hist.txt > temp1
    plot_growth.py temp1 growth.pdf

    # do core - this does not converge to a solution for n<10 samples, causes error
    pangrowth core -h hist.txt -q 0.5 > p_core

    # if [ ${params.haplotypes} -gt 8 ]; then # if the core is too small, this crashes
    plot_core.py p_core p_core.pdf     # so need at least 10 genomes
    # fi
    """
}

process PATH_FROM_GFA {
    tag {"get paths from GFA"}
    label 'vcf'
    container "pangenome/odgi:1726671973"
    containerOptions = '--entrypoint ""'

    input:
    path gfa

    output:
    path "path1"

    script:
    """
    # need to get paths from GFA
    odgi paths -i ${gfa} -L |head -n 1 > path1
    
    """
}

process VCF_FROM_GFA {
    container "chandanatpi/panalayze_env:3.0"

    tag {"get VCFs from GFA"}
    label 'vcf'

    input:
    path gfa

    output:
    path "gfavariants.vcf",emit:gfavariants 
    publishDir "results/vcf", mode: "copy"

    script:
    """
    gfautil -i ${gfa} gfa2vcf > gfavariants.vcf
    """
}

process VCF_PROCESS {
    container "chandanatpi/panalayze_env:3.0"
    tag {"Process VCF information"}
    label 'vcf'

    input:
    path pathinfo
    path vcf_file
    path gfa_file

    output:
    path "mutation_density.pdf"
    path "frequency_distribution.png"
    publishDir "results/vcf", mode: "copy"

    script:
    """
    # script to create initial input to visualise with R
    Ref=\$(cat ${pathinfo})
    for ((N=300; N<=${params.genome_length}; N+=2)) 
    # genome ends skipped, even numbers only
    do
        gfautil --quiet -t 35 -i ${gfa_file} snps --ref \$Ref  --snps \$N
    done| sort -nk 3 | grep -v \\: |grep -v path > variation_map.txt

    plot_variation_map.R || true # plot image of variation map in PDF

    # plot SNP density across genome
    plot_SNP_density.R ${vcf_file} || true

    # plot AFS (allele freq spectrum) - input VCF and number of samples
    afs.py out ${params.haplotypes} 
    """
}

process GETBASES {
    tag {"get bases"}
    label 'bases'
    container "pangenome/odgi:1726671973"
    containerOptions = '--entrypoint ""'

    input:
    path ogfile
    path (refFasta)

    output:
    path "out.flatten.fa"
    path "out.bed"
    publishDir "results/getbases", mode: "copy"

    script:
    """
    odgi flatten -i ${ogfile} -f out.flatten.fa -b out.bed -t 22
    #bases.pl ${refFasta}
    """
}

process VIZ2 {
    tag {"big viz"}
    label 'viz2'
    container "pangenome/odgi:1726671973"
    containerOptions = '--entrypoint ""'

    input:
    path ogfile

    output:
    path "out.viz.png"
    publishDir "results/viz2", mode: "copy"

    script:
    """
        # y = height of plot
        # w = step size for blocks in plot
        # c = numbers of characters for sample names
        # x = width in pixels of output
        # y = height in pixels of output
        odgi viz -i ${ogfile} -o out.viz.png -c 55 -w 50 -y 1500
    """
}

process HEAPS {
    tag {"heaps"}
    label 'heaps'
    container "pangenome/odgi:1726671973"
    containerOptions = '--entrypoint ""'

    input:
    path ogfile 

    output:
    path "heaps.txt",emit:heap_file
    publishDir "results/heaps", mode: "copy"

    script:
    """
    # visualise output, reading in heaps file, 3rd column only of interest
    for N in {1..${params.haplotypes}}
    do
       odgi heaps -i ${ogfile} -n 1000 -d \$N -t 12 | sort -nk 2 | tail -n 1
    done > heaps.txt
    """
}

process HEAPS_Visualize {
    container "chandanatpi/panalayze_env:3.0"
    tag {"heaps"}
    label 'heaps'

    input:
    path ogfile 

    output:
    path "heaps.pdf"
    publishDir "results/heaps", mode: "copy"

    script:
    """
    visualisation.R ${params.haplotypes} || true
    """
}

process PAVS {
    tag{"pavs"}
    label 'pavs'
    container "pangenome/odgi:1726671973"
    containerOptions = '--entrypoint ""'

    input:
    path ogfile

    output:
    path "out.flatten.pavs.tsv" 

    publishDir "results/pavs", mode: "copy", pattern: "flatten.pavs.count.txt"

    script:
    """
    odgi flatten -i ${ogfile} -b out.flatten.tsv -t 22

    # convert to BED
    sed '1d' out.flatten.tsv | awk -v OFS='\t' '{print(\$4,\$2,\$3,"step.rank_"\$6,".",\$5)}' > out.flatten.bed

    # get BED file
    odgi pav -i ${ogfile}  -b out.flatten.bed > out.flatten.pavs.tsv
    grep -c "" out.flatten.pavs.tsv  > flatten.pavs.count.txt # gives all PAVs

    # make a plot

    """
}

process PAVS_plot {
    container "chandanatpi/panalayze_env:3.0"
    tag{"pavs"}
    label 'pavs'

    input:
    path tsv_file

    output:
    path "out.flatten.pavs.pdf"
    publishDir "results/pavs", mode: "copy"

    script:
    """
    plot_pavs.R  ${tsv_file}
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

process BANDAGE_view {

    executor 'local' 
    cpus 1

    container "biocontainers/bandage:v0.8.1-1-deb_cv1"
    containerOptions = "--user root --env DISPLAY=$DISPLAY"
    input:
    path gfa 


    script:
    """
    # view in Bandage 
    Bandage load ${gfa} --nodewidth=22
    """
}


process COMMUNITIES {
    container 'chandanatpi/tpi:pggb'
    input:
    path communities_genome
    

    output:
    path "genomes.mapping.paf*"
    path "genomes.distances.tsv"
    path "genomes.mapping.paf",emit:paf_file  
    publishDir "results/communities", mode: "copy", pattern: "*"  

    script:
    """

    bgzip -@ 4 ${communities_genome}
    samtools faidx ${communities_genome}.gz # index

    # Community detection based on p/w alignments based on 90% ID at mash level with 6
    # mappings per segment, using k-mer of 19 and window of 67
    wfmash ${communities_genome}.gz -p 90 -n 6 -t 44 -m > genomes.mapping.paf

    # Convert PAF mappings into a network:
    paf2net.py -p genomes.mapping.paf

    # make plot
    net2communities.py -e genomes.mapping.paf.edges.list.txt -w genomes.mapping.paf.edges.weights.txt -n genomes.mapping.paf.vertices.id2name.txt --plot

    # mash-based partitioning
    mash dist ${communities_genome}.gz ${communities_genome}.gz -s 10000 -i > genomes.distances.tsv

    # get distances
    mash2net.py -m genomes.distances.tsv

    output_file="genomes.communities.mash.paths.txt"

    # Initialize a flag to check if any community files were found
    found_community_files=0

    for i in {0..3}; do
    # Check if the community file exists
        if [[ -f "genomes.mapping.paf.edges.weights.txt.community.\${i}.txt" ]]; then
    # Extract unique chromosome names from the community file
            chromosomes=\$(cut -f 3 -d '#' "genomes.mapping.paf.edges.weights.txt.community.\${i}.txt" | sort | uniq | tr '\n' ' ')
            echo "community \${i} --> \${chromosomes}"
            found_community_files=1
        else
            echo "community \${i} --> Community file not found"
        fi
    done

    # Check if any community files were found
    if [[ \$found_community_files -eq 0 ]]; then
        echo "No community files found. Exiting."
    fi  
    # Redirect the output to the specified file
    echo "Output written to \$output_file"
    # community_numbers.sh end

    pafgnostic --paf genomes.mapping.paf > genomes.mapping.paf.txt 
    #  pafgnostic outputs a tab-separated table directly to the console with the following columns:
    # - Query name, start, end, strand
    # - Target name, start, end
    # - Number of matches, mismatches, insertions, deletions
    # - Unique matches, mismatches, insertions, deletions
    # - Longest, shortest, and average lengths of matches, mismatches, insertions, deletions
    """
}

process BUSCO { 

    container "ezlabgva/busco:v5.8.0_cv1"
    containerOptions = "--user root"
    input:
    path (refFasta)

    output:
    publishDir "results/busco", mode: "copy"

    script:
    """
    busco -f -i ${refFasta} -l ${params.busco_clade} -o BUSCO -m genome
    # identifies genes using the poxvirus annotation for busco v5 process annotate {
    """
}

process PAFGNOSTIC {
    container 'chandanatpi/tpi:pggb'
    input:
    path paf_file 

    output:
    publishDir "results/pafgnostic", mode: "copy"

    script:
    """
    pafgnostic --paf ${paf_file} > pafgnostic.txt  
    """
}

process PANAROO {
    container "chandanatpi/panalayze_env:3.0"
    input:
    path (refFasta)

    output:
    publishDir "results/panaroo", mode: "copy"

    script:
    """
    # use data in PANGROWTH/SEQS/
    Rscript annotate.R gpv

    # mamba update panaroo
    ${workflow.projectDir}/bin/panaroo -i ${workflow.projectDir}/CURRENT/PANGROWTH/SEQS/*PROKKA/*.gff -o ${workflow.projectDir}/CURRENT/PANAROO --clean-mode strict 

    # generate simple plot
    Rscript ${workflow.projectDir}/bin/view_gml.R || true
    
    # visualise presence absence
    Rscript ${workflow.projectDir}/bin/panaroo_viz.R || true
    """
}

process SUMMARIZE {
    container "chandanatpi/panalayze_env:3.0"
    input:
    path template_tex 
    output:
    path "template.pdf"

    publishDir "results/summary", mode: "copy"

    script:
    """
    pdflatex --interaction=nonstopmode ${template_tex}
    """
}

process Extract_Ref{
    tag {"Extract reference"}
    label "Extract_reference"
    container "chandanatpi/panalayze_env:3.0" 
    publishDir "results/annotate", mode: "copy"

    input:
    val refid
    path refFasta
    
    output:
    path "reference.annotate.fasta"

    script:
    """
    samtools faidx ${refFasta} ${refid} > reference.annotate.fasta
    """
}

process PROKKA{
    tag {"Prokka"}
    label 'Prokka'
    container "staphb/prokka:latest"
    publishDir "results/prokka", mode: "copy"

    input:
    val refid
    path annotate_ref_fasta

    output:
    path "PROKKA"
    path "annotation.gff",emit:prokkagff


    script:
    """
    prokka --kingdom Viruses --gffver 3 --usegenus --outdir PROKKA --prefix ${refid} ${annotate_ref_fasta} --force --compliant
    #REFERENCE_ID=${refid}
    #prokka --kingdom Viruses --gffver 3 --usegenus --outdir PROKKA --prefix \${REFERENCE_ID} ${annotate_ref_fasta} --force --compliant
    cp "PROKKA/${refid}.gff" annotation.gff
    """
}

process ANNOTATE {
    tag {"annotate"}
    label "Extract_reference"
    container "pangenome/odgi:1726671973"
    containerOptions = '--entrypoint ""'
    publishDir "results/bandage", mode: "copy"

    input:
    val refid
    path prokkagff
    path ogfile

    output:
    path "PVG_annotation.csv"


    script:
    """
    INPUT_GFF="${prokkagff}"
    ANNOTATION_CSV="PVG_annotation.csv"

    string1=\$(grep "sequence-region" "\$INPUT_GFF" | awk '{print \$2}')
    sed "s/\${string1}/${refid}/g" "\$INPUT_GFF" | sed 's/^>//' > clean.gff 

    gffread -E  clean.gff -T -o clean.gtf 

    PREP_GTF="${refid}.prep.gtf" # intermediate
    # Process GTF to extract gene names or simplified gene numbers
    grep -P "transcript\\t"  clean.gtf | while IFS=\$'\\t' read -r col1 col2 col3 col4 col5 col6 col7 col8 col9; do
        # Extract gene_name if it exists
        gene_name=\$(echo "\$col9" | grep -oP 'gene_name "\\K[^"]+')

        if [ -n "\$gene_name" ]; then
            # Use gene_name if available
            annotation="\$gene_name"
        else
            # Extract gene_id and convert to simplified format
            gene_id=\$(echo "\$col9" | grep -oP 'gene_id "\\K[^"]+')
            # Extract the number part after the last underscore
            gene_num=\$(echo "\$gene_id" | sed 's/.*_//')
            annotation="gene_\$gene_num"
        fi

        # Output in GTF format with simplified annotation
        echo -e "\$col1\\t\$col2\\t\$col3\\t\$col4\\t\$col5\\t\$col6\\t\$col7\\t\$col8\\ttranscript_id \\"\$annotation\\"; gene_id \\"\$annotation\\""
    done > "\$PREP_GTF"

    # Map positions to PVG nodes
    odgi position -i ${ogfile} -E "\$PREP_GTF" > "\$ANNOTATION_CSV"

    echo "Output a PVG annotation CSV: \$ANNOTATION_CSV with"
    echo " \$(wc -l < "\$ANNOTATION_CSV") annotation entries"
    echo
    echo "To view the annotated PVG in Bandage:"
    echo "1. Launch Bandage and load the PVG in GFA format"
    echo "2. Draw the PVG in Bandage"
    echo "3. Load the annotation CSV: \$ANNOTATION_CSV"
    """
}
