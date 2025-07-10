process ANNOTATE {
    tag {"annotate"}
    label 'annotate'
    container "pangenome/odgi:1726671973"
    containerOptions = '--entrypoint ""'
    
    input:
    path ogfile
    path refFasta
    
    output:
    path "PVG_annotation.csv"
    
    publishDir "results/bandage", mode: "copy"
    
    script:
    """
    # need to ensure samtools, prokka and gffread are in the PATH
    REFERENCE_ID=\$(odgi paths -i ${ogfile} -L | head -n 1)
    samtools faidx ${refFasta} \${REFERENCE_ID} > \${REFERENCE_ID}.fasta
    REFERENCE_FASTA="\${REFERENCE_ID}.fasta"
    
    prokka --kingdom Viruses --gffver 3 --usegenus --outdir PROKKA --prefix \${REFERENCE_ID} \${REFERENCE_FASTA} --force --compliant
    
    INPUT_GFF="PROKKA/\${REFERENCE_ID}.gff"
    CLEAN_GFF="\${REFERENCE_ID}_clean.gff"
    ANNOTATION_CSV="PVG_annotation.csv"
    
    # Clean up sequence region naming in GFF
    string1=\$(grep "sequence-region" "\$INPUT_GFF" | awk '{print \$2}')
    sed "s/\${string1}/\${REFERENCE_ID}/g" "\$INPUT_GFF" | sed 's/^>//' > "\$CLEAN_GFF"
    
    # Convert GFF to GTF format
    CLEAN_GTF="\${REFERENCE_ID}_clean.gtf"
    gffread -E "\$CLEAN_GFF" -T -o "\$CLEAN_GTF"
    
    # Process GTF to create simplified annotation format
    PREP_GTF="\${REFERENCE_ID}.prep.gtf"
    
    # Extract transcript lines and clean up formatting
    grep -P "transcript\\t" "\$CLEAN_GTF" | \\
    sed 's/gene_id //g; s/"//g; s/prokka\ttranscript//g; s/_gene; gene_name//g; s/transcript_id//g; s/KCBEBCIG_0//g; s/_gene//g' | \\
    sort | uniq | \\
    awk '{
        # Print: chr, A, B, start, end, strand, score, frame, gene_name
        printf "%s\tA\tB\t%s\t%s\t%s\t%s\t%s\t", \$1, \$2, \$3, \$4, \$5, \$6
        if (\$9) {
            print \$9
        } else {
            print "gene_" \$8
        }
    }' > "\$PREP_GTF"
    
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
