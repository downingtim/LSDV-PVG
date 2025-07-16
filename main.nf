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
Modules
#==============================================
*/
include { DOWNLOAD; ALIGN; TREE; MAKE_PVG; VIZ1; ODGI; OPENNESS_PANACUS; OPENNESS_PANGROWTH; PATH_FROM_GFA;VCF_FROM_GFA;VCF_PROCESS; GETBASES; VIZ2; HEAPS; HEAPS_Visualize; PAVS; PAVS_plot; WARAGRAPH; COMMUNITIES; PANAROO; BUSCO; PAFGNOSTIC; Bandage;BANDAGE_view; GFAstat; SUMMARIZE; Extract_Ref; PROKKA; Clean_GTF; Annotate_Position  } from './modules/processes.nf'
/*
#==============================================
Modules
#==============================================
*/

workflow Main {
    // Define a single source for fasta input
    fastaCh = params.reference 
        ? channel.fromPath(params.reference, checkIfExists: true) 
        : DOWNLOAD().fasta_file
    
    // Debug print to console
    fastaCh.view { "Debug - Input file: $it" }
    
    // Proceed with downstream processes
    if (params.align) {
        alignOut = ALIGN(fastaCh)
        if (params.tree) {
            TREE(alignOut.raxml_file)
        }
    }
    
    if (params.make_pvg) {
        pvgOut = MAKE_PVG(fastaCh)
        
        if (params.viz1) {
            VIZ1(pvgOut.gfa)
        }
        
        if (params.bandage) {
            Bandage(pvgOut.gfa)
        }
        
        if (params.odgi) {
            odgiOut = ODGI(pvgOut.gfa)

	    if (params.annotate){
				Annotate(odgiOut.refid,odgiOut.ogfile,fastaCh)
                //annotation_reference = Extract_Ref(odgiOut.refid,fastaCh)
                //prokka_out = PROKKA(odgiOut.refid,annotation_reference)
                //cleaned_gtf = Clean_GTF (odgiOut.refid,prokka_out.prokkagff,odgiOut.ogfile) 
				//Annotate_Position(odgiOut.ogfile,cleaned_gtf)
	     }

	    if (params.viz2) {
                VIZ2(odgiOut.ogfile)
            }
            
            if (params.heaps) {
                HEAPS(odgiOut.ogfile) | HEAPS_Visualize
            }
            
            if (params.pavs) {
                PAVS(odgiOut.ogfile) | PAVS_plot
            }
            
            if (params.getbases) {
                GETBASES(odgiOut.ogfile, fastaCh)
            }
        }
        
        if (params.gfastat) {
            GFAstat(pvgOut.gfa, fastaCh)
        }
        
        if (params.getvcf) {
            GetVCF(pvgOut.gfa)
        }
        
        if (params.openness_panacus) {
            OPENNESS_PANACUS(pvgOut.gfa)
        }
    }
    
    if (params.openness_pangrowth) {
        pangrowthOut = OPENNESS_PANGROWTH(fastaCh)
        
        if (params.communities) {
            communitiesOut = COMMUNITIES(pangrowthOut.communities_genome)
            PAFGNOSTIC(communitiesOut.paf_file)
        }
    }
    
    if (params.busco && params.busco_clade) {
        BUSCO(fastaCh)
    }
}

workflow GetVCF(){
    take:
    gfa

    main:
    gfapath = PATH_FROM_GFA(gfa)
    gfavariants = VCF_FROM_GFA(gfa) 
    VCF_PROCESS(gfapath,gfavariants,gfa)
}

workflow Annotate(){
    take:
		refid
		ogfile
		fastaCh

	main: 
		annotation_reference = Extract_Ref(refid,fastaCh)
		prokka_out = PROKKA(refid,annotation_reference)
		cleaned_gtf = Clean_GTF (refid,prokka_out.prokkagff,ogfile) 
		Annotate_Position(ogfile,cleaned_gtf)
}

workflow Summary{
    template_ch = Channel.fromPath("template.tex") 
    SUMMARIZE(template_ch)
}

workflow {
    if(params.summary) {Summary()}
    else {Main()}
}

