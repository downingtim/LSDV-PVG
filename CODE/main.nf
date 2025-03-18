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
include { DOWNLOAD; ALIGN; TREE; MAKE_PVG; VIZ1; ODGI; OPENNESS_PANACUS; OPENNESS_PANGROWTH; PATH_FROM_GFA;VCF_FROM_GFA;VCF_PROCESS; GETBASES; VIZ2; HEAPS; HEAPS_Visualize; PAVS; PAVS_plot; WARAGRAPH; COMMUNITIES; PANAROO; BUSCO; PAFGNOSTIC; Bandage;BANDAGE_view; GFAstat; } from './modules/processes.nf'
/*
#==============================================
Modules
#==============================================
*/


workflow { 

    if (params.reference)
    {
        refFasta = channel.fromPath(params.reference, checkIfExists:true)
    }
    else 
    {
	    refFasta = DOWNLOAD() 
    }
    if (params.align) ALIGN(refFasta)
    if (params.tree) TREE(ALIGN.out.raxml_file) 
    if (params.make_pvg) PVG_out = MAKE_PVG( refFasta )
    if (params.viz1) VIZ1(PVG_out.gfa)
    if (params.gfastat) GFAstat(PVG_out.gfa,refFasta)
    if (params.bandage) Bandage(PVG_out.gfa)
	if (params.odgi) ODGI_out = ODGI(PVG_out.gfa)
    if (params.openness_panacus) OPENNESS_PANACUS(PVG_out.gfa)
    if (params.getbases) GETBASES(ODGI_out.ogfile, refFasta)
    if (params.viz2) VIZ2(ODGI_out.ogfile)
    if (params.heaps) {HEAPS(ODGI_out.ogfile)|HEAPS_Visualize}
    if (params.pavs) {PAVS(ODGI_out.ogfile)|PAVS_plot}
    if (params.getvcf) GetVCF(PVG_out.gfa)
    if (params.busco && params.busco_clade)
    {
        BUSCO(refFasta)
    }
    if (params.openness_pangrowth) Pangrowth_Out = OPENNESS_PANGROWTH(refFasta)
    if (params.communities) COMMUNITIES(Pangrowth_Out.communities_genome)|PAFGNOSTIC 
    //BANDAGE_view(PVG_out.gfa)
    //PANAROO(MAKE_PVG.out, refFasta)
}

workflow GetVCF(){
    take:
    gfa

    main:
    gfapath = PATH_FROM_GFA(gfa)
    gfavariants = VCF_FROM_GFA(gfa) 
    VCF_PROCESS(gfapath,gfavariants,gfa)
}


workflow waragraph { 
//   refFasta = channel.fromPath("goatpox_2_virus_2_2_.fasta", checkIfExists:true)
//   refFasta = channel.fromPath("lumpy_skin_disease_2_virus_2_2_.fasta", checkIfExists:true)

//      WARAGRAPH( refFasta)
//	BANDAGE( refFasta )

// add Bandage
}

