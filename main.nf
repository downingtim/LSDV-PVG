#!/usr/bin/env nextflow

/*
Optionally: Download datasets of interest, perform QC and construct a phylogeny to select suitable genomes
Main components:
[1] MAKE_PVG: nstruct a PVG using PGGB
[2] VIZ1: 
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

    refFasta = channel.fromPath(params.reference, checkIfExists:true)
    if (params.genomes)
    {
        Input_Genome = channel.fromPath(params.genomes, checkIfExists:true)
    }
    else 
    {
	    Input_Genome = DOWNLOAD() 
    }
    ALIGN(Input_Genome.first())
    TREE(ALIGN.out.raxml_file) 
    PVG_out = MAKE_PVG( refFasta )
    VIZ1(PVG_out.gfa)
    GFAstat(PVG_out.gfa,refFasta)
    Bandage(PVG_out.gfa)
	ODGI_out = ODGI(PVG_out.gfa)
    OPENNESS_PANACUS(PVG_out.gfa)
    GETBASES(ODGI_out.ogfile, refFasta)
    VIZ2(ODGI_out.ogfile)
    HEAPS(ODGI_out.ogfile)|HEAPS_Visualize
    PAVS(ODGI_out.ogfile)|PAVS_plot
    GetVCF(PVG_out.gfa)
    if (params.busco_clade)
    {
        BUSCO(refFasta)
    }
    Pangrowth_Out = OPENNESS_PANGROWTH(refFasta)
    COMMUNITIES(Pangrowth_Out.communities_genome)|PAFGNOSTIC 
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

