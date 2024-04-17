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

include { DOWNLOAD; MAKE_PVG; ODGI; OPENNESS_PANACUS; OPENNESS_PANGROWTH; GET_VCF; GETBASES; VIZ2; HEAPS; PAVS; WARAGRAPH; COMMUNITIES; PANAROO; BUSCO; BANDAGE } from './modules/processes.nf'

/*
#==============================================
Modules
#==============================================
*/

workflow { 
   faS = channel.fromPath("bin/faSplit", checkIfExists:true)
//   refFasta = channel.fromPath("goatpox_2_virus_2_2_.fasta", checkIfExists:true)

       DOWNLOAD(faS) 
       DOWNLOAD
	  .out
	  .write
	  .set { refFasta }

       MAKE_PVG( refFasta )

       COMMUNITIES(MAKE_PVG.out, refFasta)

       ODGI(MAKE_PVG.out, refFasta)
       
       OPENNESS_PANACUS(MAKE_PVG.out)

       OPENNESS_PANGROWTH(MAKE_PVG.out, refFasta)

       GETBASES(MAKE_PVG.out, refFasta)
    
       VIZ2(MAKE_PVG.out, refFasta )
       
       HEAPS(MAKE_PVG.out )

       PAVS(MAKE_PVG.out, refFasta)

       GET_VCF(MAKE_PVG.out, refFasta)

       PANAROO(MAKE_PVG.out, refFasta)

       BUSCO(MAKE_PVG.out, refFasta)

       waragraph()

}

workflow waragraph { 
//   refFasta = channel.fromPath("goatpox_2_virus_2_2_.fasta", checkIfExists:true)

//      WARAGRAPH( refFasta)
//	BANDAGE( refFasta )

// add Bandage
}
