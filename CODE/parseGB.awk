#!/bin/gawk -f
# Parse Genbank files to a tab separated table (Version 3)
# it only read several fields required by my own investigation, modify at your needs
# columns, are: 1 Accession number, 2 definition, 3 source start and end, 4 isolate, 5 host, 7 clone, 8 country, 9 collection date, 10 organism, 11 serotype, 12 strain
# as you can see is virus/patogen oriented
# it allows multiline definitions and multiple sources in one sequence.
BEGIN { 
        #Populate the columns headers and clear some variables
        ACC="Accession"; Def="Definition"; src_pos[0]="nt_pos"; hst[0]="Host"; PMID="Pubmed_ID"; is[0]="Isolate"; cl[0]="Clone"; co[0]="Country"; org[0]="Organism"; sty[0]="Serotype"; cd[0]="Collection_date"; sr[0]="Strain"; cb[0]="Collected by";op=0; src=0 };
/DEFINITION/ { 
        #write down the previous record, iterate between all the sources found
        for (f in src_pos) {
                printf ACC "\t" Def "\t" PMID "\t" src_pos[f] "\t" is[f]  "\t" hst[f] "\t" cl[f] "\t" co[f] "\t" cd[f] "\t" org[f] "\t" sty[f] "\t" sr[f] "\t" cb[f] "\n";
    };
        #clear stored info from previous recors and starts to read the posible multiline definition field.
        $1=""; gsub(/ +$/,""); Def=gensub(/%/, "pct", "g") ; delete org; delete sty; delete is; delete cl; delete co; delete cd; delete cb; delete sr; ACC=""; PMID=""; delete src_pos; src=0; op=1 };
    /^            / { 
            #check if it is reading a multiline definition
            if ( op ) { gsub(/^ +/," "); Def=Def "" gensub(/%/, "pct", "g") };};
    # read field data and store in variables. To add your own field just copy and modify host or other below that, then add the new field to the write down functions and remember to add the new variable in the headers (at the "BEGIN" block, and also to clear that new variable in "DEFINITION" block
    /ACCESSION/ { ACC=$NF ; op=0 }; 
    / +PUBMED/ { PMID=$NF };
    / {5}source +/ {src_pos[++src]=$NF;}
    /\/host=/  { hst[src] = gensub(/ +\/host="(.+)"/, "\\1", "g") };
    /\/organism=/ { org[src] = gensub(/ +\/organism="(.+)"/, "\\1", "g") };
    /\/serotype=/ { sty[src] = gensub(/ +\/serotype="(.+)"/, "\\1", "g") };
    /\/strain=/ { sr[src] = gensub(/ +\/strain="(.+)"/, "\\1", "g") };
    /\/isolate=/ { is[src] = gensub(/ +\/isolate="(.+)"/, "\\1", "g") };
    /\/clone=/ { cl[src] = gensub(/ +\/clone="(.+)"/, "\\1", "g") };
    /\/country=/ { co[src] = gensub(/ +\/country="(.+)"/, "\\1", "g") };
    /\/collection_date=/ { cd[src] = gensub(/ +\/collection_date="(.+)"/, "\\1", "g") };
    /\/collected_by=/ { cb[src] = gensub(/ +\/collected_by="(.+)"/, "\\1", "g") };
    END { 
            #write down the previous record, iterate between all the sources found
            for (f in src_pos) {
                    printf ACC "\t" Def "\t" PMID "\t" src_pos[f] "\t" is[f]  "\t" hst[f] "\t" cl[f] "\t" co[f] "\t" cd[f] "\t" org[f] "\t" sty[f] "\t" sr[f]"\t" cd[f]  "\n";
        }; };

