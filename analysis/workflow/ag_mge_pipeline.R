# ########################################################## ####
# Nonsense mutations in STAPHYLOCOCCUS AUREUS                ####
# Author:    Cara Conradsen                                  ####
# ########################################################## ####



# blast pipeline ----------------------------------------------------------
DiagrammeR::grViz("digraph{
graph [layout = dot, fontname = Arial, rankdir = TB]

node[shape = rectangle, style=\"rounded,filled\", fillcolor = white, margin = 0.25, fontname = Arial]
panX[label =<panX <br/> <i>E.coli </i>&amp; <i>S.aureus</i>>]
panX_merge[label = <Map genes onto gffs<br/><font face='Courier New'>R: data.table</font>>]
panX_faa[label = <Get consensus genes .faa <br/><font face='Courier New'>R: data.table | Biostrings</font>>]
bakta_anno[label =  <annotate proteins <br/><font face='Courier New'> bakta_proteins 1.8.1</font>>]


subgraph cluster_1 {
color = white
panX_accessions[label = 'panX accessions']
panX_pan_genes[label = 'panX genes aln.fa']
}

subgraph cluster_2 {
color = white
panX_pan_gffs[label = <panX gffs<br/><font face='Courier New'>ncbi esearch</font>>]
panX_pan_assembly[label = <scrape NCBI assembly info<br/><font face='Courier New'>ncbi esearch | wget </font>>]
}

subgraph cluster_0 {
  graph[shape = rectangle]
  bgcolor = grey90
  color = grey90
  label = \"MGE databases\";
  labeljust = l;

  node[shape = rectangle, style=\"rounded,filled\", fillcolor = white, margin = 0.25]

  // outer MGE nodes
  IS[label = <IS<br/><font face='Courier New'>ISFinder </font>>]


  // simulated subcluster
  subgraph cluster_0_IS {
    label = \"blastx\"
    style = dashed
    color = grey70
    plasmids[label = <plasmids<br/><font face='Courier New'>PLSDB<br/>v.24_05_31_v2</font>>]
    prophage[label = <prophages<br/><font face='Courier New'>PHASTEST<br/>prophage db</font>>]
    ices[label = 'ICEs\n?']
    integrons[label = <integrons<br/><font face='Courier New'>IntegronFinder</font>>]
    transposons[label = 'transposons\n?']
    mobileOG[label = <life-cycle MGEs<br/><font face='Courier New'>mobileOG-db<br/>beatrix 1.6 All</font>>]
  }
}

panX -> {panX_accessions panX_pan_genes}
panX_accessions -> {panX_pan_assembly panX_pan_gffs}
{panX_pan_gffs panX_pan_assembly panX_pan_genes} -> panX_merge 
panX_pan_genes -> panX_faa -> {plasmids prophage ices mobileOG integrons transposons} ->panX_merge
panX_faa -> bakta_anno -> IS
bakta_anno -> panX_merge
}")
 


# pipeline overview -------------------------------------------------------

DiagrammeR::grViz("digraph{
graph [layout = dot, fontname = Arial, rankdir = TB]

node[shape = rectangle, style=\"rounded,filled\", fillcolor = white, margin = 0.25, fontname = Arial]
panX[label =<panX <br/> <i>E.coli </i>&amp; <i>S.aureus</i>>]
download_genomes[label = 'download fastas']
panX_pan_gffs[label = <panX gffs<br/><font face='Courier New'>ncbi esearch</font>>]
panX_pan_core[label = 'filter out core genes']
panX_merge[label = <Map accessory genes onto gffs<br/><font face='Courier New'>R: data.table</font>>]

subgraph cluster_1 {
color = white
panX_accessions[label = 'panX accessions']
panX_pan_genes[label = 'panX genes']
}

subgraph cluster_0 {
graph[shape = rectangle]
bgcolor = grey90
color = grey90

label = \"MGE programs\";
labeljust = l;  // Set label justification to left
node[shape = rectangle, style=\"rounded,filled\", fillcolor = white, margin = 0.25]
plasmids[label = <plasmids<br/><font face='Courier New'>MOB-suite | PlasmidFinder</font>>]
prophage[label = <prophages<br/><font face='Courier New'>PHASTER</font>>]
ices[label = 'ICEs\n?']
integrons[label = <integrons<br/><font face='Courier New'>IntegronFinder</font>>]
transposons[label = 'transposons\n?']
IS[label = <IS<br/><font face='Courier New'>ISCEcan</font>>]
}

node[shape = rectangle, style=\"rounded,filled\", fillcolor = white, margin = 0.25, fontname = Arial]
anno[label = <annotate genomes<br/><font face='Courier New'> Prokka</font>>]
blast_plasmid[label = 'blast against plasmid accessions']
mges_regions[label = 'extract regions as gff']
mges[label = <assign MGEs<br/><font face='Courier New'>R: GRanges</font>>]

panX -> {panX_accessions panX_pan_genes}
panX_accessions -> download_genomes -> {plasmids prophage ices integrons transposons IS anno}
panX_accessions -> panX_pan_gffs
{prophage ices integrons transposons IS} -> mges_regions -> mges
plasmids -> blast_plasmid -> mges_regions
panX_pan_genes -> panX_pan_core
{panX_pan_gffs panX_pan_core} -> panX_merge -> mges
anno -> mges_regions
}")


