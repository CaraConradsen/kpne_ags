# ########################################################## ####
# AG introductions in KLEBSIELLA                             ####
# Author:    Cara Conradsen                                  ####
# ########################################################## ####

library(DiagrammeRsvg)
library(rsvg)
library(magick)

# AG pipeline ----------------------------------------------------------
diagram <- DiagrammeR::grViz("digraph{
graph [layout = dot, fontname = Arial, rankdir = TB, pad=0]

node[shape = rectangle, style=\"rounded,filled\", fillcolor = white, margin = 0.25, fontname = Arial]
kpne_assemblies[label =< Hybrid assemblies from E. Feil <br/> <i> K.pneumoniae</i>>]
kpne_data[label =< Detect and remove clones <br/> <i> Mash 2.3, R: stats, n = 260 </i>>]
bakta_anno[label =  <Annotate proteins <br/><font face='Courier New'> Bakta 1.11.2</font>>]
map_ags[label = <Map AGs to MGEs <br/><font face='Courier New'> R: data.table 1.17.8, GenomicRanges 1.61.1 </font>>]
remove_multiple_transfers[label = < Identify &amp; remove genes that arose from<br/>multiple horizontal transfers <br/><font face='Courier New'> R: data.table 1.17.8</font>>]

subgraph cluster_0 {
graph[shape = rectangle]
bgcolor = grey90
color = grey90

label = \"Get MGE regions\";
labeljust = l;  // Set label justification to left
node[shape = rectangle, style=\"rounded,filled\", fillcolor = white, margin = 0.25]
prophage[label = <prophages<br/><font face='Courier New'>PHASTER</font>>]
ices[label = <detect ICEs<br/><font face='Courier New'>blastn ICEberg 2.0</font>>]
ice2[label = <validate ICEs<br/><font face='Courier New'>CONJScan</font>>]
integrons[label = <integrons<br/><font face='Courier New'>IntegronFinder</font>>]
transposons[label = <transposons<br/><font face='Courier New'>TnCentral</font>>]
IS[label = <IS<br/><font face='Courier New'>ISCEcan</font>>]
}

subgraph cluster_1 {
graph[shape = rectangle]
bgcolor = darkseagreen1
color = darkseagreen1

label = \"Get syntenic regions\";
labeljust = l;  // Set label justification to left
node[shape = rectangle, style=\"rounded,filled\", fillcolor = white, margin = 0.25]
concat_fasta[label = <Concatenate gene-fastas<br/>from gffs and long read fastas<br/><font face='Courier New'>R: Biostring 2.77.2,<br/>GenomicRanges 1.61.1</font>>]
pangraph[label = <Construct gene pangraphs <br/><font face='Courier New'>Pangraph 1.2.1</font>>]
msu[label = <Detect conserved core synteny blocks<br/><font face='Courier New'>PyPanGraph 0.1.3</font>>]
junctions[label = <Assign accessory genes to syntenic<br/>sections within core blocks<br/><font face='Courier New'>R: data.table 1.17.8 </font>>]
}


subgraph cluster_2 {
graph[shape = rectangle]
bgcolor = cornsilk1
color = cornsilk1

label = \"Recombination free core phylogenetic tree\";
labeljust = l;  // Set label justification to left
node[shape = rectangle, style=\"rounded,filled\", fillcolor = white, margin = 0.25]
core_align[label = <Generate core alignment<br/><font face='Courier New'>Pangraph 1.2.1</font>>]
core_gub[label = <Detect &amp; remove recombination<br/><font face='Courier New'>Gubbins 3.4.3</font>>]
core_tree[label = <Generate core phylogentic tree<br/><font face='Courier New'>TreeTime 0.11.4</font>>]
}


subgraph cluster_3 {
graph[shape = rectangle]
bgcolor = lightblue
color = lightblue

label = \"Estimate nucleotide diversity\";
labeljust = l;  // Set label justification to left
node[shape = rectangle, style=\"rounded,filled\", fillcolor = white, margin = 0.25]
kaks[label = <Use kaks to get synonymous &pi; estimates <br/><font face='Courier New'>R: seqinr 4.2.36</font>>]
null_vs_emp[label = <Compare observed vs. neutral<br/>expectations &amp; detect outliers<br/><font face='Courier New'>R: stats 4.5.1</font>>]
slim[label = < Run forward evolutionary simulation framework<br/><font face='Courier New'> SLiM 5.1</font>>]
}

subgraph cluster_4 {
graph[shape = rectangle]
bgcolor = lightcoral
color = lightcoral

label = \"Pangenome analysis\";
labeljust = l;  // Set label justification to left
node[shape = rectangle, style=\"rounded,filled\", fillcolor = white, margin = 0.25]
PIRATE[label = <Determine accessory genes <br/><font face='Courier New'>PIRATE 1.0.5</font>>]
ag_align[label = <Generate gene alignments <br/><font face='Courier New'>PIRATE 1.0.5, MAFFT 7.310</font>>]
}

subgraph cluster_5 {
graph[shape = rectangle]
bgcolor = lightblue
color = lightblue

label = \"Neutral expectation\";
labeljust = l;  // Set label justification to left
node[shape = rectangle, style=\"rounded,filled\", fillcolor = white, margin = 0.25]
subtree[label = <Simulate neutral genealogies &amp;<br/>compute descendant length<br/><font face='Courier New'>R: ape 5.8.1</font>>]
seg_sites[label = <Compute neutral SNP expectations &amp;<br/>empirical SNP counts<br/><font face='Courier New'>R: pegas 1.3</font>>]
poissons[label = <Construct neutral Poisson distributions<br/><font face='Courier New'>R: stats 4.5.1</font>>]
}


subgraph cluster_6 {
graph[shape = rectangle]
bgcolor = peachpuff1
color = peachpuff1

label = \"Functional characterisation of outlier genes\";
labeljust = l;  // Set label justification to left
node[shape = rectangle, style=\"rounded,filled\", fillcolor = white, margin = 0.25]
cog[label = <Assign genes to COG categories<br/><font face='Courier New'>eggnog-mapper 2.1</font>>]
kegg[label = <Map genes to KEGG pathways<br/><font face='Courier New'>eggnog-mapper 2.1, bakta 1.11.2</font>>]
}

kpne_assemblies -> kpne_data -> bakta_anno -> PIRATE -> {ag_align map_ags}
{bakta_anno PIRATE} -> concat_fasta -> pangraph -> core_align -> core_gub -> core_tree-> remove_multiple_transfers
pangraph -> msu -> junctions
ag_align -> {kaks seg_sites}
kpne_data -> ices -> ice2 -> map_ags
kpne_data -> {prophage integrons transposons IS} -> map_ags
{map_ags junctions} -> remove_multiple_transfers -> kaks -> null_vs_emp
null_vs_emp -> cog -> kegg
subtree -> seg_sites -> poissons -> null_vs_emp -> slim

kaks -> cog [style=invis, weight=10];
map_ags -> subtree [style=invis, weight=-1];
concat_fasta -> core_align [style=invis, weight=10];
transposons -> ag_align [style=invis, weight=10];
}")

# Export as SVG
svg_code <- export_svg(diagram)

# Convert SVG → magick image
img <- image_read(charToRaw(svg_code))

# Trim all extra whitespace and lines
img_trimmed <- image_trim(img)

# Resize to desired dimensions first
img_resized <- image_resize(img_trimmed, "2000x2000!")  # width x height in pixels

# Then write to file
image_write(img_resized, path = paste0(outdir_fig, "/AG_pipeline_", Sys.Date(), ".png"), format = "png")

# Old code
# igraph[label = <Construct gene graphs from gffs and gfa <br/><font face='Courier New'>R: igraph 2.1.4</font>>]
# consensus_ags[label = <Extract representative AGs <br/><font face='Courier New'> PIRATE select_representative.pl </font>>]
