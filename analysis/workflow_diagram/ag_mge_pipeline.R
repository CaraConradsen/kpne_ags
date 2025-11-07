# ########################################################## ####
# AG introductions in KLEBSIELLA                             ####
# Author:    Cara Conradsen                                  ####
# ########################################################## ####



# AG pipeline ----------------------------------------------------------
DiagrammeR::grViz("digraph{
graph [layout = dot, fontname = Arial, rankdir = TB]

node[shape = rectangle, style=\"rounded,filled\", fillcolor = white, margin = 0.25, fontname = Arial]
kpne_assemblies[label =< Long read fastas from E. Feil <br/> <i> K.pneumoniae, n = 485 </i>>]
bakta_anno[label =  <Annotate proteins <br/><font face='Courier New'> bakta/1.11.2</font>>]
map_ags[label = <Map AGs to MGEs <br/><font face='Courier New'> R: data.table | GRanges </font>>]
consensus_ags[label = <Extract representative AGs <br/><font face='Courier New'> PIRATE select_representative.pl </font>>]
remove_multiple_transfers[label = < Identify genes that arose from multiple transfers <br/><font face='Courier New'> R: data.table </font>>]
slim[label = < Run forward evolutionary simulation framework <br/><font face='Courier New'> SLiM </font>>]

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
transposons[label = 'transposons\n?']
IS[label = <IS<br/><font face='Courier New'>ISCEcan</font>>]
}

subgraph cluster_1 {
graph[shape = rectangle]
bgcolor = grey90
color = grey90

label = \"Get syntenic regions\";
labeljust = l;  // Set label justification to left
node[shape = rectangle, style=\"rounded,filled\", fillcolor = white, margin = 0.25]
igraph[label = <Construct gene graphs from gffs and gfa <br/><font face='Courier New'>R: igraph 2.1.4</font>>]
}

subgraph cluster_2 {
graph[shape = rectangle]
bgcolor = lightblue
color = lightblue

label = \"Estimate nucleotide diversity\";
labeljust = l;  // Set label justification to left
node[shape = rectangle, style=\"rounded,filled\", fillcolor = white, margin = 0.25]
kaks[label = < Use kaks to get synonymous pi estimates <br/><font face='Courier New'>R: seqinr 4.2.36</font>>]
}

subgraph cluster_3 {
graph[shape = rectangle]
bgcolor = lightcoral
color = lightcoral

label = \"Pangenome analysis\";
labeljust = l;  // Set label justification to left
node[shape = rectangle, style=\"rounded,filled\", fillcolor = white, margin = 0.25]
PIRATE[label = <Determine accessory genes <br/><font face='Courier New'>PIRATE 1.0.5</font>>]
}

kpne_assemblies -> bakta_anno -> PIRATE -> map_ags
{bakta_anno PIRATE} -> igraph
kpne_assemblies -> ices -> ice2 -> map_ags
kpne_assemblies -> {prophage integrons transposons IS} -> map_ags
PIRATE -> consensus_ags
{map_ags igraph} -> remove_multiple_transfers -> kaks -> slim
}")

