# ########################################################## ####
# AG introductions in KLEBSIELLA                             ####
# Author:    Cara Conradsen                                  ####
# ########################################################## ####



# AG pipeline ----------------------------------------------------------
DiagrammeR::grViz("digraph{
graph [layout = dot, fontname = Arial, rankdir = TB]

node[shape = rectangle, style=\"rounded,filled\", fillcolor = white, margin = 0.25, fontname = Arial]
kpne_assemblies[label =< Hybrid assemblies from E. Feil <br/> <i> K.pneumoniae, n = 1,705</i>>]
kpne_QC[label = < post-QC FASTAs,<i> n = 1,695 </i><br/><font face='Courier New'>QUAST</font>>]
prokka_anno[label =  <annotate proteins <br/><font face='Courier New'> prokka 1.14.6</font>>]
PIRATE[label = < determine core/accessory genes <br/><font face='Courier New'>PIRATE 1.0.5</font>>]
map_ags[label = < map AGs to MGEs <br/><font face='Courier New'> R: data.table | GRanges </font>>]
consensus_ags[label = < extract representative AGs <br/><font face='Courier New'> PIRATE select_representative.pl </font>>]

subgraph cluster_0 {
graph[shape = rectangle]
bgcolor = grey90
color = grey90

label = \"Identify and remove plasmids contigs\";
labeljust = l;  // Set label justification to left
node[shape = rectangle, style=\"rounded,filled\", fillcolor = white, margin = 0.25]
plasmids[label = <plasmids<br/><font face='Courier New'>MOB-suite </font>>]
mlplasmids[label = <plasmids<br/><font face='Courier New'>mlplasmids</font>>]
}

subgraph cluster_1 {
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

kpne_assemblies -> kpne_QC -> prokka_anno -> PIRATE -> map_ags
kpne_QC -> ices -> ice2 -> map_ags
kpne_QC -> {plasmids mlplasmids prophage integrons transposons IS} -> map_ags
PIRATE -> consensus_ags
}")

