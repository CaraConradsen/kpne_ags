import pypangraph as pp

graph = pp.Pangraph.from_json("graph.json")

threshold_len = 120

MSU_mergers, MSU_paths, MSU_len = pp.minimal_synteny_units(graph, threshold_len)

import json

with open("MSU_mergers.json", "w") as f:
	json.dump(MSU_mergers, f)

with open("MSU_len.json", "w") as f:
	json.dump(MSU_len, f)
	
with open("MSU_paths.json", "w") as f:
	json.dump(MSU_paths, f, indent=2, default=str)

