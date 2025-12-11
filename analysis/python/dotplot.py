import os
import sys
import matplotlib
matplotlib.use('Agg')  # <-- use non-GUI backend
import matplotlib.pyplot as plt
import pypangraph as pp

# get command-line arguments
str_i = sys.argv[1]  # first genome
str_j = sys.argv[2]  # second genome
output_dir = sys.argv[3]  # output directory

# load pangraph
pan = pp.Pangraph.from_json("graph.json")

# ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# create figure
fig, ax = plt.subplots(figsize=(8, 8))
pp.dotplot(ax, str_i, str_j, pan, duplicated_color="silver", min_length=150)

# save figure with combined name
fig_name = f"{str_i}_{str_j}_dotplot.png"
fig_path = os.path.join(output_dir, fig_name)
fig.savefig(fig_path, dpi=300, bbox_inches='tight')

print(f"Figure saved to {fig_path}")
