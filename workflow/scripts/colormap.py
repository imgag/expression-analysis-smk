from itertools import cycle
import pandas as pd
import matplotlib.cm as cm
from matplotlib.colors import to_hex

metadata = pd.read_csv(snakemake.input["metadata"], sep="\t", index_col="sample")
cmap = snakemake.config["colormaps"]["group"]

cmap = cm.get_cmap(cmap)
colors = [to_hex(cmap(i), keep_alpha=True) for i in range(cmap.N)]

df = pd.DataFrame(
    zip(metadata["group"].sort_values().unique(), cycle(colors)),
    columns=["group", "color"],
)
df.to_csv(snakemake.output["colormap"], sep="\t", index=False)
