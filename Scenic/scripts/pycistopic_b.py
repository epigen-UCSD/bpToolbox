import os
import pycisTopic
from pycisTopic.cistopic_class import create_cistopic_object
pycisTopic.__version__

import pandas as pd
import numpy as np
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt
import pickle



### Read model output and add to object
wdir = sys.argv[1]
project = "ZM"

models = []
for file in os.listdir("pycistopic_outputs/model_out/"):
    if file.startswith("Topic") and file.endswith(".pkl"):
        model = pickle.load(open(os.path.join("pycistopic_outputs/model_out/", file), "rb"))
        models.append(model)


with open("pycistopic_outputs/cistopic_object.pkl", 'rb') as file:
    cistopic_object = pickle.load(file)



from pycisTopic.lda_models import evaluate_models
model = evaluate_models(models, return_model = True)


cistopic_object.add_LDA_model(model)
pickle.dump(
    cistopic_object,
    open("pycistopic_outputs/cistopic_object.pkl", "wb")
)


### Plot and binarize models 

from pycisTopic.clust_vis import (
    find_clusters,
    run_umap,
    run_tsne,
    plot_metadata,
    plot_topic,
    cell_topic_heatmap
)


find_clusters(
    cistopic_object,
    target="cell",
    k=10,
    res=[0.2, 0.4],
    prefix=f"{project}",
    scale=True
)


run_umap(cistopic_object, target="cell", scale=True)


### todo: write all plots to file  
plot_topic(
    cistopic_object,
    reduction_name="UMAP",
    target="cell",
    num_columns=5
)


from pycisTopic.topic_binarization import binarize_topics

region_bin_topics_top_3k = binarize_topics(
    cistopic_object, method='ntop', ntop=3_000,
    plot=True, num_columns=5
)

region_bin_topics_otsu = binarize_topics(
    cistopic_object, method="otsu",
    plot=True, num_columns=5
)

binarized_cell_topic = binarize_topics(
    cistopic_object,
    target="cell",
    method="li",
    plot=True,
    num_columns=5, 
    nbins=100
)

### write topics to bed 
os.makedirs(f"{wdir}/region_sets_otsu", exist_ok=True)
os.makedirs(f"{wdir}/region_sets_yen", exist_ok=True)
os.makedirs(f"{wdir}/region_sets_3k", exist_ok=True)
os.makedirs(os.path.join(f"{wdir}/region_sets_otsu", "Topics_otsu"), exist_ok=True)
os.makedirs(os.path.join(f"{wdir}/region_sets_yen", "Topics_yen"), exist_ok=True)
os.makedirs(os.path.join(f"{wdir}/region_sets_3k", "Topics_3k"), exist_ok=True)


from pycisTopic.utils import region_names_to_coordinates

for topic in region_bin_topics_otsu:
    region_names_to_coordinates(
        region_bin_topics_otsu[topic].index
    ).sort_values(["Chromosome", "Start", "End"]).to_csv(
        os.path.join("region_sets_otsu", "Topics_otsu", f"{topic}.bed"),
        sep="\t", header=False, index=False
    )


for topic in region_bin_topics_yen:
    region_names_to_coordinates(
        region_bin_topics_yen[topic].index
    ).sort_values(["Chromosome", "Start", "End"]).to_csv(
        os.path.join("region_sets_yen", "Topics_yen", f"{topic}.bed"),
        sep="\t", header=False, index=False
    )


for topic in region_bin_topics_top_3k:
    region_names_to_coordinates(
        region_bin_topics_top_3k[topic].index
    ).sort_values(["Chromosome", "Start", "End"]).to_csv(
        os.path.join("region_sets_3k", "Topics_3k", f"{topic}.bed"),
        sep="\t", header=False, index=False
    )

