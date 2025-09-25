# 📄 Turning heterogeneity of statistical epistasis networks to an advantage
This repository contains code to reproduce the results from **"Turning heterogeneity of statistical epistasis networks to an advantage"**.

The pipeline includes:

1. Selecting significant gene pairs 
2. Clustering the epistasis networks with netANOVA  
3. Visualizing protocol-specific networks  
4. Building aggregated networks and comparing them

---


## 📁 Data
The data/ folder in this repository contains the outcomes of the different epistasis detection methods.


## 🚀 Usage

Clone or download this repository, then run the following R code:

```
# Load helper functions
source("functions.r")
source("netANOVA.R")

# 1) Select significant gene pairs
result <- select_significant_pairs(
  base_dir = "path/to/data",
  out_dir  = "results/significance/Gene/",
  bonf_filtered   = 2.6e-6,
  bonf_unfiltered = 4.2e-10,
  p_cutoff        = 0.05, # for mbmdr, casmap, antEpiSeeker
  rel_cutoff      = 0.7, # for NNweights & LightGBM (as prop of max)
  epiGTBN_min     = 1, # for epiGTBN ( > epiGTBN_min)
)

# 2) Cluster the epistasis networks with netANOVA
res <- run_netanova_pipeline(
  in_dir  = "results/significance/Gene",
  out_dir = "results/netANOVA"
)

# 3) Plot the SENs
plot_protocol_networks(
  in_dir   = "results/significance/Gene/ensembl",
  pattern  = "*.txt",
  seed     = 1,
  name_map = NULL,   # optionally supply custom names for protocols
  out_dir  = "results/plots"
)

# 4) Aggregated networks (Union, Intersection, LCA)
res <- aggregate_and_compare_networks(
  clustering_file     = "results/netANOVA/netANOVAClustering_deltaCon_gigascience.txt",
  sig_folder          = "results/significance/Gene",
  pattern             = "\\.txt$",
  nclass              = 2,
  seed                = 1,
  out_plot_dir        = "results/aggregated_network",
  out_lca_symbol_file = "results/aggregated_network/lca_symbol_pairs.txt"
)
```


## 📊 Outputs
The pipeline generates the following outputs:

- results/significance/Gene/: Significant gene pairs per method

- results/netANOVA/: Clustering results

- results/plots/: Protocol-specific SEN visualizations

- results/aggregated_network/: network plots

