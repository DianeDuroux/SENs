#!/usr/bin/env Rscript

# Load main functions
source("functions.r")
source("netANOVA.R")

# Ensure results directories exist
dir.create("results/", recursive = TRUE, showWarnings = FALSE)
dir.create("results/significance/", recursive = TRUE, showWarnings = FALSE)
dir.create("results/significance/Gene", recursive = TRUE, showWarnings = FALSE)
dir.create("results/netANOVA", recursive = TRUE, showWarnings = FALSE)
dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("results/aggregated_network", recursive = TRUE, showWarnings = FALSE)
dir.create("results/significance/Gene/ensembl", recursive = TRUE, showWarnings = FALSE)

# ---- STEP 1: Select significant pairs ----
cat("Running: select_significant_pairs...\n")
res1 <- select_significant_pairs(
  base_dir       = "data",
  out_dir        = "results/significance/Gene/",
  bonf_filtered   = 2.6e-6,
  bonf_unfiltered = 4.2e-10,
  p_cutoff        = 0.05,
  rel_cutoff      = 0.7,
  epiGTBN_min     = 1
)

# ---- STEP 2: NetANOVA pipeline ----
cat("Running: run_netanova_pipeline...\n")
res2 <- run_netanova_pipeline(
  in_dir  = "results/significance/Gene",
  out_dir = "results/netANOVA"
)

# ---- STEP 3: Protocol plots ----
cat("Running: plot_protocol_networks...\n")
plot_protocol_networks(
  in_dir   = "results/significance/Gene/", ## Modification pattern.
  # in_dir   = "results/netANOVA", ## Modification pattern. 
  pattern  = "*.txt",
  seed     = 1,
  name_map = NULL,
  out_dir  = "results/plots"
)

# ---- STEP 4: Aggregated networks ----
cat("Running: aggregate_and_compare_networks...\n")
res3 <- aggregate_and_compare_networks(
  clustering_file     = "results/netANOVA/netANOVAClustering_deltaCon_SENs.txt",
  sig_folder          = "results/significance/Gene",
  pattern             = "\\.txt$",
  nclass              = 2,
  seed                = 1,
  out_plot_dir        = "results/aggregated_network",
  out_lca_symbol_file = "results/aggregated_network/lca_symbol_pairs.txt"
)

cat("Pipeline finished successfully!\n")
