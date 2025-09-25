select_significant_pairs <- function(
    base_dir,
    out_dir,
    bonf_filtered    = 2.6e-6,
    bonf_unfiltered  = 4.2e-10,
    p_cutoff         = 0.05,   # for mbmdr, casmap, antEpiSeeker
    rel_cutoff       = 0.7,    # for NNweights & LightGBM (as prop of max)
    epiGTBN_min      = 1,      # for epiGTBN (V1 > epiGTBN_min)
    do_gigascience   = TRUE,   # run the symbol->ENSEMBL conversion steps
    gigascience_symbol_dir = file.path(out_dir, "symbol"),
    counts_filename  = "pair_counts_summary.tsv"  # <- new: name of the summary table
) {
  suppressPackageStartupMessages({
    require(data.table)
  })
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  .fpath <- function(...) file.path(base_dir, ...)
  .opath <- function(...) file.path(out_dir,  ...)
  .exists_or_skip <- function(f) if (!file.exists(f)) { message("Skipping (not found): ", f); return(FALSE) } else TRUE
  
  # collector for per-file counts
  .counts <- list()
  .record_count <- function(infile, outfile, n_in, n_out, status = "ok", note = NA_character_) {
    .counts[[length(.counts) + 1L]] <<- data.table(
      input_basename  = basename(infile),
      input_path      = infile,
      output_basename = basename(outfile),
      output_path     = outfile,
      n_in            = n_in,
      n_out           = n_out,
      status          = status,
      note            = note
    )
  }
  
  filter_pairs <- function(infile, cols, predicate, outfile) {
    if (!.exists_or_skip(infile)) {
      .record_count(infile, outfile, n_in = NA_integer_, n_out = NA_integer_, status = "missing")
      return(invisible(NULL))
    }
    dt <- as.data.table(fread(infile))
    # Keep only the columns we care about (if a column is missing, error early)
    missing_cols <- setdiff(cols, names(dt))
    if (length(missing_cols)) {
      .record_count(infile, outfile, n_in = NA_integer_, n_out = NA_integer_, status = "error",
                    note = paste("Missing columns:", paste(missing_cols, collapse = ",")))
      warning("Skipping (missing columns): ", infile, " -> ", paste(missing_cols, collapse = ","))
      return(invisible(NULL))
    }
    dt <- dt[, ..cols]
    n_in <- nrow(unique(dt[, 1:2]))  # unique gene1/gene2 pairs before filtering (based on selected cols)
    keep <- predicate(dt)            # logical vector
    if (!is.logical(keep) || length(keep) != nrow(dt)) {
      .record_count(infile, outfile, n_in = n_in, n_out = NA_integer_, status = "error",
                    note = "Predicate did not return a logical vector of correct length")
      warning("Predicate returned invalid result for: ", infile)
      return(invisible(NULL))
    }
    dt_filtered <- dt[which(keep)]
    dt_filtered <- unique(dt_filtered[, 1:2])
    n_out <- nrow(dt_filtered)
    fwrite(dt_filtered, outfile, sep = "\t", col.names = FALSE)
    .record_count(infile, outfile, n_in = n_in, n_out = n_out, status = "ok")
    invisible(NULL)
  }
  
  # ------------------------------
  # plink
  # ------------------------------
  filter_pairs(
    .fpath("plink_boost_binary_filteredhla_LD_selfInteraction_gene.txt"),
    c("gene1","gene2","P"),
    function(d) d$P < bonf_filtered,
    .opath("plink_boost_binary_filteredhla_LD_selfInteraction_gene.txt")
  )
  
  filter_pairs(
    .fpath("plink_boost_binary_unfilteredhla_LD_selfInteraction_gene.txt"),
    c("gene1","gene2","P"),
    function(d) d$P < bonf_unfiltered,
    .opath("plink_boost_binary_unfilteredhla_LD_selfInteraction_gene.txt")
  )
  
  filter_pairs(
    .fpath("plink_linearRegression_filteredhla_LD_selfInteraction_gene.txt"),
    c("gene1","gene2","P"),
    function(d) d$P < bonf_filtered,
    .opath("plink_linearRegression_filteredhla_LD_selfInteraction_gene.txt")
  )
  
  filter_pairs(
    .fpath("plink_linearRegression_unfilteredhla_LD_selfInteraction_gene.txt"),
    c("gene1","gene2","P"),
    function(d) d$P < bonf_unfiltered,
    .opath("plink_linearRegression_unfilteredhla_LD_selfInteraction_gene.txt")
  )
  
  # ------------------------------
  # epiHSIC (unfiltered)
  # ------------------------------
  filter_pairs(
    .fpath("epiHSIC_unfilteredhla_LD_selfInteraction_gene.txt"),
    c("gene1","gene2","ZP"),
    function(d) d$ZP < bonf_unfiltered,
    .opath("epiHSIC_unfilteredhla_LD_selfInteraction_gene.txt")
  )
  
  # ------------------------------
  # NNweights (relative to max stat)
  # ------------------------------
  filter_pairs(
    .fpath("NNweights_corrected_unfilteredhla_LD_selfInteraction_gene.txt"),
    c("gene1","gene2","stat"),
    function(d) d$stat > rel_cutoff * max(d$stat, na.rm = TRUE),
    .opath("NNweights_corrected_unfilteredhla_LD_selfInteraction_gene_score.txt")
  )
  
  filter_pairs(
    .fpath("NNweights_noncorreted_unfilteredhla_LD_selfInteraction_gene.txt"),
    c("gene1","gene2","stat"),
    function(d) d$stat > rel_cutoff * max(d$stat, na.rm = TRUE),
    .opath("NNweights_noncorreted_unfilteredhla_LD_selfInteraction_gene_score.txt")
  )
  
  filter_pairs(
    .fpath("NNweights_top10000_corrected_filteredhla_LD_selfInteraction_gene.txt"),
    c("gene1","gene2","stat"),
    function(d) d$stat > rel_cutoff * max(d$stat, na.rm = TRUE),
    .opath("NNweights_top10000_corrected_filteredhla_LD_selfInteraction_gene_score.txt")
  )
  
  filter_pairs(
    .fpath("NNweights_top10000_notcorrected_filteredhla_LD_selfInteraction_gene.txt"),
    c("gene1","gene2","stat"),
    function(d) d$stat > rel_cutoff * max(d$stat, na.rm = TRUE),
    .opath("NNweights_top10000_notcorrected_filteredhla_LD_selfInteraction_gene_score.txt")
  )
  
  # ------------------------------
  # mbmdr (p-value < p_cutoff)
  # ------------------------------
  filter_pairs(
    .fpath("mbmdr_gammaMaxt_filteredhla_LD_selfInteraction_gene.txt"),
    c("gene1","gene2","p-value"),
    function(d) d[["p-value"]] < p_cutoff,
    .opath("mbmdr_gammaMaxt_filteredhla_LD_selfInteraction_gene_score.txt")
  )
  
  filter_pairs(
    .fpath("mbmdr_gammaMaxt_lowerOrder_unfilteredhla_LD_selfInteraction_gene.txt"),
    c("gene1","gene2","p-value"),
    function(d) d[["p-value"]] < p_cutoff,
    .opath("mbmdr_gammaMaxt_lowerOrder_unfilteredhla_LD_selfInteraction_gene_score.txt")
  )
  
  filter_pairs(
    .fpath("mbmdr_gammaMaxt_unfilteredhla_LD_selfInteraction_gene.txt"),
    c("gene1","gene2","p-value"),
    function(d) d[["p-value"]] < p_cutoff,
    .opath("mbmdr_gammaMaxt_unfilteredhla_LD_selfInteraction_gene_score.txt")
  )
  
  # ------------------------------
  # LightGBM (mean -> 1/1 + mean, then relative to max)
  # ------------------------------
  filter_pairs(
    .fpath("lightGBMfilteredCorrectedhla_LD_selfInteraction_gene.txt"),
    c("gene1","gene2","mean"),
    function(d) { m <- 1/1 + d$mean; m > rel_cutoff * max(m, na.rm = TRUE) },
    .opath("lightGBMfilteredCorrectedhla_LD_selfInteraction_gene_score.txt")
  )
  
  filter_pairs(
    .fpath("lightGBMfilteredUncorrectedhla_LD_selfInteraction_gene.txt"),
    c("gene1","gene2","mean"),
    function(d) { m <- 1/1 + d$mean; m > rel_cutoff * max(m, na.rm = TRUE) },
    .opath("lightGBMfilteredUncorrectedhla_LD_selfInteraction_gene_score.txt")
  )
  
  filter_pairs(
    .fpath("lightGBMunfilteredCorrectedhla_LD_selfInteraction_gene.txt"),
    c("gene1","gene2","mean"),
    function(d) { m <- 1/1 + d$mean; m > rel_cutoff * max(m, na.rm = TRUE) },
    .opath("lightGBMunfilteredCorrectedhla_LD_selfInteraction_gene_score.txt")
  )
  
  filter_pairs(
    .fpath("lightGBMunfilteredUncorrectedhla_LD_selfInteraction_gene.txt"),
    c("gene1","gene2","mean"),
    function(d) { m <- 1/1 + d$mean; m > rel_cutoff * max(m, na.rm = TRUE) },
    .opath("lightGBMunfilteredUncorrectedhla_LD_selfInteraction_gene_score.txt")
  )
  
  # ------------------------------
  # epiGTBN / Episcan & ReliefF (V1 > threshold)
  # ------------------------------
  for (fname in c(
    "epiGTBN_Episcan_corrected_filteredhla_LD_selfInteraction_gene.txt",
    "epiGTBN_Episcan_corrected_unfilteredhla_LD_selfInteraction_gene.txt",
    "epiGTBN_Episcan_notcorrected_filteredhla_LD_selfInteraction_gene.txt",
    "epiGTBN_Episcan_notcorrected_unfilteredhla_LD_selfInteraction_gene.txt",
    "epiGTBN_ReliefF_corrected_unfilteredhla_LD_selfInteraction_gene.txt",
    "epiGTBN_ReliefF_notcorrected_unfilteredhla_LD_selfInteraction_gene.txt"
  )) {
    filter_pairs(
      .fpath(fname),
      c("gene1","gene2","V1"),
      function(d) d$V1 > epiGTBN_min,
      .opath(sub("\\.txt$", "_score.txt", fname))
    )
  }
  
  # ------------------------------
  # epiblaster (glm_score vs Bonf)
  # ------------------------------
  filter_pairs(
    .fpath("epiblaster_corrected_filtered_removed_enhancedhla_LD_selfInteraction_gene.txt"),
    c("gene1","gene2","glm_score"),
    function(d) d$glm_score < bonf_filtered,
    .opath("epiblaster_corrected_filtered_removed_enhancedhla_LD_selfInteraction_gene.txt")
  )
  
  filter_pairs(
    .fpath("epiblaster_corrected_unfiltered_removed_enhancedhla_LD_selfInteraction_gene.txt"),
    c("gene1","gene2","glm_score"),
    function(d) d$glm_score < bonf_unfiltered,
    .opath("epiblaster_corrected_unfiltered_removed_enhancedhla_LD_selfInteraction_gene.txt")
  )
  
  filter_pairs(
    .fpath("epiblaster_uncorrected_filtered_removed_enhancedhla_LD_selfInteraction_gene.txt"),
    c("gene1","gene2","glm_score"),
    function(d) d$glm_score < bonf_filtered,
    .opath("epiblaster_uncorrected_filtered_removed_enhancedhla_LD_selfInteraction_gene.txt")
  )
  
  filter_pairs(
    .fpath("epiblaster_uncorrected_unfiltered_removed_enhancedhla_LD_selfInteraction_gene.txt"),
    c("gene1","gene2","glm_score"),
    function(d) d$glm_score < bonf_unfiltered,
    .opath("epiblaster_uncorrected_unfiltered_removed_enhancedhla_LD_selfInteraction_gene.txt")
  )
  
  # ------------------------------
  # casmap (p < p_cutoff)
  # ------------------------------
  for (fname in c(
    "casmap_region_gwas_dom_unfiltered_all_interactionshla_LD_selfInteraction_gene.txt",
    "casmap_region_gwas_dominant_filtered_all_interactionshla_LD_selfInteraction_gene.txt",
    "casmap_region_gwas_rec_unfiltered_all_interactionshla_LD_selfInteraction_gene.txt",
    "casmap_region_gwas_recessive_filtered_all_interactionshla_LD_selfInteraction_gene.txt"
  )) {
    filter_pairs(
      .fpath(fname),
      c("gene1","gene2","p"),
      function(d) d$p < p_cutoff,
      .opath(fname)
    )
  }
  
  # ------------------------------
  # antEpiSeeker (p < p_cutoff)
  # ------------------------------
  filter_pairs(
    .fpath("v1/antEpiSeeker_filteredhla_LD_selfInteraction_gene.txt"),
    c("gene1","gene2","p"),
    function(d) d$p < p_cutoff,
    .opath("antEpiSeeker_filteredhla_LD_selfInteraction_gene_score.txt")
  )
  filter_pairs(
    .fpath("v1/antEpiSeeker_unfilteredhla_LD_selfInteraction_gene.txt"),
    c("gene1","gene2","p"),
    function(d) d$p < p_cutoff,
    .opath("antEpiSeeker_unfilteredhla_LD_selfInteraction_gene_score.txt")
  )
  
  # ------------------------------
  # Bayesian (first two columns)
  # ------------------------------
  if (.exists_or_skip(.fpath("bayesian_2k_unfiltered.txt"))) {
    infile  <- .fpath("bayesian_2k_unfiltered.txt")
    outfile <- .opath("bayesian_2k_unfiltered_score.txt")
    dt <- fread(infile)
    n_in <- nrow(unique(dt[, 1:2]))
    dt <- unique(dt[, 1:2])
    fwrite(dt, outfile, sep = "\t", col.names = FALSE)
    .record_count(infile, outfile, n_in = n_in, n_out = nrow(dt), status = "ok")
  } else {
    .record_count(.fpath("bayesian_2k_unfiltered.txt"),
                  .opath("bayesian_2k_unfiltered_score.txt"),
                  n_in = NA_integer_, n_out = NA_integer_, status = "missing")
  }
  
  # ------------------------------
  # Gigascience (SYMBOL -> ENSEMBL)
  # ------------------------------
  if (isTRUE(do_gigascience)) {
    suppressPackageStartupMessages({
      require(org.Hs.eg.db)
      require(AnnotationDbi)
    })
    dir.create(gigascience_symbol_dir, recursive = TRUE, showWarnings = FALSE)
    
    convert_gigascience <- function(infile, outfile) {
      if (!.exists_or_skip(infile)) {
        .record_count(infile, outfile, n_in = NA_integer_, n_out = NA_integer_, status = "missing")
        return(invisible(NULL))
      }
      dat <- fread(infile, header = FALSE)
      n_in <- nrow(unique(dat[, 1:2]))
      keys <- unique(c(dat$V1, dat$V2))
      symbols <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                       keys     = keys,
                                       keytype  = "SYMBOL",
                                       column   = "ENSEMBL",
                                       multiVals = "first")
      symdt <- data.table(symbol = names(symbols), ensembl = unname(symbols))
      setnames(symdt, c("symbol","ensembl1"))
      dat <- merge(dat, symdt, by.x = "V1", by.y = "symbol", all.x = TRUE)
      setnames(symdt, c("symbol","ensembl2"))
      dat <- merge(dat, symdt, by.x = "V2", by.y = "symbol", all.x = TRUE)
      res <- unique(na.omit(dat[, .(ensembl1, ensembl2)]))
      fwrite(res, outfile, sep = "\t", col.names = FALSE)
      .record_count(infile, outfile, n_in = n_in, n_out = nrow(res), status = "ok")
      invisible(NULL)
    }
    
    convert_gigascience(
      .fpath("significance/Gene/symbol/gigascience_filteredhla_LD_selfInteraction_gene_score.txt.txt"),
      .opath("gigascience_filteredhla_LD_selfInteraction_gene_score.txt")
    )
    convert_gigascience(
      .fpath("significance/Gene/symbol/gigascience_unfilteredhla_LD_selfInteraction_gene_score.txt.txt"),
      .opath("gigascience_unfilteredhla_LD_selfInteraction_gene_score.txt")
    )
  }
  
  # ------------------------------
  # Write the summary table
  # ------------------------------
  if (length(.counts)) {
    counts_dt <- rbindlist(.counts, use.names = TRUE, fill = TRUE)
    #setcolorder(counts_dt, c("input_basename","n_in","n_out","status"))
    fwrite(counts_dt[,c("input_basename","n_in","n_out","status")], .opath(counts_filename), sep = "\t")
  }
  
  message("Done.")
}


run_netanova_pipeline <- function(
    in_dir,
    out_dir,
    pattern = "*.txt",
    save_prefix = "SENs",
    netanova_functions_dir = "C:/Users/Diane/Documents/2022/netANOVA/May/functions"
) {
  # === Libraries ===
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(igraph)
    library(ggplot2)
    library(reshape2)
    library(RColorBrewer)
    library(ggdendro)
    library(Matrix)
    library(SparseM)
    library(pracma)
  })
  
  # === Prepare output ===
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # === List files ===
  files <- list.files(path = in_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0L) stop("No files matched 'pattern' in 'in_dir'.")
  
  # === Read & merge only NON-EMPTY files ===
  data <- NULL
  dimTools <- integer(0)
  used_files <- character(0)
  
  for (f in files) {
    message("Reading: ", basename(f))
    tmp <- tryCatch(fread(f, header = FALSE), error = function(e) data.table())
    if (nrow(tmp) == 0L) {
      message("  -> Skipping (empty): ", basename(f))
      next
    }
    tmp$V3 <- 1L
    colnames(tmp) <- c("V1", "V2", basename(f))
    tmp <- unique(tmp)
    
    if (is.null(data)) {
      data <- tmp
    } else {
      data <- unique(merge(data, tmp, all = TRUE, by = c("V1", "V2"), allow.cartesian = TRUE))
    }
    dimTools <- c(dimTools, nrow(tmp))
    used_files <- c(used_files, f)
  }
  
  if (is.null(data)) {
    stop("All input files are empty. Nothing to process.")
  }
  if (ncol(data) < 3L) {
    stop("Only one non-empty file found (no tool columns to compare). Need at least two non-empty files.")
  }
  if (length(used_files) < 2L) {
    stop("Need at least two non-empty files to compute distances and clustering.")
  }
  
  dimTools <- data.frame(files = basename(used_files), dimTools, row.names = NULL)
  
  # === Pair occurrences ===
  occurencePairs <- rowSums(data[, 3:ncol(data)], na.rm = TRUE)
  message("Pair occurrence summary:"); print(summary(occurencePairs))
  
  # most frequent
  mostRec <- data[which(occurencePairs == max(occurencePairs)), ]
  secondMostRec <- data[which(occurencePairs == (max(occurencePairs) - 1)), ]
  
  # === Build adjacency matrices ===
  nodes <- unique(c(data$V1, data$V2))
  adj <- matrix(0, nrow = length(nodes), ncol = length(nodes),
                dimnames = list(nodes, nodes))
  
  G <- list()
  nameG <- character(0)
  a <- 1L
  for (i in colnames(data[, 3:ncol(data)])) {
    tmp_data <- as.matrix(data %>% dplyr::select("V1", "V2", i))
    tmp_data <- na.omit(tmp_data)
    if (nrow(tmp_data) > 1L) {  # keep >1 to avoid trivial single-edge graphs
      tmp <- adj
      tmp[tmp_data[, 1:2]] <- 1
      tmp[tmp_data[, 2:1]] <- 1
      G[[a]] <- tmp
      nameG <- c(nameG, i)
      a <- a + 1L
    }
  }
  names(G) <- nameG
  
  if (length(G) < 2L) {
    stop("After removing empty/degenerate inputs, fewer than two networks remain—cannot compute distances.")
  }
  
  # === Compute distances ===
  init <- initialization(G, meth = "deltaCon")
  save(init, file = file.path(out_dir, paste0("init_0Pairs_DeltaCon_", save_prefix, ".RData")))
  
  # === Heatmap ===
  heatmap <- init[[2]]
  colnames(heatmap) <- nameG
  rownames(heatmap) <- nameG
  heatmap <- melt(heatmap)
  p <- ggplot(heatmap, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() + theme(text = element_text(size = 10))
  ggsave(file.path(out_dir, paste0("heatmap_", save_prefix, ".png")), plot = p, width = 10, height = 8)
  
  # === Run netANOVA clustering ===
  output <- netANOVA(Dist = init[[3]], t = 2, method_clust = "complete",
                     perturbation = 0.2, permutations = 9909)
  save(output, file = file.path(out_dir, paste0("output_0Pairs_DeltaCon_", save_prefix, ".RData")))
  
  netclustering <- output[[1]]
  netclustering$net <- nameG
  netclustering <- netclustering[order(netclustering$group_id), ]
  fwrite(netclustering, file.path(out_dir, paste0("netANOVAClustering_deltaCon_", save_prefix, ".txt")))
  
  # === Dendrogram plot ===
  disthc <- init[[3]]
  colnames(disthc) <- nameG
  rownames(disthc) <- nameG
  hc <- hclust(as.dist(disthc), method = "complete")
  dendr <- dendro_data(hc, type = "rectangle")
  clust <- output[[1]]$group_id
  names(clust) <- nameG
  clust.df <- data.frame(label = names(clust), cluster = factor(clust))
  dendr[["labels"]] <- merge(dendr[["labels"]], clust.df, by = "label")
  
  getPalette <- colorRampPalette(brewer.pal(8, "Dark2"))
  dend_plot <- ggplot() +
    geom_segment(data = segment(dendr), aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_text(data = label(dendr),
              aes(x, y, label = label, hjust = 0, color = cluster), size = 3.5) +
    coord_flip() + scale_y_reverse(expand = c(0.8, 0)) +
    theme_minimal()
  ggsave(file.path(out_dir, paste0("dendrogram_", save_prefix, ".png")), plot = dend_plot, width = 12, height = 8)
  
  # === Graph visualizations ===
  for (i in seq_along(G)) {
    a <- graph_from_adjacency_matrix(G[[i]], mode = "undirected")
    Isolated <- which(degree(a) == 0)
    a <- delete.vertices(a, Isolated)
    png(file.path(out_dir, paste0("graph_", i, "_", save_prefix, ".png")),
        width = 800, height = 800)
    plot(a, vertex.label.cex = 1, vertex.size = 10, edge.width = 3)
    dev.off()
  }
  
  message("Pipeline complete. Results written to: ", out_dir)
  return(list(data = data,
              mostRec = mostRec,
              secondMostRec = secondMostRec,
              clustering = netclustering,
              used_files = basename(used_files),
              dimTools = dimTools))
}



plot_protocol_networks <- function(
    in_dir,
    pattern = "*.txt",
    seed = 29,
    # Optional labels for plots (one per file/column). If NULL, uses cleaned filenames.
    name_map = NULL,
    # Layout + visuals
    layout_fn = igraph::layout_in_circle,   # function taking (graph, order=...) or (graph) signature
    weight_divisor = 10,                    # edge weight scaling (your code used /10)
    label_degree_threshold = 20,            # hide labels if > this many non-isolated nodes
    vertex_size = 1,
    vertex_label_cex = 0.8,
    # Output: if NULL, plots to current device; if a directory, saves PNGs there
    out_dir = NULL,
    png_width = 1200, png_height = 1200, png_res = 160
) {
  suppressPackageStartupMessages({
    library(data.table)
    library(igraph)
  })
  
  stopifnot(dir.exists(in_dir))
  files <- list.files(path = in_dir, pattern = pattern, full.names = TRUE, recursive = FALSE)
  if (length(files) == 0L) stop("No files found matching '", pattern, "' in: ", in_dir)
  
  # Friendly names
  clean_names <- basename(files)
  clean_names <- gsub("hla.*", "", clean_names)  # keep your cleaning step
  if (is.null(name_map)) {
    plot_names <- clean_names
  } else {
    if (is.character(name_map) && is.null(names(name_map))) {
      if (length(name_map) != length(files)) stop("name_map length != number of files.")
      plot_names <- name_map
    } else if (!is.null(names(name_map))) {
      # named vector or list: map by basename(files)
      plot_names <- name_map[basename(files)]
      if (any(is.na(plot_names))) stop("Some file basenames not present in names(name_map).")
    } else stop("Unsupported name_map format.")
  }
  
  # Load & merge all as presence matrices (1 for present, NA/0 for absent)
  message("Reading and merging files...")
  dt <- fread(files[[1]], header = FALSE)
  dt$V3 <- 1
  colnames(dt) <- c("V1", "V2", basename(files[[1]]))
  
  if (length(files) > 1) {
    for (f in files[-1]) {
      tmp <- fread(f, header = FALSE)
      tmp$V3 <- 1
      colnames(tmp) <- c("V1", "V2", basename(f))
      tmp <- unique(tmp)
      dt <- unique(merge(dt, tmp, all = TRUE, by = c("V1", "V2"), allow.cartesian = TRUE))
    }
  }
  
  # Convert NAs to 0 for presence/absence
  for (j in 3:ncol(dt)) set(dt, which(is.na(dt[[j]])), j, 0L)
  
  # Prepare layout (fixed across all plots)
  set.seed(seed)
  genes <- sample(unique(c(dt$V1, dt$V2)))
  # Use a “skeleton” graph to compute node order once
  skel <- graph.data.frame(dt[, .(V1, V2)], directed = FALSE)
  # Ensure all genes are present (even if isolated in some layers)
  skel <- igraph::add_vertices(skel, nv = sum(!(genes %in% V(skel)$name)),
                               name = genes[!(genes %in% V(skel)$name)])
  # Circle layout with stable order
  coords <- layout_fn(skel, order = as.numeric(factor(genes)))
  
  # Make sure output directory exists if saving
  if (!is.null(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Plot each protocol
  message("Plotting ", length(files), " graphs...")
  for (i in seq_len(ncol(dt) - 2L)) {
    col_idx <- i + 2L
    tmp <- dt[, .(gene1 = V1, gene2 = V2, weight = get(colnames(dt)[col_idx]))]
    g <- graph.data.frame(tmp, directed = FALSE)
    # Drop zero-weight edges
    z <- which(E(g)$weight == 0)
    if (length(z)) g <- delete.edges(g, z)
    
    # Add any missing vertices so coords align
    missing_vertices <- setdiff(V(skel)$name, V(g)$name)
    if (length(missing_vertices)) g <- add_vertices(g, nv = length(missing_vertices), name = missing_vertices)
    
    # Reorder layout coords to match current graph's vertex order
    coords_ord <- coords[match(V(g)$name, V(skel)$name), , drop = FALSE]
    
    # Styling
    E(g)$weight <- E(g)$weight / weight_divisor
    V(g)$size <- vertex_size
    # Label only if few non-isolated nodes
    active_nodes <- names(which(degree(g) > 0))
    if (length(active_nodes) > label_degree_threshold) {
      V(g)$label <- ""
    } else {
      V(g)$label <- ifelse(V(g)$name %in% active_nodes, V(g)$name, NA)
    }
    V(g)$label.cex <- vertex_label_cex
    
    main_title <- plot_names[[i]]
    
    if (is.null(out_dir)) {
      plot(g, layout = coords_ord, main = main_title)
    } else {
      png(filename = file.path(out_dir, sprintf("network_%02d.png", i)),
          width = png_width, height = png_height, res = png_res)
      plot(g, layout = coords_ord, main = main_title)
      dev.off()
    }
  }
  
  invisible(list(merged = as.data.frame(dt), file_names = basename(files), plot_names = plot_names))
}




aggregate_and_compare_networks <- function(
    clustering_file,
    sig_folder,
    pattern = "\\.txt$",
    # LCA options
    nclass = 2,
    seed = 2022,
    # Output (optional)
    out_plot_dir = NULL,      # if set, saves union.png, intersection.png, lca.png
    out_lca_symbol_file = NULL  # if set, writes SYMBOL pairs for LCA-selected edges
) {
  suppressPackageStartupMessages({
    library(data.table)
    library(igraph)
    library(poLCA)
    # for optional SYMBOL export
    if (!is.null(out_lca_symbol_file)) {
      library(org.Hs.eg.db)
      library(AnnotationDbi)
    }
  })
  set.seed(seed)
  
  # ---------- 1) Take largest group_id ----------
  clus <- fread(clustering_file)
  max_group <- names(which.max(table(clus$group_id)))
  nets_vec  <- as.character(clus[group_id == max_group, net])
  if (!length(nets_vec)) stop("No nets for largest group_id: ", max_group)
  
  # ---------- 2) Load matching files (no header) ----------
  files <- list.files(sig_folder, pattern = pattern, full.names = TRUE)
  if (!length(files)) stop("No files found in: ", sig_folder)
  
  basenames <- basename(files)
  file_keys <- sub(pattern, "", basenames)       # remove ".txt"
  net_keys  <- sub("\\.txt$", "", nets_vec)      # normalize names from clustering
  
  keep_idx   <- file_keys %in% net_keys
  files_keep <- files[keep_idx]
  names_keep <- file_keys[keep_idx]
  if (!length(files_keep)) stop("No matching files for nets in largest group.")
  
  out_list <- lapply(files_keep, function(f) fread(f, header = FALSE))
  names(out_list) <- names_keep
  
  # ---------- 3) Presence matrix on (V1,V2) ----------
  combine_gene_pairs <- function(lst) {
    std_list <- lapply(seq_along(lst), function(i) {
      dt <- unique(lst[[i]][, .(V1 = as.character(V1), V2 = as.character(V2))])
      dt[, val := 1L]
      setnames(dt, "val", names(lst)[i])
      dt
    })
    Reduce(function(x, y) merge(x, y, by = c("V1","V2"), all = TRUE), std_list)
  }
  pres <- combine_gene_pairs(out_list)
  meth <- colnames(pres)[-(1:2)]
  if (length(meth) == 0) stop("No method columns found after merging.")
  
  # We'll use 0/1 presence
  pres[is.na(pres)] <- 0L
  
  # Helper to build graph from edge list (V1,V2)
  build_graph <- function(edges_dt) {
    if (!nrow(edges_dt)) return(make_empty_graph())  # empty
    g <- graph_from_data_frame(d = data.frame(id1 = edges_dt$V1,
                                              id2 = edges_dt$V2),
                               directed = FALSE)
    g
  }
  # Helper to compute metrics
  graph_metrics <- function(g) {
    if (gorder(g) == 0) {
      return(list(nodes = 0L, edges = 0L, clusters = 0L))
    }
    comp <- igraph::components(g)
    list(nodes = gorder(g), edges = gsize(g), clusters = length(comp$csize))
  }
  
  # ---------- 4) UNION network ----------
  union_edges <- pres[rowSums(pres[, ..meth]) > 0, .(V1, V2)]
  g_union <- build_graph(union_edges)
  m_union <- graph_metrics(g_union)
  
  # ---------- 5) INTERSECTION network ----------
  intersection_edges <- pres[rowSums(pres[, ..meth]) == length(meth), .(V1, V2)]
  g_intersection <- build_graph(intersection_edges)
  m_intersection <- graph_metrics(g_intersection)
  
  # ---------- 6) LCA network (minority class = selected) ----------
  # Prepare LCA input: categorical 1..K
  X <- as.data.frame(pres[, ..meth]) + 1L
  colnames(X) <- paste0("V", seq_len(ncol(X)))
  f <- as.formula(paste0("cbind(", paste(colnames(X), collapse = ","), ") ~ 1"))
  
  M1 <- poLCA(f, data = X, nclass = nclass, na.rm = TRUE, verbose = FALSE)
  pred <- M1$predclass
  minority <- as.integer(names(which.min(table(pred))))
  predicted01 <- as.integer(pred == minority)  # 1 = selected, 0 = not
  lca_edges <- data.table(V1 = pres$V1, V2 = pres$V2, sel = predicted01)[sel == 1, .(V1, V2)]
  g_lca <- build_graph(lca_edges)
  m_lca <- graph_metrics(g_lca)
  
  # ---------- 7) Optional: plots ----------
  if (!is.null(out_plot_dir)) {
    dir.create(out_plot_dir, recursive = TRUE, showWarnings = FALSE)
    
    plot_and_save <- function(g, file, main) {
      if (gorder(g) == 0) {
        png(file, width = 1200, height = 1200, res = 150)
        plot.new(); title(main = paste0(main, " (empty)"))
        dev.off()
      } else {
        png(file, width = 1200, height = 1200, res = 150)
        V(g)$size <- 5
        V(g)$label.cex <- 0.8
        V(g)$frame.color <- "white"
        V(g)$color <- "lightblue"
        E(g)$label <- ""
        plot(g, main = main, margin = c(0,0,0,0))
        dev.off()
      }
    }
    
    plot_and_save(g_union,        file.path(out_plot_dir, "union.png"),        "Union network")
    plot_and_save(g_intersection, file.path(out_plot_dir, "intersection.png"), "Intersection network")
    plot_and_save(g_lca,          file.path(out_plot_dir, "lca.png"),          "LCA-selected network")
  }
  
  # ---------- 8) Optional: export LCA SYMBOL pairs ----------
  if (!is.null(out_lca_symbol_file)) {
    keys <- unique(c(lca_edges$V1, lca_edges$V2))
    sym_map <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = keys, keytype = "ENSEMBL",
                                     column = "SYMBOL", multiVals = "first")
    sym_dt <- data.table(ensembl = names(sym_map), symbol = unname(sym_map))
    lca_sym <- merge(lca_edges, setNames(sym_dt, c("V1","symbol1")), by = "V1", all.x = TRUE)
    lca_sym <- merge(lca_sym, setNames(sym_dt, c("V2","symbol2")), by = "V2", all.x = TRUE)
    out <- unique(na.omit(lca_sym[, .(symbol1, symbol2)]))
    setnames(out, c("gene1","gene2"))
    dir.create(dirname(out_lca_symbol_file), recursive = TRUE, showWarnings = FALSE)
    fwrite(out, out_lca_symbol_file, sep = "\t", col.names = FALSE)
  }
  
  # ---------- 9) Comparison table ----------
  comparison <- data.table(
    method   = c("Union", "Intersection", "LCA"),
    nodes    = c(m_union$nodes, m_intersection$nodes, m_lca$nodes),
    edges    = c(m_union$edges, m_intersection$edges, m_lca$edges),
    clusters = c(m_union$clusters, m_intersection$clusters, m_lca$clusters)
  )
  
  return(list(
    presence_matrix   = pres,
    edges = list(union = union_edges, intersection = intersection_edges, lca = lca_edges),
    graphs = list(union = g_union, intersection = g_intersection, lca = g_lca),
    metrics = comparison,
    lca_model = M1
  ))
}
