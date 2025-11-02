devtools::load_all()
data("celiac_gwas_example", package = "CeDExplorer")

# Extract gene symbols from example dataset
genes <- unique(na.omit(celiac_gwas_example$mapped_gene))

# Run enrichment
res <- mapToPathways(genes)

# View top hits
head(res$results)

# Plot
res$plot
