#' Celiac disease GWAS example data
#'
#' @name celiac_gwas_example
#' @docType data
#' @usage data(celiac_gwas_example)
#' @format A data frame with columns:
#' \describe{
#'   \item{SNP}{SNP identifier}
#'   \item{chromosome}{Chromosome location}
#'   \item{mapped_gene}{Mapped gene symbol}
#'   \item{odds_ratio}{Odds ratio}
#'   \item{p_value}{P-value}
#'   \item{position}{Genomic position}
#'   \item{source}{Data source}
#' }
#' @keywords datasets
#' 
#' @details
#' **Data Source and Curation:**
#' This dataset is simulated based on established celiac disease genetic 
#' associations from major GWAS studies. The variants and effect sizes are 
#' modeled after findings from:
#' 
#' - The NHGRI-EBI GWAS Catalog (https://www.ebi.ac.uk/gwas/)
#' - Trynka et al. (2011) Nature Genetics 43:1193-1201
#' - Dubois et al. (2010) Nature Genetics 42:295-302
#'
#' **Key Features:**
#' - Includes top HLA region associations (HLA-DQA1, HLA-DQB1, HLA-DRB1)
#' - Contains established non-HLA risk loci (IL2-IL21, CCR3, etc.)
#' - Effect sizes reflect realistic odds ratios from celiac disease GWAS
#' - P-values span the range from highly significant to suggestive associations
#'
#' **Usage:**
#' This dataset is intended for demonstration and testing purposes in the 
#' CeDExplorer package. For actual research, users should obtain current 
#' GWAS data from the GWAS Catalog or conduct their own association studies.
#'
#' @references
#' Jabri, B., & Sollid, L. M. (2009). T-cell and dendritic cell responses 
#' in celiac disease. Nature Reviews Gastroenterology & Hepatology, 6(4), 
#' 220-227. https://doi.org/10.1038/nrgastro.2009.11
#'
#' Di Niro, R., D'Angelo, S., Secca, P., & De Re, V. (2011). Gene expression 
#' profiling of the intestinal mucosa of celiac patients. Cellular and 
#' Molecular Immunology, 8(4), 301-309. https://doi.org/10.1038/cmi.2011.10
#'
#' Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold 
#' change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15(12), 
#' 550. https://doi.org/10.1186/s13059-014-0550-8
#'
#' @source 
#' Simulated data based on established celiac disease associations from:
#' - NHGRI-EBI GWAS Catalog (accession ranges: GCST003044, GCST003045)
#' - ImmunoBase celiac disease portal
#' - Literature-curated associations from major celiac disease GWAS
#'
#' @examples
#' \dontrun{
#' # Load the dataset
#' data(celiac_gwas_example)
#' 
#' # View the top associations
#' head(celiac_gwas_example)
#' 
#' # Filter for genome-wide significant hits (p < 5e-8)
#' significant_hits <- subset(celiac_gwas_example, p_value < 5e-8)
#' 
#' # Plot the top hits
#' plotTopGwasHits(celiac_gwas_example, top_n = 20)
#' }
"celiac_gwas_example"

#' Example gene expression data
#'
#' @name example_expression
#' @docType data
#' @usage data(example_expression)
#' @format A data frame with columns:
#' \describe{
#'   \item{gene}{Gene symbol}
#'   \item{sample}{Sample identifier}
#'   \item{condition}{Condition (CeD or Control)}
#'   \item{expression}{Expression value}
#' }
#' @keywords datasets
#' 
#' @references
#' Data simulation methodology based on:
#' 
#' Jabri, B., & Sollid, L. M. (2009). T-cell and dendritic cell responses 
#' in celiac disease. *Nature Reviews Gastroenterology & Hepatology*, *6*(4), 
#' 220–227. https://doi.org/10.1038/nrgastro.2009.11
#' 
#' Expression patterns modeled after:
#' Duggan, S. P., Gallagher, W. M., Fox, E. J., Abdel-Latif, M. M., 
#' Reynolds, J. V., & Kelleher, D. (2006). Low pH conditions mimic 
#' duodenal exposure to non-steroidal anti-inflammatory drugs (NSAIDs) 
#' and activate expression of genes involved in tissue homeostasis. 
#' *Journal of Biological Chemistry*, *281*(42), 31268–31278. 
#' https://doi.org/10.1074/jbc.M604986200
#' 
#' Statistical distribution based on:
#' Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation 
#' of fold change and dispersion for RNA-seq data with DESeq2. 
#' *Genome Biology*, *15*(12), 550. https://doi.org/10.1186/s13059-014-0550-8
#' 
#' HLA gene expression patterns referenced from:
#' Trynka, G., Hunt, K. A., Bockett, N. A., Romanos, J., Mistry, V., 
#' Szperl, A., Bakker, S. F., Bardella, M. T., Bhaw-Rosun, L., Castillejo, G., 
#' de la Concha, E. G., de Almeida, R. C., Dias, K. R., van Diemen, C. C., 
#' Dubois, P. C., Duerr, R. H., Edkins, S., Franke, L., Fransen, K., ... 
#' van Heel, D. A. (2011). Dense genotyping identifies and localizes multiple 
#' common and rare variant association signals in celiac disease. 
#' *Nature Genetics*, *43*(12), 1193–1201. https://doi.org/10.1038/ng.998
#' 
#' Immune gene signatures modeled after:
#' Hänzelmann, S., Castelo, R., & Guinney, J. (2013). GSVA: Gene set 
#' variation analysis for microarray and RNA-seq data. *BMC Bioinformatics*, 
#' *14*, 7. https://doi.org/10.1186/1471-2105-14-7
#' 
#' @source
#' Simulated expression data modeled after:
#' - Celiac disease transcriptomic studies
#' - Typical RNA-seq count distributions
#' - Established expression differences in immune genes
"example_expression"

#' Example Expression Matrix
#'
#' A simulated gene expression matrix in typical RNA-seq analysis format
#' with genes as rows and samples as columns.
#'
#' @format A matrix with 200 rows (genes) and 50 columns (samples)
#' 
#' @references
#' Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2010). edgeR: A 
#' Bioconductor package for differential expression analysis of digital gene 
#' expression data. Bioinformatics, 26(1), 139-140. 
#' https://doi.org/10.1093/bioinformatics/btp616
#'
#' Ritchie, M. E., Phipson, B., Wu, D., Hu, Y., Law, C. W., Shi, W., & Smyth, 
#' G. K. (2015). limma powers differential expression analyses for 
#' RNA-sequencing and microarray studies. Nucleic Acids Research, 43(7), e47.
#' https://doi.org/10.1093/nar/gkv007
#' @source  
#' Simulated data following typical RNA-seq count distributions
#' and celiac disease expression patterns
"example_expression_matrix"

#' Example metadata
#'
#' @name example_metadata
#' @docType data
#' @usage data(example_metadata)
#' @format A data frame with columns:
#' \describe{
#'   \item{sample}{Sample identifier}
#'   \item{condition}{Condition (CeD or Control)}
#' }
#' @keywords datasets
#' 
#' @references
#' Leek, J. T., Scharpf, R. B., Bravo, H. C., Simcha, D., Langmead, B., 
#' Johnson, W. E., ... & Irizarry, R. A. (2010). Tackling the widespread and 
#' critical impact of batch effects in high-throughput data. Nature Reviews 
#' Genetics, 11(10), 733-739. https://doi.org/10.1038/nrg2825
#' 
#' @source
#' Simulated sample information for demonstration purposes
"example_metadata"

#' Curated Protein-Protein Interaction Data
#'
#' A curated dataset containing protein-protein interactions for
#' celiac disease-related genes. Includes interactions from immune signaling
#' pathways, antigen presentation complexes, and inflammatory responses.
#'
#' @format A list with two components:
#' \describe{
#'   \item{edges}{A data frame with interaction data:
#'     \describe{
#'       \item{gene1}{First interacting gene}
#'       \item{gene2}{Second interacting gene}
#'       \item{score}{Interaction confidence score (0-1)}
#'       \item{evidence}{Type of evidence: complex, signaling, costimulation, predicted}
#'     }
#'   }
#'   \item{nodes}{A data frame with gene information:
#'     \describe{
#'       \item{gene}{Gene symbol}
#'       \item{type}{Gene type: MHC, Cytokine/Chemokine, Transcription Factor, Cell Surface, Other}
#'       \item{importance}{Relative importance score (0-1)}
#'     }
#'   }
#' }
#'
#' @references
#' Szklarczyk, D., Gable, A. L., Lyon, D., Junge, A., Wyder, S., 
#' Huerta-Cepas, J., ... & Mering, C. V. (2019). STRING v11: Protein-protein 
#' association networks with increased coverage, supporting functional 
#' discovery in genome-wide experimental datasets. Nucleic Acids Research, 
#' 47(D1), D607-D613. https://doi.org/10.1093/nar/gky1131
#'
#' Oughtred, R., Rust, J., Chang, C., Breitkreutz, B. J., Stark, C., 
#' Willems, A., ... & Tyers, M. (2021). The BioGRID database: A comprehensive 
#' biomedical resource of curated protein, genetic, and chemical interactions. 
#' Protein Science, 30(1), 187-200. https://doi.org/10.1002/pro.3978
#'
#' Sollid, L. M., Qiao, S. W., Anderson, R. P., Gianfrani, C., & Koning, F. 
#' (2012). Nomenclature and listing of celiac disease relevant gluten T-cell 
#' epitopes restricted by HLA-DQ molecules. Immunogenetics, 64(6), 455-460.
#' https://doi.org/10.1007/s00251-012-0599-z
#'
#' @source
#' Curated from STRING database, BioGRID, and literature on celiac disease 
#' immune pathways. Includes interactions relevant to antigen presentation, 
#' T-cell activation, and inflammatory signaling.
#'
#' @examples
#' data(ppi_data)
#' head(ppi_data$edges)
#' head(ppi_data$nodes)
"ppi_data"

# [END] 