#' Example Celiac Disease GWAS Data
#'
#' A curated dataset containing simulated genome-wide association study (GWAS) 
#' results for celiac disease. This dataset mimics real GWAS findings with 
#' strong associations in the HLA region and other established risk loci.
#'
#' @format A data frame with 50 rows and 8 variables:
#' \describe{
#'   \item{association_id}{Unique identifier for each association}
#'   \item{SNP}{Variant identifier (rsID)}
#'   \item{mapped_gene}{Nearest gene to the variant}
#'   \item{p_value}{Association p-value}
#'   \item{odds_ratio}{Effect size as odds ratio}
#'   \item{beta}{Effect size as beta coefficient (NA for some variants)}
#'   \item{standard_error}{Standard error of effect size estimate}
#'   \item{risk_allele}{Risk allele for the association}
#'   \item{source}{Data source identifier}
#' }
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

#' Example Gene Expression Data
#'
#' A simulated dataset containing gene expression values for celiac disease 
#' and control samples. Data mimics RNA-seq expression patterns with 
#' biologically realistic distributions and condition-specific differences.
#'
#' @format A data frame with 1000 rows and 6 variables:
#' \describe{
#'   \item{sample}{Sample identifier (S1-S1000)}
#'   \item{condition}{Condition, either "CeD" or "Control"}
#'   \item{HLA-DQA1}{HLA-DQA1 expression values}
#'   \item{HLA-DQB1}{HLA-DQB1 expression values} 
#'   \item{IL2}{IL2 expression values}
#'   \item{IFNG}{IFNG expression values}
#' }
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

#' Example Sample Metadata
#'
#' Simulated sample metadata matching the example expression datasets.
#'
#' @format A data frame with 50 rows and 3 variables:
#' \describe{
#'   \item{sample_id}{Sample identifier}
#'   \item{condition}{Disease condition}
#'   \item{batch}{Processing batch}
#' }
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
