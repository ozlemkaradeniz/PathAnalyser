#' ER gene signature
#'
#' A list of up-regulated and down-regulated genes constituting the gene
#' signature of ER pathway activation. The ER signature was obtained from the
#' sensitivity to endocrine therapy (SET) genomic index defined using a list of
#' genes co-expressed with the estrogen receptor (ESR1) from microarray data
#' (437 profiles) originating from newly diagnosed breast cancer independent of
#' outcome and treatment (Symmans et al., 2010).
#' @format A data frame with 160 rows and 2 variables:
#' \describe{
#'    \item{gene}{The name of a highly differentially expressed gene when ER
#'    pathway is active.}
#'    \item{expression}{The change in expression of the gene when the pathway
#'    is active i.e. -1 for down-regulated genes and 1 for up-regulated genes.}
#' }
#' @source SET index defined by Symmans et al. 2010 (see
#' \url{https://pubmed.ncbi.nlm.nih.gov/20697068/}).
#' (PMID: 20697068
#' PMCID: PMC2953969 DOI: 10.1200/JCO.2010.28.4273)
#' @references Symmans, W.F., Hatzis, C., Sotiriou, C., Andre, F., Peintinger,
#' F., Regitnig, P., Daxenbichler, G., Desmedt, C., Domont, J., Marth, C. and
#' Delaloge, S., 2010. Genomic index of sensitivity to endocrine therapy for
#' breast cancer. Journal of clinical oncology, 28(27), p.4111.
"ER_sig_df"


#' HER2 gene signature
#'
#' A list of up-regulated and down-regulated genes constituting the gene
#' signature of ERBB2 pathway activation. The list of differentially expressed
#' genes were acquired from MCF-7 cells (breast cancer cell line)
#' which were ESR1 (estrogen receptor) positive and engineered to express
#' ligand-activatable ERBB2.
#'
#'
#' @format A data frame with 387 rows and 2 variables:
#' \describe{
#'    \item{gene}{The name of a highly differentially expressed gene when HER2
#'    pathway is active.}
#'    \item{expression}{The change in expression of the gene when the pathway
#'    is active i.e. -1 for down-regulated genes and 1 for up-regulated genes.}
#' }
#'
#' @source MSigDB (For up-regulated component of the:
#' \url{https://www.gsea-msigdb.org/gsea/msigdb/geneset_page.jsp?geneSetName=ERBB2_UP.V1_UP})
#' @source MSigDB (For down-regulated component of the signature:
#'  \url{https://www.gsea-msigdb.org/gsea/msigdb/cards/ERBB2_UP.V1_DN.html}
#' )
#' @references Creighton CJ, Hilger AM, Murthy S, Rae JM, Chinnaiyan AM,
#' El-Ashry D. Activation of mitogen-activated protein kinase in estrogen
#' receptor alpha-positive breast cancer cells in vitro induces an in vivo
#' molecular phenotype of estrogen receptor alpha-negative human breast tumors.
#' Cancer Res. 2006 Apr 1;66(7):3903-11. doi: 10.1158/0008-5472.CAN-05-4363.
#' PMID: 16585219.
"HER2_sig_df"

#' ER gene expression data set 1
#'
#' A gene expression matrix containing raw read counts from RNA-seq for 20
#' samples (10 ER positive and 10 ER negative).
#' @format A data frame with 20,124 rows and 20 variables:
#' \describe{
#'    \item{gene}{The name of the gene for which expression data is provided.}
#'    \item{sample}{The sample or case ID.}
#' }
#' @source \url{}).
"ER_data_se1"

#' HER gene expression data set 1
#'
#' A gene expression matrix for 20 samples (10 ER positive and 10 ER negative).
#' @format A data frame with 20,124 rows and 20 variables:
#' \describe{
#'    \item{gene}{The name of the gene for which expression data is provided.}
#'    \item{sample}{The sample or case ID.}
#' }
#' @source TCGA \url{}).
"HER2_data_se1"

#' ER gene expression data set 2
#'
#' A gene expression matrix for 20 samples (10 ER positive and 10 ER negative).
#' @format A data frame with 23,113 rows and 20 variables:
#' \describe{
#'    \item{gene}{The name of the gene for which expression data is provided.}
#'    \item{sample}{The sample or case ID.}
#' }
#' @source cBioPortal \url{}).
"ER_data_se2"

#' HER gene expression data set 2
#'
#' A gene expression matrix for 20 samples (10 HER2 positive and 10 HER2 negative).
#' @format A data frame with 23,113 rows and 20 variables:
#' \describe{
#'    \item{gene}{The name of the gene for which expression data is provided.}
#'    \item{sample}{The sample or case ID.}
#' }
#' @source cBioPortal \url{}).
"HER2_data_se2"
