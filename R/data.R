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
#' A list of up-regulated and down-regulated genes of the ERBB2 (HER2) sub type
#' (which is characterized by high expression of ERBB2) of 344 primary breast
#' tumors from lymph node-negative patients (Smid et al. 2008).
#'
#'
#' @format A data frame with 156 rows and 2 variables:
#' \describe{
#'    \item{gene}{The name of a highly differentially expressed gene when HER2
#'    pathway is active.}
#'    \item{expression}{The change in expression of the gene when the pathway
#'    is active i.e. -1 for down-regulated genes and 1 for up-regulated genes.}
#' }
#'
#' @source MSigDB (For up-regulated component of the:
#' \url{https://www.gsea-msigdb.org/gsea/msigdb/cards/SMID_BREAST_CANCER_ERBB2_UP.html})
#' @source MSigDB (For down-regulated component of the signature:
#'  \url{https://www.gsea-msigdb.org/gsea/msigdb/cards/SMID_BREAST_CANCER_ERBB2_DN.html}
#' )
#' @references
#' Smid, M., Wang, Y., Zhang, Y., Sieuwerts, A.M., Yu, J., Klijn, J.G., Foekens,
#' J.A. and Martens, J.W., 2008. Subtypes of breast cancer show preferential
#' site of relapse. Cancer research, 68(9), pp.3108-3114.
"HER2_sig_df"

#' Gene signature for ER and HER2 extracted from PAM50
#'
#' A list of two data frames containing the gene names and their corresponding
#' expression values constituting ER and HER gene signatures obtained from
#' PAM50.
#' @format List of two of two data frames: one for ER signature, the other for
#' HER2 signature:
#' \describe{
#'    \item{ER}{data frame containing gene names and their relative expression
#'    value (-1 for down-regulated genes and 1 for up-regulated genes) in the
#'    ER signature}
#'    \item{HER2}{data frame containing gene names their relative expression
#'    value (-1 for down-regulated genes and 1 for up-regulated genes) in the
#'    HER2 signature}
#' }
"pam50"

#' ER gene expression data set 1
#'
#' A gene expression matrix containing raw read counts from RNA-seq for 60
#' samples (30 ER positive and 30 ER negative). The Raw RNA-seq counts data
#' were obtained from TCGA. The data set is a subset of a larger data set of
#' filtered for breast cancer-primary tumors from 1,101 samples collected from
#' TCGA.
#' @format A data frame with 20,124 rows and 60 variables:
#' \describe{
#'    \item{TCGA.5L.AAT0}{The sample or case ID.}
#'    \item{TCGA.A1.A0SD}{The sample or case ID.}
#'     .
#' }
#' @source \url{https://portal.gdc.cancer.gov}).
"ER_data_se1"

#' HER2 gene expression data set 1
#'
#' A gene expression matrix containing raw read counts from RNA-seq for 60
#' samples (30 HER2 positive and 30 HER2 negative). The Raw RNA-seq counts data
#' were obtained from TCGA. The data set is a subset of a larger data set of
#' filtered for breast cancer-primary tumors from 1,101 samples collected from
#' TCGA.
#' @format A data frame with 20,124 rows and 60 variables:
#' \describe{
##'   \item{TCGA.A1.A0SK}{The sample or case ID.}
#'    \item{TCGA.A2.A04X}{The sample or case ID.}
#' }
#' @source TCGA \url{https://portal.gdc.cancer.gov}).
"HER2_data_se1"



