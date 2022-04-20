#' ER gene signature
#'
#' A list of up-regulated and down-regulated genes constituting the gene
#' signature of ER pathway activation. The ER signature was obtained from the
#' sensitivity to endocrine therapy (SET) genomic index defined using a list of
#' genes co-expressed with the estrogen receptor (ESR1 or ER) from microarray
#' data (437 profiles) originating from newly diagnosed breast cancer
#' independent of outcome and treatment (Symmans et al., 2010).
#' @format A data frame with 160 rows (gene names) and 2 variables:
#' \describe{
#'    \item{gene}{The HUGO gene symbols of a highly differentially expressed gene
#'    when ER pathway is active.}
#'    \item{expression}{The change in expression of the gene when the pathway
#'    is active i.e. -1 for down-regulated genes and 1 for up-regulated genes.}
#' }
#' @source SET index defined by Symmans et al. 2010 (see
#' \url{https://pubmed.ncbi.nlm.nih.gov/20697068/}).
#' @references Symmans, W.F., Hatzis, C., Sotiriou, C., Andre, F., Peintinger,
#' F., Regitnig, P., Daxenbichler, G., Desmedt, C., Domont, J., Marth, C. and
#' Delaloge, S., 2010. Genomic index of sensitivity to endocrine therapy for
#' breast cancer. Journal of clinical oncology, 28(27), p.4111. doi:
#' \url{https://doi.org/10.1200/jco.2010.28.4273}
"ER_sig_df"


#' HER2 gene signature
#'
#' A list of up-regulated and down-regulated genes of the ERBB2 (HER2) sub type
#' (which is characterized by high expression of ERBB2) of 344 primary breast
#' tumors from lymph node-negative patients (Smid et al. 2008).
#'
#'
#' @format A data frame with 156 rows (gene names) and 2 variables:
#' \describe{
#'    \item{gene}{The HUGO gene symbols of a highly differentially expressed gene
#'    when HER2 pathway is active.}
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
#' site of relapse. Cancer research, 68(9), pp.3108-3114. doi:
#' \url{https://doi.org/10.1158/0008-5472.CAN-07-5644}
"HER2_sig_df"

#' Gene signature for ER and HER2 extracted from PAM50
#'
#' A list of two data frames containing the gene names and their corresponding
#' expression values constituting ER and HER gene signatures obtained from
#' PAM50 gene signature (Parker et al. 2009). The ER signature was obtained by
#' combining the list of genes for the Luminal A and B intrinsic subtypes
#' (centroids) predicted by PAM50 classifier in "pam50" list featured in the
#' genefu R BiocManager package. Genes that were the most informative as
#' demonstrated by a highly deviant centroid index for either Luminal A or B
#' (estrogen receptor positive breast cancer subtypes) were selected as part
#' of the ER signature for PAM50. Similarly, genes that were highly deviant in
#' their centroid index in the HER2 intrinsic subtype (centroid) from PAM50 were
#' selected as the HER2 signature for PAM50 using pam50 from genefu.
#' @format A list holding two data frames: one for ER signature, the other for
#' HER2 signature extracted from the PAM50 gene signature:
#' \describe{
#'    \item{ER}{data frame containing HUGO gene symbols and their relative
#'    expression value (-1 for down-regulated genes and 1 for up-regulated
#'    genes) in the ER signature}
#'    \item{HER2}{data frame containing gene names their relative expression
#'    value (-1 for down-regulated genes and 1 for up-regulated genes) in the
#'    HER2 signature}
#' }
#' @references
#' Parker, J.S., Mullins, M., Cheang, M.C., Leung, S., Voduc, D., Vickery, T.,
#' Davies, S., Fauron, C., He, X., Hu, Z. and Quackenbush, J.F., 2009.
#' Supervised risk predictor of breast cancer based on intrinsic subtypes.
#' Journal of clinical oncology, 27(8), p.1160. doi:
#' \url{https://doi.org/10.1200/JCO.2008.18.1370}
#' @source PAM50 publication: \url{https://pubmed.ncbi.nlm.nih.gov/19204204/}).
#' @source genefu: \url{https://bioconductor.org/packages/release/bioc/html/genefu.html}
"pam50"

#' ER RNA-seq gene expression data set from TCGA
#'
#' A gene expression matrix containing RNA-seq raw read counts for 60 human
#' primary breast tumour samples (30 estrogen receptor (ER) positive and 30 ER
#' negative samples were selected at random). This data set is a subset of a
#' much larger data set containing 1,101 primary breast tumour samples collected
#' from The Cancer Genome Atlas (TCGA).
#' @format  A matrix containing 20,124 HUGO gene symbols (row names) and 60 breast
#' cancer tumour samples IDs (columns) given in the form of TCGA barcodes for
#' each sample. For further information on TCGA bar code semantics, please see
#' the NIH GDC documentation
#' \url{https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/}.
#' @source \url{https://portal.gdc.cancer.gov}).
"ER_TCGA_RNAseq"

#' HER2 RNA-seq gene expression data set from TCGA
#'
#' A gene expression matrix containing RNA-seq raw read counts for 60 primary
#' human breast tumour samples (30 human epidermal growth receptor (HER2)
#' positive and 30 HER2 negative samples were selected at random). This data set
#' is a subset of a much larger data set containing 1,101 primary breast tumour
#' samples collected from The Cancer Genome Atlas (TCGA).
#' @format A matrix containing 20,124 HUGO gene symbols (row names) and 60 breast
#' cancer tumour samples IDs (columns) given in the form of TCGA barcodes for
#' each sample. For further information on TCGA bar code semantics, please see
#' the NIH GDC documentation
#' \url{https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/}.
#' @source TCGA \url{https://portal.gdc.cancer.gov}).
"HER2_TCGA_RNAseq"

#' ER microarray data set obtained from GEO
#'
#' A matrix containing microarray data for 60 human breast tumour samples (30
#' estrogen receptor (ER) positive and 30 ER negative samples were selected at
#' random) and were extracted from GSE31448 dataset (Sabatier et al. 2011) using
#' the GEO query library. Multi-gene probes were excluded and only cancer samples
#' from human breast cancer tumours were selected for constructing this matrix.
#' @format  A matrix containing 21,656 HUGO gene symbols (row names) and 60
#' breast cancer tumour samples IDs (columns) which are given as sample
#' accession numbers.
#' @references
#' Sabatier, R., Finetti, P., Adelaide, J., Guille, A., Borg, J.P.,
#' Chaffanet, M., Lane, L., Birnbaum, D. and Bertucci, F., 2011. Down-regulation
#' of ECRG4, a candidate tumor suppressor gene, in human breast cancer.
#' PloS one, 6(11), p.e27656. doi: https://doi.org/10.1371/journal.pone.0027656
#' @source GSE31448 series link on GEO:
#' \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31448}).
"ER_GEO_microarr"

#' HER2 microarray data set obtained from GEO
#'
#' A matrix containing microarray data for 60 human breast tumour samples (30
#' ERBB2 (HER2) positive and 30 HER2 negative samples were selected at random)
#' and were extracted from GSE31448 dataset (Sabatier et al. 2011) using the GEO
#' query library.Multi-gene probes were excluded and only cancer samples from
#' human breast cancer tumours were selected for constructing this matrix.
#' @format A matrix containing 21,656 HUGO gene symbols (row names) and 60
#' breast cancer tumour samples IDs (columns) which are given as sample
#' accession numbers.
#' @references
#' Sabatier, R., Finetti, P., Adelaide, J., Guille, A., Borg, J.P.,
#' Chaffanet, M., Lane, L., Birnbaum, D. and Bertucci, F., 2011. Down-regulation
#' of ECRG4, a candidate tumor suppressor gene, in human breast cancer.
#' PloS one, 6(11), p.e27656. doi: https://doi.org/10.1371/journal.pone.0027656
#' @source GSE31448 series link on GEO:
#' \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31448}).
"HER2_GEO_microarr"

