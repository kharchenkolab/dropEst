#' Results of processing of SRR3879617 dataset, neccesary for low-quality cells filtration.
#'
#' @format A list with the next fields:
#' \describe{
#'   \item{price}{price, in US dollars}
#'   \item{carat}{weight of the diamond, in carats}
#'   ...
#' }
#' @source \url{https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR3879617}
"lq_cells_info"

#' Results of processing of 10x Frozen BMMCs (Healthy Control 1) dataset, neccesary for UMI errors correction.
#'
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/frozen_bmmc_healthy_donor1}
"reads_per_umi_per_cb"

#' Results of saturation analysis of SRR1784312 dataset (mouse ES cells).
#'
#' @format A list with the next fields:
#' \describe{
#'   \item{sat}{saturation curve}
#'   \item{cur}{observed values}
#'   ...
#' }
#' @source \url{https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1784312}
"saturation_srr1784312"
