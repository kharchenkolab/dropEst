# #' Results of processing of SRR3879617 dataset, neccesary for low-quality cells filtration.
# #'
# #' @format A list with the next fields:
# #' \describe{
# #'   \item{price}{price, in US dollars}
# #'   \item{carat}{weight of the diamond, in carats}
# #'   ...
# #' }
# #' @source \url{https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR3879617}
# "lq_cells_info"

#' Results of processing of mouse BC dataset (SCG71), neccesary for UMI errors correction for the 100 largest cells.
#'
"reads_per_umi_per_cell"

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
