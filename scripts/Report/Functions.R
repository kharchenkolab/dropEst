library(org.Mm.eg.db)
library(GO.db)
library(EMCluster)

get_otsu_threshold <- function(probs) {
    probs_sorted <- sort(probs)

    total <- length(probs_sorted)
    w_sum <- sum(probs_sorted)

    sum_b <- 0
    w_b <- 0
    w_f <- 0
    var_max <- 0
    threshold <- 0

    t <- probs_sorted[1]
    for (t in probs_sorted) {
        w_b <- w_b + 1
        w_f <- total - w_b
        if (w_f == 0) break
        sum_b <- sum_b + t

        m_b <- sum_b / w_b
        m_f <- (w_sum - sum_b) / w_f

        var_between <- w_b * w_f * (m_b - m_f) ** 2
        if (var_between > var_max) {
          var_max <- var_between
          threshold <- t
        }
    }

    return(threshold)
}


get_genesets <- function(genes) {
    # translate gene names to ids
    ids <- unlist(lapply(mget(genes, org.Mm.egALIAS2EG, ifnotfound = NA), function(x) x[1]))
    rids <- names(ids); names(rids) <- ids
    # convert GO lists from ids to gene names
    go.env <- lapply(mget(ls(org.Mm.egGO2ALLEGS), org.Mm.egGO2ALLEGS), function(x) as.character(na.omit(rids[x])));
    names(go.env) <- paste(names(go.env),unlist(lapply(mget(names(go.env),GOTERM),function(x) x@Term)))

    genesets <- list(
    apoptotic=unique(go.env$`GO:0006915 apoptotic process`),
    cytoplasm=unique(go.env$`GO:0005737 cytoplasm`),
    extracellular=unique(go.env$`GO:0005576 extracellular region`),
    membrane=unique(go.env$`GO:0016020 membrane`),
    metabolic=unique(go.env$`GO:0008152 metabolic process`),
    ribosome=unique(go.env$`GO:0005840 ribosome`),
    mitochondrion=unique(go.env$`GO:0005739 mitochondrion`)
    )
    return(genesets)
}

get_em <- function(fraction, good_number) {
    labels <- c(rep(1, good_number), rep(0, length(fraction) - good_number))
    em_1_class_res <- init.EM(matrix(fraction), nclass = 1, method = "Rnd.EM")
    em_res <- init.EM(matrix(fraction), nclass = 2, lab = labels, method = "Rnd.EM")

    if (em_res$nc[2] / em_res$nc[1] < 0.1) return(NULL)
    if (em_res$LTSigma[1] ** 0.5 * 2 > abs(em_res$Mu[1] - em_res$Mu[2])) return(NULL)

    return(em_res)
}

arr_deriv <- function(x, y, lag=1) {diff(y, lag) / diff(x, lag)}

logspace <- function(low, up, len) { round(exp(log(10)*seq(log10(low),log10(up),by=(log10(up) - log10(low)) / (len - 1)))) }

get_cells_number <- function(umis_counts, lag) {
    log_umis_counts <- log(umis_counts)
    log_rank <- log(1:length(umis_counts))

    x <- log_rank[(1+lag):length(log_rank)]; y <- arr_deriv(log_rank, log_umis_counts, lag)
    x2 <- x[(1+lag):length(x)]; y2 <- arr_deriv(x, y, lag)

    lens <- rle(as.vector(y2 >= 0))$lengths
    start_inds <- cumsum(lens)
    cells_number=start_inds[which(lens == max(lens[y2[start_inds] > 0])) - 1] + lag
    colors <- c(rep('Good', cells_number), rep('Bad', length(umis_counts) - cells_number))

    return(list(cells_number=cells_number, colors=colors, inds=(1+lag):(length(umis_counts)-lag)))
}
