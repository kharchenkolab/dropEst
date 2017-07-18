
##' Draw a panel of diagnostic plots for an indrop dataset
##'
##' @title Diagnostic plots for indrop dataset
##' @param x dataset read in using read.indropest()
##' @param n.cells max number of cells to show in the cell number estimate analysis
##' @param merge.threshold maximum merge probability threshold after which cell barcodes are considered false (default p=0.08)
basic.plots <- function(x,n.cells=min(ncol(x$cm)*2,length(x$umig.cov)),merge.threshold=0.08) {
  # start with the raw output - NOTE: this should probably be made into a matrix in the C code
  #cm <- t(matrix(x$cm,nrow=length(x$cell.names)))
  #rownames(cm) <- x$gene.names;
  #colnames(cm) <- x$cell.names;
  cm <- x$cm;
  
  df <- data.frame(rank=c(1:min(length(x$merge.n),n.cells)),s=(sign(rev(x$merge.n))==-1)[1:min(length(x$merge.n),n.cells)]);
  m <- glm(cbind(s==1,s==0)~rank,family=binomial(logit),data=df)
  if(max(m$fitted<merge.threshold)) {
    ti <- length(m$fitted)
  } else {
    ti <- which(m$fitted>merge.threshold)[1];
  }

  #l <- layout(matrix(c(seq(1,5),5),nrow=3,byrow=T));
  #par(mar = c(3.5,3.5,1.0,1.0), mgp = c(2,0.65,0),cex=1)
  par(mfrow=c(3,2), mar = c(3.5,3.5,1.0,1.0), mgp = c(2,0.65,0),cex=1)
  hist(log10(colSums(x$cm)),xlab="log10( UMIs/cell )",main="",col="wheat")
  hist(log10(rowSums(x$cm)),xlab="log10( UMIs/gene )",main="",col="wheat")
  smoothScatter(colSums(x$cm),x$rpu,xlab="UMIs/cell",ylab="reads/UMI")
  plot((cumsum(x$umig.cov))[1:n.cells],type='l',xlab="cell rank",ylab="number of unique UMI+g")
  abline(v=ncol(x$cm),col=2,lty=2)
  abline(v=ti,lty=2,col=8)
  #barplot(x$exonic.chr,las=2,main="exonic reads")
  #barplot(x$nonexonic.chr,las=2,main="non-exonic reads")
  
  barplot(rbind(colSums(x$nonex_cells_chr_counts),colSums(x$ex_cells_chr_counts)),col=c("gray50","blue"),las=3,ylab="reads")
  legend(x="top",fill=c("gray50","blue"),legend=c("non-exonic","exonic"),horiz=T,bty="n")
  
  smoothScatter(df$s,xlab="cell rank",ylab="merge p",bandwidth=0.01,ylim=c(0,1))
  lines(df$rank,m$fitted,col=2)
  abline(v=ti,lty=2,col=8)
  abline(v=ncol(x$cm),lty=2,col=2)
  legend(x="bottomright",bty='n',col=8,lty=2,legend=paste("cell",ti))

}

get.preseq <- function(reads_by_class,max.size=NULL,steps=100) {
  if(is.null(max.size)) { max.size=10*sum(counts) }
  counts <- table(reads_by_class)
  freqs <- as.integer(names(counts));
  reads_count <- sum(reads_by_class)
  if(is.null(max.size)) { max.size <- read_count*10 }

  x <- preseqR.pf.mincount(ss=round(max.size / steps),n=cbind(freqs, as.vector(counts)),max.extrapolation=max.size,mt=200)
  list(sat=x$yield.estimates[[1]],cd=reads_count,cy=length(reads_by_class))
}

preseqr.plot.resl <- function(resl,max.size=max(unlist(lapply(resl,function(d) max(d$sat[,1])))),cols=rainbow(length(resl)),lwd=2,legendx="topleft") {
  maxy <- max(unlist(lapply(resl,function(d) max(d$sat[d$sat[,1]<=max.size,2]))))
  #par(mar = c(3.5,3.5,2.0,1.5), mgp = c(2,0.65,0), cex = 0.95);
  plot(c(),c(),xlab="sequencing depth",ylab="unique molecules",xlim=c(0,max.size),ylim=c(0,maxy),xaxs="i",yaxs="i")

  lapply(1:length(resl),function(i) {
    d <- resl[[i]];
    vi <- d$sat[,1] <= d$cd;
    lines(d$sat[vi,],col=cols[i],lty=1,lwd=lwd);
    vi <- d$sat[,1] > d$cd;
    if(any(vi)) {
      lines(d$sat[vi,],col=cols[i],lty=2,lwd=lwd);
    }
    points(d$cd,d$cy,pch=19,col=cols[i])
  })
  legend(x=legendx,lty=rep(1,length(resl)),col=cols,lwd=rep(lwd,length(resl)),legend=names(resl),bty='n')
}


new.plots <- function(x,top.cells=1000,n.breaks=50,max.extrapolation=1e8) {
  # size curve
  arr_deriv <- function(x, y, lag=1) {diff(y, lag) / diff(x, lag)}
  get.cell.number <- function(umis_counts, lag) {
    log_umis_counts <- log(umis_counts)
    log_rank <- log(1:length(umis_counts))
    
    x <- log_rank[(1+lag):length(log_rank)]; y <- arr_deriv(log_rank, log_umis_counts, lag)
    x2 <- x[(1+lag):length(x)]; y2 <- arr_deriv(x, y, lag)
    
    lens <- rle(as.vector(y2 >= 0))$lengths
    start_inds <- cumsum(lens)
    return(start_inds[which(lens == max(lens[y2[start_inds] > 0])) - 1] + lag)
  }
  umi.counts <- sort(table(x$reads_by_umig_cbs),decreasing=T);
  
  n.cells <- get.cell.number(umi.counts, min(100, as.integer(0.1 * length(umi.counts))))[1]
  
  par(mfrow=c(3,2), mar = c(3.5,3.5,1.0,1.0), mgp = c(2,0.65,0),cex=1)
  plot(log10(1:length(umi.counts)),log10(as.integer(umi.counts)),type='l',lwd=2,xlab="log10[ cell rank ]",ylab="log10[ UMIs ]",panel.first=grid());
  abline(v=log10(n.cells),col=2,lty=2); legend(x='topright',legend=paste("N =",n.cells),bty='n')

  # histogram
  h <- hist(log10(umi.counts),breaks=n.breaks,plot=F)
  y <- h$counts*(10^h$mids); y[y<0] <- 0;
  plot(c(),c(),xlab="log10[ UMIs ]",ylab="UMIs * # of cells",ylim=c(0,max(y)),xlim=range(h$mids))
  rect(h$breaks[-length(h$breaks)],0,h$breaks[-1],y,col='wheat')
  abline(v=log10(umi.counts[n.cells]),col=2,lty=2)
         
  
  # saturation
  require(preseqR)

  vb <- names(umi.counts)[1:top.cells]
  sat <- get.preseq(x$reads_by_umig[x$reads_by_umig_cbs %in% vb],max.size=max.extrapolation)
  asat <- readRDS("/home/pkharchenko/drop/cp/scripts/SRR1784321.sat.rds")
  
  preseqr.plot.resl(list("this"=sat,"mES"=asat),cols=c(1,8))
  
  # chromsome distribution
  barplot(rbind(colSums(x$nonex_cells_chr_counts),colSums(x$ex_cells_chr_counts)),col=c("gray50","blue"),las=3,ylab="reads")
  legend(x="top",fill=c("gray50","blue"),legend=c("non-exonic","exonic"),horiz=T,bty="n")
  
  # UMI scatter
  reads.per.cell <- tapply(x$reads_by_umig,as.factor(x$reads_by_umig_cbs),sum)
  #smoothScatter(log10(1:length(umi.counts)),log10((reads.per.cell[names(umi.counts)]/as.integer(umi.counts))))
  smoothScatter(log10(1:length(umi.counts)),log10((reads.per.cell[names(umi.counts)]/as.integer(umi.counts))),nrpoints=1e3,xlab="log10[ cell rank ]",ylab="log10[ reads/UMI ]")
  
  # genes per (selected) cell
  hist(log10(colSums(x$cm>0)),xlab="log10[ number of genes ]",col='wheat',main="")
  
}

test.int <- function() {
  # diagnostic plots

  
  CairoPNG(file=paste(fname,"info.png",sep="."),width=700,height=800);
  par(mfrow = c(3,2), mar = c(3.5,3.5,1.0,1.0), mgp = c(2,0.65,0),cex=1)
  hist(log10(colSums(x$cm)),xlab="log10( UMIs/cell )",main="",col="wheat")
  hist(log10(rowSums(x$cm)),xlab="log10( UMIs/gene )",main="",col="wheat")
  smoothScatter(colSums(x$cm),x$rpu,xlab="UMIs/cell",ylab="reads/UMI")
  plot((cumsum(x$umig.cov))[1:min(ncol(x$cm)*2,length(x$umig.cov))],type='l',xlab="cell rank",ylab="number of unique UMI+g")
  abline(v=ncol(x$cm),col=2,lty=2)
  #barplot(x$exonic.chr,las=2,main="exonic reads")
  #barplot(x$nonexonic.chr,las=2,main="non-exonic reads")

  par(mar = c(4,5,1.0,1.0), mgp = c(2,0.65,0),cex=1)
  barplot(rbind(x$nonexonic.chr,x$exonic.chr[names(x$nonexonic.chr)]),col=c("gray50","blue"),las=2,ylab="reads")
  legend(x="top",fill=c("gray50","blue"),legend=c("non-exonic","exonic"),horiz=T,bty="n")
  dev.off();
  
  
  
  plot(colSums(cm>0),x$rpu)
  barplot(ec,las=2)
  barplot(nonec,las=2)
  str(cm)
  
  plot(cumsum(x$umig.cov))
}
