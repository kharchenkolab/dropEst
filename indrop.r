
##' Draw a panel of diagnostic plots for an indrop dataset
##'
##' @title Diagnostic plots for indrop dataset
##' @param x dataset read in using read.indropest()
##' @param n.cells max number of cells to show in the cell number estimate analysis
##' @param merge.threshold maximum merge probability threshold after which cell barcodes are considered false (default p=0.08)
basic.plots <- function(x,n.cells=min(ncol(x$cm)*2,length(x$umig.cov)),merge.threshold=0.08) {
  # start with the raw output - NOTE: this should probably be made into a matrix in the C code
  cm <- t(matrix(x$cm,nrow=length(x$cell.names)))
  rownames(cm) <- x$gene.names;
  colnames(cm) <- x$cell.names;
  
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
