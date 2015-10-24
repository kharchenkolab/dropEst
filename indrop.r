require(Rcpp);
require(inline)

# routines for reading in binary output of indropest
indrop.compile.plug <- Rcpp:::Rcpp.plugin.maker(include.before = '',
                                         libs = paste(
                                                      "-lz -lm -lbam -lboost_iostreams -lboost_serialization"))

registerPlugin("indropCompilePlugin",indrop.compile.plug)

settings=getPlugin("indropCompilePlugin")
settings$env$PKG_CXXFLAGS='-I /home/pkharchenko/drop/cp'

internal.read.indropest <- cxxfunction(signature(Fname="character"),
                              includes='#include <string>
                                        #include <boost/archive/binary_iarchive.hpp>
                                        #include <fstream>
                                        #include "indrop_results.h"',
                              body='
      using namespace std;
      string fname=Rcpp::as<string>(Fname);
      ifstream ifs(fname.c_str());
      if(!ifs) { stop("cannot open specified file"); }
      boost::archive::binary_iarchive ia(ifs);
      indrop_results res;
      ia >> res;
      return List::create(Named("cell.names") = wrap(res.cm.cell_names),
                          Named("gene.names") = wrap(res.cm.gene_names),
                          Named("cm") = wrap(res.cm.counts),
                          Named("none_c") = wrap(res.non_exon_counts),
                          Named("none_chr") = wrap(res.non_exon_count_names),
                          Named("e_c") = wrap(res.exon_counts),
                          Named("e_chr") = wrap(res.exon_count_names),
                          Named("rpu") = wrap(res.reads_per_umi),
                          Named("umig.cov") = wrap(res.umig_covered),
                          Named("merge.n") = wrap(res.merge_n),
                          Named("fname") = wrap(fname));
   ',plugin="indropCompilePlugin",settings=settings)

read.indropest <- function(fname) {
  x <- internal.read.indropest(fname)
  cm <- t(matrix(x$cm,nrow=length(x$cell.names)))
  rownames(cm) <- x$gene.names;
  colnames(cm) <- x$cell.names;
  nonec <- x$none_c; names(nonec) <- x$none_chr;
  ec <- x$e_c; names(ec) <- x$e_chr;
  
  return(list(cm=cm,umig.coverage=x$umig.cov,rpu=x$rpu,exonic.chr=ec,nonexonic.chr=nonec,mergen=x$merge.n,name=x$fname))
}


##' Draw a panel of diagnostic plots for an indrop dataset
##'
##' @title Diagnostic plots for indrop dataset
##' @param x dataset read in using read.indropest()
##' @param n.cells max number of cells to show in the cell number estimate analysis
##' @param merge.threshold maximum merge probability threshold after which cell barcodes are considered false (default p=0.08)
basic.plots <- function(x,n.cells=min(ncol(x$cm)*2,length(x$umig.cov)),merge.threshold=0.08) {
  df <- data.frame(rank=c(1:n.cells),s=(sign(rev(x$mergen))==-1)[1:n.cells]);
  m <- glm(cbind(s==1,s==0)~rank,family=binomial(logit),data=df)
  ti <- which(m$fitted>merge.threshold)[1]

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
  
  #par(mar = c(4,5,1.0,1.0), mgp = c(2,0.65,0),cex=1)
  barplot(rbind(x$nonexonic.chr,x$exonic.chr[names(x$nonexonic.chr)]),col=c("gray50","blue"),las=3,ylab="reads")
  legend(x="top",fill=c("gray50","blue"),legend=c("non-exonic","exonic"),horiz=T,bty="n")

  smoothScatter(df$s,xlab="cell rank",ylab="merge p")
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
