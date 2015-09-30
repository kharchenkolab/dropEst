require(Rcpp);
require(inline)

# routines for reading in binary output of indropest
indrop.compile.plug <- Rcpp:::Rcpp.plugin.maker(include.before = '',
                                         libs = paste("-lz -lm -lbam -lboost_iostreams -lboost_serialization"))

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
                          Named("fname") = wrap(fname));
   ',plugin="indropCompilePlugin",settings=settings)

read.indropest <- function(fname) {
  x <- internal.read.indropest(fname)
  cm <- t(matrix(x$cm,nrow=length(x$cell.names)))
  rownames(cm) <- x$gene.names;
  colnames(cm) <- x$cell.names;
  nonec <- x$none_c; names(nonec) <- x$none_chr;
  ec <- x$e_c; names(ec) <- x$e_chr;

  return(list(cm=cm,umig.coverage=x$umig.cov,rpu=x$rpu,exonic.chr=ec,nonexonic.chr=nonec))
}

test.int <- function() {
  fname <- "cell.counts.txt.bin"
  x <- read.indropest(fname)

  # diagnostic plots

  require(Cairo)
  
  par(mfrow = c(2,2), mar = c(3.5,3.5,1.0,1.0), mgp = c(2,0.65,0))
  smoothScatter(colSums(cm),x$rpu,xlab="UMIs/cell",ylab="reads/UMI")
  plot((cumsum(x$umig.cov))[1:3000],type='l',xlab="cell rank",ylab="number of unique UMI+g")
  abline(v=ncol(cm),col=2,lty=2)
  barplot(x$exonic.chr,las=2,main="exonic reads")
  barplot(x$nonexonic.chr,las=2,main="non-exonic reads")
  
  
  
  plot(colSums(cm>0),x$rpu)
  barplot(ec,las=2)
  barplot(nonec,las=2)
  str(cm)
  
  plot(cumsum(x$umig.cov))
}
