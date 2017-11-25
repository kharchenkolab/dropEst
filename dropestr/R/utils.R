.onUnload <- function (libpath) {
  library.dynam.unload("dropestr", libpath)
}

Rcpp::loadModule("CppMapModule", TRUE)
