.onUnload <- function (libpath) {
  library.dynam.unload("dropestr", libpath)
}
