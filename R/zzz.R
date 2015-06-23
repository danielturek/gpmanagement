
.onAttach <- function(libname, pkgname) {
  nfFile1 <- file.path(path.package('gpmanagement'), 'inst', 'include', 'nfs.R')
  nfFile2 <- file.path(path.package('gpmanagement'),         'include', 'nfs.R')
  if(file.exists(nfFile1)) {
      source(nfFile1)
      return(invisible(NULL))
  } else if(file.exists(nfFile2)) {
      source(nfFile2)
      return(invisible(NULL))
  } else stop('cannot find nfs.R')
}


