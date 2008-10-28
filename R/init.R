.onLoad <- function(libname, pkgname) {
  s <- search() 
  library.dynam("flagme",pkgname,libname,now=FALSE)
}

