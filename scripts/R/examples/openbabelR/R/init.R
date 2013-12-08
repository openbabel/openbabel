
.onLoad <- function(libname,pkgname) {

	library.dynam(pkgname,package=pkgname,lib.loc=libname,local=FALSE)

}
