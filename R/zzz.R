.onLoad <- function(libname, pkgname) {
  op <- options()
  op.devtools <- list(
    devtools.path = "~/git/R-dev",
    devtools.install.args = "",
    devtools.name = "Erik Schutte",
    devtools.desc.author = '"Erik Schutte <schutte.erik@hotmail.com> [aut, cre]"',
    devtools.desc.license = "GPL-2",
    devtools.desc.suggests = NULL,
    devtools.desc = list()
  )
  toset <- !(names(op.devtools) %in% names(op))
  if(any(toset)) options(op.devtools[toset])

  invisible()
}
