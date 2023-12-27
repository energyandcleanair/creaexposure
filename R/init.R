init <- function(){
  readRenviron(".Renviron")
  dir.create("cache", showWarnings = F)
  dir.create("results", showWarnings = F)
}
