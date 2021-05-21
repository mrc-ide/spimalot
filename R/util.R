`%||%` <- function(x, y) { # nolint
  if (is.null(x)) y else x
}


data_frame <- function(...) {
  data.frame(..., stringsAsFactors = FALSE)
}


read_csv <- function(...) {
  read.csv(..., stringsAsFactors = FALSE, check.names = FALSE)
}


last <- function(x) {
  x[[length(x)]]
}


squote <- function(x) {
  sprintf("'%s'", x)
}


spimalot_file <- function(...) {
  system.file(..., package = "spimalot", mustWork = TRUE)
}
