`%||%` <- function(x, y) { # nolint
  if (is.null(x)) y else x
}


data_frame <- function(...) {
  data.frame(..., stringsAsFactors = FALSE, check.names = FALSE)
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


## See mcstate:::abind3 for inspiration
abind1 <- function(a, b) {
  na <- dim(a)[1]
  nb <- dim(b)[1]
  nab <- dim(a)[2:3]
  ret <- array(NA_real_, c(na + nb, nab))
  ret[seq_len(na), , ] <- a
  ret[seq_len(nb) + na, , ] <- b
  rownames(ret) <- c(rownames(a), rownames(b))
  ret
}
