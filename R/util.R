`%||%` <- function(x, y) { # nolint
  if (is.null(x)) y else x
}


data_frame <- function(...) {
  data.frame(..., stringsAsFactors = FALSE, check.names = FALSE)
}


read_csv <- function(...) {
  read.csv(..., stringsAsFactors = FALSE, check.names = FALSE)
}


write_csv <- function(...) {
  write.csv(..., row.names = FALSE)
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


list_transpose <- function(x) {
  nms <- lapply(x, names)
  stopifnot(length(unique(nms)) == 1L)
  ret <- lapply(nms[[1]], function(el) lapply(x, "[[", el))
  names(ret) <- nms[[1]]
  ret
}


vnapply <- function(X, FUN, ...) {
  vapply(X, FUN, numeric(1), ...)
}


vlapply <- function(X, FUN, ...) {
  vapply(X, FUN, logical(1), ...)
}


abind_quiet <- function(...) {
  suppressWarnings(abind::abind(...))
}


set_names <- function(x, nms) {
  names(x) <- nms
  x
}


nlayers <- function(x) {
  dim(x)[[3L]]
}


str_collapse <- function(str) {
  sprintf("{%s}", paste0(str, sep = ", "))
}


seq_rows <- function(x) {
  seq_len(nrow(x))
}


mean_ci <- function(x) c(mean = mean(x), quantile(x, c(0.025, 0.975)))


immutable_error <- function(x) {
  stop(sprintf(
    "Objects of class '%s' are immutable and you may not alter them",
    class(x)[[1]]))
}

##' @export
`[[<-.immutable` <- function(x, i, ..., value) {
  immutable_error(x)
}


##' @export
`[<-.immutable` <- function(x, i, ..., value) {
  immutable_error(x)
}


##' @export
`$<-.immutable` <- function(x, name, value) {
  immutable_error(x)
}
