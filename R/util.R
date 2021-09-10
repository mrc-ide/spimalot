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


add_dimension <- function(x, d, nm = list(NULL)) {
  old_dim <- dim(x)
  old_nms <- dimnames(x)
  if (is.null(old_nms)) {
    old_nms <- rep(list(NULL), length(old_dim))
  }

  if (d < 1) {
    stop("d < 1")
  } else if (d > length(old_dim) + 1) {
    stop("d > dim(x) + 1")
  } else if (d > length(old_dim)) {
    new_dim <- c(old_dim, 1)
    new_nms <- c(old_nms, nm)
  } else if (d == 1) {
    new_dim <- c(1, old_dim)
    new_nms <- c(nm, old_nms)
  } else {
    new_dim <- c(old_dim[seq(d - 1)], 1, old_dim[seq(d, length(old_dim))])
    new_nms <- c(old_nms[seq(d - 1)], nm, old_nms[seq(d, length(old_nms))])
  }

  if (all(vlapply(new_nms, is.null))) {
    new_nms <- NULL
  }

  array(x, new_dim, new_nms)
}


drop_dimension <- function(x, i) {
  array(x, dim(x)[-i], dimnames(x)[-i])
}