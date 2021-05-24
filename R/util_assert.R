assert_has_names <- function(x, required, name = deparse(substitute(x))) {
  msg <- setdiff(required, names(x))
  if (length(msg) > 0) {
    stop(sprintf("Required names missing from '%s': %s",
                 name, paste(squote(msg), collapse = ", ")))
  }
  invisible(x)
}


assert_setequal <- function(a, b,
                            name_a = deparse(substitute(a)),
                            name_b = deparse(substitute(b)),
                            unique = TRUE) {
  if (unique) {
    assert_unique(a, name_a)
    assert_unique(b, name_b)
  }
  if (!setequal(a, b)) {
    ## TODO: describe difference here
    stop(sprintf("'%s' and '%s' are not setequal", name_a, name_b))
  }
}


assert_unique <- function(x, name = deparse(substitute(x))) {
  if (anyDuplicated(x)) {
    ## TODO: describe duplicates here
    stop(sprintf("'%s' contains duplicate elements", name))
  }
}


assert_is <- function(x, what, name = deparse(substitute(x))) {
  if (!inherits(x, what)) {
    stop(sprintf("'%s' must be a %s", name, paste(what, collapse = " / ")),
         call. = FALSE)
  }
  invisible(x)
}


assert_file_exists <- function(path, dir = NULL, name = "File") {
  if (is.null(dir)) {
    if (!file.exists(path)) {
      stop(sprintf("%s '%s' does not exist", name, path))
    }
  } else {
    if (!file.exists(file.path(dir, path))) {
      stop(sprintf("%s '%s' does not exist (relative to '%s')",
                   name, path, dir))
    }
  }
}
