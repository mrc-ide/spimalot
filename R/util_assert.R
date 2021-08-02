assert_has_names <- function(x, required, name = deparse(substitute(x))) {
  msg <- setdiff(required, names(x))
  if (length(msg) > 0) {
    stop(sprintf("Required names missing from '%s': %s",
                 name, paste(squote(msg), collapse = ", ")))
  }
  invisible(x)
}


assert_length <- function(x, len, name = deparse(substitute(x))) {
  if (length(x) != len) {
    stop(sprintf("'%s' must have length %d", name, len))
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


assert_scalar <- function(x, name = deparse(substitute(x))) {
  if (length(x) != 1L) {
    stop(sprintf("'%s' must be a scalar", name), call. = FALSE)
  }
  invisible(x)
}


assert_scalar_character <- function(x, name = deparse(substitute(x))) {
  assert_scalar(x, name)
  assert_character(x, name)
}


assert_scalar_logical <- function(x, name = deparse(substitute(x))) {
  force(name)
  assert_scalar(x, name)
  assert_logical(x, name)
  invisible(x)
}


assert_scalar_numeric <- function(x, name = deparse(substitute(x))) {
  assert_scalar(x, name)
  assert_numeric(x, name)
}


assert_numeric <- function(x, name = deparse(substitute(x))) {
  if (!is.numeric(x)) {
    stop(sprintf("'%s' must be numeric", name), call. = FALSE)
  }
}


assert_logical <- function(x, name = deparse(substitute(x))) {
  if (!(is.logical(x))) {
    stop(sprintf("'%s' must be a logical", name), call. = FALSE)
  }
  invisible(x)
}


assert_character <- function(x, name = deparse(substitute(x))) {
  if (!is.character(x)) {
    stop(sprintf("'%s' must be character", name), call. = FALSE)
  }
}


assert_scalar_positive_integer <- function(x, name = deparse(substitute(x))) {
  force(name)
  assert_scalar(x, name)
  x <- assert_integer(x, name)
  if (x < 1L) {
    stop(sprintf("'%s' must be at least 1", name), call. = FALSE)
  }
  invisible(x)
}


assert_integer <- function(x, name = deparse(substitute(x)),
                           what = "integer") {
  if (!(is.integer(x))) {
    eps <- sqrt(.Machine$double.eps)
    usable_as_integer <- is.numeric(x) && (max(abs(round(x) - x)) < eps)
    if (!usable_as_integer) {
      stop(sprintf("'%s' must be an %s", name, what), call. = FALSE)
    }
    x <- as.integer(round(x))
  }
  invisible(x)
}


match_value <- function(arg, choices, name = deparse(substitute(arg))) {
  assert_scalar_character(arg)
  if (!(arg %in% choices)) {
    stop(sprintf("%s must be one of %s", name,
                 paste(squote(choices),  collapse = ", ")), call. = FALSE)
  }
  arg
}
