##' Create an mcstate "initial" function from a previous final state
##' matrix
##'
##' @title Create restart initial conditions
##'
##' @param state Final state at the point of the restart (matrix with
##'   rows corresponding to model state and columns corresponding to
##'   samples)
##'
##' @param date The sircovid date of the restart; this will be
##'   converted into the initial step.
##'
##' @export
spim_restart_initial <- function(state, date) {
  n_state <- ncol(state)
  function(info, n_particles, pars) {
    i <- sample.int(n_state, n_particles, replace = TRUE)
    list(state = state[, i], step = date * pars$steps_per_day)
  }
}
