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
##' @param multistrain Logical, whether the restart fits are multistrain or not
##'
##' @export
spim_restart_initial <- function(state, date, multistrain) {
  e <- new.env(parent = topenv())
  e$state <- state
  e$date <- sircovid::as_sircovid_date(date)
  e$multistrain <- multistrain
  with(e, {
    n_state <- ncol(state)
    function(info, n_particles, pars) {
      if (multistrain) {
        ## Then seed strain 2 at some fraction
        compartment <- c("E", "I_A", "I_P", "I_C_1")
        for (i in seq_along(compartment)) {
          dim <- info$dim[[compartment[[i]]]]
          idx <- info$index[[compartment[[i]]]]
          new_state <- array(state[idx, ], c(dim, ncol(state)))

          if (length(dim) != 4) {
            stop(sprintf(
              "Unexpected dimensions (%d) in move_strain_compartment",
              length(dim)))
          }

          n <- new_state[, 1, , , ]
          new_state[, 2, , , ] <- rbinom(length(n), n, pars$prop_strain_2)
          new_state[, 1, , , ] <- n - new_state[, 2, , , ]

          state[idx, ] <- new_state
        }
      }

      i <- sample.int(n_state, n_particles, replace = TRUE)
      list(state = state[, i], step = date * pars$steps_per_day)
    }
  })
}
