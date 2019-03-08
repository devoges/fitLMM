#' @importFrom GMMAT glmmkin
fit_gmmat <- function(Y,
                      X,
                      Slist) {
  data <- as.data.frame(X)
  data$Y <- Y

  glmmkin(fixed = Y ~ .,
          data = data,
          family = gaussian(),
          kins = Slist)
}
