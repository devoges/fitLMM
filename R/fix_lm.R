det <- function(x) {
  r <- determinant(x)$modulus
  attr(r, "logarithm") <- NULL
  r
}

#' @export
fit_lm <- function(Y,
                   X,
                   REML = TRUE) {
  model <- lm(Y ~ -1 + .,
              data = as.data.frame(X))
  X <- as.matrix(X)
  beta <- model$coefficients
  r <- Y - X %*% beta
  rss <- {t(r) %*% r}[1, 1]
  if (!REML) {
    vcs <- c("Error" = rss / nrow(X))
  } else {
    vcs <- c("Error" = rss / (nrow(X) - ncol(X)))
  }

  lnum <- log(2 * pi * rss)
  d <- 0
  dd <- det(t(X) %*% X)

  ll_REML <- -1 * (d + dd + (nrow(X) - ncol(X)) * (1 + lnum - log(nrow(X) - ncol(X)))) / 2
  ll_ML <-  -1 * (d + nrow(X) * (1 + lnum - log(nrow(X)))) / 2

  beta <- beta[order(names(beta))]

  list("type" = "lm",
       "model" = model,
       "beta" = beta,
       "vcs" = vcs,
       "rss" = rss,
       "ll_REML" = ll_REML,
       "ll_ML" = ll_ML,
       "method" = ifelse(REML, "REML", "ML"))
}
