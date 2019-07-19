det <- function(x) {
  r <- determinant(x)$modulus
  attr(r, "logarithm") <- NULL
  r
}

#' @export
#' @importFrom GMMAT glmmkin
fit_gmmat <- function(Y,
                      X,
                      Slist,
                      REML = TRUE) {
  n <- nrow(X)
  data <- as.data.frame(X)
  data$Y <- Y
  data$id <- 1:n
  Slist_new <- lapply(1:length(Slist), function(i) {
    out <- Slist[[i]]
    colnames(out) <- 1:n
    rownames(out) <- 1:n
    out
  })
  names(Slist_new) <- names(Slist)

  method <- ifelse(REML, "REML", "ML")

  model <- glmmkin(fixed  = Y ~ -1 + . - id,
                   data   = data,
                   family = gaussian(),
                   kins   = Slist_new,
                   tol    = .Machine$double.eps^0.5,
                   id     = "id",
                   method = method)
  vcs <- setNames(model$theta, c("Error", names(Slist)))
  #vcs <- vcs[order(names(vcs))]
  beta <- setNames(model$coefficients, colnames(X))
  #beta <- beta[order(names(beta))]

  v <- lapply(names(vcs) %>% remove_these("Error"),
              function(l) (vcs[l]/vcs["Error"]) * Slist[[l]]) %>%
    Reduce("+", .)
  diag(v) <- diag(v) + 1
  vinv <- solve(v)
  r <- Y - as.matrix(X) %*% beta
  rss <- {t(r) %*% vinv %*% r}[1, 1]
  lnum <- log(2 * pi * rss)
  d <- det(v)
  dd <- det(t(X) %*% vinv %*% X)
  ll_REML <- -1 * (d + dd + (nrow(X) - ncol(X)) * (1 + lnum - log(nrow(X) - ncol(X)))) / 2
  ll_ML <-  -1 * (d + nrow(X) * (1 + lnum - log(nrow(X)))) / 2

  list("type" = "GMMAT",
       "model" = model,
       "beta" = beta,
       "vcs" = vcs,
       "rss" = rss,
       "ll_REML" = ll_REML,
       "ll_ML" = ll_ML,
       "method" = method)
}
