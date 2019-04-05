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
  vcs <- vcs[order(names(vcs))]
  beta <- setNames(model$coefficients, colnames(X))

  list("type" = "GMMAT",
       "model" = model,
       "beta" = beta,
       "vcs" = vcs,
       "method" = method)
}
