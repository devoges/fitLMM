#' @export
fit_lmm <- function(Y,
                    X,
                    Zlist = NULL,
                    Slist = NULL,
                    REML) {
  if (!is.null(Zlist)) {
    fit_lmer(Y = Y,
             X = X,
             Zlist = Zlist,
             REML = REML)
  } else if (!is.null(Slist)) {
    fit_gmmat(Y = Y,
              X = X,
              Slist = Slist,
              REML = REML)
  } else {
    fit_lm(Y = Y,
           X = X,
           REML = REML)
  }
}
