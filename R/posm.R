#' Ordinal data regression using the Partial Ordered Stereotype Model (POSM).
#'
#' Fit a regression model to an ordered factor response. The model is NOT a
#' logistic or probit model because the link function is not the logit, but the
#' link function is log-based.
#'
#'This function should be used in a very similar way to \code{clustord::osm}, and
#' some of the arguments are the same as \code{osm}, but the ordinal model used
#' here is less restrictive as it allaw multiple score parameters sets within the same model.
#'
#' @param formula an object of class formula with description of the model to be fitted.
#' @param data a data frame o tibble containing the variables in the model.
#' @param phi1.iszero logicals.  If \code{TRUE} means we fixed \eqn{\phi_1=0} and \eqn{\phi_q = 1}, else \eqn{\phi_1=1} and \eqn{\phi_q = 0}
#' @param repar logicals.  If \code{TRUE} means the response variable is ordinal and the reparametrization to ensure the score constraint is applied.
#' @param Approx.start.values logicals.  If \code{TRUE} means that the initial values for the parameters to be optimized over are
#' approximated, else are assumed to be 0.1.
#' @param method The method to be used to optimize the log likelihood. By default is \code{BFGS} (quasi-Newton method).
#' @param control a list to control some procedures in \code{optim()}.
#' @param HES  logicals.  If \code{TRUE} means that the output will include the Hessian matrix.
#' @param SE  logicals.  If \code{TRUE} means that the output will include the standard errors of the coefficients.
#' @param ... For \code{OSM()}: additional arguments to be passed to the low level in \code{optim()}.
#'
#' @return \code{OSM} returns a list containing at least the following components:
#'
#'  \code{alpha}:  Matrix with the estimated parameters vector $\\{ \\alpha_k \\}$ and their standard error if \code{SE = TRUE}.
#'
#'  \code{phi}: Matrix with the estimated parameters vector $\\{\\phi_k\\}$ and their standard error if \code{SE = TRUE}.
#'
#'  \code{beta}: Matrix with the estimated parameters vector $\\{\\beta_S \\}$ and their standard error if \code{SE = TRUE}.
#'
#'  \code{logLike}: Estimated optimal value of the menus log-likelihood.
#'
#'  \code{AIC}, \code{BIC}: the calculed values of Akaike information criterion (AIC) and Bayesian information criterion (BIC)
#'
#'  There are some optional outputs:
#'
#'  \code{phi.beforeRep}: If \code{repar = TRUE}. Matrix with the estimated parameters vector $\\{\\phi_k\\}$ before to invert reparametrization.
#'
#'  \code{hessian}: If \code{HES = TRUE}. Hessian Matrix.
#'
#' @references
#' Fernandez, D., Arnold, R., & Pledger, S. (2016). Mixture-based clustering for the ordered stereotype model. *Computational Statistics & Data Analysis*, 93, 46-75.
#' Anderson, J. A. (1984). Regression and ordered categorical variables. *Journal of the Royal Statistical Society: Series B (Methodological)*, 46(1), 1-22.
#' Agresti, A. (2010). *Analysis of ordinal categorical data* (Vol. 656). John Wiley & Sons.
#'
#' @export
posm <- function(formula, grouping, data, weights, start, ..., subset,
                na.action, Hess = FALSE, model = TRUE, grouping.offsets = NULL)
{

  ## Create an object containing the original function call, which will
  ## later be used to obtain the model frame
  ## "expand.dots=FALSE" means those parts are left as "..." instead of
  ## converted to named arguments
  m <- match.call(expand.dots = FALSE)

  ## Convert the data from a matrix to a data frame if required
  if(is.matrix(eval.parent(m$data))) m$data <- as.data.frame(data)

  ## Delete any input arguments not needed to create the model frame,
  ## including the ... arguments
  m$start <- m$Hess <- m$model <- m$grouping <- m$... <- NULL


  ## Keep the relevant inner parts of the function call, but change it to
  ## a call to "model.frame()"
  m[[1L]] <- quote(stats::model.frame)

  ## Evaluate this call to model.frame() in the parent frame i.e. where
  ## osm() was called from. This gets the model frame that contains only
  ## the data rows in the subset, and only the variables required by the
  ## formula
  m <- eval.parent(m)

  ## Also attach the Terms object of the formula to the model frame object
  Terms <- attr(m, "terms")

  ## Now create the model matrix, i.e. just the predictors, not the response,
  ## and with each categorical predictor changed to multiple dummy variables
  ## in the manner required by the "contrasts" option
  x <- model.matrix(Terms, m)

  ## Get the column index of the (Intercept) column that was created by
  ## the call to model.matrix()
  xint <- match("(Intercept)", colnames(x), nomatch = 0L)

  ## Count rows and columns in the model matrix
  n <- nrow(x)
  num_beta <- ncol(x)

  ## Drop the (Intercept) column from the model matrix
  if(xint > 0L) {
    x <- x[, -xint, drop = FALSE]
    num_beta <- num_beta - 1L
  } else warning("an intercept is needed and assumed")

  ## Fetch the weights, or generate ones
  wt <- model.weights(m)
  if(!length(wt)) wt <- rep(1, n)

  ## Fetch the offsets, or generate zeros
  offset <- model.offset(m)
  if(length(offset) <= 1L) offset <- rep(0, n)
  # if(length(offset) <= 1L) offset <- matrix(1, ncol = H, nrow = n)

  ## Fetch the response variable, and its levels, and count them
  y <- model.response(m)
  if(!is.factor(y)) stop("response must be a factor")
  lev <- levels(y); llev <- length(lev)
  if(llev <= 2L) stop("response must have 3 or more levels")
  y <- unclass(y)
  qminus <- llev - 1L

  ## Fetch the grouping of covariates
  if(length(grouping) != length(attr(Terms, "term.labels"))) stop("grouping must have the same length as the covariates")

  # We create a new grouping according to the new data matrix with the dummy variables
  names(grouping) <- attr(Terms, "term.labels")
  newgrouping <- numeric(num_beta)
    for(i in attr(Terms, "term.labels")){
    newgrouping[grep(i, colnames(x))] <- grouping[i]
  }


  # newgrouping <- as.numeric(factor(newgrouping))
  if(!is.numeric(newgrouping)) stop("grouping must be numeric")
  if(newgrouping[1L] != 1) stop("grouping must have 1 as its first element")
  groups <- sort.int(unique(newgrouping))
  H <- length(groups)
  if(!all(seq_len(H) == groups)) stop("grouping must be consecutive numbers")
  ind_H <- seq_len(H)

  ## Generate starting values for optimization
  if(missing(start)) {
    # try logistic/probit regression on 'middle' cut to find starting
    # values for the coefficients of the predictors
    # q1 is the level at, or just before, halfway through the levels of y
    q1 <- llev %/% 2L

    ## y1 is a binary variable with y1=0 if y <= q1 and y1=1 if y > q1
    y1 <- (y > q1)

    ## Construct a new model matrix and add an intercept column to it
    X <- cbind(Intercept = rep(1, n), x)

    ## Now attempt to fit logistic regression to the binary response y1
    fit <- glm.fit(X, y1, wt, family = binomial(), offset = offset)
    if(!fit$converged)
      stop("attempt to find suitable starting values failed")
    coefs <- fit$coefficients
    if(any(is.na(coefs))) {
      warning("design appears to be rank-deficient, so dropping some coefs")
      keep <- names(coefs)[!is.na(coefs)]
      coefs <- coefs[keep]
      x <- x[, keep[-1L], drop = FALSE]
      num_beta <- ncol(x)
    }

    ## The other parameters are labelled as alphas in Agresti's definition
    ## of the proportional odds model. They are the base probabilities
    ## for each level of the response variable, and must be strictly increasing.
    ## Generate them initially assuming they're evenly spaced across the
    ## range, convert to the linear predictor space using the logit link
    ## and adjust them to incorporate the fact that the logistic
    ## regression fitting produced an intercept term coefs[1L] for the
    ## q1 level
    logit <- function(p) log(p/(1 - p))
    spacing <- logit((1L:(qminus))/(qminus+1)) # just a guess
    gammas <- -coefs[1L] + spacing - spacing[q1]

    ## Also generate starting values for phi, assuming equal spacing in
    ## the space of phi and converting using the logit link to the space
    ## of the auxiliary variable u
    startingphi <- (1:(qminus-1))/(qminus)
    u2 <- logit(startingphi)[1]
    us <- log(diff(logit(startingphi)))

    ## Construct the full starting values vector, using the fact that
    ## coefs[1L] has already been incorporated into the gammas object
    start <- c(coefs[-1L], gammas, rep(c(u2, us), H))
  } else if(length(start) != num_beta + (llev-1) + H*(llev-2))
    stop("'start' is not of the correct length")

  ## Now run the fitting
  ans <- posm.fit(x, y, newgrouping, wt, start, offset, hessian = Hess, ...)


  ## Extract parts of the fitted model, to use when calculating the fitted
  ## values for each observation
  ## "res" is the output object from optim(), which contains the hessian
  ## object if requested
  beta <- ans$beta
  mu <- c(0,ans$mu)
  phi <- ans$phi
  res <- ans$res
  deviance <- ans$deviance
  u <- ans$u

  ## Calculate the fitted values of each observation, which are the probabilities
  ## of getting each of the levels of the response
  eta <- if(num_beta) offset + sapply(ind_H, function(h) drop(x[,which(newgrouping == h)] %*% beta[which(newgrouping == h)]))
 else offset + rep(0, n)
  fitted <- matrix(1,n,llev)
  for (k in 2:(llev)) {
    fitted[,k] <- exp(pmax(pmin(50,mu[k]+drop(eta %*% phi[k,])),-100))
  }
  fitted <- fitted/rowSums(fitted)
  dimnames(fitted) <- list(row.names(m), lev)

  ## Count the number of calls to the function, some of which were used
  ## to numerically calculate the gradient
  niter <- c(f.evals = res$counts[1L], g.evals = res$counts[2L])

  groups <- sapply(ind_H, function(h) paste("Group", h, ":", paste(colnames(x)[which(newgrouping == h)], collapse = ", ")))

  ## Construct the output object
  fit <- list(beta = beta, mu = mu, phi = phi, u = u, deviance = deviance,
              groups = groups, newgrouping = newgrouping,
              fitted.values = fitted, lev = lev, terms = Terms,
              df.residual = sum(wt) - num_beta - qminus - (qminus-1),
              edf = num_beta + qminus + (qminus-1), n = sum(wt),
              nobs = sum(wt), call = match.call(),
              convergence = res$convergence, niter = niter, eta = eta)

  if(Hess) {
    dn <- c(names(beta), names(ans$mu), paste0(rep(rownames(ans$u), times = ncol(ans$u)), rep(colnames(ans$u), each = nrow(ans$u))))
    H <- res$hessian
    dimnames(H) <- list(dn, dn)
    fit$Hessian <- H
  }
  if(model) fit$model <- m
  fit$na.action <- attr(m, "na.action")
  fit$xlevels <- .getXlevels(Terms, m)
  class(fit) <- "posm"
  fit
}

#' @importFrom stats optim
posm.fit <- function(x, y, newgrouping, wt, start, offset, ...)
{
  ## Set up the function call to use in optim(), which extracts the parameters
  ## from the parameter vector and calculates the negative of the log-likelihood
  fmin <- function(coefficients) {
    mu <- c(0, coefficients[num_beta + ind_mu_k])

    ## u are the auxiliary parameters for phi, which are used because they
    ## are free to take any value between -Inf and Inf, and don't have to be
    ## ordered, whereas the phi values, from their construction, will be
    ## increasing, and the end values will be 0 and 1
    u <- matrix(coefficients[num_beta + num_mu_k + ind_phi_totalk], ncol = H, nrow = num_phi_k)
    phi <- apply(u, 2, function(v) c(0 , expit(cumsum(c(v[1L], exp(v[-1L])))), 1))
    eta <- offset
    if (num_beta) eta <- eta + sapply(ind_H, function(h) drop(x[,which(newgrouping == h)] %*% coefficients[which(newgrouping == h)]))


    ## Construct the probabilities of getting each level of the response for
    ## each observation
    theta <- matrix(1, nrow=n, ncol=num_mu_k+1)
    for (k in 2:(num_mu_k+1)) {
      theta[,k] <- exp(pmax(pmin(50,mu[k]+drop(eta %*% phi[k,])),-100))
    }
    theta <- theta/rowSums(theta)

    ## Now calculate the components of the likelihood for each observation
    pr <- vapply(1:n, function(i) theta[i,y[i]], 1)

    ## Construct the negative log-likelihood
    if (all(pr > 0)) -sum(wt * log(pr)) else Inf
  }

  ## Count the number of rows and columns in the model matrix, and the number
  ## of predictors (including dummy variables, not the original categorical variables)
  n <- nrow(x)
  num_beta <- ncol(x)
  ind_beta <- seq_len(num_beta)

  ## Count the number of groups
  H <- length(unique(newgrouping))
  ind_H <- seq_len(H)

  ## Count the number of levels of y, and calculate q and q2 (there will be
  ## q independent values of the alpha parameters, and q2 = q-1 independent
  ## values of phi)
  lev <- levels(y)
  if(length(lev) <= 2L) stop("response must have 3 or more levels")
  y <- unclass(y)
  num_mu_k <- length(lev) - 1L
  ind_mu_k <- seq_len(num_mu_k)

  num_phi_totalk <- (length(lev) - 2L)*H
  ind_phi_totalk <- seq_len(num_phi_totalk)
  num_phi_k <- length(lev) - 2L
  ind_phi_k <- seq_len(num_phi_k)

  ## Run optim, and extract the results
  res <- optim(start, fmin, method="L-BFGS-B", ...)
  beta <- res$par[ind_beta]
  mu <- res$par[num_beta + ind_mu_k]
  u <- matrix(res$par[num_beta + num_mu_k + ind_phi_totalk], ncol = H, nrow = num_phi_k)
  phi <- phi <- apply(u, 2, function(v) c(0 , expit(cumsum(c(v[1L], exp(v[-1L])))), 1))
  deviance <- 2 * res$value
  names(mu) <- paste(lev[1L], lev[-1L], sep="|")
  rownames(phi) <- lev
  colnames(phi) <- paste("Group", ind_H)
  rownames(u) <- paste0("phiAux",lev[-c(1L,length(lev))])
  colnames(u) <- paste("Group", ind_H)
  if(num_beta) names(beta) <- colnames(x)
  list(beta = beta, mu = mu, phi = phi, u=u, deviance = deviance, res = res)
}

#' @importFrom stats naprint
#' @export
print.posm <- function(x, ...)
{
  if(!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl, control=NULL)
  }
  if(length(x$beta)) {
    cat("\nCoefficients beta:\n")
    print(x$beta, ...)
  } else {
    cat("\nNo coefficients\n")
  }
  cat("\nIntercepts mu:\n")
  print(x$mu, ...)
  cat("\nScore parameters phi:\n")
  print(x$phi, ...)
  cat("\nGrouping:\n")
  print(x$groups)
  cat("\nResidual Deviance:", format(x$deviance, nsmall=2L), "\n")
  cat("AIC:", format(x$deviance + 2*x$edf, nsmall=2L), "\n")
  cat("BIC:", format(x$deviance + x$edf*log(x$n), nsmall=2L), "\n")
  if(nzchar(mess <- naprint(x$na.action))) cat("(", mess, ")\n", sep="")
  if(x$convergence > 0)
    cat("Warning: did not converge as iteration limit reached\n")
  invisible(x)
}

#' @export
vcov.posm <- function(object, ...){

  pc <- length(object$beta)
  llev <- length(object$lev)
  num_mu_k <- llev - 1L
  ind_mu_k <- seq_len(num_mu_k)
  H <- ncol(object$u)
  ind_H <- seq_len(H)
  num_u_k <- llev-2
  ind_u_k <- seq_len(num_u_k)

  num_u_total <- num_u_k*H

  if(is.null(object$Hessian)) {
    message("\nRe-fitting to get Hessian\n")
    utils::flush.console()
    object <- update(object, Hess = TRUE,
                     start = c(object$beta, object$mu[ind_mu_k + 1L], c(object$phi[ind_u_k + 1L, ])))
  }

  vc <- ginv(object$Hessian)

  # delta method
  u <- object$u
  u.ind <- pc + num_mu_k + seq_len(num_u_total)

  # Define reparametrization in formula format
  create_formula <- function(indices) {
    paste0("~ 1 / (1 + exp(-((", paste(indices, collapse = ") + exp("), "))))")
  }

  syms <- paste0("x", ind_u_k)
  formulas <- lapply(seq_along(syms), function(i) { create_formula(syms[1:i])})
  trans.formula <- lapply(formulas, as.formula)

  J <- lapply(ind_H, function(h) {
    u.aux <- u[, h]
    for (i in ind_u_k) {
      assign(syms[i], u.aux[i])
    }
    sapply(trans.formula, function(form) {
      as.numeric(attr(eval(deriv(form, syms)), "gradient"))
    })
  })


  A <- diag(pc + num_mu_k + num_u_total)
  for (h in ind_H) {
    idx <- u.ind[(h - 1) * num_u_k + ind_u_k]
    A[idx, idx] <- J[[h]]
  }

  V <- t(A) %*% vc %*% A

  structure(V, dimnames = lapply(dimnames(object$Hessian), function(x) gsub("phiAux", "phi", x)))
}

#' @export
summary.posm <- function(object, digits = max(3, .Options$digits - 3), correlation = FALSE,
                         signif.stars = getOption("show.signif.stars"),...){
  pc <- length(object$beta)
  q <- nrow(object$phi)
  H <- ncol(object$phi)
  hind <- seq_len(H)
  cc <- c(object$beta, object$mu[-1L])
  coef <- matrix(0, pc+q-1L, 4L, dimnames=list(names(cc),
                                               c("Value", "Std. Error", "t value", "Pr(>|t|)")))


  # coef <- matrix(0, pc+(q-1)+(q-2)*H, 3L, dimnames=list(names(cc),
  #                                               c("Value", "Std. Error", "t value")))
  coef[, 1L] <- cc
  vc <- vcov(object)
  sd <- sqrt(diag(vc))
  coef[, 2L] <- sd[seq_len(pc+q-1)]
  coef[, 3L] <- coef[, 1L]/coef[, 2L]
  coef[, 4L] <- 2 * pnorm(-abs(coef[, 3L]))

  phi <- object$phi
  rownames(phi) <- paste0("phi", rownames(phi))
  namesphi <- matrix(paste("phi", rep(rownames(object$phi), times = ncol(object$phi)), "|",
                           rep(colnames(object$phi), each = nrow(object$phi)), sep = ""),
                     ncol = H, nrow = q)

  coef.phi <- matrix(0, (q-2L)*H, 2L, dimnames=list(c(namesphi[c(-1L, -q),]),c("Value", "Std. Error")))
  coef.phi[, 1L] <- c(phi[c(-1L, -q),])
  coef.phi[, 2L] <- sd[pc+q-1L + seq_len((q-2L)*H)]


  vc.phi <- lapply(hind, function(h) cbind(0, rbind(0,vc[pc+q-1L + (h-1) + seq_len(q-2L),pc+q-1L + (h-1) + seq_len(q-2L)],0),0))
  v.phi <- lapply(hind, function(h) diag(vc.phi[[h]]))
  k <- 2:(q-1)
  test.phi <- replicate(H, matrix(0, q-1L, 2L, dimnames=list(c(paste(rownames(phi)[1L], "vs", rownames(phi)[2L]),
                                                  paste(rownames(phi)[k], "vs", rownames(phi)[k+1])),
                                                c("t value", "Pr(>|t|)"))), simplify = FALSE)
  for(h in hind){
    rownames(test.phi[[h]]) <- paste0("Group ", h, ": ", rownames(test.phi[[h]]))
    test.phi[[h]][1L, 1L] <- (phi[2L,h]-phi[1L,h])/sqrt(v.phi[[h]][2L] + v.phi[[h]][1L] - 2*vc.phi[[h]][2L,1L])
    test.phi[[h]][k, 1L] <- (phi[k,h]-phi[k+1L,h])/sqrt(v.phi[[h]][k] + v.phi[[h]][k+1L] - 2*vc.phi[[h]][cbind(k, k + 1L)])
    test.phi[[h]][, 2L] <- 2 * pnorm(-abs(test.phi[[h]][, 1L]))
  }

  object$coefficients <- coef
  object$coefficients.phi <- coef.phi
  object$test.phi <- do.call("rbind", test.phi)
  object$pc <- pc
  object$q <- q
  object$digits <- digits
  object$signif.stars <- signif.stars
  if(correlation)
    object$correlation <- (vc/sd)/rep(sd, rep(pc+q-1+(q-2)*H, pc+q-1+(q-2)*H))
  class(object) <- "summary.posm"
  object
}

#' @export
print.summary.posm <- function(x, digits = x$digits, signif.stars = x$signif.stars, ...){
  if(!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl, control=NULL)
  }
  coef <- x$coefficients
  coef.phi <- x$coefficients.phi
  test.phi <- x$test.phi
  pc <- x$pc
  q <- x$q
  if(pc > 0) {
    cat("\nCoefficients (beta):\n")
    printCoefmat(coef[seq_len(pc), , drop=FALSE], digits = digits, quote = FALSE,
                 signif.stars = signif.stars, signif.legend = FALSE, na.print = "NA", ...)
  } else {
    cat("\nNo coefficients\n")
  }
  cat("\nIntercepts (mu):\n")
  printCoefmat(coef[(pc+1L):nrow(coef), , drop=FALSE], digits = digits, quote = FALSE,
               signif.stars = signif.stars, signif.legend = FALSE, na.print = "NA", ...)

  cat("\nScore parameters (phi):\n")
  print(coef.phi[,,drop=FALSE], quote = FALSE,
        digits = digits, ...)

  cat("\nTest Score parameters (phi):\n")
  printCoefmat(test.phi, digits = digits, quote = FALSE, P.values = TRUE, has.Pvalue = TRUE,
               signif.stars = signif.stars, signif.legend = FALSE, na.print = "NA", ...)


  cat("\nGrouping:\n")
  print(x$groups)

  cat("\nResidual Deviance:", format(x$deviance, nsmall=2L), "\n")
  cat("AIC:", format(x$deviance + 2*x$edf, nsmall=2L), "\n")
  cat("BIC:", format(x$deviance + x$edf*log(x$n), nsmall=2L), "\n")
  if(nzchar(mess <- naprint(x$na.action))) cat("(", mess, ")\n", sep="")
  if(!is.null(correl <- x$correlation)) {
    cat("\nCorrelation of Coefficients:\n")
    ll <- lower.tri(correl)
    correl[ll] <- format(round(correl[ll], digits))
    correl[!ll] <- ""
    print(correl[-1L, -ncol(correl)], quote = FALSE, ...)
  }
  invisible(x)
}
