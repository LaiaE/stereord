#' A exhaustive search model selection procedure for the
#' Partial Ordered Stereotype Model (POSM).
#'
#' A model selection strategy for the POSM that
#' builds the model incrementally performs an exhaustive search of the overall
#' best model. The overall best model is the model with the lowest AIC among
#' all the candidates. It fits all possible models for each number of possible
#' covariates in each step until no further improvement in AIC is observed,
#' that is, the AIC stops decreasing.
#'
#'
#'@param Outcome name of the outcome (response) variable.
#'
#'@param Covariates a character vector containing the names of the potential
#'covariates.
#'
#'@param data a data frame, list, or environment in which to evaluate the
#'variables used in the model selection procedure.
#'
#' @param split_categorical logical for whether different levels of a categorical
#' variable should be forced into the same group. If TRUE, dummy variables derived
#' from the same categorical covariate are treated as a independent covariates.
#'
#' @export
exhaustive.search.POSM <- function(Outcome, Covariates, data, split_categorical = FALSE, keep = FALSE){
  start_time <- Sys.time()
  if(!is.factor(data[,Outcome])) stop("response must be a factor")
  if (any(!(Covariates %in% names(data))))
    warning(paste0("Variables ",
      paste0(Covariates[!(Covariates %in% names(data))], collapse = ", "),
      " are not in the data"))

  Covariates <- names(data)[names(data) %in% Covariates]

  if (!(Outcome %in% names(data)))
    stop("The outcome variable is not found in the data.")
  if (length(Covariates) <= 1)
    stop("At least two covariates are required to use this method.")

  if (split_categorical) {
    matrixdata <- model.matrix(as.formula(paste0("~ ", paste(
      Covariates, collapse = " + "
    ))),
    model.frame(as.formula(paste0(
      "~ ", paste(Covariates, collapse = " + ")
    )), data, na.action = na.pass))
    matrixdata <- as.data.frame(matrixdata[, -which(colnames(matrixdata) == "(Intercept)")])

    data <- cbind(data[, Outcome], matrixdata)
    names(data) <- c(Outcome, Covariates)
    Covariates <- names(matrixdata)
  }

  formula_aux <- paste0(Outcome, " ~ ")
  S <- length(Covariates)

  # Step 1: Fit univariate OSM with each potential covariate and keep the AIC
  AIC_aux <- sapply(Covariates, function(var)
    as.numeric(try(posm(formula = as.formula(paste0(formula_aux, var)),
                        grouping = 1, data = data)$aic,
                   silent = TRUE)
    ))
  MinimumAIC <- min(AIC_aux, na.rm = TRUE)
  if (is.infinite(MinimumAIC)) stop("No model configuration led to convergence.")
  Chosen.Variable <- names(AIC_aux[which.min(AIC_aux)])
  formula_minAIC <- paste0(formula_aux, Chosen.Variable)
  grouping_aux <- 1
  if (keep) history <- paste0(formula_minAIC, " (grouping: ", paste(grouping_aux, collapse = ", "),
      "; AIC: ", round(MinimumAIC, 3), ") Time: ",
      round(difftime(time1 = Sys.time(),time2 = start_time,units = "secs"),3), " secs")

  CONT <- TRUE
  step <- 2

  while(CONT){

    # Step 3: Fit possible OSM and POSM adding 1 covariate
    Posible.cov <- combn(Covariates, m = step)

    Possible.SetCov <- expand.grid(rep(list(1:step), step))
    ind <- apply(Possible.SetCov, 1,
                 function(l){l <- unique(as.numeric(l));
                 ifelse(!identical(l, as.numeric(1:max(l))), 2, 1)})
    Possible.SetCov <- Possible.SetCov[ind == 1,]
    Possible.SetCov <- split(as.matrix(Possible.SetCov), seq(nrow(Possible.SetCov)))
    Possible.SetCov <- lapply(Possible.SetCov, as.numeric)

    AIC_aux <- lapply(Possible.SetCov, function(SC){
      apply(Posible.cov, 2, function(cov){
        as.numeric(try(posm(as.formula(paste0(formula_aux,  paste(cov, collapse = " + "))),
                                 data = data, grouping = SC )$aic, silent = TRUE))
      })
    })

    if(!any(unlist(AIC_aux) < MinimumAIC, na.rm = TRUE)){
      CONT <- FALSE
      if (keep) history <- history[-length(history)]
    } else if(step == S){
      CONT <- FALSE
      MinimumAIC <- min(unlist(AIC_aux), na.rm = TRUE)
      ind1 <- which(sapply(AIC_aux, function(x) any(x == MinimumAIC)))
      grouping_aux <- Possible.SetCov[[ind1]]
      Chosen.Variable <- Posible.cov[,which.min(AIC_aux[[ind1]])]
      formula_minAIC <- paste0(formula_aux,  paste(Chosen.Variable, collapse = " + "))
    } else {
      MinimumAIC <- min(unlist(AIC_aux), na.rm = TRUE)
      ind1 <- which(sapply(AIC_aux, function(x) any(x == MinimumAIC)))
      grouping_aux <- Possible.SetCov[[ind1]]
      Chosen.Variable <- Posible.cov[,which.min(AIC_aux[[ind1]])]
      formula_minAIC <- paste0(formula_aux,  paste(Chosen.Variable, collapse = " + "))
      if (keep) history <- c(history, paste0(formula_minAIC, " (grouping: ", paste(grouping_aux, collapse = ", "),
                                  "; AIC: ", round(MinimumAIC, 3), ") Time: ",
                                  round(difftime(time1 = Sys.time(),time2 = start_time,units = "secs"),3), " secs"))

      step <- step + 1
    }
  }
  end_time <- Sys.time()
  OUT <- paste0(formula_minAIC, " (SetCov: ", paste(grouping_aux, collapse = ", "), "; AIC: ", round(MinimumAIC,3), ") Time: ",
                round(difftime(time1 = end_time, time2 = start_time, units = "secs"),3), " secs")
  if (keep) OUT <- c(history, OUT)
  if (keep) names(OUT) <- paste("Step", seq_len(length(OUT)), ":")
  return(OUT)
}
