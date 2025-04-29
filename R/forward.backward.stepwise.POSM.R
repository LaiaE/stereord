#' A forward-backward stepwise model selection procedure specifically adapted to the
#' Partial Ordered Stereotype Model (POSM).
#'
#' A model selection strategy for the POSM that builds the model incrementally
#' using a forward stepwise approach. At each step, the method selects both
#' the covariate and its group assignment that achieve the greatest improvement
#' in the Akaike Information Criterion (AIC). Therefore, this algorithm not only
#' selects the relevant predictors, but also determines the appropriate number
#' of covariate groups and assigns predictors to groups accordingly.
#'
#'@param Outcome name of the outcome (response) variable.
#'
#'@param Covariates a character vector containing the names of the potential
#'covariates.
#'
#'@param data a data frame, list, or environment in which to evaluate the
#'variables used in the model selection procedure.
#'
#' @param IC Information criterion to be used as the stopping rule. Must be
#' either 'AIC' or 'BIC'.
#'
#' @param NominalVarSameGroup logical for whether different levels of a categorical
#' variable should be forced into the same group.
#'
#' @param keep logical for whether to keep a record of all selection steps.
#'
#' @param maxSteps the maximum number of steps allowed during the model selection process.
#'
#' @export
forward.backward.stepwise.POSM <- function(Outcome,
                                           Covariates,
                                           data,
                                           NominalVarSameGroup = FALSE,
                                           keep = FALSE,
                                           maxSteps = 100){
  start_time <- Sys.time()

  if (any(!(Covariates %in% names(data))))
    warning(paste0(
      "Variables ",
      paste0(Covariates[!(Covariates %in% names(data))], collapse = ", "),
      " are not in the data"
    ))

  Covariates <- names(data)[names(data) %in% Covariates]

  if (!(Outcome %in% names(data)))
    stop("The outcome variable is not found in the data.")
  if (length(Covariates) <= 1)
    stop("At least two covariates are required to use this method.")
  if (!(IC %in% c("AIC", "BIC")))
    stop("The information criterion (IC) must be either 'AIC' or 'BIC'.")

  if (NominalVarSameGroup) {
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


  AIC_aux <- sapply(Covariates, function(var)
    as.numeric(try(posm(formula = as.formula(paste0(formula_aux, var)),
                        grouping = 1,
                        data = data)$aic,
                   silent = TRUE)
    ))
  MinimumAIC <- min(AIC_aux, na.rm = TRUE)
  if (is.infinite(MinimumAIC)) stop("No model configuration led to convergence.")
  Chosen.Variable <- names(AIC_aux[which.min(AIC_aux)])
  formula_minAIC <- paste0(formula_aux, Chosen.Variable)
  Model.Vars <- Chosen.Variable
  Covariates <- Covariates[!(Covariates %in% Chosen.Variable)]
  grouping_aux <- 1
  if (keep)
    history <- paste0(
      formula_minAIC, " (grouping: ", paste(grouping_aux, collapse = ", "), "; AIC: ", round(MinimumAIC, 3), ") Time: ",
      round(difftime(time1 = Sys.time(), time2 = start_time, units = "secs"), 3)," secs"
    )

  CONT <- TRUE
  step <- 2
  add <- TRUE

  while(CONT){
    # Add
    Posible.grouping <- lapply(c(unique(grouping_aux), max(grouping_aux) + 1L), function(x) c(grouping_aux, x))

    AIC_aux <- lapply(Posible.grouping, function(SC) {
      sapply(Covariates, function(var)
        as.numeric(try(posm(formula = as.formula(paste0(formula_minAIC, "+", var)),
          data = data, grouping = SC)$aic, silent = TRUE))
        )
    })
    names(AIC_aux) <- rep("ADD", length(AIC_aux))

    # Remove
    if(length(Model.Vars)>2){
      Posible.cov_rem <- combn(Model.Vars, m = (length(Model.Vars)-1))
      if(add){Posible.cov_rem <- Posible.cov_rem[,apply(Posible.cov_rem, 2, function(x) any(x == Chosen.Variable))]}

      Posible.SetCov_rem <- apply(Posible.cov_rem, 2, function(x) grouping_aux[which(Model.Vars %in% x)])

      # Asses if there is some SetCov not addecuated (Not start by 1 and have all the numbers)
      ind <- apply(Posible.SetCov_rem, 2, function(l){l <- unique(as.numeric(l)); !identical(l, as.numeric(1:max(l)))})
      if(any(ind)){
        if(sum(ind) == 1){
          Posible.SetCov_rem[,ind] <- rep(1:length(table(Posible.SetCov_rem[,ind])), table(Posible.SetCov_rem[,ind]))
        } else{
          Posible.SetCov_rem[,ind] <- apply(Posible.SetCov_rem[,ind], 2, function(x) rep(1:length(table(x)), table(x)))
        }
      }


      AIC_aux_rem <- sapply(1:ncol(Posible.cov_rem), function(i){
        as.numeric(try(posm(as.formula(paste0(formula_aux,  paste(Posible.cov_rem[,i], collapse = " + "))),
                                 data = data, grouping = Posible.SetCov_rem[,i])$aic, silent = TRUE))
      })


      AIC_aux <- c(AIC_aux, REMOVE = list(AIC_aux_rem))
    }


    if(!any(unlist(AIC_aux) < MinimumAIC, na.rm = TRUE)){
      CONT <- FALSE
    } else if(length(Covariates) == 1 | step == maxSteps){

      CONT <- FALSE
      MinimumAIC <- min(unlist(AIC_aux), na.rm = TRUE)
      ind1 <- which(sapply(AIC_aux, function(x) any(x == MinimumAIC)))

      if(names(ind1) == "ADD"){
        grouping_aux <- Posible.grouping[[ind1]]
        Chosen.Variable <- names(which.min(AIC_aux[[ind1]]))
        Model.Vars <- c(Model.Vars, Chosen.Variable)

      } else {
        ind2 <- which(AIC_aux[[ind1]] == MinimumAIC)
        grouping_aux <- Posible.SetCov_rem[,ind2]
        Model.Vars <- Posible.cov_rem[,ind2]

      }

      formula_minAIC <- paste0(formula_aux, paste0(Model.Vars, collapse = " + "))
      if (keep) history <- c(history, paste0(formula_aux," (grouping: ", paste(grouping_aux, collapse = ", "), "; AIC: ",
                                             round(MinimumAIC, 3), ") Time: ",
                                             round(difftime(time1 = Sys.time(), time2 = start_time, units = "secs"), 3), " secs"))

      # Last step remove
      if(length(Covariates) == 1 & length(Model.Vars)>2){
        Posible.cov_rem <- combn(Model.Vars, m = (length(Model.Vars)-1))
        if(add){Posible.cov_rem <- Posible.cov_rem[,apply(Posible.cov_rem, 2, function(x) any(x == Chosen.Variable))]}

        Posible.SetCov_rem <- apply(Posible.cov_rem, 2, function(x) grouping_aux[which(Model.Vars %in% x)])

        # Asses if there is some SetCov not addecuated (Not start by 1 and have all the numbers)
        ind <- apply(Posible.SetCov_rem, 2, function(l){l <- unique(as.numeric(l)); !identical(l, as.numeric(1:max(l)))})
        if(any(ind)){
          if(sum(ind) == 1){
            Posible.SetCov_rem[,ind] <- rep(1:length(table(Posible.SetCov_rem[,ind])), table(Posible.SetCov_rem[,ind]))
          } else{
            Posible.SetCov_rem[,ind] <- apply(Posible.SetCov_rem[,ind], 2, function(x) rep(1:length(table(x)), table(x)))
          }
        }


        AIC_aux_rem <- sapply(1:ncol(Posible.cov_rem), function(i){
          as.numeric(try(posm(as.formula(paste0(formula_aux,  paste(Posible.cov_rem[,i], collapse = " + "))),
                                   data = data, grouping = Posible.SetCov_rem[,i])$aic, silent = TRUE))
        })

        if(any(AIC_aux_rem < MinimumAIC, na.rm = TRUE)){
          if (keep) history <- c(history, paste0(formula_aux," (grouping: ", paste(grouping_aux, collapse = ", "), "; AIC: ",
                                                 round(MinimumAIC, 3), ") Time: ",
                                                 round(difftime(time1 = Sys.time(), time2 = start_time, units = "secs"), 3), " secs"))
          MinimumAIC <- min(AIC_aux_rem, na.rm = TRUE)
          ind2 <- which(AIC_aux_rem == MinimumAIC)
          grouping_aux <- Posible.SetCov_rem[,ind2]
          Model.Vars <- Posible.cov_rem[,ind2]
          formula_minAIC <- paste0(formula_aux, paste0(Model.Vars, collapse = " + "))
        }

      }


    } else {

      MinimumAIC <- min(unlist(AIC_aux), na.rm = TRUE)
      ind1 <- which(sapply(AIC_aux, function(x) any(x == MinimumAIC, na.rm = TRUE)))

      if(names(ind1) == "ADD"){
        grouping_aux <- Posible.grouping[[ind1]]

        Chosen.Variable <- names(which.min(AIC_aux[[ind1]]))
        Model.Vars <- c(Model.Vars, Chosen.Variable)
        Covariates <- Covariates[!(Covariates %in% Chosen.Variable)]

        add <- TRUE

      } else {
        ind2 <- which(AIC_aux[[ind1]] == MinimumAIC)
        grouping_aux <- Posible.SetCov_rem[,ind2]
        Model.Vars <- Posible.cov_rem[,ind2]
        add <- FALSE

      }

      formula_minAIC <- paste0(formula_aux, paste0(Model.Vars, collapse = " + "))
      step <- step + 1

    }
  }

  end_time <- Sys.time()
  OUT <- paste0(formula_minAIC, " (SetCov: ", paste(grouping_aux, collapse = ", "),
                "; AIC: ", round(as.numeric(MinimumAIC),3), ") Time: ", difftime(time1 = end_time, time2 = start_time, units = "secs"), " secs")
  if (keep) OUT <- c(history, OUT)
  if (keep) names(OUT) <- paste("Step", seq_len(length(OUT)), ":")
  return(OUT)
}
