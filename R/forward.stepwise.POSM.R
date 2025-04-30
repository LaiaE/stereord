#' A forward stepwise model selection procedure specifically adapted to the
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
#' @export
forward.stepwise.POSM <- function(Outcome, Covariates, data, IC = "AIC",
                                  NominalVarSameGroup = FALSE, keep = FALSE) {
  start_time <- Sys.time()
  if(!is.factor(data[,Outcome])) stop("response must be a factor")

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
  if(IC == "AIC"){
  # Step 1: Fit univariate OSM with each potential covariate and keep the AIC
  AIC_aux <- sapply(Covariates, function(var)
    as.numeric(try(posm(formula = as.formula(paste0(formula_aux, var)),
                        grouping = 1,
                        data = data)$aic,
                   silent = TRUE)
    ))
  MinimumAIC <- min(AIC_aux, na.rm = TRUE)
  if (is.infinite(MinimumAIC)) stop("No model configuration led to convergence.")
  Chosen.Variable <- names(AIC_aux[which.min(AIC_aux)])
  formula_aux <- paste0(formula_aux, Chosen.Variable)
  Covariates <- Covariates[!(Covariates %in% Chosen.Variable)]
  grouping_aux <- 1
  if (keep)
    history <- paste0(
      formula_aux,
      " (grouping: ",
      paste(grouping_aux, collapse = ", "),
      "; AIC: ",
      round(MinimumAIC, 3),
      ") Time: ",
      round(
        difftime(
          time1 = Sys.time(),
          time2 = start_time,
          units = "secs"
        ),
        3
      ),
      " secs"
    )

  CONT <- TRUE
  while (CONT) {
    # Step 2: Fit possible OSM and POSM adding 1 covariate
    Posible.grouping <- lapply(c(unique(grouping_aux), max(grouping_aux) +
                                   1), function(x)
                                     c(grouping_aux, x))

    AIC_aux <- lapply(Posible.grouping, function(SC) {
      sapply(Covariates, function(var)
        as.numeric(try(posm(
          formula = as.formula(paste0(formula_aux, "+", var)),
          data = data,
          grouping = SC
        )$aic,
        silent = TRUE)
        ))
    })

    if (!any(unlist(AIC_aux) < MinimumAIC, na.rm = TRUE)) {
      CONT <- FALSE
      if (keep) history <- history[-length(history)]
    } else if (length(Covariates) == 1) {
      CONT <- FALSE
      MinimumAIC <- min(unlist(AIC_aux), na.rm = TRUE)
      ind1 <- which(sapply(AIC_aux, function(x)
        any(x == MinimumAIC, na.rm = TRUE)))
      grouping_aux <- Posible.grouping[[ind1]]
      Chosen.Variable <- names(which.min(AIC_aux[[ind1]]))
      formula_aux <- paste0(formula_aux, " + ", Chosen.Variable)
    } else {
      MinimumAIC <- min(unlist(AIC_aux), na.rm = TRUE)
      ind1 <- which(sapply(AIC_aux, function(x)
        any(x == MinimumAIC, na.rm = TRUE)))
      grouping_aux <- Posible.grouping[[ind1]]
      Chosen.Variable <- names(which.min(AIC_aux[[ind1]]))
      formula_aux <- paste0(formula_aux, " + ", Chosen.Variable)
      Covariates <- Covariates[!(Covariates %in% Chosen.Variable)]
      if (keep) history <- c(history, paste0(formula_aux," (grouping: ", paste(grouping_aux, collapse = ", "), "; AIC: ",
            round(MinimumAIC, 3), ") Time: ", round(difftime(time1 = Sys.time(), time2 = start_time, units = "secs"), 3), " secs")
        )
    }
  }

  OUT <- paste0(formula_aux, " (grouping: ", paste(grouping_aux, collapse = ", "), "; AIC: ",
    round(MinimumAIC, 3),") Time: ", round(difftime(time1 = Sys.time(),time2 = start_time,units = "secs"), 3)," secs")
  } else {

    # Step 1: Fit univariate OSM with each potential covariate and keep the BIC
    BIC_aux <- sapply(Covariates, function(var)
      as.numeric(try(posm(formula = as.formula(paste0(formula_aux, var)),
                          grouping = 1,
                          data = data)$bic,
                     silent = TRUE)
      ))
    MinimumBIC <- min(BIC_aux, na.rm = TRUE)
    if (is.infinite(MinimumBIC))
      stop("No model configuration led to convergence.")
    Chosen.Variable <- names(BIC_aux[which.min(BIC_aux)])
    formula_aux <- paste0(formula_aux, Chosen.Variable)
    Covariates <- Covariates[!(Covariates %in% Chosen.Variable)]
    grouping_aux <- 1
    if (keep)
      history <- paste0(
        formula_aux,
        " (grouping: ",
        paste(grouping_aux, collapse = ", "),
        "; BIC: ",
        round(MinimumBIC, 3),
        ") Time: ",
        round(
          difftime(
            time1 = Sys.time(),
            time2 = start_time,
            units = "secs"
          ),
          3
        ),
        " secs"
      )

    CONT <- TRUE
    while (CONT) {
      # Step 2: Fit possible OSM and POSM adding 1 covariate
      Posible.grouping <- lapply(c(unique(grouping_aux), max(grouping_aux) +
                                     1), function(x)
                                       c(grouping_aux, x))

      BIC_aux <- lapply(Posible.grouping, function(SC) {
        sapply(Covariates, function(var)
          as.numeric(try(posm(
            formula = as.formula(paste0(formula_aux, "+", var)),
            data = data,
            grouping = SC
          )$bic,
          silent = TRUE)
          ))
      })

      if (!any(unlist(BIC_aux) < MinimumBIC, na.rm = TRUE)) {
        CONT <- FALSE
        if (keep) history <- history[-length(history)]
      } else if (length(Covariates) == 1) {
        CONT <- FALSE
        MinimumBIC <- min(unlist(BIC_aux), na.rm = TRUE)
        ind1 <- which(sapply(BIC_aux, function(x)
          any(x == MinimumBIC, na.rm = TRUE)))
        grouping_aux <- Posible.grouping[[ind1]]
        Chosen.Variable <- names(which.min(BIC_aux[[ind1]]))
        formula_aux <- paste0(formula_aux, " + ", Chosen.Variable)
      } else {
        MinimumBIC <- min(unlist(BIC_aux), na.rm = TRUE)
        ind1 <- which(sapply(BIC_aux, function(x)
          any(x == MinimumBIC, na.rm = TRUE)))
        grouping_aux <- Posible.grouping[[ind1]]
        Chosen.Variable <- names(which.min(BIC_aux[[ind1]]))
        formula_aux <- paste0(formula_aux, " + ", Chosen.Variable)
        Covariates <- Covariates[!(Covariates %in% Chosen.Variable)]
        if (keep)
          history <- c(
            history,
            paste0(
              formula_aux,
              " (grouping: ",
              paste(grouping_aux, collapse = ", "),
              "; BIC: ",
              round(MinimumBIC, 3),
              ") Time: ",
              round(
                difftime(
                  time1 = Sys.time(),
                  time2 = start_time,
                  units = "secs"
                ),
                3
              ),
              " secs"
            )
          )
      }
    }

    OUT <- paste0(formula_aux, " (SetCov: ", paste(grouping_aux, collapse = ", "), "; BIC: ", round(MinimumBIC, 3), ") Time: ",
      round(difftime(time1 = Sys.time(), time2 = start_time, units = "secs"), 3)," secs")
  }
  if (keep)
    OUT <- c(history, OUT)
  if (keep)
    names(OUT) <- paste("Step", seq_len(length(OUT)), ":")
  return(OUT)
}
