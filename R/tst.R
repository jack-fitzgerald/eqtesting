### INPUTS ###
#estimate: The estimate of interest; numeric scalar
#se: The estimate of interest's standard error; numeric scalar greater than zero
#ROPE: Can either be a strictly positive numeric scalar (interpreted as the width of a symmetric ROPE around zero), or a vector of two different numeric scalars
#df: If added, must be a positive integer. If left blank, asymptotic normal approximations are reported. If provided, exact results are reported
#alpha: Defaults to 0.05. If provided, must be a numeric scalar strictly between 0 and 0.5
#power: Defaults to 0.8. If provided, must be a numeric scalar strictly between 0.5 and 1
### OUTPUTS ###
#bounds: data.frame consisting of CI boundaries (asymptotic or exact, depending on whether degees of freedom are offered)
#test: Only generated if ROPE is provided; data.frame consisting of the t-statistic and TOST p-value for an equivalence test within the provided ROPE

tst = function(estimate, se, ROPE, df = NA, alpha = 0.05, power = 0.8) {

  ##################
  ##### ERRORS #####
  ##################

  #If any argument of the function is a list, data.frame, or matrix...
  if (is.list(estimate) | is.list(se) | is.list(ROPE) | is.list(df) | is.list(alpha) | is.list(power) | is.data.frame(estimate) | is.data.frame(se) | is.data.frame(ROPE) | is.data.frame(df) | is.data.frame(alpha) | is.data.frame(power) | is.matrix(estimate) | is.matrix(se) | is.matrix(ROPE) | is.matrix(df) | is.matrix(alpha) | is.matrix(power)) {

    #... then stop the function
    stop("'tst' does not accept lists, dataframes, or matrices")

  }

  #If estimate is not a numeric scalar...
  if (!(is.numeric(estimate) & length(estimate) == 1)) {

    #... then stop the function
    stop("'estimate' must be a numeric scalar")

  }

  #If se is not a numeric scalar...
  if (!(is.numeric(se) & length(se) == 1)) {

    #... then stop the function
    stop("'se' must be a numeric scalar")

  }
  #If se is not greater than zero...
  if (se <= 0) {

    #... then stop the function
    stop("'se' must be strictly greater than 0")

  }

  #If alpha is not a numeric scalar...
  if (!(is.numeric(alpha) & length(alpha) == 1)) {

    #... then stop the function
    stop("'alpha' must be a numeric scalar")

  }
  #If alpha is not between 0 and 0.5...
  if (alpha <= 0 | alpha >= 0.5) {

    #... then stop the function
    stop("'alpha' must be strictly between 0 and 0.5")

  }

  #If power is not a numeric scalar...
  if (!(is.numeric(power) & length(power) == 1)) {

    #... then stop the function
    stop("'power' must be a numeric scalar")

  }
  #If power is not between 0.5 and 1...
  if (power <= 0.5 | power >= 1) {

    #... then stop the function
    stop("'power' must be strictly between 0.5 and 1")

  }

  #If ROPE is not of length 1 or 2...
  if (!(length(ROPE) %in% c(1, 2))) {

    #... then stop the function
    stop("'ROPE' must be of length 1 or 2")

  }
  #If ROPE is of length 1...
  if (length(ROPE) == 1) {

    #If ROPE is not a numeric scalar...
    if (!(is.numeric(ROPE) & length(ROPE) == 1) | ROPE <= 0) {

      #... then stop the function
      stop("If 'ROPE' is of length 1, then it must be a numeric scalar")

    }
    #If ROPE is not greater than zero...
    if (ROPE <= 0) {

      #... then stop the function
      stop("If 'ROPE' is of length 1, then it must be strictly greater than zero")

    }

  }
  #If ROPE is of length 2...
  if (length(ROPE) == 2) {

    #If one of the two elements of ROPE is not of length 1...
    if (length(ROPE[1]) != 1 | length(ROPE[2]) != 1) {

      #... then stop the function
      stop("If 'ROPE' is of length 2, then its two elements must each be of length 1")

    }
    #If the first and second element of ROPE are equal to one another...
    if (ROPE[1] == ROPE[2]) {

      #... then stop the function
      stop("If 'ROPE' is of length 2, then its two elements must not be equal to one another")

    }
    #If one of the two elements of ROPE are non-numeric...

    if (!is.numeric(ROPE[1]) | !is.numeric(ROPE[2])) {

      #... then stop the function
      stop("All elements of 'ROPE' must be numeric")

    }

  }

  #If ROPE is a positive numeric scalar...
  if (length(ROPE) == 1) {

    #Generate a symmetric region of practical equivalence of length ROPE around zero
    ROPE = c(-ROPE, ROPE)

  }
  #If ROPE is a list of two different numeric scalars...
  if (length(ROPE) == 2) {

    #Generate a region of practical equivalence of [min(ROPE), max(ROPE)]
    ROPE = c(min(ROPE), max(ROPE))

  }

  #Generate bounds dataframe
  bounds = as.data.frame(matrix(nrow = 4, ncol = 2))
  colnames(bounds) = c("Lower Bound", "Upper Bound")
  rownames(bounds) = c(paste0(round((1 - alpha)*100, 3), "% equivalence confidence interval (ECI)"),
                       paste0(round((1 - alpha)*100, 3), "% confidence interval (CI)"),
                       paste0(round((1 - alpha)*100, 3), "% ECI adjusted to ", round(power*100, 3), "% power"),
                       paste0(round((1 - alpha)*100, 3), "% CI adjusted to ", round(power*100, 3), "% power"))

  #If the estimate is exactly midway between the lower and upper bounds of the ROPE...
  if (estimate == (ROPE[1] + ROPE[2])/2) {

    #... then select the upper boundary as the relevant bound
    bound = ROPE[2]

  }
  #... otherwise...
  else {

    #... the closer bound to the estimate is the relevant TOST bound
    bound = ROPE[which((c(abs(estimate - ROPE[1]), abs(estimate - ROPE[2])) == min(c(abs(estimate - ROPE[1]), abs(estimate - ROPE[2])))))]

  }

  ########################
  ##### SUB-ROUTINES #####
  ########################

  #If df is not provided...
  if ((is.na(df) | is.null(df))) {

    #Generate test dataframe
    test = as.data.frame(matrix(nrow = 3, ncol = 5))
    colnames(test) = c("ROPE Lower Bound", "ROPE Upper Bound", "z-statistic", "p-value", "Relevant")
    rownames(test) = c("Test: Estimate Bounded Above ROPE (Two-Sided)",
                       "Test: Estimate Bounded Within ROPE (TOST)",
                       "Test: Estimate Bounded Below ROPE (Two-Sided)")

    #Generate the bounds of the ECI
    bounds[1, 1] = estimate - qnorm(p = 1 - alpha)*se
    bounds[1, 2] = estimate + qnorm(p = 1 - alpha)*se

    #Generate the bounds of the CI
    bounds[2, 1] = estimate - qnorm(p = 1 - alpha/2)*se
    bounds[2, 2] = estimate + qnorm(p = 1 - alpha/2)*se

    #Generate the bounds of the power-adjusted ECI
    bounds[3, 1] = ROSE(estimate, se, alpha, power)$ROSE["Lower bound"]
    bounds[3, 2] = ROSE(estimate, se, alpha, power)$ROSE["Upper bound"]

    #Generate the bounds of the power-adjusted CI
    bounds[4, 1] = ROSE(estimate, se, alpha/2, power)$ROSE["Lower bound"]
    bounds[4, 2] = ROSE(estimate, se, alpha/2, power)$ROSE["Upper bound"]

    #Store the ROPE
    test[, 1] = rep(ROPE[1], 3)
    test[, 2] = rep(ROPE[2], 3)

    #Store the z-statistic and p-value of the two-sided test for bounding above the ROPE
    test[1, 3] = (estimate - ROPE[2])/se
    test[1, 4] = min(pnorm(test[1, 3], lower.tail = FALSE)*2, 1)

    #If the lower bound of the ROPE is the relevant TOST bound...
    if (bound == ROPE[1]) {

      #Store the z-statistic as estimate - min(ROPE) in standard error units
      test[2, 3] = (estimate - ROPE[1])/se
      #Store the p-value of the one-sided test in the upper tail
      test[2, 4] = pnorm(test[2, 3], lower.tail = FALSE)

    }
    #If the upper bound of the ROPE is the relevant TOST bound...
    if (bound == ROPE[2]) {

      #Store the z-statistic as estimate - max(ROPE) in standard error units
      test[2, 3] = (estimate - ROPE[2])/se
      #Store the p-value of the one-sided test in the lower tail
      test[2, 4] = pnorm(test[2, 3], lower.tail = TRUE)

    }

    #Store the z-statistic and p-value of the two-sided test for bounding below the ROPE
    test[3, 3] = (estimate - ROPE[1])/se
    test[3, 4] = min(pnorm(test[3, 3], lower.tail = TRUE)*2, 1)

    #Determine relevant test
    if (estimate > ROPE[2]) {
      relevant = 1
    }
    if (estimate >= ROPE[1] & estimate <= ROPE[2]) {
      relevant = 2
    }
    if (estimate < ROPE[1]) {
      relevant = 3
    }

    #Mark test relevance
    test[relevant, 5] = "Y"
    test[setdiff(c(1:3), relevant), 5] = "N"

    #If no p-value is below alpha...
    if (min(test[, 4]) >= alpha) {

      #Report that the result is inconclusive
      conclusion = "The significance of the estimate is inconclusive."

    }
    #If the test bounding the estimate above the ROPE is significant...
    if (test[1, 4] < alpha) {

      #Report the result
      conclusion = "The estimate is significantly bounded above the ROPE."

    }
    #If the TOST procedure is significant...
    if (test[2, 4] < alpha) {

      #Report the result
      conclusion = "The estimate is significantly bounded within the ROPE."

    }
    #If the test bounding the estimate below the ROPE is significant...
    if (test[3, 4] < alpha) {

      #Report the result
      conclusion = "The estimate is significantly bounded below the ROPE."

    }

  #Print citation disclaimer
  print(noquote("Asymptotically approximate equivalence confidence intervals and three-sided testing (TST) results reported"))
  print(noquote("If using for academic/research purposes, please cite the paper underlying this program:"))
  print(noquote("Fitzgerald, Jack (2024). The Need for Equivalence Testing in Economics. Institute for Replication Discussion Paper Series No. 125. https://www.econstor.eu/handle/10419/296190."))
  #Store output
  output = list(bounds, test, conclusion)
  names(output) = c("bounds", "test", "conclusion")
  #Return bounds
  return(output)

  }

  #If df is provided...
  if (!(is.na(df) | is.null(df))) {

    #If df is not a positive numeric scalar...
    if (!(is.numeric(df) & length(df) == 1) | df <= 0) {

      #... then stop the function
      stop("'df' must be a positive numeric scalar")

    }

    #Generate test dataframe
    test = as.data.frame(matrix(nrow = 3, ncol = 5))
    colnames(test) = c("ROPE Lower Bound", "ROPE Upper Bound", "t-statistic", "p-value", "Relevant")
    rownames(test) = c("Test: Estimate Bounded Above ROPE (Two-Sided)",
                       "Test: Estimate Bounded Within ROPE (TOST)",
                       "Test: Estimate Bounded Below ROPE (Two-Sided)")

    #Generate the bounds of the ECI
    bounds[1, 1] = estimate - qt(p = 1 - alpha, df = df)*se
    bounds[1, 2] = estimate + qt(p = 1 - alpha, df = df)*se

    #Generate the bounds of the ECI
    bounds[2, 1] = estimate - qt(p = 1 - alpha/2, df = df)*se
    bounds[2, 2] = estimate + qt(p = 1 - alpha/2, df = df)*se

    #Generate the bounds of the power-adjusted ECI
    bounds[3, 1] = ROSE(estimate, se, alpha, power, df)$ROSE["Lower bound"]
    bounds[3, 2] = ROSE(estimate, se, alpha, power, df)$ROSE["Upper bound"]

    #Generate the bounds of the power-adjusted CI
    bounds[4, 1] = ROSE(estimate, se, alpha/2, power, df)$ROSE["Lower bound"]
    bounds[4, 2] = ROSE(estimate, se, alpha/2, power, df)$ROSE["Upper bound"]

    #Store the ROPE
    test[, 1] = rep(ROPE[1], 3)
    test[, 2] = rep(ROPE[2], 3)

    #Store the t-statistic and p-value of the two-sided test for bounding above the ROPE
    test[1, 3] = (estimate - ROPE[2])/se
    test[1, 4] = min(pt(test[1, 3], df = df, lower.tail = FALSE)*2, 1)

    #If the lower bound of the ROPE is the  TOST bound...
    if (bound == ROPE[1]) {

      #Store the t-statistic as estimate - min(ROPE) in standard error units
      test[2, 3] = (estimate - ROPE[1])/se
      #Store the p-value of the one-sided test in the upper tail
      test[2, 4] = pt(test[2, 3], df = df, lower.tail = FALSE)

    }
    #If the upper bound of the ROPE is the  TOST bound...
    if (bound == ROPE[2]) {

      #Store the t-statistic as estimate - max(ROPE) in standard error units
      test[2, 3] = (estimate - ROPE[2])/se
      #Store the p-value of the one-sided test in the lower tail
      test[2, 4] = pt(test[2, 3], df = df, lower.tail = TRUE)

    }

    #Store the t-statistic and p-value of the two-sided test for bounding below the ROPE
    test[3, 3] = (estimate - ROPE[1])/se
    test[3, 4] = min(pt(test[3, 3], df = df, lower.tail = TRUE)*2, 1)

    #Determine relevant test
    if (estimate > ROPE[2]) {
      relevant = 1
    }
    if (estimate >= ROPE[1] & estimate <= ROPE[2]) {
      relevant = 2
    }
    if (estimate < ROPE[1]) {
      relevant = 3
    }
    
    #Mark test relevance
    test[relevant, 5] = "Y"
    test[setdiff(c(1:3), relevant), 5] = "N"

    #If no p-value is below alpha...
    if (min(test[, 4]) >= alpha) {

      #Report that the result is inconclusive
      conclusion = "The significance of the estimate is inconclusive."

    }
    #If the test bounding the estimate above the ROPE is significant...
    if (test[1, 4] < alpha) {

      #Report the result
      conclusion = "The estimate is significantly bounded above the ROPE."

    }
    #If the TOST procedure is significant...
    if (test[2, 4] < alpha) {

      #Report the result
      conclusion = "The estimate is significantly bounded within the ROPE."

    }
    #If the test bounding the estimate below the ROPE is significant...
    if (test[3, 4] < alpha) {

      #Report the result
      conclusion = "The estimate is significantly bounded below the ROPE."

    }


    #Print citation disclaimer
    print(noquote("Exact (equivalence) confidence intervals (ECIs) and three-sided testing (TST) results reported"))
    print(noquote("If using for academic/research purposes, please cite the paper underlying this program:"))
    print(noquote("Fitzgerald, Jack (2024). The Need for Equivalence Testing in Economics. Institute for Replication Discussion Paper Series No. 125. https://www.econstor.eu/handle/10419/296190."))
    #Store output
    output = list(bounds, test, conclusion)
    names(output) = c("bounds", "test", "conclusion")
    #Return bounds
    return(output)

    }

  }
