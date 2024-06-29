lddtest = function(runvar, cutpoint, epsilon, alpha = 0.05, power = 0.8, bin = NULL, bw = NULL, verbose = FALSE, plot = TRUE) {
  
  #If epsilon is not a numeric scalar...
  if (!(is.numeric(epsilon) & length(epsilon) == 1)) {
    
    #... then stop the function
    stop("'epsilon' must be a numeric scalar")
    
  }
  #If epsilon is not greater than one...
  if (epsilon <= 1) {
    
    #... then stop the function
    stop("'epsilon' must be strictly greater than 1")
    
  }
  
  #Run DCdensity
  mccrary = DCdensity(runvar, cutpoint, bin = bin, bw = bw, verbose = verbose, plot = plot, ext.out = TRUE, htest = FALSE)
  estimate = mccrary$theta
  se = mccrary$se
  
  #Organize ROPE
  ROPE = c(-log(epsilon), log(epsilon))
  
  #Generate bounds dataframe
  bounds = as.data.frame(matrix(nrow = 2, ncol = 2))
  colnames(bounds) = c("Lower Bound", "Upper Bound")
  rownames(bounds) = c(paste0(round((1 - alpha)*100, 3), "% equivalence confidence interal (ECI)"),
                       paste0(round((1 - alpha)*100, 3), "% region of statistical equivalence (ROSE) with ", round(power*100, 3), "% power"))
  
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
  
  #Generate the bounds of the ECI
  bounds[1, 1] = estimate - qnorm(p = 1 - alpha)*se
  bounds[1, 2] = estimate + qnorm(p = 1 - alpha)*se
  
  #Generate the bounds of the ROSE
  bounds[2, 1] = ROSE(estimate, se, alpha, power)$ROSE["Lower bound"]
  bounds[2, 2] = ROSE(estimate, se, alpha, power)$ROSE["Upper bound"]
  
  #Generate test dataframe
  test = as.data.frame(matrix(nrow = 1, ncol = 6))
  colnames(test) = c("Epsilon Lower Bound", "Epsilon Upper Bound", "Theta", "SE", "Equivalence z-statistic", "p-value")
  
  #Store the epsilon boundaries
  test[1, 1] = ROPE[1]
  test[1, 2] = ROPE[2]
  
  #Store the estimate and SE
  test[1, 3] = estimate
  test[1, 4] = se
  
  #If the lower bound of the ROPE is the relevant TOST bound...
  if (bound == ROPE[1]) {
    
    #Store the z-statistic as estimate - min(ROPE) in standard error units
    test[1, 5] = (estimate - ROPE[1])/se
    #Store the p-value of the one-sided test in the upper tail
    test[1, 6] = pnorm(test[1, 5], lower.tail = FALSE)
    
  }
  #If the upper bound of the ROPE is the relevant TOST bound...
  if (bound == ROPE[2]) {
    
    #Store the z-statistic as estimate - max(ROPE) in standard error units
    test[1, 5] = (estimate - ROPE[2])/se
    #Store the p-value of the one-sided test in the lower tail
    test[1, 6] = pnorm(test[1, 5], lower.tail = TRUE)
    
  }
  
  #If the p-value is below alpha...
  if (test[1, 6] <= alpha) {
    
    #... then conclude the LDD is significantly bounded
    conclusion = paste0("The running variable's density discontinuity at the cutoff is significantly bounded beneath a ratio of ",
                        epsilon,
                        " at a ",
                        round(alpha*100, 3),
                        "% significance level.")
    
  }
  
  #Otherwise...
  if (test[1, 6] > alpha) {
    
    #... then conclude the LDD is significantly bounded
    conclusion = paste0("The running variable's density discontinuity at the cutoff is NOT significantly bounded beneath a ratio of, ",
                        epsilon,
                        " at a ",
                        round(alpha*100, 3),
                        "% significance level.")
    
  }
  
  #Print citation disclaimer
  print(noquote("Please cite the paper underlying this program:"))
  print(noquote("Fitzgerald, Jack (2024). Manipulation Tests in Regression Discontinuity: The Need for Equivalence Testing. Working paper. https://jack-fitzgerald.github.io/files/RDD_equivalence.pdf."))
  #Store output
  output = list(bounds, test, conclusion)
  names(output) = c("bounds", "test", "conclusion")
  #Return bounds
  return(output)
  
}
