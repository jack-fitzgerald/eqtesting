llddtest = function(runvar, data, cutpoint, epsilon, alpha = 0.05, cluster = "", bootstrap = FALSE, breps = 1000, bin = NULL, bw = NULL, verbose = FALSE, plot = TRUE) {
  
  #If runvar is not a single string...
  if (!is.character(runvar) | length(runvar) != 1) {
    
    #... then stop the function
    stop("'runvar' must be a single character string")
    
  }
  #If cluster is not a single string...
  if (!is.character(runvar) | length(runvar) != 1) {
    
    #... then stop the function
    stop("'cluster' must be a single character string")
    
  }
  #If data is not a data.frame...
  if (!is.data.frame(data)) {
    
    #... then stop the function
    stop("'data' must be a data.frame")
    
  }
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
  #If bootstrap repetitions are changed but bootstrap itself is not specified as TRUE...
  if (breps != 1000 & bootstrap == FALSE) {
    
    #... then stop the function
    stop("'breps' is changed from its default value but 'bootstrap' is specified as FALSE; change one or the other")
    
  }
  
  #Run DCdensity
  mccrary = DCdensity(data[, runvar], cutpoint, bin = bin, bw = bw, verbose = verbose, plot = plot, ext.out = TRUE, htest = FALSE)
  estimate = mccrary$theta
  se = mccrary$se
  
  #Organize ROPE
  ROPE = c(-log(epsilon), log(epsilon))
  
  #Generate bounds dataframe
  bounds = as.data.frame(matrix(nrow = 1, ncol = 2))
  colnames(bounds) = c("Lower Bound", "Upper Bound")
  rownames(bounds) = c(paste0(round((1 - alpha)*100, 3), "% equivalence confidence interal (ECI)"))
  
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
  
  #Generate test dataframe
  test = as.data.frame(matrix(nrow = 1, ncol = 8))
  colnames(test) = c("Epsilon Lower Bound", "Epsilon Upper Bound", "Theta", "SE", "ECI Lower Bound", "ECI Upper Bound", "Equivalence z-statistic", "p-value")
  
  #Store the epsilon boundaries
  test[1, "Epsilon Lower Bound"] = ROPE[1]
  test[1, "Epsilon Upper Bound"] = ROPE[2]
  
  #Store the estimate
  test[1, "Theta"] = estimate
  
  #If bootstrap is unnecessary...
  if (bootstrap == FALSE & breps == 1000 & is.na(cluster)) {
    
    #Store the SE
    test[1, "SE"] = se
    
    #Store the ECI bounds
    test[1, "ECI Lower Bound"] = estimate - qnorm(1 - alpha)*se
    test[1, "ECI Upper Bound"] = estimate + qnorm(1 - alpha)*se
    
    #If the lower bound of the ROPE is the relevant TOST bound...
    if (bound == ROPE[1]) {
      
      #Store the z-statistic as estimate - min(ROPE) in standard error units
      test[1, "Equivalence z-statistic"] = (estimate - ROPE[1])/se
      #Store the p-value of the one-sided test in the upper tail
      test[1, "p-value"] = pnorm(test[1, "Equivalence z-statistic"], lower.tail = FALSE)
      
    }
    #If the upper bound of the ROPE is the relevant TOST bound...
    if (bound == ROPE[2]) {
      
      #Store the z-statistic as estimate - max(ROPE) in standard error units
      test[1, "Equivalence z-statistic"] = (estimate - ROPE[2])/se
      #Store the p-value of the one-sided test in the lower tail
      test[1, "p-value"] = pnorm(test[1, "Equivalence z-statistic"], lower.tail = TRUE)
      
    }
    
  }
  
  #If bootstrap is necessary...
  if (bootstrap == TRUE | !is.na(cluster)) {
    
    #Initialize bootstrap LDD list and count
    ldd_list = rep(NA, breps)
    i = 1
    
    #If conventional bootstrap is necessary...
    if (cluster == "") {
      
      #Define bootstrap data
      bdata = data[which(!is.na(data[, runvar])), ]
      
    }
    
    #If cluster bootstrap is necessary
    if (cluster != "") {
      
      #Define bootstrap data
      bdata = data[which(!is.na(data[, runvar] & !is.na(data[, cluster]))), ]
      
    }
    
    #Set data.table
    setDT(bdata)
    
    #Until breps defined logarithmic density discontinuities are obtained...
    while (i <= breps) {
      
      #If conventional bootstrap is necessary...
      if (cluster == "") {
        
        #Resample rows of the dataset with replacement
        bsample <- bdata[sample(.N, replace = TRUE)]
        
        #Obtain the LDD estimate for the bootstrap sample, holding bin sizes and bandwidths constant
        bDCdensity = DCdensity(bsample[[runvar]], bw = mccrary$bw, bin = mccrary$binsize, ext.out = TRUE, plot = FALSE)
        
        #If the estimate is nonmissing...
        if (!is.nan(bDCdensity$theta)) {
          
          #Store the estimate
          ldd_list[i] = bDCdensity$theta
          
          #Up the count
          i = i + 1
          
        }
        
      }
      
      #If cluster bootstrap is necessary...
      if(cluster != "") {
        
        #Identify unique clusters
        unique_clusters <- unique(bdata[[cluster]])
        
        #Obtain resampling clusters
        sampled_clusters <- sample(unique_clusters, replace = TRUE)
        
        #Count occurrences of each cluster in the sample
        cluster_counts <- table(sampled_clusters)
        
        #Obtain bootstrap sample
        bsample <- do.call(rbind, lapply(names(cluster_counts), function(cl) {
          # Select rows for this cluster and repeat them as needed
          subset <- bdata[bdata[[cluster]] == cl]
          subset[rep(1:.N, times = cluster_counts[cl]), ]
        }))
        
        #Obtain the LDD estimate for the bootstrap sample, holding bin sizes and bandwidths constant
        bDCdensity = DCdensity(bsample[[runvar]], bw = mccrary$bw, bin = mccrary$binsize, ext.out = TRUE, plot = FALSE)
        
        #If the estimate is nonmissing...
        if (!is.nan(bDCdensity$theta)) {
          
          #Store the estimate
          ldd_list[i] = bDCdensity$theta
          
          #Up the count
          i = i + 1
          
        }
        
      }
      
    }
    
    #Compute the bootstrap standard error
    test[1, "SE"] = sd(ldd_list)
    
    #Compute the bootstrap ECI bounds
    test[1, "ECI Lower Bound"] = sort(ldd_list)[floor(alpha*breps)]
    test[1, "ECI Upper Bound"] = sort(ldd_list)[ceiling((1 - alpha)*breps)]
    
    #Compute the bootstrap p-value
    test[1, "p-value"] = sum(ifelse(ldd_list >= -log(epsilon) & ldd_list <= log(epsilon), 0, 1))/breps
    
  }
  
  #If the p-value is below alpha...
  if (test[1, "p-value"] <= alpha) {
    
    #... then conclude the LDD is significantly bounded
    conclusion = paste0("The running variable's density discontinuity at the cutoff is significantly bounded beneath a ratio of ",
                        epsilon,
                        " at the ",
                        round(alpha*100, 3),
                        "% significance level.")
    
  }
  
  #Otherwise...
  if (test[1, "p-value"] > alpha) {
    
    #... then conclude the LDD is NOT significantly bounded
    conclusion = paste0("The running variable's density discontinuity at the cutoff is NOT significantly bounded beneath a ratio of, ",
                        epsilon,
                        " at the ",
                        round(alpha*100, 3),
                        "% significance level.")
    
  }
  
  #Print citation disclaimer
  print(noquote("Please cite the paper underlying this program:"))
  print(noquote("Fitzgerald, Jack (2024). Manipulation Tests in Regression Discontinuity Design: The Need for Equivalence Testing. Institute for Replication Discussion Paper Series, No. 125. https://hdl.handle.net/10419/300277."))
  #Store output
  output = list(test, conclusion)
  names(output) = c("test", "conclusion")
  #Return bounds
  return(output)
  
}
