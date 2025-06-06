
construct_intervals <- function(T, alpha = sqrt(1/2), m = 4) {
  layers <- ceiling(log(T, 1/alpha))
  
  # Since we flatten at the end anyway, no need to store as a list
  intervals <- matrix(NA, nrow = 0, ncol = 2) 
  
  for (k in 1:layers) {
    # Slight variation in these calculations
    # Since I_1 is (1, T] not (0, T] (1 indexing in R)
    nk <- 2 * ceiling((1/alpha)^(k-1)) - 1
    lk <- (T) * alpha^(k - 1)
    sk <- ifelse(nk == 1, 0, ((T) - lk) / (nk - 1))
    
    # We can always end in this case
    if(ceiling(lk) < m) {break}
    
    for (i in 1:nk) {
      start <- floor((i - 1) * sk)
      end <- min(ceiling((i - 1) * sk + lk), T)
      
      # In case of lk < m, but ceiling(lk) = m
      # some intervals may still have length = m 
      if ((end - start) >= m) {
        intervals <- rbind(intervals, c(start, end))
      }
    }
    
  }
  # Eliminate any overlapping or duplicate intervals
  return(unique(intervals))
}

plot_intervals <- function(intervals) {
  # Prepare the plot area
  plot(1, type = "n", xlab = "Time", ylab = "Layer", 
       xlim = c(0, max(intervals[, 2])), 
       ylim = c(0, nrow(intervals)), yaxt = 'n', 
       main = "Interval Plot")
  
  # Set the initial y position
  y_pos <- 0
  height_increment <- 0.8  # Controls how much the y position increases for each interval
  
  # Iterate over each row of the intervals matrix (reversed order)
  for (i in nrow(intervals):1) {
    start <- intervals[i, 1]
    end <- intervals[i, 2]
    # Plot each interval as a horizontal line
    segments(x0 = start, x1 = end, y0 = y_pos, y1 = y_pos, col = "blue", lwd = 2)
    # Increase the y position for the next interval
    y_pos <- y_pos + height_increment
  }
}

cusum_on_intervals <- function(CUSUM, obj, intervals, verbose = TRUE, ...) {
  # Fit CUSUM statistic on intervals
  # Intervals should be a list, or a single vector (which is converted)
  
  # ... caries arguments for CUSUM(obj, s, e, t, ...) 
  # pass by name when calling
  # e.g. seeded_binary_segmentation(my_CUSUM, my_obj, 40, rank = hat.rank)
  
  if (is.vector(intervals)) {
    intervals <- matrix(intervals, nrow = 1)
  }
  
  results <- matrix(nrow = nrow(intervals), ncol = 4)
  colnames(results) = c("Candidates", "Gain", "Start", "End")
  
  for (i in 1:nrow(intervals)) {
    # Calculate candidate and maximal gain, keep s and e
    max_gain = 0
    candidate = NA
    
    epsilon <- ceiling((intervals[i, 2]-intervals[i, 1])/64)
    s = intervals[i, 1] + epsilon
    e = intervals[i, 2] - epsilon

    if ((e-s) > 1) {
      for (t in (s+1):(e-1)) {
        # Calculate candidate and maximal gain, keep s and e
        gain <- CUSUM(obj, s, e, t, verbose, ...)  # Pass additional arguments to CUSUM 
        if (gain > max_gain) {
          max_gain <- gain
          candidate <- t
        }
      }
      if (verbose) {
        print(paste0("candidate = ", candidate, ", max_gain = ", max_gain, ", s = ", s, ", e = ", e, ", interval = ", i, "."))
      }
    }
    
    results[i, ] <- c(candidate, max_gain, s, e)
  }
  
  return(results[!is.na(results[, 1]), , drop = FALSE])
  
}

seeded_binary_seg <- function(CUSUM, obj, T, threshold = NULL, alpha = sqrt(1/2), m = 4, method = "Greedy", 
                              CUSUM_res = NULL, verbose = TRUE, ...) {
  # Method can be "Greedy" or "Narrowest"
    # Greedy selects the highest Gain over the threshold first
    # Narrowest selects shorted interval with Gain over threshold first
  # Thresholds should be numeric vector, or left unspecified
    # If unspecified, no selection is done
    # Otherwise, solution path computed for each threshold. 
  # if CUSUM results on intervals already exist, they can be passed in to skip computation
  # ... caries arguments for CUSUM(obj, s, e, t, ...) 
    # pass by name when calling
    # e.g. seeded_binary_segmentation(my_CUSUM, my_obj, 40, rank = hat.rank)
  
  
  if (is.null(CUSUM_res)) {
    intervals <- construct_intervals(T, alpha, m)
    results <- cusum_on_intervals(CUSUM, obj, intervals, verbose, ...)
  } else {
    results <- CUSUM_res
  }
  
  # Implement selection
  solution_path = list(list(threshold = -1, results = results))
  if (is.null(threshold)) {return(solution_path)}
  
  threshold <- sort(threshold, decreasing = TRUE)
  if (method == "Greedy") {
    results <- results[order(results[, 2], decreasing = TRUE), ]
    
    # For greedy, index list is outside band loop
    # Since decreasing threshold will not eliminate past candidates
    index_list <- numeric()
    i = 1
    cur_thres_idx = 1
    
    if (verbose) {
      cat("Number of rows ", nrow(results), "\n")
      cat("Number of thresholds ", length(threshold), "\n")
    }

    while(i <= nrow(results) && cur_thres_idx <= length(threshold)){
      if (verbose) {cat("Interval", results[i, 3], "-", results[i, 4], " with threshold index", cur_thres_idx, "\n")}
      # Go to next threshold (if threshold too large)
      # Select candidate (if above threshold and not covered yet)
      # Next row (if above threshold but covered)
      
      if (results[i, 2] < threshold[cur_thres_idx]) {
        subset = results[index_list, , drop = FALSE]
        solution_path <- c(solution_path, 
                           list(list(threshold = threshold[cur_thres_idx], results = subset)))
        cur_thres_idx <- cur_thres_idx + 1
      } 
      else if (all(results[index_list, 1] < results[i, 3] | results[index_list, 1] > results[i, 4])) {
        index_list <- c(index_list, i)
        i <- i+1
      } 
      else {
        i <- i+1
      }
    }
    
    while (cur_thres_idx <= length(threshold)) {
      # This means all rows were above the current threshold
      # So we fill in the remaining threshold entries with all 
      subset = results[index_list, , drop = FALSE]
      solution_path <- c(solution_path, 
                         list(list(threshold = threshold[cur_thres_idx], results = subset)))
      cur_thres_idx <- cur_thres_idx + 1
    }

    return(solution_path)
  }
  
  else if (method == "Narrowest") {
    if (verbose) {
      cat("Number of rows ", nrow(results), "\n")
      cat("Number of thresholds ", length(threshold), "\n")
    }
    
    results <- results[order(results[, 4] - results[, 3], -results[, 2]), ]
    for (band in threshold) {
      if (verbose) {cat("Current threshold: ", band, "\n")}
      index_list <- numeric()
      i = 1
      while(i <= nrow(results)){
        if (verbose) {cat("\tInterval", results[i, 3], "-", results[i, 4], " has cp ", results[i, 1], " with gain ", results[i, 2], "\n")}
        if (results[i, 2] >= band && all(results[index_list, 1] < results[i, 3] | 
                                        results[index_list, 1] > results[i, 4])) {
          index_list <- c(index_list, i)
        }
        i <- i+1
      }
      subset = results[index_list, , drop = FALSE]
      solution_path <- c(solution_path, list(list(threshold = band, results = subset)))
    }
    return(solution_path)
  }
  else {
    stop("Method must be \"Greedy\" or \"Narrowest\"")
  }
  
}


# Example usage
#T <- 150
#alpha <- sqrt(1/2)
#intervals <- construct_intervals(T, alpha, m = 4)
#print(intervals)
#plot_intervals(intervals)

