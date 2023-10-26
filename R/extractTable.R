# Function to screen for Compounds in the Exported Waters file

extractTable <- function(lcms) {

  compoundsindex <- str_which(lcms[,1], pattern = "Compound .*:", negate = FALSE)  # Check that your compounds match the position

  firstcomp <- as.numeric(compoundsindex[1])  #Position of the first compound
  lastcomp <- as.numeric(compoundsindex[length(compoundsindex)]) # Position of the last compound
  row_betweencomp <- as.numeric(compoundsindex[2]) - as.numeric(compoundsindex[1]) # Rows between compounds
  injections <- row_betweencomp-4 #Number of samples is calculated by subtraction

  results_set <- as.data.frame(matrix(data = NA, ncol = ncol(lcms), nrow = injections*length(compoundsindex))) #Initializate a matrix
  colnames(results_set) <- colnames(lcms)

  row_n <-  1

  for (i in compoundsindex) {
    comp_id <- strsplit(lcms[i,1],split=':  ')[[1]][2]   # Select the line of the compound and split /Extract name from list
    print(comp_id)
    jinit <- i+3  # Position of the first sample
    jend <- i+3+injections-1   # Position of the last sample

    for (j in jinit:jend) {
      d <- lcms[j,] %>% mutate(ID = comp_id)

      results_set[row_n,] <- d
      row_n <- row_n +1

    # results_set <- rbind(results_set, d)
    }
  }

  return(results_set)

}





