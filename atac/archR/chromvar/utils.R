
.customDeviations <- function(
  countsMatrix = NULL,
  annotationsMatrix = NULL,
  backgroudPeaks = NULL,
  expectation = NULL,
  prefix = "",
  out = c("deviations", "z"),
  threads = 1,
  verbose = TRUE
  ){

  # sanity checks
  stopifnot(nrow(countsMatrix) == nrow(backgroudPeaks))
  stopifnot(length(expectation) == nrow(countsMatrix))

  colData <- DataFrame(seq_len(ncol(countsMatrix)), row.names = colnames(countsMatrix))[,FALSE]
  norm_expectation <- expectation / sum(expectation) #Double check this sums to 1!
  countsPerSample <- Matrix::colSums(countsMatrix)

  d <- max(floor(ncol(annotationsMatrix)/20), 1)
  m <- 0
  results <- .safelapply(seq_len(ncol(annotationsMatrix)), function(x){
    if(x %% d == 0){
      m <- 1 #Print to console
    }
    if(x %% max(floor(d/5), 2) == 0){
      if(m != 1){
      }else{
        m <- 0 #Reset
      }
    }
    if(x %% max(c(d, 10)) == 0){
      gc()
    }
    .customDeviationsSingle(
      annotationsVector = annotationsMatrix[, x, drop=FALSE],
      countsMatrix = countsMatrix,
      backgroudPeaks = backgroudPeaks,
      countsPerSample = countsPerSample,
      expectation = norm_expectation,
      out = out,
      prefix = prefix
    )
  }, threads = threads)
  cn <- colnames(countsMatrix)
  rm(countsMatrix)
  gc()

  # parse output
  if("z" %in% tolower(out)){
    z <- t(vapply(results, function(x) x[["z"]], rep(0, length(cn))))
  }else{
    z <- matrix(0, nrow = ncol(annotationsMatrix), ncol = length(cn))
  }
  if("deviations" %in% tolower(out)){
    dev <- t(vapply(results, function(x) x[["dev"]], rep(0, length(cn))))
  }else{
    dev <- matrix(0, nrow = ncol(annotationsMatrix), ncol = length(cn))
  }
  colnames(z) <- cn
  colnames(dev) <- cn

  #Check First
  nullOverlap <- is.null(results[[1]]$overlap)
  rowData <- lapply(seq_along(results), function(x){
      resx <- results[[x]]
      if(nullOverlap){
        data.frame(fractionMatches = resx$matches)
      }else{
        data.frame(fractionMatches = resx$matches, fractionBackgroundOverlap = resx$overlap)
      }
    }) %>% Reduce("rbind",.)
  rownames(rowData) <- colnames(annotationsMatrix)
  
  # craete output summarized experiment
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(deviations = dev, z = z),
    colData = colData,
    rowData = rowData
  )
  SummarizedExperiment::assays(se) <- SummarizedExperiment::assays(se)[tolower(out)]
  
  return(se)

}

.customDeviationsSingle <- function(
  annotationsVector = NULL,
  countsMatrix = NULL,
  countsPerSample = NULL,
  backgroudPeaks = NULL,
  out = c("deviations", "z"),
  expectation = NULL,
  threshold = 1,
  prefix = ""
  ){

  .binarizeMat <- function(mat = NULL){
    mat@x[mat@x > 0] <- 1
    mat
  }

  if (length(annotationsVector@x) == 0) {
      out <- list(
        z = rep(NA, ncol(countsMatrix)), 
        dev = rep(NA, ncol(countsMatrix)), 
        expFG = NA,
        expBG = NA,
        matches = 0, 
        overlap = NA
        )
    return(out)
  }

  outList <- tryCatch({

    ################################
    # Fore Ground Deviations
    ################################
    observed <- as.vector(Matrix::t(annotationsVector) %*% countsMatrix)
    expected <- as.vector(Matrix::t(annotationsVector) %*% expectation %*% countsPerSample)
    observed_deviation <- (observed - expected)/expected

    #Filter those with no matches at all
    fail_filter <- which(expected == 0)
    
    ################################
    # Back Ground Deviations
    ################################
    if("z" %in% tolower(out)){

      #Compute Background Null Per Iteration
      niterations <- ncol(backgroudPeaks)
      sampleMat <- Matrix::sparseMatrix(
          j = as.vector(backgroudPeaks[annotationsVector@i + 1, seq_len(niterations)]),
          i = rep(seq_len(niterations), each = length(annotationsVector@x)),
          x = rep(annotationsVector@x, niterations),
          dims = c(niterations, nrow(countsMatrix))
      )  
      sampled <- as.matrix(sampleMat %*% countsMatrix)
      sampledExpected <- sampleMat %*% expectation %*% countsPerSample
      sampledDeviation <- (sampled - sampledExpected)/sampledExpected
      bgOverlap <- Matrix::mean(.binarizeMat(sampleMat) %*% .binarizeMat(annotationsVector)) / length(annotationsVector@x)
      
      #Summary
      meanSampledDeviation <- Matrix::colMeans(sampledDeviation)
      sdSampledDeviation <- apply(as.matrix(sampledDeviation), 2, sd)

      #Norm Deviation
      normdev <- (observed_deviation - meanSampledDeviation)
      z <- normdev/sdSampledDeviation
      if (length(fail_filter) > 0) {
        z[fail_filter] <- NA
        normdev[fail_filter] <- NA
      }

    }else{

      #Compute Background Null Per Iteration
      niterations <- ncol(backgroudPeaks)
      sampleMat2 <- Matrix::sparseMatrix(
          j = as.vector(backgroudPeaks[annotationsVector@i + 1, seq_len(niterations)]),
          i = rep(1, niterations * length(annotationsVector@x)),
          x = rep(annotationsVector@x, niterations),
          dims = c(1, nrow(countsMatrix))
      )
      sampled2 <- (sampleMat2 %*% countsMatrix)[1,]
      sampledExpected2 <- (sampleMat2 %*% expectation %*% countsPerSample)[1,]
      ######################
      # Equivalent to above
      # colMeans(sampled) - colMeans(sampledExpected))/colMeans(sampledExpected)
      ######################
      sampledDeviation2 <- (sampled2 - sampledExpected2)/sampledExpected2
      bgOverlap <- NA

      #Norm Deviation
      normdev <- (observed_deviation - sampledDeviation2)
      z <- NULL
      if (length(fail_filter) > 0) {
        normdev[fail_filter] <- NA
      }

    }

    outList <- list(
      z = z, 
      dev = normdev, 
      matches = length(annotationsVector@x) / nrow(countsMatrix), 
      overlap = bgOverlap
    )

    outList

  }, error = function(e){

    errorList <- list(
      annotationsVector = annotationsVector,
      observed = if(exists("observed", inherits = FALSE)) observed else "observed",
      expected = if(exists("expected", inherits = FALSE)) expected else "expected",
      sampleMat = if(exists("sampleMat", inherits = FALSE)) sampleMat else "sampleMat",
      sampleMat2 = if(exists("sampleMat", inherits = FALSE)) sampleMat2 else "sampleMat2",
      sampledDeviation = if(exists("sampledDeviation", inherits = FALSE)) sampledDeviation else "sampledDeviation",
      sampledDeviation2 = if(exists("sampledDeviation2", inherits = FALSE)) sampledDeviation2 else "sampledDeviation2",
      normdev = if(exists("normdev", inherits = FALSE)) normdev else "normdev",
      z = if(exists("z", inherits = FALSE)) z else "z"
    )

  })

  return(outList)

}


.safelapply <- function(..., threads = 1, preschedule = FALSE){

  if(threads > 1){

    o <- mclapply(..., mc.cores = threads, mc.preschedule = preschedule)

    errorMsg <- list()

    for(i in seq_along(o)){ #Make Sure this doesnt explode!
      if(inherits(o[[i]], "try-error")){
        capOut <- utils::capture.output(o[[i]])
        capOut <- capOut[!grepl("attr\\(\\,|try-error", capOut)]
        capOut <- head(capOut, 10)
        capOut <- unlist(lapply(capOut, function(x) substr(x, 1, 250)))
        capOut <- paste0("\t", capOut)
        errorMsg[[length(errorMsg) + 1]] <- paste0(c(paste0("Error Found Iteration ", i, " : "), capOut), "\n")
      }
    }

    if(length(errorMsg) != 0){

      errorMsg <- unlist(errorMsg)
      errorMsg <- head(errorMsg, 50)
      errorMsg[1] <- paste0("\n", errorMsg[1])
      stop(errorMsg)

    }

  } else{

    o <- lapply(...)

  }

  o

}