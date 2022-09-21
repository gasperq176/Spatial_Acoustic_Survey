#=====================================#
#      Main simulation function       #
#=====================================#
#' @export
scr.sim = function(lambda0, sigma, traps,
                   density = 50,
                   distr = "pois",
                   limits = list(xlim = NULL, ylim = NULL),
                   draw = FALSE,
                   binom.n = 1,
                   acoustic = FALSE,
                   lambda_c = NULL,
                   toa = FALSE,
                   sigma_toa = 0.002,
                   speed_sound = 330,
                   ...) {
  ## Setting up the total survey area
  ##  - Survey area (vs. trap area) is based on extreme trap co-ordinates
  ##    - Ranges are taken from each x and y column
  ##    - Co-ordinates are extended by (5 * sigma); stored as survey area box
  ##  - Then some error handling (i.e. only x or only y coordinates missing)
  if(is.null(limits$xlim) & is.null(limits$ylim)) {
    limits = list(xlim = range(traps[, 1], na.rm = TRUE),
                  ylim = range(traps[, 2], na.rm = TRUE))
    limits = lapply(limits, "+", c(-1, 1) * (5 * sigma))
  } else if(is.null(limits$xlim) & !is.null(limits$ylim)) {
    stop("X co-ordinates missing")
  } else if(!is.null(limits$xlim) & is.null(limits$ylim)) {
    stop("Y coordinates missing")
  }

  ## Generating random points on the area
  ##  - Area calculated from xlim and ylim
  ##  - Total number of animals N ~ Pois(DA)
  ##    - Density (D) is given per hectare, so divide by 10,000 to get per metre
  area = (limits$xlim[2] - limits$xlim[1]) * (limits$ylim[2] - limits$ylim[1])
  n = rpois(1, (density / 10000) * area)
  coords = pointgen(n, xlim = limits$xlim, ylim = limits$ylim)

  ## Plotting the traps and activity centres, provided draw = TRUE
  if(draw) {
    ## Setting up the survey area
    plot.new()
    plot.window(xlim = limits$xlim,
                ylim = limits$ylim,
                xaxs = "i", yaxs = "i")
    box()

    ## Setting up the traps
    points(traps, pch = 3, col = "red")

    ## Plotting the activity centres
    points(coords, ...)
  }

  ## Acoustic captures
  ## - Each animal has a matrix of captures - 1 per call
  ##    - Obviously, if none of the "traps" captured a call, it isn't recorded
  ## - Have some average # of calls (lambda_c)
  ##    - Randomly generate (lambda_c) capture histories for each call
  if(acoustic && is.null(lambda_c)) {
    stop("Acoustic calls need a mean number of calls (lambda_c)")
  }

  ## Setting up the random count generation - depending on the distribution
  if(distr == "pois" & acoustic == FALSE) {
    rDistr = paste0("r", distr, "(length(d), lambda0 * exp(-d^2 / (2 * sigma^2)))")
  } else if(distr == "bernoulli" | distr == "binom" | acoustic) {
    rDistr = paste0("r", "binom", "(length(d),", binom.n, ", lambda0 * exp(-d^2 / (2 * sigma^2)))")
  } else if(distr == "negbin" | distr == "nbinom") {
    ifelse(!is.null(list(...)$size), size <- list(...)$size, size <- 2)
    rDistr = paste0("r", "nbinom", "(length(d), mu = lambda0 * exp(-d^2 / (2 * sigma^2)), size = size)")
  }

  ## Calculating the distances between each activity centre and every trap
  ## - Distances are only calculated once; more efficient than calculating it in the C++ files.
  distances = eucdist_nll(coords, traps)

  ## Filling the omega matrix row-by-row
  ##  - Counts are randomly generated based on the specified distribution
  ##  - A row will only be added if its sum > 0; i.e. if at least 1 of the traps had a detection.
  ##  - The simulated counts are removed at the end of each loop, just to keep things tidy.
  ##
  ## Acoustic and regular SCR are differentiated here
  ## - If it's acoustic, then simCounts is repeated rpois(1, lambda_c) times and bound
  ##    - Note that if rpois(1, lambda_c) == 0, then you need to generate it again
  omega = id = toa.mat = NULL
  if(acoustic) {
    idIterator = 1
    for(i in 1:nrow(coords)) {
      d = distances[i, ]

      ## Generating a random number of calls
      nCalls = rpois(1, lambda_c)
      while(nCalls == 0) {
        nCalls = rpois(1, lambda_c)
      }

      ## Generating [nCalls] detection vectors for a given animal's calls
      simCounts = t(replicate(nCalls, eval(parse(text = rDistr))))
      #parseEval(rDistr)))

      ## Dropping out the 0 counts
      ## Also storing animal labels in a separate vector
      ## - Label only added if counts > 1
      simCounts = simCounts[as.logical(rowSums(simCounts)), ]
      if(length(simCounts) / nrow(traps) > 0) {
        ## Adding label
        id = c(id, rep(idIterator, length(simCounts) / nrow(traps)))

        ## Keeping track of labels
        idIterator = idIterator + 1
      }

      ## Generating a "time of arrival" matrix
      ## - [nrow(simCounts)] rows of times of arrival are generated
      ## - TOA = (1 / [speed of sound]) * distance + error (sigma_toa)
      if(toa && (length(simCounts) / nrow(traps) > 0)) {
        toa.mat = rbind(toa.mat,
                        t(replicate(length(simCounts) / nrow(traps),
                                    (((1/speed_sound) * d) + rnorm(length(d), 0, sigma_toa)))) * simCounts)
      }

      ## Binding the matrices to the "grand matrix"
      omega = rbind(omega,
                    simCounts)

      ## Cleaning up
      rm(simCounts)
    }
    #omega = cbind(omega, id)
  } else {
    for(i in 1:nrow(coords)) {
      d = distances[i, ]
      simCounts = eval(parse(text = rDistr))

      ## Dropping out the 0 counts
      if(sum(simCounts) != 0) {
        omega = rbind(omega,
                      simCounts)
      }

      ## Cleaning up
      rm(simCounts)
    }
  }

  ## Removing the row names given as a result of rbind()
  rownames(omega) = NULL

  ## Converting the count data to binary, if the count type = "binary"
  ## Acoustic counts are also returned as binary
  if(acoustic) {
    omega = cbind(omega, id)
  } #else if(distr == "binom") {
## THIS IS REDUNDANT
   # omega = ifelse(omega > 0, 1, 0)
  #}

  ## Returning the result
  ## - Returns either omega matrix or list of omega matrix and times of arrival
  if(toa) {
    list("omega" = omega,
         "toa" = toa.mat)
  } else {
    omega
  }
}

#==========================================================================#
#==========================================================================#

#=====================================#
#    make.capthist format function    #
#=====================================#
#' @export
toCapthist = function(captures) {
  formatted = NULL
  for(rowNum in 1:nrow(captures)) {
    ## Control flow for matrices with/without names for traps
    ##  - If captures DON'T have names for traps, will expand the trap/individual combo
    ##    - Specifically, the trap given is the trap name
    ##  - If captures have names for traps, then it'll just give the trap (column) number
    if(is.null(colnames(captures))) {
      call = "rep(which(captures[rowNum, ] != 0), captures[rowNum, which(captures[rowNum, ] != 0)])"
    } else {
      call = "rep(as.numeric(names(which(captures[rowNum, ] != 0))), captures[rowNum, which(captures[rowNum, ] != 0)])"
    }
    trapNums = eval(parse(text = call))

    ## Putting together all of the columns in the final matrix
    ## Also converting it to a data frame, with column names
    formatted = rbind(formatted, matrix(c(rep(1, length(trapNums)),
                                          rep(rowNum, length(trapNums)),
                                          rep(1, length(trapNums)),
                                          trapNums),
                                        ncol = 4))
  }
  formatted = as.data.frame(formatted)
  colnames(formatted) = c("session", "ID", "occasion", "trap")
  formatted
}

#==========================================================================#
#==========================================================================#

#=====================================#
#      read.traps format function     #
#=====================================#
#' @export
toTraps = function(traps) {
  colnames(traps) = c("x", "y")
  data.frame(testID = 1:nrow(traps), traps)
}

#==========================================================================#
#==========================================================================#

#=====================================#
#          Package postamble          #
#=====================================#

#' @import Rcpp
#' @useDynLib scr
NULL
