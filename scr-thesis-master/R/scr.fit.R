#=====================================#
#       Main fitting function         #
#=====================================#
#' @export
scr.fit = function(capthist, traps, mask,
                   start = NULL, acoustic = FALSE, binom = FALSE,
                   toa = NULL, fix.g0 = FALSE, speed_sound = 330, trace = FALSE, method = "Nelder-Mead") {
    ## General error/exception handling
    if(is.null(start)) {
        warning("Initial pararameter values defaulting to c(D = 50, lambda0 = 5, sigma = 15)")
        start = c(50, 5, 15)
    } else if(is.null(start) && acoustic) {
        warning("Initial pararameter values defaulting to c(D = 50, lambda0 = 5, sigma = 15), lambda_c = 10")
        start = c(50, 5, 15, 10)
    } else if(is.null(start) && !is.null(toa)) {
        warning("Initial pararameter values defaulting to c(D = 50, lambda0 = 5, sigma = 15), lambda_c = 10, sigma_toa = 0.002")
        start = c(50, 5, 15, 10, 0.002)
    }
    if(length(start) == 4 && acoustic == FALSE) {
        warning("Data treated as acoustic captures (4 start parameters)")
    }
    if(!is.null(toa) && acoustic == FALSE) {
        warning("Data treated as acoustic captures (!is.null(toa))")
        acoustic = TRUE
    }
    if(is.null(toa) && length(start) == 5) {
        warning("Give time of arrival matrix")
    } else if(!is.null(toa) && length(start) != 5) {
        warning("Check time of arrival matrix has corresponding start parameter")
    }
    
    ## Checking to see if things need unpacking
    if(class(capthist) == "capthist") {
        capthist = capthist[, 1, ]
        traps = capthist$traps
    }
    if (fix.g0){
        g0.fixed <- start[2]
    } else {
        g0.fixed <- 1
    }
    ## Transforming the start values
    ## - Note that:
    ##      - acoustic captures have one extra parameter (lambda_c)
    ##      - time of arrival (toa) has parameter sigma_toa
    ## - Machine minimum subtracted so as to avoid log errors
    start = start - .Machine$double.xmin
    if(acoustic | binom) {
        start = c(log(start[1]),
                  qlogis(start[2]),
                  log(start[3:length(start)]))
    } else {
        ## 3 parameters, POISSON: all logged
        ## - May have 4th parameter (sigma_toa)
        start = log(start)
    }
    if (fix.g0){
        start <- start[-2]
    }
    ## Calculating mask distances before giving to optim
    ## - More efficient
    ## - Note: for some reason,the basic scr.nll requires eucdist_nll to have (traps, mask)
    maskDists = eucdist_nll(mask, traps)
    
    ## Setting `use_toa` for the likelihood
    ## - If the TOA matrix hasn't been provided:
    ##    - use_toa = FALSE
    ##    - toa (matrix in likelihood) set as dummy-matrix
    ## - Otherwise:
    ##    - use_toa = TRUE
    ##    - toa given
    if(is.null(toa)){
        use_toa = FALSE
        toa_ssq = toa = matrix()
    } else {
        use_toa = TRUE
        
        ## Calculating the sums of squares of TOA
        ## - Note: call to eucdist_nll has arguments reversed c.f. maskDists.
        toa_ssq = make_toa_ssq(toa, eucdist_nll(traps, mask), speed_sound)
    }
    ## Taking the capthist, traps, and mask and maximising likelihood
    ## - Note that likelihood function changes for acoustic
    ## - Requires start values
    if(acoustic) {
        fit = optim(start, scr.nll.acoustic,
                    caps = capthist,
                    traps = traps,
                    mask = mask,
                    maskDists = maskDists,
                    toa = toa,
                    toa_ssq = toa_ssq,
                    use_toa = use_toa,
                    is_g0_fixed = fix.g0,
                    g0_fixed = g0.fixed,
                    trace = trace,
                    method = method,
                    hessian = TRUE)
    } else {
        fit = optim(start, scr.nll,
                    caps = capthist,
                    traps = traps,
                    mask = mask,
                    maskDists = maskDists,
                    binom = binom,
                    method = method,
                    hessian = TRUE)
    }
    
    ## Calculating confidence intervals
    ## - Using the (sqrt of) diagonals of (-ve) Hessian obtained from optim
    ##    - i.e. Information matrix
    ## - Wald CIs calculated by sapply() loop
    ##    - Loops through each of fitted parameters and calculates lower/upper bounds
    ## Note: fitted pars must be on LINK scale
    ##     : if matrix is singular, none of the SEs or CIs are calculated (inherits/try statement)
    fittedPars = fit$par
    if(inherits(try(solve(fit$hessian), silent = TRUE), "try-error")) {
        ## Hessian is singular
        warning("Warning: singular hessian")
        ## SE and Wald CIs not calculated
        se = NA
        waldCI = matrix(NA, nrow = length(fittedPars), ncol = 2)
        ## But columns still need to be returned (if/when simulations are run)
        cnames = c("Estimate", "SE", "Lower", "Upper")
    } else {
        ## Calculating basic Wald CIs
        se = sqrt(diag(solve(fit$hess)))
        waldCI = t(sapply(1:length(fittedPars),
                          function(i) fittedPars[i] + (c(-1, 1) * (qnorm(0.975) * se[i]))))
        ## Back-transforming the confidence limits, depending on whether we're using lambda0 or g0
        if(acoustic || binom) {
            waldCI = rbind(exp(waldCI[1, ]),
                           plogis(waldCI[2, ]),
                           exp(waldCI[3:length(fittedPars), ]))
        } else {
            waldCI = exp(waldCI)
        }
        
        ## Using the delta method to get the standard errors
        ## - G = jacobian matrix of partial derivatives of back-transformed
        ##    - i.e. log(D) -> exp(D) -- deriv. --> exp(D)
        ##    - Note: 1st deriv of plogis (CDF) = dlogis (PDF)
        G = diag(length(fittedPars)) * c(exp(fittedPars[1]),
                                         dlogis(fittedPars[2]),
                                         exp(fittedPars[3:length(fittedPars)]))
        se = sqrt(diag(G %*% solve(fit$hessian) %*% t(G)))
    }
    
    
    ## Returning the fitted parameters in a named vector
    ## - First checks to see if TOA is being used,
    ##    then inserts par names in front of "sigma_toa"
    parNames = NULL
    if(!(is.null(toa) || is.na(toa))) {
        parNames = "sigma_toa"
    }
    if(acoustic) {
        parNames = c("lambda_c", parNames)
    }
    if(binom || acoustic) {
        if (fix.g0){
            parNames = c("D", "sigma", parNames)
            fittedPars = c(exp(fittedPars[1]),
                           exp(fittedPars[2:length(fittedPars)]))
        } else {
            parNames = c("D", "g0", "sigma", parNames)
            fittedPars = c(exp(fittedPars[1]),
                           plogis(fittedPars[2]),
                           exp(fittedPars[3:length(fittedPars)]))
        }
        
    } else {
        parNames = c("D", "lambda0", "sigma", parNames)
        
        fittedPars = exp(fittedPars)
    }
    results = cbind(fittedPars, se, waldCI)
    dimnames(results) = list(parNames, c("Estimate", "SE", "Lower", "Upper"))
    results
}

#==========================================================================#
#==========================================================================#


#=====================================#
#          Package postamble          #
#=====================================#

#' @import Rcpp
#' @useDynLib scr
NULL
