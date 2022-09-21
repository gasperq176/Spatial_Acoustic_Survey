#=====================================#
#       Log-likelihood function       #
#=====================================#
#' @export
scr.nll = function(pars, caps, traps, mask, maskDists, binom) {
  scr_nll(pars, caps, traps, mask, maskDists, binom)
}

#==========================================================================#
#==========================================================================#

#=====================================#
#  Log-likelihood function (acoustic) #
#=====================================#
#' @export
scr.nll.acoustic = function(pars, caps, traps, mask, maskDists, toa, toa_ssq, use_toa, is_g0_fixed, g0_fixed, trace) {
  nCalls = table(caps[, ncol(caps)])
  scr_nll_acoustic(pars = pars,
                   caps = caps[, -ncol(caps)],
                   traps = traps,
                   mask = mask,
                   maskDists = maskDists,
                   nCalls = nCalls,
                   toa = toa,
                   toa_ssq = toa_ssq,
                   use_toa = use_toa,
                   is_g0_fixed = is_g0_fixed,
                   g0_fixed = g0_fixed,
                   trace = trace)
}

#==========================================================================#
#==========================================================================#


#=====================================#
#          Package postamble          #
#=====================================#

#' @import Rcpp
#' @useDynLib scr
NULL
