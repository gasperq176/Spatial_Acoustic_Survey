#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;


// ============================== //
//        eucdist (matrix)        //
// ============================== //
/*
 * Calculating the Euclidean distance between a point and each trap.
 * Returns a vector of distances.
 * - Not exported.
 */
// [[Rcpp::export]]
NumericMatrix eucdist_nll(NumericMatrix points,
                          NumericMatrix traplocations) {
  NumericMatrix dists(points.nrow(), traplocations.nrow());
  for(int i = 0;  i < points.nrow(); i++) {
    for(int j = 0; j < traplocations.nrow(); j++) {
      dists(i, j) = sqrt(pow(points(i, 0) - traplocations(j, 0), 2.0)
                         + pow(points(i, 1) - traplocations(j, 1), 2.0));
    }
  }
  return dists;
}

// =================================================================================== //
// =================================================================================== //

// ============================== //
//            scr_nll             //
// ============================== //
/*
* Calculating the Euclidean distance between a point and each trap.
* Returns a vector of distances.
*/
// [[Rcpp::export]]
double scr_nll(NumericVector pars,
               NumericMatrix caps,
               NumericMatrix traps,
               NumericMatrix mask,
               NumericMatrix maskDists,
               bool binom) {
  // Storing/initialising (starting) parameter values.
  double D = exp(pars[0]);
  double g0;
  double sigma = exp(pars[2]);
  if(binom) {
    g0 = R::plogis(pars[1], 0, 1, 1, 0);//exp(pars[1]);
  } else {
    g0 = exp(pars[1]);
  }

  // Number of animals detected.
  int n = caps.nrow();
  // Number of traps. NB: NOT USED
  //int nTraps = traps.nrow();
  // Number of mask points.
  int nMask = mask.nrow();
  // Area of a single mask pixel.
  double area = mask.attr("area");

  /*
   * Constructing distance matrix.
   * - Element (i, j) gives dist. b/w ith mask pint and jth trap.
   */
  //NumericMatrix maskDists = eucdist_nll(mask, traps);


  /*
   * Constructing distance matrix.
   * - Element (i, j) gives dist. b/w ith mask pint and jth trap.
   */
  //NumericMatrix maskDists = eucdist_nll(mask, traps);

  /*
  * Constructing a detection probability matrix.
  * - Element (i, j) gives prob. of animal @ ith mask pt. being detected @ jth trap.
  * - Line that fills in maskProbs(i, j) is the detection function
  * - Need to include machine minimum so that we don't get any Infs
  *
  * Detection function for the non-binomial/binary data is the HHN
  * Binomial data have HN detection function
  */
  NumericMatrix maskProbs(maskDists.nrow(), maskDists.ncol());
  NumericMatrix maskER(maskDists.nrow(), maskDists.ncol());

  if(binom) {
    for(int i = 0; i < maskDists.nrow(); i++) {
      for(int j = 0; j < maskDists.ncol(); j++) {
        maskProbs(i, j) = g0 * exp(-pow(maskDists(i, j), 2.0) / (2 * pow(sigma, 2.0))) + DBL_MIN;
      }
    }
  } else {
    for(int i = 0; i < maskDists.nrow(); i++) {
      for(int j = 0; j < maskDists.ncol(); j++) {
        maskER(i, j) = g0 * exp(-pow(maskDists(i, j), 2.0) / (2 * pow(sigma, 2.0))) + DBL_MIN;
        maskProbs(i, j) = 1 - exp(-maskER(i, j));
      }
    }
  }

  /*
   * Constructing a detection probability vector
   * - ith element = P(animal @ ith mask pt. is detected by >= 1 trap)
   */
  NumericVector pAvoid(maskProbs.nrow());
  for(int i = 0; i < maskProbs.nrow(); i++) {
    pAvoid[i] = 1 - maskProbs(i, 0);
    for(int j = 1; j < maskProbs.ncol(); j++) {
      pAvoid[i] *= 1 - maskProbs(i, j);
    }
  }

  /*
   * Vector of probability of detection at each mask point.
   * - Pc(s) in the acoustic form
   */
  NumericVector pDetected = 1 - pAvoid;

  /*
   * Calculating likelihood contribution
   * - Contribution from each detected animal's capt. hist.
   */
  NumericVector fCapt(n);
  NumericVector logfCapt_givenS(nMask);
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < nMask; j++) {
      /*
       * Calculating log-probability of animal's capt. hist,
       * conditional on being at the jth mask point.
       * - Note that 'caps' and 'maskProbs' both have the same ncol();
       *    - i.e. # traps
       *  - Also note: R::dbinom(double x, double n, double p, int lg)
       *    - Where 'int lg' is 0 = F, 1 = T.
       */
      if(binom) {
        logfCapt_givenS[j] = R::dbinom(caps(i, 0), 1, maskProbs(j, 0), 1);
        for(int k = 1; k < caps.ncol(); k++) {
          logfCapt_givenS[j] += R::dbinom(caps(i, k), 1, maskProbs(j, k), 1);
        }
      } else {
        logfCapt_givenS[j] = R::dpois(caps(i, 0), maskProbs(j, 0), 1);
        for(int k = 1; k < caps.ncol(); k++) {
          logfCapt_givenS[j] += R::dpois(caps(i, k), maskProbs(j, k), 1);
        }
      }
    }
    // Summing probabilities over all mask points.
    fCapt[i] = sum(exp(logfCapt_givenS));
  }

  /*
   * Log-likelihood contribution from all capture histories
   * - Calculated by log of sum of individual likelihood contributions.
   */
  double logfCapt = sum(log(fCapt + DBL_MIN));

  // Calculating effective survey area (unused in likelihood).
  double esa = area * sum(pDetected);

  // Log-likelihood contribution from number of animals detected.
  double logf_n = R::dpois(n, D * esa, 1);

  // Overall log-likelihood.
  double logLik = logf_n + logfCapt - n * log(sum(pDetected));

  // Printing log-likelihood.
  std::cout << "LL: " << logLik << std::endl;
  
  return -logLik;
}

// =================================================================================== //
// =================================================================================== //
