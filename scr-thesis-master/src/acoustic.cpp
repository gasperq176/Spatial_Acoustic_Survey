#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// =================================================================================== //
// =================================================================================== //

// ============================== //
//            numNonZero          //
// ============================== //
/*
* Counts the number of non-zero elements in a vector
* - Not exported.
*/
// [[Rcpp::export]]
int numNonZero(NumericVector x) {
  int count = 0;
  for(int i = 0; i < x.length(); i++) {
    if(x[i] > 0) {
      count++;
    }
  }

  return count;
}

// =================================================================================== //
// =================================================================================== //

// ============================== //
//          Acoustic NLL          //
// ============================== //
/*
 * Calculates the log-likelihood for cue-based (acoustic) SCR data
 * - Arguments:
 *    - pars:       Vector of (density, g0, sigma, lambda_c [, sigma_toa])
 *    - caps:       Matrix of captures; last column is animal ID
 *    - toa_ssq:    Matrix of times of arrival for each mask point
 *    - mask:       Matrix of mask points (coordinates)
 *    - maskDists:  Matrix of distances between each trap and mask point
 *    - nCalls:     Vector of counts; number of calls by each animal ID
 *    - use_toa:    TRUE/FALSE; whether TOA matrix is used
 */
// [[Rcpp::export]]
double scr_nll_acoustic(NumericVector pars,
                        NumericMatrix caps,
                        NumericMatrix traps,
                        NumericMatrix mask,
                        NumericMatrix maskDists,
                        NumericVector nCalls,
                        NumericMatrix toa,
                        NumericMatrix toa_ssq,
                        bool use_toa,
			bool is_g0_fixed,
			double g0_fixed,
			bool trace) {
  /*
   *  Storing/initialising (starting) parameter values.
   *  - Note that parameters are back-transformed
   *  - Also note that if use_toa = FALSE, there will be no sigma_toa to estimate
   */
  double D = exp(pars[0]);
  double g0;
  double sigma;
  double lambda_c;
  double sigma_toa;
  if (is_g0_fixed){
    g0 = g0_fixed;
    sigma = exp(pars[1]);
    lambda_c = exp(pars[2]);
  } else {
    g0 = R::plogis(pars[1], 0, 1, 1, 0);//exp(pars[1]) / (1 + exp(pars[1]));
    sigma = exp(pars[2]);
    lambda_c = exp(pars[3]);
  }
  // Number of animals
  int nAnimals = nCalls.size();
  // Number of calls detected in total.
  //int n = caps.nrow();
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
  * Constructing a detection probability matrix.
  * - Element (i, j) gives prob. of animal @ ith mask pt. being detected @ jth trap.
  * - Line that fills in maskProbs(i, j) is the Hazard Half-Normal function (HHN)
  */
  NumericMatrix maskProbs(maskDists.nrow(), maskDists.ncol());
  for(int i = 0; i < maskDists.nrow(); i++) {
    for(int j = 0; j < maskDists.ncol(); j++) {
      maskProbs(i, j) = g0 * exp(-pow(maskDists(i, j), 2.0) / (2 * pow(sigma, 2.0))) + DBL_MIN;
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
  *  Probability of detecting one specific call emitted from s
  *
  */
  NumericVector pDetected = 1 - pAvoid;

  /*
  * Probability of detecting at least one call on at least one microphone
  *  from an individual located at s
  */
  NumericVector pAnimal = 1 - exp(-lambda_c * pDetected);

  // ========================================= //
  // ========================================= //
  NumericVector fCapt(nAnimals);
  // Probability for an animal being at each mask point
  NumericVector logfCapt_givenNS(nMask);
  NumericVector logfn_givenS(nMask);
  // Row index for matrix subset
  double subRow = 0;
  // Creating subTOAs (regardless of whether use_toa = T/F)
  NumericMatrix subTOAs_raw;
  NumericMatrix subTOAs_ssq;
  // Number of traps that detected a call (i.e. non-zero elements in a given row of the sub-matrix)
  //int trapsHeard; - NOW FOUND VIA numNonZero()

  // Looping through all animals
  for (int i = 0; i < nAnimals; i++) {
    /*
     * Subsetting the capture matrix and TOA matrix
     * - Sub-matrices: all cols; first row of sub-mat --- nCalls - 1
     * Note: TOA matrix is subset IFF use_toa = TRUE
     */
    if(use_toa) {
      subTOAs_raw = toa(Range(subRow, subRow + nCalls[i] - 1), _);
      subTOAs_ssq = toa_ssq(Range(subRow, subRow + nCalls[i] - 1), _);
    }
    NumericMatrix subCaps = caps(Range(subRow, subRow + nCalls[i] - 1), _);
    subRow += nCalls[i];

    // Looping through each mask point
    for (int j = 0; j < nMask; j++) {
      logfn_givenS[j] = log(R::dpois(nCalls[i], lambda_c * pDetected[j], 0) + DBL_MIN);

      // Looping through the calls (each sub-matrix)
      for (int k = 0; k < nCalls[i]; k++) {
        logfCapt_givenNS[j] = -log(pDetected[j] + DBL_MIN);

        /*
         * Checking to see whether Time Of Arrival (TOA) has been specified
         * - TRUE: adds to log-likelihood
         * - FALSE: ignored
         */
        if (use_toa){
	  if (is_g0_fixed){
	    sigma_toa = exp(pars[3]);
	  } else {
	    sigma_toa = exp(pars[4]);
	  }
          logfCapt_givenNS[j] += (1 - numNonZero(subCaps(k, _))) * log(sigma_toa) - (subTOAs_ssq(k, j) / (2 * pow(sigma_toa, 2)));
        }

        // Looping through each trap
        for (int m = 0; m < traps.nrow(); m++) {
          logfCapt_givenNS[j] += R::dbinom(subCaps(k, m), 1, maskProbs(j, m), 1);
        }
      }
    }
    fCapt[i] = sum(exp(logfn_givenS + logfCapt_givenNS));
  }

  // ========================================= //
  // ========================================= //

  /*
  * Log-likelihood contribution from all capture histories
  * - Calculated by log of sum of individual likelihood contributions.
  */
  double logfCapt = sum(log(fCapt + DBL_MIN));

  // Calculating effective survey area (unused in likelihood).
  double esa = area * sum(pAnimal);

  // Log-likelihood contribution from number of animals detected.
  double logf_n = R::dpois(nAnimals, D * esa, 1);

  // Overall log-likelihood.
  double logLik = logf_n + logfCapt - nAnimals * log(sum(pAnimal));
  if (trace){
    std::cout << "D: " << D << ", g0: " << g0 << ", sigma: " << sigma << ", lambda_c: " << lambda_c << ", sigma_toa: " << sigma_toa << std::endl;
  }
  // Returning log-likelihood
  return -logLik;
}
// =================================================================================== //
// =================================================================================== //
