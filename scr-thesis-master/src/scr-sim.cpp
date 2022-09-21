#include <Rcpp.h>
using namespace Rcpp;

// ============================== //
//            eucdist             //
// ============================== //
/*
 * Calculating the Euclidean distance between a point and each trap.
 * Returns a vector of distances.
 */
/*
// [[Rcpp::export]]
NumericVector eucdist(NumericVector point,
                      NumericMatrix traplocations) {
  NumericVector dists(traplocations.nrow());
  for(int i = 0;  i < traplocations.nrow(); i++) {
    dists[i] = sqrt(pow(point[0] - traplocations(i, 0), 2.0) + pow(point[1] - traplocations(i, 1), 2.0));
  }
  return dists;
}
*/
// =================================================================================== //
// =================================================================================== //

// ============================== //
//            pointgen            //
// ============================== //
/*
 * Randomly generates points on the plotting/survey area.
 *  - Generates 50 random coordinates
 *  - Then combines them into an n x 2 matrix and returns
 */
// [[Rcpp::export]]
NumericMatrix pointgen(int n = 50,
                       NumericVector xlim = NumericVector::create(0, 100),
                       NumericVector ylim = NumericVector::create(0, 100)) {
  NumericVector xcoords = runif(n, xlim[0], xlim[1]);
  NumericVector ycoords = runif(n, ylim[0], ylim[1]);

  NumericMatrix coords(xcoords.length(), 2);
  coords(_, 0) = xcoords;
  coords(_, 1) = ycoords;
  return coords;
}


// =================================================================================== //
// =================================================================================== //

/*
 * Formatting the matrix returned by scr_sim,
 * so that it can be read by secr.fit
 */
/*
// [[Rcpp::export]]
NumericMatrix toCapthist_matrix(NumericMatrix captures) {
  NumericMatrix formatted(sum(captures), 2);
  for(int rowNum = 0; rowNum < captures.nrow(); rowNum++) {
    NumericVector trapNums = rep(whichCpp(captures(rowNum, _)), 3);//captures(rowNum, whichCpp(captures(rowNum, _))));
    formatted(rowNum, 0) = 1;
    formatted(rowNum, 1) = rowNum;
    formatted(rowNum, 2) = trapNums[1];
  }
  return formatted;
}
*/

// =================================================================================== //
// =================================================================================== //

// ============================== //
//          omega_acoustic        //
// ============================== //
/*
 *
 */


// =================================================================================== //
// =================================================================================== //

// ============================== //
//              rdistr            //
// ============================== //
/*
*
*/
/*
// [[Rcpp::export]]
NumericVector rdistr(String distribution,
              double lambda0,
              int n = 1,
              double prob = NULL,
              double mu = NULL,
              double size = NULL) {
  // Initialising result
  NumericVector result;

  // String matching
  if(distribution == "pois") {
    result = rpois(n, lambda0);
  } else if(distribution == "binom") {
    result = rbinom(n, 1, prob);
  } else if(distribution == "nbinom") {
    result = rnbinom(n, mu = lambda0, size = size);
  } else {
    return 0;
  }

  return result;

}
*/

// =================================================================================== //
// =================================================================================== //
