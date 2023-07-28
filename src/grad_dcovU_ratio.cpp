#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double grad_dcovU_ratio(NumericMatrix a, NumericMatrix a2, NumericMatrix b, NumericVector row_sum_b, double sum_b,int n) {
  int nrow = a.nrow(), ncol = a.ncol();
  double total = 0;

  for (int i = 0; i < nrow; i++){
    for (int j = 0; j < ncol; j++){
      total += a(i,j)*a2(i,j)*(b(i,j)/n/(n-3) - 2*row_sum_b[i]/n/(n-2)/(n-3) + sum_b/n/(n-1)/(n-2)/(n-3));
    }
  }
  return total;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

// # /*** R
// # timesTwo(42)
// # */
