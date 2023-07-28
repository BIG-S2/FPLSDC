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
NumericVector grad_dcov_U(NumericMatrix sgn_a, NumericMatrix b, NumericMatrix  a_u,NumericVector row_sum_b, double sum_b,int n) {
  int nrow = a_u.nrow(), ncol = a_u.ncol();
  NumericVector out(ncol);
  NumericMatrix out_u(nrow,nrow);
  NumericVector total(ncol);

  for (int u = 0; u < ncol; u++){

    for (int i_u = 0; i_u < nrow; i_u++) {
      for (int j_u = 0; j_u < nrow; j_u++) {
        out_u(i_u,j_u)= a_u(i_u,u) - a_u(j_u,u);
      }
    }

    for (int i = 0; i < nrow; i++){
      for (int j = 0; j < nrow; j++){
        total[u] += sgn_a(i,j)*out_u(i,j)*(b(i,j)/n/(n-3) - 2*row_sum_b[i]/n/(n-2)/(n-3) + sum_b/n/(n-1)/(n-2)/(n-3));
      }
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
