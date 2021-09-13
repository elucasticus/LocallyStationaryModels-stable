
#include <RcppEigen.h>
#include "traits.hpp"



// [[Rcpp::export]]
Rcpp::List rcppeigen_add22(const Eigen::MatrixXd & x, const Eigen::MatrixXd & y) {
  Eigen::MatrixXd mult = x*y*2;
  Eigen::MatrixXd add = x+y*2;
  return Rcpp::List::create(Rcpp::Named("Multiplication")=mult,
                            Rcpp::Named("Addition")=add);
}
