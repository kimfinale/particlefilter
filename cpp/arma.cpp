#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
// So the question is whether you can extract matrix using an underscore from
// an array

// //[[Rcpp::export]]
// double R_C_sum(arma::cube c) {
//   double d;
//   for (int i = 0; i < c.n_rows; i++) {
//     for (int j = 0; j < c.n_cols; j++) {
//       for (int k = 0; k < c.n_slices; k++) {
//         d += c(i, j, k);
//       }
//     }
//   }
//   return d;
// }
//
// //[[Rcpp::export]]
// NumericMatrix export_mat(arma::cube c){
//   NumericMatrix m(c.n_rows, c.n_cols);
//   for (int i = 0; i < c.n_rows; i++) {
//     for (int j = 0; j < c.n_cols; j++) {
//       m(i, j) = c(i, j, 0);
//     }
//   }
//   return m;
// }
//

//[[Rcpp::export]]
cube export_array(cube Q, int slice){
  cube Q2 = Q(span(), span(), span(slice));
  return Q2;
}

// //[[Rcpp::export]]
// NumericMatrix export_mat(cube Q, int slice){
//   NumericMatrix m = as<NumericMatrix>(Q(span(), span(), span(slice)));
//   return m;
// }
//[[Rcpp::export]]
NumericVector export_vec(cube Q, int col, int slice){
  NumericVector v = wrap(Q(span(), span(col), span(slice)));
  return v;
}

//[[Rcpp::export]]
IntegerVector export_vec_from_mat(NumericMatrix m, int col){
  NumericVector v = m(_, col);
  return as<IntegerVector>(v);
}

// [[Rcpp::export]]
IntegerVector to_int(NumericVector x) {
  IntegerVector iv = as<IntegerVector>(x);
  return iv;
}

