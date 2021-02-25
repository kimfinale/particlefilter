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

// [[Rcpp::export]]
CharacterVector add_x(CharacterVector x) {
  CharacterVector y = x;
  return y;
}

// [[Rcpp::export]]
CharacterVector all_x(CharacterVector x) {
  CharacterVector y = x;
  CharacterVector test = {"x"};
  if (x[0] == test[0]){
  }
  else{
    y = test;
  };
  return y;
}
// [[Rcpp::export]]
IntegerVector int_sample(int n) {
  IntegerVector sam = Rcpp::sample(n, n, true);
  return sam;
}

// [[Rcpp::export]]
IntegerVector int_wt_sample(int n, NumericVector wt) {
  IntegerVector sam = Rcpp::sample(n, n, true, wt);
  return sam;
}

//[[Rcpp::export]]
void Rprintf_test(int n){
  for (int i = n; i > 1; i--) {
    Rprintf("i = %i\n", i);
  }
}

//[[Rcpp::export]]
IntegerVector sample_test(int n, NumericVector wt){
  IntegerVector npart_seq = seq(0, n - 1);
  IntegerVector ids_resampled = sample(npart_seq, n, true, wt);
  return ids_resampled;
}

