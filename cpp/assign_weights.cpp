#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;



// [[Rcpp::export]]
NumericVector assign_weights(NumericMatrix state_now,
                             NumericMatrix state_before,
                            int t,
                            NumericMatrix data,
                            CharacterVector data_type) {

  NumericVector inc_model(state_now.nrow()); // incidence predicted by the model
  NumericVector inc_data(1); //incidence from data
  //type <- match.arg(data_type)
  // c("S", "E1", "E2", "I", "R", "CE1", "CE2", "CI")
  CharacterVector std_data_types = {"infection", "symptom onset", "confirmation"};
  if (data_type[0] == std_data_types[0]) {
    inc_model[Rcpp::Range(0, inc_model.size()-1)] = state_now(_, 5) - state_before(_, 5); // CE1
    inc_data[0] = data(t, 3); // daily infection
  }
  else if (data_type[0] == std_data_types[1]) {
    inc_model[Rcpp::Range(0, inc_model.size()-1)] =  state_now(_, 7) - state_before(_, 7); // CI
    inc_data[0] = data(t, 4); // daily onset
  }
  else if (data_type[0] == std_data_types[2]) {
    inc_model[Rcpp::Range(0, inc_model.size()-1)] = state_now(_, 4) - state_before(_, 4); //R
    inc_data[0] = data(t, 5);// daily confirmation
  }

  NumericVector expected(inc_model.size());
  NumericVector log_lik(inc_model.size());
  //allow 0 and positive numbers, same as pmax(o, inc_model)
  for (int i = 0; i < inc_model.size(); i++) {
    if (inc_model[i] < 0){
      inc_model[i] = 0;
    }
  }
  //IntegerVector inc_data_int = as<IntegerVector>(inc_data);
  for (int i = 0; i < inc_model.size(); i++) {
    if (!R_IsNA(inc_data[0])) {
      log_lik[i] = R::dpois(inc_data[0], inc_model[i], true); // vector input
    }
    else {
    log_lik[i] = R_NegInf;
    }
  }

  return exp(log_lik); //# convert to normal probability
}

