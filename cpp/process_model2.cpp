#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix process_model(NumericVector params,
                            NumericMatrix y,
                            double tbegin,
                            double tend,
                            double dt,
                            NumericVector beta) {

  for (int p = 0; p < npart; p++) {
    for (int s = 0; s < y.length(); s++) {
      lat_var(p, s) = latent_var(p, t - 1, s);
    }
  }

  NumericVector S = y(_, 0);
  NumericVector E1 = y(_, 1);
  NumericVector E2 = y(_, 2);
  NumericVector I = y(_, 3);
  NumericVector R = y(_, 4);
  NumericVector CE1 = y(_, 5);
  NumericVector CE2 = y(_, 6);
  NumericVector CI = y(_, 7);

  NumericVector pop = S + E1 + E2 + I + R;
  NumericVector infect_rate = beta / pop * dt;
  double incub_rate = params["sigma"] * 2 * dt; //# as two compartments
  double removal_rate = params["gamma"] * dt;

  // for (ii in seq((tbegin + dt), tend, dt)) {
  int reps = ((tend - tbegin) / dt);
  for (int i = 0; i < reps; i++) {

    NumericVector foi = I * infect_rate;
    NumericVector S_to_E1 = S * foi;
    NumericVector E1_to_E2 = E1 * incub_rate;
    NumericVector E2_to_I = E2 * incub_rate;
    NumericVector I_to_R = I * removal_rate;

//Process model for SEIR
    S = S - S_to_E1;
    E1 = E1 + S_to_E1 - E1_to_E2;
    E2 = E2 + E1_to_E2 - E2_to_I;
    I = I + E2_to_I - I_to_R;
    R = R + I_to_R;
    CE1 = CE1 + S_to_E1;
    CE2 = CE2 + E1_to_E2;
    CI = CI + E2_to_I;
  }

  y(_, 0) = S;
  y(_, 1) = E1;
  y(_, 2) = E2;
  y(_, 3) = I;
  y(_, 4) = R;
  y(_, 5) = CE1;
  y(_, 6) = CE2;
  y(_, 7) = CI;

  return y;
}

