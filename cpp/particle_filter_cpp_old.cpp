#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;
using namespace arma;

//
// // [[Rcpp::export]]
// NumericMatrix process_model_cpp(NumericVector params,
//                                 NumericMatrix y,
//                                 double tbegin,
//                                 double tend,
//                                 double dt,
//                                 NumericVector beta) {
//
//   NumericVector S = y(_, 0);
//   NumericVector E1 = y(_, 1);
//   NumericVector E2 = y(_, 2);
//   NumericVector I = y(_, 3);
//   NumericVector R = y(_, 4);
//   NumericVector CE1 = y(_, 5);
//   NumericVector CE2 = y(_, 6);
//   NumericVector CI = y(_, 7);
//
//   NumericVector pop = S + E1 + E2 + I + R;
//   NumericVector infect_rate = beta / pop * dt;
//   double incub_rate = params["sigma"] * 2 * dt; //# as two compartments
//   double removal_rate = params["gamma"] * dt;
//
//   // for (ii in seq((tbegin + dt), tend, dt)) {
//   int reps = ((tend - tbegin) / dt);
//   for (int i = 0; i < reps; i++) {
//
//     NumericVector foi = I * infect_rate;
//     NumericVector S_to_E1 = S * foi;
//     NumericVector E1_to_E2 = E1 * incub_rate;
//     NumericVector E2_to_I = E2 * incub_rate;
//     NumericVector I_to_R = I * removal_rate;
//
//     //Process model for SEIR
//     S = S - S_to_E1;
//     E1 = E1 + S_to_E1 - E1_to_E2;
//     E2 = E2 + E1_to_E2 - E2_to_I;
//     I = I + E2_to_I - I_to_R;
//     R = R + I_to_R;
//     CE1 = CE1 + S_to_E1;
//     CE2 = CE2 + E1_to_E2;
//     CI = CI + E2_to_I;
//   }
//
//   y(_, 0) = S;
//   y(_, 1) = E1;
//   y(_, 2) = E2;
//   y(_, 3) = I;
//   y(_, 4) = R;
//   y(_, 5) = CE1;
//   y(_, 6) = CE2;
//   y(_, 7) = CI;
//
//   return y;
// }
// [[Rcpp::export]]
NumericMatrix process_model_cpp(NumericVector params,
                                NumericMatrix y,
                                double tbegin,
                                double tend,
                                double dt,
                                NumericVector beta) {

  // You could just modify the input matrix y without creating another matrix out (commented out above)
  // They have different implications


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

  NumericMatrix out(y.nrow(), y.ncol());
  out(_, 0) = S;
  out(_, 1) = E1;
  out(_, 2) = E2;
  out(_, 3) = I;
  out(_, 4) = R;
  out(_, 5) = CE1;
  out(_, 6) = CE2;
  out(_, 7) = CI;

  return out;
}


// [[Rcpp::export]]
NumericVector assign_weights_cpp(NumericMatrix state_now,
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
    inc_data[0] = data(t - 1, 2); // daily infection
  }
  else if (data_type[0] == std_data_types[1]) {
    inc_model[Rcpp::Range(0, inc_model.size()-1)] =  state_now(_, 7) - state_before(_, 7); // CI
    inc_data[0] = data(t - 1, 3); // daily symptom onset
  }
  else if (data_type[0] == std_data_types[2]) {
    inc_model[Rcpp::Range(0, inc_model.size()-1)] = state_now(_, 4) - state_before(_, 4); //R
    inc_data[0] = data(t - 1, 4);// daily confirmation
  }

  IntegerVector inc_data_int = as<IntegerVector>(inc_data);
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
    if (!R_IsNA(inc_data_int[0])) {
      log_lik(i) = R::dpois(inc_data_int[0], inc_model[i], true); // vector input
    }
    else {
      log_lik(i) = R_NegInf;
    }
  }
  return exp(log_lik); //# convert to normal probability
}


// [[Rcpp::export]]
List particle_filter_cpp(NumericVector params,
                     NumericVector y,
                     NumericMatrix data,
                     CharacterVector data_type,
                     int npart,
                     int tend,
                     double dt) {


  int nstate = y.size();
  CharacterVector type = data_type;
  arma::cube latent_var(npart, tend, nstate); //RcppArmadillo 3D array
  latent_var.zeros(); // initialize the array with ones

  for (int p = 0; p < npart; p++) {
    for (int s = 0; s < nstate; s++) {
      latent_var(p, 0, s) = y(s);
      // Rprintf("latent_var(%i,0,%i) = %.1f\n", p, s, latent_var(p, 0, s));
    }
  }
  double beta = params["R0"] * params["gamma"];
  NumericMatrix beta_vol(tend, npart);
  for (int t = 0; t < tend; t++) {
    for (int p = 0; p < npart; p++) {
      beta_vol(t, p) = R::rnorm(0, params["betavol"]);
      // Rprintf("beta_vol(%i,%i) = %.3f\n", t, p, beta_vol(t, p));
    }
  }

  for(int p = 0; p < npart; p++) {
    beta_vol(0, p) = beta * exp(beta_vol(0, p)); // initial beta
  }
  // beta_vol(0, _) = rep(beta, npart) * exp(beta_vol(0, _)); // initial beta

  // for (int p = 0; p < npart; p++) {
  //   Rprintf("beta_vol(0,%i) = %.3f\n", p, beta_vol(0, p));
  // }
  NumericMatrix wt(npart, tend);
  wt(_, 0) = rep(1.0/npart, npart); // initial weights
  // for (int p = 0; p < npart; p++) {
  //   Rprintf("wt(%i, 0) = %.3f\n", p, wt(p, 0));
  // }
  NumericMatrix W(npart, tend); // normalized weights
  IntegerMatrix A(npart, tend); // Resample according to the normalized weight


  //begin particle loop
  for (int t = 1; t < tend; t++) {  //DEBUG  t=2
    beta_vol(t, _) = beta_vol(t - 1, _) * exp(beta_vol(t, _));

    NumericMatrix lat_var_t_neg_one(npart, nstate); // temporary state variables
    for (int p = 0; p < npart; p++) {
       for (int s = 0; s < nstate; s++) {
         lat_var_t_neg_one(p, s) = latent_var(p, t - 1, s);
       }
    }

    //run process model
    NumericVector beta_t = beta_vol(t, _);
    NumericMatrix lat_var_t =
      process_model_cpp(params,
                        lat_var_t_neg_one,
                        t - 1,
                        t,
                        dt,
                        beta_t);

    // for (int p = 0; p < 3; p++) {
    //   for (int s = 0; s < nstate; s++) {
    //     Rprintf("t = %i, lat_var_t_neg_one(%i, %i) = %.6f\n", t, p, s, lat_var_t_neg_one(p, s));
    //     Rprintf("t = %i, lat_var_t(%i, %i) = %.6f\n", t, p, s, lat_var_t(p, s));
    //   }
    // }
    //updatea master variables latent_vat
    for (int p = 0; p < npart; p++) {
      for (int s = 0; s < nstate; s++) {
        latent_var(p, t, s) = lat_var_t(p, s);
      }
    }

  // calculate weights (likelihood)
    wt(_, t) =
      assign_weights_cpp(lat_var_t,
                         lat_var_t_neg_one,
                         t,
                         data,
                         type);

// normalize particle weights
    W(_, t) = wt(_, t) / sum(wt(_, t));
// resample particles by sampling parent particles according to weights
    NumericVector Wvec = W(_, t);
    IntegerVector npart_seq = seq(0, npart - 1);
    IntegerVector ids_resampled = sample(npart_seq, npart, true, Wvec);

    // for (int p = 0; p < 10; p++) {
    //   Rprintf("t = %i, ids_resampled(%i) = %i\n", t, p, ids_resampled(p));
    // }

    A(_, t) = ids_resampled;
    // Resample particles for corresponding variables
    for (int p = 0; p < npart; p++) {
      beta_vol(t, p) = beta_vol(t, ids_resampled(p)); //#- needed for random walk on beta
      for (int s = 0; s < nstate; s++) {
        latent_var(p, t, s) = latent_var(ids_resampled(p), t, s);
      }
    }
  } //end particle loop

// Marginal likelihoods:
  NumericVector lik_values(tend, 0.0);
  for (int t = 0; t < tend; t++) {
    lik_values(t) = log(sum(wt(_, t))); // log-likelihoods
  }
  double loglik = - tend * log(npart) + sum(lik_values);// # averaged log likelihoods log(L/(npart^tend))

  // =====================BACKWARD SAMPLING ====================================
  // Latent variables sampled from the filtered distribution (backward sampling)
  NumericMatrix traj(tend, nstate + 1);
  CharacterVector nm = {"S", "E1", "E2", "I", "R", "CE1", "CE2", "CI", "beta"};
  colnames(traj) = nm;

  IntegerVector loc(tend, NA_INTEGER);
  NumericVector Wvec_tend = W(_, tend - 1);
  IntegerVector sample_at_tend = sample(npart, 1, true, Wvec_tend);
  loc(tend - 1) = sample_at_tend(0);
  // particle for the last time step
  traj(tend - 1, 8) = beta_vol(tend, loc(tend - 1)); //beta
  for (int s = 0; s < 8; s++) {
    traj(tend - 1, s) = latent_var(loc[tend - 1], tend - 1, s);
  }
  // arma::cube latent_var(npart, tend, nstate);
  // update backward
  for (int t = (tend - 1); t > 0; t--) {
    loc(t - 1) = A(loc(t), t);
    traj(t - 1, 8) = beta_vol(t-1, loc[t - 1]);
    for (int s = 0; s < 8; s++) {
      traj(t - 1, s) = latent_var(loc[t - 1], t - 1, s);
    }
  }


  List list = List::create(Named("trace") = traj,
                        Named("lik_marginal") = lik_values,
                        Named("lik_overall_average") = loglik,
                        Named("latent_var_filtered") = latent_var,
                        Named("beta_filtered") = beta_vol);

//# DEBUG  plot(Rep_traj[,1]-C_traj[,1])
  return list;
};

