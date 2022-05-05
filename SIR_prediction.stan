functions {
  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {
      
      real beta = theta[1];
      real gamma = theta[2];
      real i0 = theta[3];

      int N = x_i[1];

      real init[3] = {N - i0, i0, 0};
      real S = y[1] + init[1];
      real I = y[2] + init[2];
      real R = y[3] + init[3];

      
      real dS_dt = -beta * I * S / N;
      real dI_dt = beta * I * S / N - gamma * I;
      real dR_dt =  gamma * I;
      
      return {dS_dt, dI_dt, dR_dt};
  }
}
data {
  int<lower=1> n_days;
  int<lower=1> n_cases;
  real t0;
  real ts[n_days];
  //int <lower=1> n_forecast_days;
  //real forecastts[n_forecast_days];
  int N;
  int cases[n_cases];
}
transformed data {
  real x_r[0];
  int x_i[1] = { N };
}
parameters {
  real<lower=0> gamma;
  real<lower=0> beta;
  real<lower=0> phi_inv;
  real<lower=0> i0; 
  real<lower=0, upper=1> p_reported; 
  //proportion of infected people reported
}
transformed parameters{
  real y[n_days, 3];
  real incidence[n_days - 1];
  real phi = 1. / phi_inv;
  //real forecast[n_forecast_days,3];
  ///real f_incidence[n_forecast_days-1];
  
  real theta[3];
  theta = {beta, gamma, i0};
  y = integrate_ode_rk45(sir, rep_array(0.0, 3), t0, ts, theta, x_r, x_i);
  
  
  //forecast= integrate_ode_rk45(sir, rep_array(0.0, 3), t0, forecastts, theta, x_r, x_i);
  
  for (i in 1:n_days-1){
    incidence[i] =  (y[i, 1] - y[i+1, 1])* p_reported; 
    //incidence = (S(t) - S(t-1))*p_reported
  }
  
  //for (j in 1:n_forecast_days-1){
    //f_incidence[j] = (forecast[j,1] - forecast[j+1,1])*p_reported;
  //}
}
model {
  //priors
  beta ~ normal(2, 1);
  gamma ~ normal(0.1, 0.3);
  i0 ~ normal(10*cases[1],2*cases[1]);
  phi_inv ~ exponential(5);
  p_reported ~ beta(1,2); 
  //similar to exponential, between 0 and 1, shouldn't close to 1
  
  //sampling distribution
  //col(matrix x, int n) - The n-th column of matrix x. Here the number of infected people 
  cases ~ neg_binomial_2(incidence[1:n_cases], phi);
}
generated quantities {
  real R0 = beta / gamma;
  real recovery_time = 1 / gamma;
  real pred_cases[n_days-1];
  //real forecast_cases[n_forecast_days-1];
 
  pred_cases = neg_binomial_2_rng(incidence, phi);
  //forecast_cases =  neg_binomial_2_rng(f_incidence, phi);
}

