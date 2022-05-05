functions {
int real_to_int(real x, int min_val, int max_val){
    // This assumes that min_val >= 0 is the minimum integer in range, 
    //  max_val > min_val,
    // and that x has already been rounded. 
    //  It should find the integer equivalent to x.
    int range = (max_val - min_val+1)/2; // We add 1 to make sure that truncation doesn't exclude a number
    int mid_pt = min_val + range;
    int out;
    while(range > 0) {
      if(x == mid_pt){
        out = mid_pt;
        range = 0;
      } else {
        // figure out if range == 1
        range =  (range+1)/2; 
        mid_pt = x > mid_pt ? mid_pt + range: mid_pt - range; 
        }
    }
    return out;
  }  
  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {
      real lambda = theta[1];
      real sigma = theta[2];
      real mu = theta[3];
      real alpha = theta[4];
      real E0 = theta[5];
      real I0 = theta[6];
      
      real N = x_i[1];
      real n_days = x_r[1];
      
      real init[4] = {N-I0-E0, E0, I0, 0};
      real S = y[1] + init[1];
      real E = y[2] + init[2];
      real I = y[3] + init[3];
      real R = y[4] + init[4];
      
      int tindex = real_to_int(ceil(t),0,size(x_r)-1);
      real dS_dt = -lambda*S*I/N - alpha*x_r[tindex+1];
      real dE_dt = lambda*S*I/N - sigma*E ;
      real dI_dt = sigma*E - mu*I;
      real dR_dt = mu*I + alpha*x_r[tindex+1];
      
      return {dS_dt, dE_dt, dI_dt, dR_dt};

               
             }
}

data {
  int<lower=1> n_days;
  int <lower=1> n_cases;
  real t0;
  real ts[n_days];
  real x_r[n_days +1];
  int cases[n_cases];
  int N;
}



transformed data {
  int x_i[1] = { N };
}

parameters {
  real<lower=0> lambda;
  real<lower=0> sigma;
  real<lower=0> mu;
  real<lower=0> alpha;

  real<lower=0> phi_inv;
  real<lower=0> I0; 
  real<lower=0> E0;
  real<lower=0, upper=1> p_reported; 
  //proportion of infected people reported
}

transformed parameters{
  real y[n_days,4];
  real incidence[n_days - 1];
  real phi = 1. / phi_inv;
  //real forecast[n_forecast_days,4];
 // real f_incidence[n_forecast_days-1];

  real theta[6];
  theta = {lambda, sigma,mu, E0, I0, alpha};
  y = integrate_ode_rk45(sir, rep_array(0.0, 4), t0, ts, theta, x_r, x_i);
  
  //forecast = integrate_ode_rk45(sir, rep_array(0.0, 4), t0, forecastts, theta, x_r, x_i);
  
  for (i in 1:n_days-1){
    incidence[i] =  (y[i, 1] - y[i+1, 1] + y[i,2] - y[i+1,2])* p_reported; 
    //incidence = (S(t)-S(t+1)+E(t)-E(t+1))*p_reported
  }
  
  //for (j in 1:n_forecast_days-1){
   // f_incidence[j] = (forecast[j,1] - forecast[j+1,1] + forecast[j,2] - forecast[j+1,2])* p_reported;
  //}
}


model {
  //priors
  lambda ~ normal(0.21, 0.5);
  sigma ~ normal(0.2, 0.5);
  mu ~normal(0.07,0.3);
  alpha ~ normal(0.95,0.01);
  I0 ~ normal(6*cases[1],2*cases[1]);
  E0 ~ normal(4*cases[1],2*cases[1]); 
  phi_inv ~ exponential(5);
  p_reported ~ beta(1,2); 
 
  cases ~ neg_binomial_2(incidence[1:n_cases], phi);
  
}
generated quantities {
  real R0 = lambda / mu;
  real recovery_time = 1 / mu;
  real incubation_time = 1 / sigma;
  real pred_cases[n_days-1];
  //real forecast_cases[n_forecast_days-1];
  
  pred_cases = neg_binomial_2_rng(incidence, phi);
  //forecast_cases =  neg_binomial_2_rng(f_incidence, phi);
}
