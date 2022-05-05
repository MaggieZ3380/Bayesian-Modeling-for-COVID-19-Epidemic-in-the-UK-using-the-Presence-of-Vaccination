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
  
      real lambda1 = theta[1];
      real lambda2 = theta[2];
      real sigma = theta[3];
      real mu = theta[4];
      real alpha = theta[5];
      real E0 = theta[6];
      real I0 = theta[7];
      real delta = theta[8];
      //real epsilon = theta[8];
      
      
      
      
      real N = x_i[1];
      real n_days = x_r[1];
      
      real lambdaP; 
      
      real init[5] = {(N*0.9)-I0-E0, (N*0.1)-I0-E0, E0, I0, 0};
      real P = y[1] + init[1];
      real S = y[2] + init[2];
      real E = y[3] + init[3];
      real I = y[4] + init[4];
      real R = y[5] + init[5];
      
      int tindex = real_to_int(ceil(t),0,size(x_r)-1);
      if (tindex < 58 + delta){
        lambdaP = 0; 
      }else{
        lambdaP=lambda2;
      }  
      
      real dP_dt = -lambdaP*S*I/N - 0.9*alpha*x_r[tindex+1];
      real dS_dt = lambdaP*S*I/N -lambda1*S*I/(N) - 0.1*alpha*x_r[tindex+1];
      real dE_dt = lambda1*S*I/N - sigma*E ;
      real dI_dt = sigma*E - mu*I;
      real dR_dt = mu*I + alpha*x_r[tindex+1];
        
      return {dP_dt, dS_dt, dE_dt, dI_dt, dR_dt};

            
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
  real<lower=0> lambda1;
  real<lower=0> lambda2;
  real<lower=0> sigma;
  real<lower=0> mu;
  real<lower=0> alpha;
  
  real<lower=0> delta;
  //real epsilon;

  real<lower=0> phi_inv;
  real<lower=0> I0; 
  real<lower=0> E0;
  real<lower=0, upper=1> p_reported; 


}

transformed parameters{
  real y[n_days,5];
  real incidence[n_days - 1];
  real phi = 1. / phi_inv;

  real theta[8];
  theta = {lambda1, lambda2, sigma, mu, alpha, E0, I0,  delta};
  y = integrate_ode_rk45(sir, rep_array(0.0, 5), t0, ts, theta, x_r, x_i);
  
  
  for (i in 1:n_days-1){
    incidence[i] =  (y[i,1] + y[i,2] - y[i+1,1] - y[i+1, 2] + y[i,3] - y[i+1,3])* p_reported; 
    //incidence = (P(t) + S(t)-P(t+1) -S(t+1)+E(t)-E(t+1))*p_reported
  }
  
 
}


model {
  //priors
  delta ~ normal(11, 3);
  lambda1 ~ normal(0.21, 0.5);
  //epsilon ~ normal(0.01, 0.5);
  lambda2 ~ normal(0.1,0.5);
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
  //real R0 = lambda / mu;
  real recovery_time = 1 / mu;
  real incubation_time = 1 / sigma;
  real pred_cases[n_days-1];
  //real forecast_cases[n_forecast_days-1];
  
  pred_cases = neg_binomial_2_rng(incidence, phi);
  //forecast_cases =  neg_binomial_2_rng(f_incidence, phi);
}

