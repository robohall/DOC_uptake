
data {
  int<lower=0> N;
  int<lower=0> N_groups;
  vector[N] DIC_flux;
  vector[N] R_g;
   vector <lower=0> [N]   MinFrom0;
  int<lower=0> group_id [N] ;
  vector[N] removed_c;
  vector[N] q;

  
}



parameters {
  vector <lower=0,upper=1> [N_groups]   fi;
  real <lower=0, upper=1> phi;
  real <lower=0.1> lambda;
  
  vector <lower=0> [N_groups]   kc;
  real <lower=0> kc_mean;
  real <lower=0> kc_sd;
  
  real<lower=0> sigma;
}

transformed parameters{
  vector[N] pred_imm;
  vector[N] pred_delayed;
  real <lower=0> alpha;
  real <lower=0> beta;



for (i in 1:N){
   pred_imm[i] = fi[group_id[i]] * removed_c[i]* q[i];
   
  pred_delayed[i]= (1-fi[group_id[i]])*R_g[i]*kc[group_id[i]]/1440 * exp(-kc[group_id[i]]*MinFrom0[i]/1440);
}
  
alpha = lambda * phi;
beta = lambda * (1 - phi);
  
}

model {
  DIC_flux ~ normal(pred_imm + pred_delayed, sigma);
  
  fi~beta (alpha, beta);
  kc~normal (kc_mean, kc_sd);
  
  phi~beta (1,1);
  kc_mean~ gamma(1,0.5) ;
  
  lambda~gamma(1,5);
  kc_sd~normal(0,1);
  
  sigma~ normal(0,1);
}

