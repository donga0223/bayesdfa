functions {

  // this function subsets a matrix by dropping the row/column labeled 'drop'. P represents dimensions

  matrix block_diag_matrix(matrix A, int H, int P) {
    //matrix[H,H] O;
    matrix[H*P,H*P] rtn;
    //O = rep_matrix(0, H,H);
    
    for(i in 1:P){
      rtn[((i-1)*H+1):(i*H),] = append_col(append_col(rep_matrix(0,H,H*(i-1)),A), rep_matrix(0,H,H*(P-i)));
    }  
    return rtn;
  }

  matrix block_diag_matrix_1(matrix A, int rep) {
    int a[2] = dims(A);
    //matrix[a[1], a[2]] O;
    matrix[a[1]*rep, a[2]*rep] rtn;
    //O = rep_matrix(0, a[1],a[2]);

    for(i in 1:rep){
      rtn[((i-1)*a[1]+1):(i*a[1]),] = append_col(append_col(rep_matrix(0,a[1],a[2]*(i-1)),A), rep_matrix(0,a[1],a[2]*(rep-i)));
    }
    return rtn;
  } 

  vector repeat_vector(vector input, int K) {
    int N = rows(input);
    vector[N*K] repvec; // stack N-vector K times

    for (k in 1:K) {
      for (i in 1:N) {
        repvec[i+(k-1)*N] = input[i]; // assign i-th value of input to i+(k-1)*N -th value of repvec
      }
    }
    return repvec;
  }
}


data{
  int<lower=0> N; // number of data points
  int<lower=0> P; // number of time series of data
  int<lower=0> K; // number of trends
  int<lower=0> H; // number of horizons  *****horizon dimension include*****
  int<lower=0> nZ; // number of unique z elements
  int<lower=0> n_pos; // number of non-missing observations
  real z_bound[2];
  int<lower=0> var_horizonIndx[H];
  int<lower=0> row_indx_pos[n_pos]; // row indices of non-missing obs
  int<lower=0> col_indx_pos[n_pos]; // col indices of non-missing obs
  int<lower=0> hor_indx_pos[n_pos]; // hor indices of non-missing obs   *****horizon dimension include*****
  real y[n_pos]; // vectorized matrix of observations
  int<lower=0> row_indx[nZ];
  int<lower=0> col_indx[nZ];
  int<lower=0> nZero;
  int<lower=0> row_indx_z[nZero];
  int<lower=0> col_indx_z[nZero];
  int<lower=0> num_obs_covar; // number of unique observation covariates, dimension of matrix
  int obs_covar_index[num_obs_covar,4];// indexed by time, trend, covariate   , covariate value. +1 because of indexing issues
  int<lower=0> est_rw; 
  int<lower=0> n_obs_covar; // number of unique covariates included
  real obs_covar_value[num_obs_covar];
  int estimate_nu; // Estimate degrees of freedom?
  int<lower=0> est_sigma_process; // optional, 0 == not estimate sigma_pro (default), 1 == estimate
  int<lower=0> n_sigma_process; // single value, or equal number of tre
  int use_normal; // flag, for large values of nu > 100, use normal instead
  int est_cor; // whether to estimate correlation in obs error (=1) or not (=0)
  real<lower=1> nu_fixed; // df on student-t
  int<lower=0, upper=1> use_expansion_prior;
  int<lower=0> proportional_model;
}

transformed data{
  int n_pcor; // dimension for cov matrix
  matrix[H, 1] Jone; 
  matrix[P*H, P] Aone;
  real lower_bound_z;
  int n_loglik; // dimension for loglik calculation
  vector[K*proportional_model] alpha_vec;

  lower_bound_z = -100;
  if(use_expansion_prior==1) lower_bound_z = 0;
  n_pcor = P;
  
  if(H == 1){
    Aone = identity_matrix(P);
  }else{
    for(i in 1:H){
      Jone[i,1]=1;
    }
    Aone = block_diag_matrix_1(Jone, P);
  }
  

  

  // for compositional model
  if(proportional_model==1) {
    for(k in 1:K) alpha_vec[k] = 1;
  }

  
}

parameters{
  vector<lower=z_bound[1],upper=z_bound[2]>[nZ*(1-proportional_model)] z; // estimated loadings in vec form
  vector<lower=lower_bound_z>[K*(1-proportional_model)] zpos; // constrained positive values
  simplex[K] p_z[P*proportional_model]; // alternative for proportional Z
  matrix[K * est_rw,(N-1) * est_rw] devs; // random deviations of trends
  vector[K] x0; // initial state
  vector<lower=0>[K*(1-proportional_model)*use_expansion_prior] psi; // expansion parameters
  vector[n_obs_covar] b_obs; // coefficients on observation model
  real<lower=-1,upper=1> phi[K]; // AR(1) coefficients specific to each trend
  real<lower=0> sigma_h[H];
  real<lower=2> nu[estimate_nu]; // df on student-t
  real<lower=0> sigma_process[est_sigma_process * n_sigma_process]; // process variances, potentially unique
}



transformed parameters{
  matrix[P*H,N] pred; //vector[P] pred[N];
  matrix[P,K] Z;
  matrix[P*H,K] AZ;
  matrix[P*H,N] yall;  //*****horizon dimension include*****
  vector<lower=0>[P*H] sigma_vec;
  vector<lower=0>[H] sigma_horizon;
  matrix[K,N] x; //vector[N] x[P]; // random walk-trends
  vector[K] indicator; // indicates whether diagonal is neg or pos
  vector<lower=0>[K*use_expansion_prior] psi_root; // derived sqrt(expansion parameter psi)
  vector[K] phi_vec; // for AR(1) part
  vector[K] sigma_pro;
  vector[K] theta_vec; // for MA(1) part



  for(k in 1:K) {
    sigma_pro[k] = 1; // default constraint of all DFAs
    if(est_sigma_process==1) {
      if(n_sigma_process==1) {
        sigma_pro[k] = sigma_process[1];
      } else {
        sigma_pro[k] = sigma_process[k];
      }
    }
  }

  for(k in 1:K) {phi_vec[k] = phi[k];}

  for(k in 1:K) {theta_vec[k] = 0;}

  
  if(H == 1){
    sigma_horizon[1] = sigma_h[var_horizonIndx[1]];
    sigma_vec = rep_vector(sigma_horizon[1], P);
  }else{
    for(h in 1:H) {sigma_horizon[h] = sigma_h[var_horizonIndx[h]];}
    sigma_vec = repeat_vector(sigma_horizon, P);
  }
  


  if(proportional_model == 0) {

    for(i in 1:nZ) {
      Z[row_indx[i],col_indx[i]] = z[i]; // convert z to from vec to matrix
    }
    
    // fill in zero elements in upper diagonal
    if(nZero > 2) {
      for(i in 1:(nZero-2)) {
        Z[row_indx_z[i],col_indx_z[i]] = 0;
      }
    }
    for(k in 1:K) {
      Z[k,k] = zpos[k];// add constraint for Z diagonal
    }
    // this block is for the expansion prior
    if(use_expansion_prior==1) {
      for(k in 1:K) {
        if(zpos[k] < 0) {
          indicator[k] = -1;
        } else {
          indicator[k] = 1;
        }
        // psi_root here should really be named inv_psi_root, because it's the inv
        psi_root[k] = sqrt(psi[k]);
        for(p in 1:P) {
          // see Ghosh & Dunson 2009 eq 3
          Z[p,k] = Z[p,k] * indicator[k] * (1/psi_root[k]);
        }
      }
    }

    // initial state for each trend
    if(est_rw == 1) {
      for(k in 1:K) {
        x[k,1] = x0[k];
        // trend is modeled as random walk, with optional
        // AR(1) component = phi, and optional MA(1) component
        // theta. Theta is included in the model block below.
        for(t in 2:N) {
          x[k,t] = phi_vec[k]*x[k,t-1] + devs[k,t-1];
        }
      }
    }

    // this block also for the expansion prior, used to convert trends
    if(use_expansion_prior==1) {
      for(k in 1:K) {
        //x[k,1:N] = x[k,1:N] * indicator[k] * psi_root[k];
        // see Ghosh and Dunson 2009 eq 3. psi_root[k] here is really psi^(-1/2)
        x[k] = x[k] * indicator[k] * psi_root[k];
      }
    }
  }

  if(proportional_model == 1) {
    // initial state for each trend
    if(est_rw == 1) {
      for(k in 1:K) {
        x[k,1] = x0[k];
        // trend is modeled as random walk, with optional
        // AR(1) component = phi, and optional MA(1) component
        // theta. Theta is included in the model block below.
        for(t in 2:N) {
          x[k,t] = phi_vec[k]*x[k,t-1] + devs[k,t-1];
        }
      }
    }

    // proportional model
    for(p in 1:P) {
      //for(k in 1:K) {
      //  Z[p,k] = p_z[p,k]; // compositions sum to 1 for a time series
      //}
      Z[p] = to_row_vector(p_z[p]);
    }
  }

  
  
  

  AZ = Aone * Z;
  pred = AZ * x;
  

  if(num_obs_covar > 0) {
    for(i in 1:num_obs_covar) {
      pred[((obs_covar_index[i,3]-1)*H+obs_covar_index[i,2]),obs_covar_index[i,1]] += b_obs[obs_covar_index[i,4]] * obs_covar_value[i];
    }
  }

}

model{
  // initial state for each trend
  x0 ~ normal(0, 1); // initial state estimate at t=1

  if(use_expansion_prior==1) {
    psi ~ gamma(2, 1); // expansion parameter for par-expanded priors
  }

  // prior on AR(1) component if included
  phi ~ normal(0,1); // K elements
  //phi ~ uniform(-1,1);
  
  if(est_sigma_process) {
    sigma_process ~ normal(0,1);
  }

  // observation variance, which depend on family
  //sigma_l ~ student_t(3, 0, 1);
  sigma_h ~ student_t(3, 0, 1);


  if(est_rw == 1) {
    for(k in 1:K) {
      if(use_normal == 0) {
        for(t in 1:1) {
          if (estimate_nu == 1) {
            devs[k,t] ~ student_t(nu[1], 0, sigma_pro[k]); // random walk
          } else {
            devs[k,t] ~ student_t(nu_fixed, 0, sigma_pro[k]); // random walk
          }
        }
        for(t in 2:(N-1)) {
          // if MA is not included, theta_vec = 0
          if (estimate_nu == 1) {
            devs[k,t] ~ student_t(nu[1], theta_vec[k]*devs[k,t-1], sigma_pro[k]); // random walk
          } else {
            devs[k,t] ~ student_t(nu_fixed, theta_vec[k]*devs[k,t-1], sigma_pro[k]); // random walk
          }
        }
      } else {
        devs[k,1] ~ normal(0, 1);
        for(t in 2:(N-1)) {
          // if MA is not included, theta_vec = 0
          devs[k,t] ~ normal(theta_vec[k]*devs[k,t-1], sigma_pro[k]);
        }
      }
    }
  }

  if(proportional_model == 0) {
    // prior on loadings
    z ~ std_normal(); // off-diagonal
    zpos ~ std_normal();// diagonal
  } else {
    for(p in 1:P) {
      p_z[p] ~ dirichlet(alpha_vec);
    }
  }


  

  // likelihood for independent
      
  for(i in 1:P*H) target += normal_lpdf(yall[i] | pred[i], sigma_vec[i]);
}



generated quantities{
  vector[n_loglik] log_lik;
  
  int<lower=0> j;
  j = 0;

  
  for(n in 1:N) {
    for(p in 1:P*H) {
      j = j + 1;
      log_lik[j] = normal_lpdf(yall[p,n] | pred[p,n], sigma_vec[p]);
    }
  }
   
}

