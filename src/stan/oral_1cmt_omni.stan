functions {
  #include "utils/utils.stan"
  #include "utils/mpi_utils.stan"
  #include "utils/model_lib.stan"
  #include "utils/model_lib_2.stan"
  
  //linear allometric scaling for log of clearance or volume given
  //the log of reference value75 (lref), power (alpha) and size
  //of the patient (wi). If alpha is not learned from the data,
  //principled default value for it is 3/4. This function doesn't
  //take into account the personal effect of the patient.
  vector log_lin_allo(real lref, real alpha, vector wi){
    return(lref + alpha *( log(wi) - log(75.0)));
  }

  // log-version for 1-cmt model
  vector pk_mr(vector theta, vector eta, real[] x_r, int[] x_i) {
    // x_r[1] = tad; x_r[2] = log(amt); x_r[3] = tau
    // x_i[1] = cmt; x_i[2] = addl
    // relative bioavailability
    real lfrel = (theta[4] < -99) ? 0.0: log(1E3) + theta[4] * (x_r[2] - log(5.0));
    //real log(1E3) + theta[4] * (x_r[2] - log(5.0)); #x_r[2];
    return lfrel + lpk_1cmt_oral_mr(theta[1:3], rep_vector(0.0, 0), x_r, x_i);
  }

#include "utils/model_lib_2_post.stan"

  vector subject_obs_mr(vector theta, vector xi, real[] x_r, int[] x_i) {
    int num_obs = x_i[1];
    int offset_obs_int  = x_i[2];
    int offset_obs_real = x_i[3];
    int ioffset = 3;
    int obs_cmt[num_obs] = extract_1d_int(x_i, ioffset, num_obs, 1);
    real lV = xi[2]; // theta[4]; TODO
    real lka = theta[1];
    real lfrel = theta[2];
    // log(ke) = log(CL) - log(V)
    real lke = xi[1] - lV;
    // ka, ke, log(ka/(ka-ke)), lfrel
    vector[4] pk_theta = [exp(theta[1]), exp(lke), lka - log_diff_exp_abs(lka, lke), lfrel ]';
    //print(theta);
    //print(xi);
    //real tmp = sum(lpk_subject_mr(pk_theta, rep_vector(0, 0), x_r[offset_obs_real+1 : size(x_r)], x_i[3+offset_obs_int+1 : size(x_i)]) - lV);
    //if( tmp != tmp ) { #catches nan
    //  print(pk_theta);
    //  print(x_r);
    //  print(x_i);
    //}
    return lpk_subject_mr(pk_theta, rep_vector(0, 0), x_r[offset_obs_real+1 : size(x_r)], x_i[3+offset_obs_int+1 : size(x_i)]) - lV;
  }

  vector subject_lpdf_mr(vector theta, vector xi, real[] x_r, int[] x_i) {
    int num_obs = x_i[1];
    real sigma_p = theta[5];
    real sigma_y = theta[6];
    vector[num_obs] lmu = subject_obs_mr(theta, xi, x_r, x_i);
    vector[num_obs] sigma = sqrt( square(sigma_y) + square(sigma_p * exp(-lmu)));
    return [ normal_lpdf(extract_vector(x_r, 0, num_obs, 1)| lmu, sigma) ]';
  }

  // grr... needed due to parsers constraints on map_rect calls inside
  // functions
  vector sim_obs_map_rect(int[] sidx, vector theta, vector[] xi, real[,] x_r, int[,] x_i) {
    int N = num_elements(sidx) - 1;
    int num_outputs = sidx[N+1]-1;
    vector[num_outputs] sim;

    for(i in 1:N) {
      int start = sidx[i];
      int end = sidx[i+1]-1;
      sim[start:end] = subject_obs_mr(theta, xi[i], x_r[i], x_i[i]);
    }
    return sim;
  }

  // simulates for a posterior draw the full data set based on the same input as the stan model itself
  vector simulate_model(int N, int[] id, vector lndv, vector time, vector amt, int[] cmt, int[] mdv, int[] evid, int[] addl, vector tau, vector phi, vector[] Eta) {
    vector[N] dv = lndv;
    int use_dose_log_amt = 1;
    #include "utils/nm_defs.stan"
    #include "utils/nm_defs_2.stan"
    vector[O] lmu;
    //vector[1] Eta_internal[J];
    // #include "utils/nm_checks.stan" avoid print outputs

    //for(j in 1:J)
    //  Eta_internal[j,1] = Eta[j,1];
    lmu = sim_obs_map_rect(obs_sidx, phi, Eta, rect_nm_real, rect_nm_int);
    return lmu;
  }

  // simulates for a posterior draw the full data set based on the same input as the stan model itself
  vector simulate_model_rng(int N, int[] id, vector lndv, vector time, vector amt, int[] cmt, int[] mdv, int[] evid, int[] addl, vector tau, vector phi, vector[] Eta) {
    int O = N-sum(mdv);
    vector[O] lmu = simulate_model(N, id, lndv, time, amt, cmt, mdv, evid, addl, tau, phi, Eta);
    real sigma_p = phi[5];
    real sigma_y = phi[6];
    vector[O] sigma = sqrt( square(sigma_y) + square(sigma_p * exp(-lmu)));
    vector[O] sim_pred;

    for (i in 1:O) {
      sim_pred[i] = normal_rng(lmu[i], sigma[i]);
    }
    return sim_pred;
  }

  // defines pk_system and more (DO NOT USE IT RIGHT NOW; HAS BUGS)
  
  // we fit a 1-cmt oral dosing situation
  matrix pk_system(vector lref, vector Dt, vector theta, real[] x_r, int[] x_i) {
    // as we fitting a 1-cmt oral dosing situation such that k1=k12 (all
    // mass out of 1 goes to 2)
    //reject("Please do not use this pk_system facility right now!");
    return(pk_1cmt_metabolite(lref, Dt, theta[1], theta[1], theta[2], 0, 0));
  }

  // note that n are the additional doses to be added such that in total
  // n+1 are added
  matrix pk_system_addl(vector lref, vector Dt, int cmt, real lamt, real tau, int n, vector theta, real[] x_r, int[] x_i) {
    matrix[num_elements(Dt), num_elements(lref)] lstate;
    matrix[num_elements(Dt), num_elements(lref)] lstate_mdose;
    vector[num_elements(lref)] lref_mdose;
    int S;
    //reject("Please do not use this pk_system_addl facility right now!");
    
    // evolve reference state freely...
    lstate = pk_system(lref, Dt, theta, x_r, x_i);
    
    // ... and add the extra doses correctly time-shifted
    S = num_elements(lref);
    lref_mdose = rep_vector(-35, S);
    lref_mdose[cmt] = lamt;
    //if(prod(Dt - tau * n) < 0) reject("All requested times must be past the last dosing addl event.");
    /*
      for(i in 1:num_elements(Dt))
      if((Dt[i] - tau * n) < 0)
    reject();
    ("All requested times must be past the last dosing addl event.");
    */
    lstate_mdose = pk_1cmt_metabolite(lref_mdose, Dt - tau * n, theta[1], theta[1], theta[2], tau, n+1);
    for(s in 1:S)
      for(t in 1:num_elements(Dt))
        lstate[t,s] = log_sum_exp(lstate_mdose[t,s], lstate[t,s]);
    return(lstate);
  }

  matrix distance(real[] x1,
  real[] x2,
  int N1,
  int N2,
  real length_scale) {
    matrix[N1, N2] dist;
    for (n1 in 1:N1) {
      for (n2 in 1:N2) {
        dist[n1,n2] = (x1[n1]- x2[n2])/length_scale/length_scale;
      }
    }       
    return dist;
  }
  matrix cov_lin(real[] x1, real[] x2, int N1, int N2, real c, real b, real v) {
    matrix[N1, N2] K;
    for (n1 in 1:N1) {
      for (n2 in 1:N2) {
        K[n1,n2] = b + v*(x1[n1]-c)*(x2[n2]-c);
      }
    }  
    return K;
  }
  matrix dx1_cov_lin(real[] x1, real[] x2, int N1, int N2, real c, real b, real v) {
    matrix[N1, N2] K;
    for (n1 in 1:N1) {
      for (n2 in 1:N2) {
        K[n1,n2] =v*(x2[n2]-c);
      }
    }  
    return K;
  }
  matrix dx2_cov_lin(real[] x1, real[] x2, int N1, int N2, real c, real b, real v) {
    matrix[N1, N2] K;
    for (n1 in 1:N1) {
      for (n2 in 1:N2) {
        K[n1,n2] =v*(x1[n1]-c);
      }
    }  
    return K;
  }
  matrix dx1dx2_cov_lin(real[] x1, real[] x2, int N1, int N2, real c, real b, real v) {
    matrix[N1, N2] K;
    for (n1 in 1:N1) {
      for (n2 in 1:N2) {
        K[n1,n2] = v;
      }
    }  
    return K;
  }
  matrix dx1_cov_exp_quad(real[] x1d,
  real[] x2,
  int N1d,
  int N2,
  real alpha,
  real length_scale) {
    matrix[N1d, N2] dist;
    matrix[N1d, N2] Sigma_reg;
    matrix[N1d, N2] Sigma;
    dist = distance(x1d,x2, N1d, N2, length_scale);
    Sigma_reg = cov_exp_quad(x1d, x2, alpha, length_scale);
    Sigma = -1. * dist .* Sigma_reg;
    return Sigma;
  }
  matrix dx1dx2_cov_exp_quad(real[] x1d,
  real[] x2d,
  int N1d,
  int N2d,
  real alpha,
  real length_scale) {
    matrix[N1d, N2d] dist;
    matrix[N1d, N2d] Sigma_reg;
    matrix[N1d, N2d] Sigma;
    dist = distance(x1d,x2d, N1d, N2d, length_scale);
    Sigma_reg = cov_exp_quad(x1d, x2d, alpha, length_scale);
    Sigma = -1*Sigma_reg .* dist .* dist + Sigma_reg ./ length_scale ./ length_scale;
    return Sigma;
  }
  
  matrix cov_exp_quad_full(real[] x,
  real[] xd,
  int N,
  int Nd,
  real alpha,
  real length_scale, real c, real b, real v) {
    matrix[N + Nd, N + Nd] Sigma;
    vector[N + Nd] eig;
    Sigma[1:N,1:N] = cov_exp_quad(x,alpha,length_scale) .* cov_lin(x,x,N,N,c,b,v);
    Sigma[N+1:N+Nd, 1:N] = dx1_cov_exp_quad(xd,x, Nd,N,alpha,length_scale) .* cov_lin(xd,x,Nd,N,c,b,v) + cov_exp_quad(xd,x,alpha,length_scale) .* dx1_cov_lin(xd,x,Nd,N,c,b,v);
    Sigma[1:N, N+1:N+Nd] = -1.0*dx1_cov_exp_quad(x,xd,N,Nd,alpha,length_scale) .* cov_lin(x,xd,N,Nd,c,b,v) + cov_exp_quad(x,xd,alpha,length_scale) .* dx2_cov_lin(x,xd,N,Nd,c,b,v);
    Sigma[N+1:N+Nd, N+1:N+Nd] = dx1dx2_cov_exp_quad(xd,xd,Nd,Nd,alpha, length_scale) .* cov_lin(xd,xd,Nd,Nd,c,b,v) + cov_exp_quad(xd,alpha,length_scale) .* dx1dx2_cov_lin(xd,xd,Nd,Nd,c,b,v) + dx1_cov_exp_quad(xd,xd,Nd,Nd,alpha,length_scale) .* dx2_cov_lin(xd,xd,Nd,Nd,c,b,v) -1.0*dx1_cov_exp_quad(xd,xd,Nd,Nd,alpha,length_scale) .* dx1_cov_lin(xd,xd,Nd,Nd,c,b,v);
    eig = eigenvalues_sym(Sigma);
    #for (i in 1:N+Nd) {
    #  if (eig[i] < 0) {
    #    print("Info: eig = ", eig[i]);
    #    print("Info: length_scale = ", length_scale);
    #    print("Info: alpha = ", alpha);
    #    print("Info: c = ", c);
    #    print("Info: b = ", b);
    #    print("Info: v = ", v);
    #    break;
    #  }
    #}
    return Sigma;
  }

  matrix cov_exp_quad_full2(real[] x1,
  real[] xd1,
  real[] x2,
  real[] xd2,
  int N1,
  int Nd1,
  int N2,
  int Nd2,
  real alpha,
  real length_scale, real c, real b, real v) {
    matrix[N1 + Nd1, N2 + Nd2] Sigma;
    real det;
    Sigma[1:N1,1:N2] = cov_exp_quad(x1,x2,alpha,length_scale) .* cov_lin(x1,x2,N1,N2,c,b,v);
    Sigma[N1+1:N1+Nd1, 1:N2] = dx1_cov_exp_quad(xd1,x2, Nd1,N2,alpha,length_scale) .* cov_lin(xd1,x2,Nd1,N2,c,b,v) + cov_exp_quad(xd1,x2,alpha,length_scale) .* dx1_cov_lin(xd1,x2,Nd1,N2,c,b,v);
    Sigma[1:N1, N2+1:N2+Nd2] = -1.0*dx1_cov_exp_quad(x1,xd2,N1,Nd2,alpha,length_scale) .* cov_lin(x1,xd2,N1,Nd2,c,b,v) + cov_exp_quad(x1,xd2,alpha,length_scale) .* dx2_cov_lin(x1,xd2,N1,Nd2,c,b,v);
    Sigma[N1+1:N1+Nd1, N2+1:N2+Nd2] = dx1dx2_cov_exp_quad(xd1,xd2,Nd1,Nd2,alpha, length_scale) .* cov_lin(xd1,xd2,Nd1,Nd2,c,b,v) + cov_exp_quad(xd1,xd2,alpha,length_scale) .* dx1dx2_cov_lin(xd1,xd2,Nd1,Nd2,c,b,v) + dx1_cov_exp_quad(xd1,xd2,Nd1,Nd2,alpha,length_scale) .* dx2_cov_lin(xd1,xd2,Nd1,Nd2,c,b,v) -1.0*dx1_cov_exp_quad(xd1,xd2,Nd1,Nd2,alpha,length_scale) .* dx1_cov_lin(xd1,xd2,Nd1,Nd2,c,b,v);
    #det = determinant(Sigma);
    #if (det < 0) {
    #  print("Info: length_scale = ", length_scale);
    #  print("Info: alpha = ", alpha);
    #  print("Info: c = ", c);
    #  print("Info: b = ", b);
    #  print("Info: v = ", v);
    #}
    return Sigma;
  }
  // recommended by Rob in Stan manual as prior for length-scale
  real generalized_inverse_gaussian_lpdf(real x, int p, real a, real b) {
    return p * 0.5 * log(a / b)
      - log(2 * modified_bessel_second_kind(p, sqrt(a * b)))
      + (p - 1) * log(x)
      - (a * x + b / x) * 0.5;
  }
  // Gives latent f given all real + fake data. x contains the true values, xd contains the derivative values.
  matrix latent_L_K(real[] x, real[] xd, int n, int nd, real rho, real alpha, real delta, real c, real b, real v) {
    // can we avoid the rho < 0 if statement here?
    matrix[n+nd, n+nd] K = cov_exp_quad_full(x, xd, n, nd, alpha, rho, c,b,v); // < 0 ? rho_fixed : rho
    matrix[n+nd, n+nd] L_K;
    for (j in 1:n+nd)
    	K[j, j] += delta;
    L_K = cholesky_decompose(K);
    return L_K;
  }

  vector gp_pred_rng(real[] x2,
                     real[] xd2,
                       real[] x1,
                       real[] xd1,
                       vector f1,
                       matrix L_K,
                       real rho,
                       real alpha,
                       real delta,
                       real c,
                       real b,
                       real v) {
      int n1 = size(x1);
      int nd1 = size(xd1);
      int n2 = size(x2);
      int nd2 = size(xd2);
      vector[n2 + nd2] f2;
      {
        vector[n1+nd1] L_div_f1;
        matrix[n1+nd1, n2+nd2] k_x1_x2;
        matrix[n1+nd1, n2+nd2] L_div_k_x1_x2;
        matrix[n2+nd2, n2+nd2] diag_delta;
        vector[n2+nd2] f2_mu;
        matrix[n2+nd2, n2+nd2] f2_cov;
        L_div_f1 = mdivide_left_tri_low(L_K, f1);
        k_x1_x2 = cov_exp_quad_full2(x1, xd1, x2, xd2, n1, nd1, n2, nd2, alpha, rho, c,b,v); //< 0 ? rho_fixed : rho
        L_div_k_x1_x2 = mdivide_left_tri_low(L_K, k_x1_x2);
        f2_mu = L_div_k_x1_x2' * L_div_f1;
        f2_cov = cov_exp_quad_full(x2, xd2, n2, nd2, alpha, rho, c, b, v) - L_div_k_x1_x2' * L_div_k_x1_x2; //< 0 ? rho_fixed : rho
        diag_delta = diag_matrix(rep_vector(delta, n2+nd2));
        f2 = multi_normal_rng(f2_mu, f2_cov+diag_delta);
      }
    return f2;
  }

}
data {
  #include "utils/nm_data.stan"
  //#include "../utils/allo_data.stan"
  vector[N] age;
  vector[N] weight;
  int<lower=0> stud[N];
  vector[N] lndv;
  int num_omega; //Number of individual effects if 1, only clearance has individual effect, if 2, also volume has individual effect.
  int<lower=0,upper=1> known[5+2+num_omega];
  vector[5+2+num_omega] prior_mean;
  matrix[5+2+num_omega,5+2+num_omega] prior_Sigma;
  vector[2] mat_mean; //h, lEC50
  vector[2] mat_omega; // omega h, omega lEC50
  int<lower=0,upper=2> maturation_type;  // 0 = no maturation, 1 = parametric maturation, 2 = GP maturation
  real alpha;   //exponent for clearance log linear scaling
  real beta;    //exponent for volume log linear scaling
  int<lower=0> Nt;            // #number of test observations
  real aget[Nt];              // #ages of test observations
  vector[Nt] weightt;             // #weghts of test observations
  int<lower=0> Nf;              // #number of fixed discrepancy values
  int<lower=0> Nds;             // #number of derivative sign values
  int<lower=0> Ndf;             // #number of fixed derivative values
  vector[Nf] agef;                  // #ages where we know the error exactly
  vector[Nds] ageds;                // #ages of derivative sign observations
  vector[Ndf] agedf;                // #ages of fixed derivative observations
  vector[Nf] ef;                  // #values of fixed error where we know it
  vector[Nds] eds;                // #Values of derivative sign observations
  vector[Ndf] edf;                // #Values of fixed derivative observations
  real delta;
  int rho_p;
  real<lower=0> rho_a;
  real<lower=0> rho_b;
  int rho_p2;
  real<lower=0> rho_a2;
  real<lower=0> rho_b2;
  real nu_gp;                   // #Scaling parameter for virtual derivative observations
  real sigma_ef;                // #noise variance for exact error observations
  real sigma_edf;               // #noise variance for exact derivative observations
  real lin_min_mean;
  real lin_min_std;
  real ref_weight_kid;
  real<lower=0> delta_min;
}
transformed data {
  vector[N] dv = lndv;
  int use_dose_log_amt = 1;
  #include "utils/nm_defs.stan"
  #include "utils/nm_defs_2.stan"
  
  int<lower=0> stud_0[J] = stud[first_ind];
  vector[J] ls_weight_0 = log(weight[first_ind]/ref_weight_kid);
  vector[J] age_s = age[first_ind];
  vector[J] weight_s = weight[first_ind];
  int num_param = 5+2+num_omega;
  int<lower=0> num_unknown = num_param-sum(known);
  int<lower=0> num_known = sum(known);
  int<lower=1,upper=num_param> unknown_idx[num_unknown] = which_elem(known, 0);
  int<lower=1,upper=num_param> known_idx[num_known] = which_elem(known, 1);
  cholesky_factor_cov[num_unknown] prior_Sigma_L = cholesky_decompose(prior_Sigma[unknown_idx,unknown_idx]);
  int Nr = J+Nf;
  int Nd = Nds + Ndf;
  real agetd[0];
  real ager[Nr];
  real aged[Nd];
  int ones[Nds];
  real c = log(0.7);
  real b = 1.0;
  
  for (n in 1:J) // we cannot subindex because ages are vectors, but we need reals
    ager[n] = age_s[n];
  for (n in 1:Nf)
    ager[J+n] = agef[n];
  for (n in 1:Nds) {
    aged[n] = ageds[n];
    ones[n] = 1;
  }
  for (n in 1:Ndf)
    aged[Nds+n] = agedf[n];
  #include "utils/nm_checks.stan"

  print("Info: Unknown parameters = ", unknown_idx);
  print("Info: Known   parameters = ", known_idx);
  for(i in 1:num_known) {
    int idx = known_idx[i];
    print("nu[", idx, "] = ", prior_mean[idx]);
  }
}
parameters {
  // log(ka)
  real theta_lka;
  // log(frel)
  real theta_bioavail;
  // log(CL) = log(ke) + theta_lV
  real theta_lCL;
  // log(V)
  real theta_lV; // Real lv75
  // allo_peds
  real theta_allo_peds;
  real lsigma_p;
  real lsigma_y;
  vector[num_omega] Xi[J]; // 1: Personal effect for clearance 2: Personal effect for volume
  vector[num_omega] lomega; // log of standard deviation of personal effect from the group effect
  real<lower=0> h;
  real lEC50;
  real<lower=delta_min> gp_ls;
  real<lower=delta_min> gp_s;
  vector[Nr+Nd] eta;
  real<lower=0.7, upper=1.0> L;
}
transformed parameters {
  vector[5] theta;
  vector[num_unknown] nu_raw = append_row([theta_lka, theta_bioavail, theta_lCL, theta_lV, theta_allo_peds, lsigma_p, lsigma_y]', lomega)[unknown_idx];//vector[5+2+num_omega] nu_raw = append_row([theta_lka, theta_bioavail, theta_lCL, theta_lV, theta_allo_peds, lsigma_p, lsigma_y]', lomega);
  vector[5+2+num_omega] nu;
  real sigma_p;
  real sigma_y;
  vector[num_omega] omega;
  vector[4+2] phi;
  vector[2] Eta[J]; //clearance and volume
  vector[J] cl_clean;
  vector[Nr+Nd] f = rep_vector(0, Nr+Nd);
  matrix[Nr+Nd, Nr+Nd] L_K = rep_matrix(0.0, Nr+Nd, Nr+Nd);
  vector[J] lVi;
  vector[J] lCLi;
  real v = (L-1.0)/((log(5.7)-c)^2);

  // fix the known parameters to their mean values
  nu[known_idx]   = prior_mean[known_idx];
  nu[unknown_idx] = nu_raw;
  
  lVi = log_lin_allo(nu[4], beta, weight_s);
  lCLi = log_lin_allo(nu[3], alpha, weight_s);
  
  theta = nu[1:5];
  sigma_p = exp(nu[6]);
  sigma_y = exp(nu[7]);
  omega = exp(nu[8:]);

  phi = append_row(theta[1:4], [sigma_p, sigma_y]');
  if(maturation_type==1) {
    f[1:J] = log( inv_logit( h * ( log(age_s) - lEC50 ) ));
  } else if(maturation_type==2) {
    L_K = latent_L_K(log(ager), log(aged), Nr, Nd, gp_ls, gp_s, delta, c, b, v);
    f = L_K*eta;
  }
  
  cl_clean = lCLi + f[1:J];
  
  if(num_omega > 0) {
    for(j in 1:J) {
      Eta[j,1] = omega[1] * Xi[j,1] + lCLi[j] + f[j];
      Eta[j,2] = lVi[j];
    }
    if(num_omega > 1) {
      for(j in 1:J) {
        Eta[j,2] = Eta[j,2] + omega[2] * Xi[j,2];
      }
    }
  }
  for(j in 1:J) {
    if(stud_0[j] == 2305) {
      Eta[j,1] += nu[5] * ls_weight_0[j];
    }
  }

}
model {
  // Parametric maturation function
  h ~ normal(mat_mean[1], mat_omega[1]);
  lEC50 ~ normal(mat_mean[2], mat_omega[2] );
  
  // GP priors
  target += generalized_inverse_gaussian_lpdf(gp_ls | rho_p, rho_a, rho_b);   // length scale
  target += generalized_inverse_gaussian_lpdf(gp_s | rho_p2, rho_a2, rho_b2); // magnitude
  L ~ normal(lin_min_mean, lin_min_std);
  
  // gp posterior 
  eta ~ normal(0, 1);
  if(maturation_type==2) {
    // approximates a normal cdf function
    ones ~ bernoulli(inv_logit(nu_gp * eds .* f[Nr+1:Nr+Nds]));      // #Likelihood of virtual derivative observations
      
      // encodes that we know that the maturation is zero for old ones and the derivative is 0 as well by continuity consideration
    edf ~ normal(f[Nr+Nds+1:Nr+Nds+Ndf], sigma_edf);      // #Likelihood of direct derivative observations (for old ones the derivative is 0)
    ef ~ normal(f[J+1:J+Nf], sigma_ef);      // #Likelihood of direct discrepancy observations (for old ones the function value is 1 in natural space / 0 in log space)
  }
  
  // PK model
  nu_raw ~ multi_normal_cholesky(prior_mean[unknown_idx], prior_Sigma_L);
  Xi ~ multi_normal_cholesky(rep_vector(0.0, num_omega), diag_matrix(rep_vector(1.0, num_omega)));

  //stan_next: target += sum(map_rect(subject_lpdf_mr, phi, Eta, rect_nm_real, rect_nm_int));
}
generated quantities {
  real ka = exp(theta[1]);
  real frel = exp(theta[2]);
  real CL = exp(theta[3]);
  real V = exp(theta[4]);
  real CL_alpha_peds = theta[5];
  real ke = exp(theta[3] - theta[4]);
  real Thalf = log(2.0)/ke;
  vector[Nt] ft = rep_vector(0, Nt);
  vector[Nt] clt;
  vector[J] cl;
  vector[J] vol;
  vector[J] y_cl;
  vector[J] lin;
  vector[J] y_lin;
  vector[O] y_pred_clean;
  vector[O] y_pred;
  vector[O] log_lik;
  y_pred_clean = simulate_model(N,id,lndv,time,amt,cmt,mdv,evid,addl,tau,phi,Eta);
  for (i in 1:O) {
    log_lik[i] = normal_lpdf( lndv[obs_ind[i]] | y_pred_clean[i], sqrt( square(sigma_y) + square(sigma_p * exp(-y_pred_clean[i]))));
    y_pred[i] = normal_rng(y_pred_clean[i], sqrt( square(sigma_y) + square(sigma_p * exp(-y_pred_clean[i]))));
  }
  if(maturation_type==1) {
    ft = log(inv_logit( h * ( log(to_vector(aget)) - lEC50 ) ));
  } else if(maturation_type==2){
    ft = gp_pred_rng(log(aget), agetd, log(ager), log(aged), f, L_K, gp_ls, gp_s, delta, c, b, v);
  }
  if(J>0 && stud_0[1] == 2305) {
    ft += nu[5] * log(weightt/ref_weight_kid);
  }
  clt = log_lin_allo(nu[3], alpha, weightt) + ft;
  lin = lCLi;
  cl = lCLi + f[1:J];
  vol = lVi;
  for(j in 1:J) {
    y_cl[j] = Eta[j,1];
    y_lin[j] = Eta[j,1]- f[j];
  }
  y_cl -= to_vector(lin);
  // print("map rect");
  // print(map_rect(subject_lpdf_mr, phi, Eta, rect_nm_real, rect_nm_int));
  // print("obs ind");
  // print(obs_ind);
  //  we need to do the allometric scaling for Theta[2] and/or Lscale before passing to evaluate_model_fast??
}
