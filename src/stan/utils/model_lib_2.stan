  // return the index vectors for the dose and according observation id for each term
  int[,] linear_superimpose_munge_idx(vector obs_time, vector dose_time) {
   int num_obs = num_elements(obs_time);
   int num_dose = num_elements(dose_time);
   int obs_id[num_obs] = seq_int(1, num_obs);
   int dose_id[num_dose] = seq_int(1, num_dose);
   int last_active_dose[num_obs] = find_interval(obs_time, dose_time);
   int num_dose_events = sum(last_active_dose);
   int terms_oid[num_dose_events] = rep_each(obs_id, last_active_dose);
   int terms_did[num_dose_events];
   int terms_sidx[num_obs+1] = make_slice_index(count_elems(terms_oid, obs_id));

    for(i in 1:num_obs) {
      int start = terms_sidx[i];
      int end = terms_sidx[i+1] - 1;
      int num_active = last_active_dose[i];
      terms_did[start:end] = seq_int(1, num_active);
   }

   return {terms_oid, terms_did};
  }

  // create slice index for observation of terms vector
  int[] linear_superimpose_map2obs(vector obs_time, vector dose_time) {
    int num_obs = num_elements(obs_time);
    int obs_id[num_obs] = seq_int(1, num_obs);
    return make_slice_index(count_elems(linear_superimpose_munge_idx(obs_time, dose_time)[1], obs_id));
  }

  // return the number of terms needed
  int linear_superimpose_num_terms(vector obs_time, vector dose_time) {
    return size(linear_superimpose_munge_idx(obs_time, dose_time)[1]);
  }

  real[,] linear_superimpose_munge_real(vector obs_time, vector dose_time, vector dose_amt, vector dose_tau) {
    int num_dose_events = size(linear_superimpose_munge_idx(obs_time, dose_time)[1]);
    int terms_idx[2,num_dose_events] = linear_superimpose_munge_idx(obs_time, dose_time);
    int terms_idx_oid[num_dose_events] = terms_idx[1];
    int terms_idx_did[num_dose_events] = terms_idx[2];
    int num_obs = num_elements(obs_time);
    int num_dose = num_elements(dose_time);
    int obs_id[num_obs] = seq_int(1, num_obs);
    int dose_id[num_dose] = seq_int(1, num_dose);
    vector[num_dose_events] terms_dose_time;
    vector[num_dose_events] terms_dose_amt;
    vector[num_dose_events] terms_dose_tau;
    vector[num_dose_events] terms_obs_time;
    vector[num_dose_events] terms_tad;

    terms_dose_time = left_join_vector(terms_idx_did, dose_time, dose_id);
    terms_dose_amt = left_join_vector(terms_idx_did, dose_amt, dose_id);
    terms_dose_tau = left_join_vector(terms_idx_did, dose_tau, dose_id);
    terms_obs_time = left_join_vector(terms_idx_oid, obs_time, obs_id);
    terms_tad = terms_obs_time - terms_dose_time;

    return {to_array_1d(terms_tad), to_array_1d(terms_dose_amt), to_array_1d(terms_dose_tau)};
  }

  int[,] linear_superimpose_munge_int(vector obs_time, vector dose_time, vector dose_tau, int[] dose_cmt, int[] dose_addl) {
    int num_dose_events = size(linear_superimpose_munge_idx(obs_time, dose_time)[1]);
    int terms_idx[2,num_dose_events] = linear_superimpose_munge_idx(obs_time, dose_time);
    int terms_idx_oid[num_dose_events] = terms_idx[1];
    int terms_idx_did[num_dose_events] = terms_idx[2];
    int num_obs = num_elements(obs_time);
    int num_dose = num_elements(dose_time);
    int obs_id[num_obs] = seq_int(1, num_obs);
    int dose_id[num_dose] = seq_int(1, num_dose);
    int terms_dose_cmt[num_dose_events];
    int terms_dose_addl[num_dose_events];
    vector[num_dose_events] terms_dose_time;
    vector[num_dose_events] terms_obs_time;
    vector[num_dose_events] terms_dose_tau;
    vector[num_dose_events] terms_tad;

    terms_dose_cmt = left_join_int(terms_idx_did, dose_cmt, dose_id);
    terms_dose_addl = left_join_int(terms_idx_did, dose_addl, dose_id);

    terms_dose_time = left_join_vector(terms_idx_did, dose_time, dose_id);
    terms_dose_tau = left_join_vector(terms_idx_did, dose_tau, dose_id);
    terms_obs_time = left_join_vector(terms_idx_oid, obs_time, obs_id);
    terms_tad = terms_obs_time - terms_dose_time;

    // correct the addl records such that these do not extend past the actual observed time
    for(i in 1:num_dose_events) {
      while(terms_dose_tau[i] * terms_dose_addl[i] > terms_tad[i]) {
        terms_dose_addl[i] = terms_dose_addl[i] - 1;
      }
    }

    return { terms_dose_cmt, terms_dose_addl };
  }

// grouped variants of the above
int[] linear_superimpose_num_terms_grouped(int[] obs_sidx,  vector obs_time,
                                         int[] dose_sidx, vector dose_time) {
  int num_groups = size(obs_sidx)-1;
  int num_terms[num_groups] = rep_array(0, num_groups);
  for(i in 1:num_groups) {
    int obs_start = obs_sidx[i];
    int obs_end = obs_sidx[i+1] - 1;
    int dose_start = dose_sidx[i];
    int dose_end = dose_sidx[i+1] - 1;
    num_terms[i] = linear_superimpose_num_terms(obs_time[obs_start:obs_end], dose_time[dose_start:dose_end]);
  }
  return num_terms;
}


// the grouped variant returns the data in a form as needed for map_rect.
// Thus, we return per group (subject) one row with all real data concatenated
real[,] linear_superimpose_munge_grouped_real(int[] obs_sidx, vector obs_time,
                                              int[] dose_sidx, vector dose_time,
                                              vector dose_amt, vector dose_tau) {
  int num_groups = size(obs_sidx)-1;
  int num_terms[num_groups] = linear_superimpose_num_terms_grouped(obs_sidx, obs_time, dose_sidx, dose_time);
  int num_cols = max(num_terms) * 3;
  real rdata[num_groups, num_cols] = rep_array(0.0, num_groups, num_cols);
  for(i in 1:num_groups) {
    int obs_start = obs_sidx[i];
    int obs_end = obs_sidx[i+1] - 1;
    int dose_start = dose_sidx[i];
    int dose_end = dose_sidx[i+1] - 1;
    real gdata[3,num_terms[i]]
      = linear_superimpose_munge_real(obs_time[obs_start:obs_end],
                                      dose_time[dose_start:dose_end],
                                      dose_amt[dose_start:dose_end],
                                      dose_tau[dose_start:dose_end]
                                      );
    rdata[i,1:3*num_terms[i]] = concatenate_real_array(gdata);
  }
  return rdata;
}

int[,] linear_superimpose_munge_grouped_int(int[] obs_sidx, vector obs_time,
                                            int[] dose_sidx, vector dose_time,
                                            vector dose_tau, int[] dose_cmt, int[] dose_addl) {
  int num_groups = size(obs_sidx)-1;
  int num_terms[num_groups] = linear_superimpose_num_terms_grouped(obs_sidx, obs_time, dose_sidx, dose_time);
  int num_obs[num_groups] = adjacent_difference_int(obs_sidx);
  // # of terms, # of observations, int terms, obs sidx wrt to terms
  int num_cols = 1 + 1 + max(num_terms) * 2 + max(num_obs) + 1;
  int idata[num_groups, num_cols] = rep_array(0, num_groups, num_cols);
  for(i in 1:num_groups) {
    int obs_start = obs_sidx[i];
    int obs_end = obs_sidx[i+1] - 1;
    int dose_start = dose_sidx[i];
    int dose_end = dose_sidx[i+1] - 1;
    int gdata[2,num_terms[i]]
      = linear_superimpose_munge_int(obs_time[obs_start:obs_end],
                                     dose_time[dose_start:dose_end],
                                     dose_tau[dose_start:dose_end],
                                     dose_cmt[dose_start:dose_end],
                                     dose_addl[dose_start:dose_end]
                                     );
    idata[i,1] = num_terms[i];
    idata[i,2] = num_obs[i];
    idata[i,3:2+2*num_terms[i]] = concatenate_int_array(gdata);
    idata[i,2+2*num_terms[i]+1:2+2*num_terms[i]+num_obs[i]+1]
      = linear_superimpose_map2obs(obs_time[obs_start:obs_end],
                                   dose_time[dose_start:dose_end]);
  }
  return idata;
}

// log-pk version for oral 1cmt
// expects as input: ka, ke, log(abs(ka/(ka-ke)))
vector lpk_1cmt_oral_mr(vector theta, vector eta, real[] x_r, int[] x_i) {
  // x_r[1] = tad; x_r[2] = log(amt); x_r[3] = tau
  // x_i[1] = cmt; x_i[2] = addl
  real ka = theta[1];
  real ke = theta[2];
  real lk_scale = theta[3]; ## log( abs(ka/(ka-ke)) ) 
  //real lscale = log(x_r[2]) + lk_scale;
  real lscale = x_r[2] + lk_scale;
  int num_dose = x_i[2] + 1;
  return [ log_diff_exp_abs( log_scaled_geometric_series(lscale - ke * x_r[1], ke * x_r[3], num_dose),
                             log_scaled_geometric_series(lscale - ka * x_r[1], ka * x_r[3], num_dose) ) ]';
}


// log-version for oral 2cmt
// inputs: lka, lalpha, lbeta, lA, lB, ka, alpha, beta, A, B
vector lpk_2cmt_oral_mr(vector theta, vector eta, real[] x_r, int[] x_i) {
  // x_r[1] = tad; x_r[2] = log(amt); x_r[3] = tau
  // x_i[1] = cmt; x_i[2] = addl
  real ka = theta[6];
  real alpha = theta[7];
  real beta = theta[8];
  real lka = theta[1];
  real lscale = x_r[2] + lka;
  real lAscaled = theta[4] + lscale;
  real lBscaled = theta[5] + lscale;
  //real lA = log(theta[4]) - log_diff_exp_abs(lka, log(alpha));
  //real lB = log(theta[5]) - log_diff_exp_abs(lka, log(beta));
  real lscaled_term1 = log_scaled_geometric_series(lAscaled - alpha * x_r[1], x_r[3] * alpha, x_i[2] + 1);
  real lscaled_term2 = log_scaled_geometric_series(lBscaled - beta * x_r[1], x_r[3] * beta, x_i[2] + 1);
  real lscaled_term_ka = log_scaled_geometric_series(- ka * x_r[1], x_r[3] * ka, x_i[2] + 1);
  return [ log_sum_exp( log_diff_exp_abs(lscaled_term1, lAscaled + lscaled_term_ka),
                        log_diff_exp_abs(lscaled_term2, lBscaled + lscaled_term_ka) ) ]';
}
