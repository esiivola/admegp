
// this function must have been defined prior to this include
// forward declare pk_mr function which handles a dose-stratum with the data input
// x_r[1] = tad; x_r[2] = amt; x_r[3] = tau
// x_i[1] = cmt; x_i[2] = addl
// vector pk_mr(vector theta, vector eta, real[] x_r, int[] x_i);

// pk per subject in the map_rect version for PK formulated on the
// linear scale
vector pk_subject_mr(vector theta, vector eta, real[] x_r, int[] x_i) {
  int num_dose_terms = x_i[1];
  int num_obs = x_i[2];
  int ioffset = 2;
  int terms_dose_cmt[num_dose_terms] = extract_1d_int(x_i, ioffset, num_dose_terms, 1);
  int terms_dose_addl[num_dose_terms] = extract_1d_int(x_i, ioffset, num_dose_terms, 2);
  int terms_oid_sidx[num_obs+1] = x_i[ioffset + 2*num_dose_terms + 1 : ioffset + 2*num_dose_terms + num_obs + 1];
  vector[num_dose_terms] terms_conc;
  vector[num_obs] obs_conc;

  //stan_next: terms_conc = map_rect(pk_mr,
  //stan_next:                       append_row(theta, eta),
  //stan_next:                       rep_array(rep_vector(0, 0), num_dose_terms),
  //stan_next:                       to_array_2d_colwise_real(x_r[1:3*num_dose_terms], num_dose_terms, 3),
  //stan_next:                       swap_index_int({terms_dose_cmt, terms_dose_addl}));
  //to_array_2d_colwise_int(x_i[ioffset+1 : ioffset+2*num_dose_terms], 2, num_dose_terms ));

  for(i in 1:num_obs) {
    int start = terms_oid_sidx[i];
    int end = terms_oid_sidx[i+1] - 1;
    if(end < start) {
      obs_conc[i] = 0.0;
    } else {
      obs_conc[i] = sum(terms_conc[start:end]);
    }
  }

  return obs_conc;
}

// pk per subject in the map_rect version for PK formulated on the
// log scale
vector lpk_subject_mr(vector theta, vector eta, real[] x_r, int[] x_i) {
  int num_dose_terms = x_i[1];
  int num_obs = x_i[2];
  int ioffset = 2;
  int terms_dose_cmt[num_dose_terms] = extract_1d_int(x_i, ioffset, num_dose_terms, 1);
  int terms_dose_addl[num_dose_terms] = extract_1d_int(x_i, ioffset, num_dose_terms, 2);
  int terms_oid_sidx[num_obs+1] = x_i[ioffset + 2*num_dose_terms + 1 : ioffset + 2*num_dose_terms + num_obs + 1];
  vector[num_dose_terms] terms_lconc;
  vector[num_obs] obs_lconc;

  //stan_next: terms_lconc = map_rect(pk_mr,
  //stan_next:                       append_row(theta, eta),
  //stan_next:                       rep_array(rep_vector(0, 0), num_dose_terms),
  //stan_next:                       to_array_2d_colwise_real(x_r[1:3*num_dose_terms], num_dose_terms, 3),
  //stan_next:                       swap_index_int({terms_dose_cmt, terms_dose_addl}));
  //to_array_2d_colwise_int(x_i[ioffset+1 : ioffset+2*num_dose_terms], 2, num_dose_terms ));

  for(i in 1:num_obs) {
    int start = terms_oid_sidx[i];
    int end = terms_oid_sidx[i+1] - 1;
    if(end < start) {
      obs_lconc[i] = -25.0;
    } else {
      obs_lconc[i] = log_sum_exp(terms_lconc[start:end]);
    }
  }

  return obs_lconc;
}
