
int num_dose_terms[J] = linear_superimpose_num_terms_grouped(obs_sidx, obs_time, dose_sidx, dose_time);

real rect_dose_real[J,max(num_dose_terms)*3] = linear_superimpose_munge_grouped_real(obs_sidx, obs_time,
                                                                                     dose_sidx, dose_time,
                                                                                     use_dose_log_amt ? dose_lamt : dose_amt,
                                                                                     dose_tau);

int rect_dose_int[J,1+1+max(num_dose_terms)*2 + max(obs_M) + 1] = linear_superimpose_munge_grouped_int(obs_sidx, obs_time,
                                                                                                       dose_sidx, dose_time,
                                                                                                       dose_tau, dose_cmt, dose_addl);

real rect_obs_real[J,2*max(obs_M)] = ragged2rect_real(obs_M, append_col(obs_dv, obs_time));
int rect_obs_int[J,max(obs_M)] = ragged2rect_int(obs_M, swap_index_int({obs_cmt}));

int num_rect_cols_real = 3*max(num_dose_terms) + 2*max(obs_M);
int num_rect_cols_int = 1+1+max(num_dose_terms)*2 + max(obs_M) + 1 + max(obs_M);

real rect_nm_real[J, num_rect_cols_real] = append_real_array_cols(rect_obs_real, rect_dose_real);
int rect_nm_int[J, 1+2+num_rect_cols_int] = append_int_array_cols_4(swap_index_int({obs_M}), rep_array({max(obs_M), 2*max(obs_M)}, J), rect_obs_int, rect_dose_int);
