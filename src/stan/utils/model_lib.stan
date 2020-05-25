int[] count_dose_given(vector time, vector dose_time, vector dose_tau, int[] dose_addl) {
  int dose_count[num_elements(time)];
  int time_rank[num_elements(time)];
  int o;
  int O;
  o = 1;
  O = num_elements(time);
  time_rank = find_interval(time, dose_time);
  dose_count = rep_array(0, O);
  //print("time_rank = ", time_rank);
  // first skip all dosings before the first time
  while(o < O && time_rank[o] == 0) { o = o + 1; }
  //print("o = ", o);
  for(i in o:O) {
    int d;
    d = time_rank[i];
    if(dose_tau[d] > 0)
      dose_count[i] = min(floor_div_int(time[i] - dose_time[d], dose_tau[d]), dose_addl[d]);
  }
  return dose_count;
}

int[] count_dose_given_blocked(int[] M, vector time, int[] M_dose, vector dose_time, vector dose_tau, int[] dose_addl) {
  int dose_count[num_elements(time)];
  int B;
  int tl;
  int dl;
  B = num_elements(M);
  tl = 1;
  dl = 1;
  for(b in 1:B) {
    int tu;
    int du;
    tu = tl + M[b] - 1;
    du = dl + M_dose[b] - 1;
    dose_count[tl:tu] = count_dose_given(time[tl:tu], dose_time[dl:du], dose_tau[dl:du], dose_addl[dl:du]);
    tl = tu + 1;
    dl = du + 1;
  }
  return dose_count;
}

/*
 * Need to define an operator which develops a state forward in time
 * by an amount Dt. Input is the state at t, the log-coefficients of
 * the exponentials and the exponential exponents and Dt. Output is
 * the new state after Dt time units for one of the cmts (depends on
 * given coefs).
 *
 * lstate_ref is the initial
 *
 * coefsN is a 2-row matrix with the first row being the
 * log-coefficients and the second row the exponents.
 */
matrix evolve_lsystem(int S, vector Dt, matrix coefs, int[] coefs_map) {
  matrix[num_elements(Dt),S] lsystem;
  int T;
  T = num_elements(Dt);

  // initialize to zero
  lsystem = rep_matrix(-500, T, S);

  for(o in 1:cols(coefs)) {
    int s;
    vector[T] term;
    s    = coefs_map[o];
    term = coefs[1,o] + Dt * coefs[2,o];
    for(t in 1:T)
      lsystem[t,s] = log_sum_exp(lsystem[t,s], term[t]);
  }
  return(lsystem);
}

real geometric_series(real r, int n) {
  if(n == 1) return 1.0;
  return (1.0 - pow(r, n)) / (1.0 - r);
}


// log_s = log scaling factor
// log_r = log of summed factor raised to n
// N = number of terms of geometric series
real log_scaled_geometric_series(real log_s, real log_r, int N) {
  // in natural space: s * sum_{n=0}^{N-1} r^n = s * (1-r^N)/(1-r)
  // if s is very small, we need to be careful
  if(N == 1) return log_s;
  if(log_r == 0) {
    return log_s + log(N);
  } else if(log_r > 0) {
    real log_r_diff = log_diff_exp(log_r, 0.0);
    return log_diff_exp( log_s + N * log_r - log_r_diff, log_s - log_r_diff );
  } else { // log_r < 0
    real log_r_diff = log1m_exp(log_r);
    return log_diff_exp( log_s - log_r_diff, log_s + N * log_r - log_r_diff );
  }
  return not_a_number();
}

real lgeometric_series(real log_r, int N) {
  return log_scaled_geometric_series(0.0, log_r, N);
}

/*
 * models 2 cmts 1 and 2 which each have an elimination rate k1 / k20
 * and we have a flow from 1 to 2, k12. This models a 1cmt oral dosing
 * (k1=k12) and/or a building block of a metabolite. Finally all mass
 * exiting cmt 2 is recorded in cmt 3.
 *
 * cmt1 -- k12 --> cmt2
 *  |               |
 *  k1              k20
 *  |               |
 * None            cmt3 / Depot
 */
matrix pk_1cmt_metabolite_depot(vector lref, vector Dt, real lk1, real lk12, real lk20, real tau, int n) {
  matrix[2,8] coefs1;
  int coefs1_map[8];
  int coefs1_zero[8];
  matrix[2,6] coefs2;
  int coefs2_map[6];
  int coefs2_zero[6];
  // positive terms
  matrix[num_elements(Dt),3] lsystem1;
  // negative terms (due to Bateman function)
  matrix[num_elements(Dt),2] lsystem2;
  real ldeltaSq;
  real nk1;
  real nk20;

  coefs1_zero = rep_array(0, 8);
  coefs2_zero = rep_array(0, 6);

  ldeltaSq = 2*log_diff_exp_abs(lk1, lk20);

  nk1  = -exp(lk1);
  nk20 = -exp(lk20);

  // setup coefficient matrix and coef index vectors
  coefs1_map[1] = 1;
  coefs1[1,1] = lref[1];
  coefs1[2,1] = nk1;

  coefs1_map[2] = 2;
  coefs1[1,2] = lref[1] + lk12 - ldeltaSq + lk1;
  coefs1[2,2] = nk20;
  coefs1_map[3] = 2;
  coefs1[1,3] = lref[1] + lk12 - ldeltaSq + lk20;
  coefs1[2,3] = nk1;

  coefs1_map[4] = 2;
  coefs1[1,4] = lref[2];
  coefs1[2,4] = nk20;

  // whatever is in the depot cmt doesn´t go away
  coefs1_map[5] = 3;
  coefs1[1,5] = log_sum_exp(lref[3], lref[2]);
  coefs1[2,5] = 0;
  coefs1_zero[5] = 1;

  coefs1_map[6] = 3;
  coefs1[1,6] = lref[1] + lk12 + lk20 - ldeltaSq + log_sum_exp(lk1 - lk20, lk20 - lk1);
  coefs1[2,6] = 0;
  coefs1_zero[6] = 1;

  coefs1_map[7] = 3;
  coefs1[1,7] = lref[1] + lk12 + lk20 - ldeltaSq;
  coefs1[2,7] = nk1;

  coefs1_map[8] = 3;
  coefs1[1,8] = lref[1] + lk12 + lk20 - ldeltaSq;
  coefs1[2,8] = nk20;

  // for the negative terms we only use a two cmts; hence 2 is
  // relabeled to 1, and 3 to 2
  coefs2_map[1] = 1;
  coefs2[1,1] = lref[1] + lk12 - ldeltaSq + lk1;
  coefs2[2,1] = nk1;
  coefs2_map[2] = 1;
  coefs2[1,2] = lref[1] + lk12 - ldeltaSq + lk20;
  coefs2[2,2] = nk20;

  coefs2_map[3] = 2;
  coefs2[1,3] = lref[2];
  coefs2[2,3] = nk20;

  coefs2_map[4] = 2;
  coefs2[1,4] = lref[1] + lk12 - ldeltaSq + lk20 + log(2);
  coefs2[2,4] = 0;
  coefs2_zero[4] = 1;

  coefs2_map[5] = 2;
  coefs2[1,5] = lref[1] + lk12 - ldeltaSq + lk20 + lk1 - lk20;
  coefs2[2,5] = nk20;

  coefs2_map[6] = 2;
  coefs2[1,6] = lref[1] + lk12 - ldeltaSq + lk20 + lk20 - lk1;
  coefs2[2,6] = nk1;

  // in case the initial state is dosed in a regular pattern, we can
  // take advantage of the geometric series here by modifing the
  // coefficients
  if(n>1) {
    real logn;
    logn = log(n);
    for(i in 1:8) {
      if(coefs1_zero[i]) {
        coefs1[1,i] = coefs1[1,i] + logn;
      } else {
        coefs1[1,i] = coefs1[1,i] + lgeometric_series(coefs1[2,i] * tau, n);
      }
    }
    for(i in 1:6) {
      if(coefs2_zero[i]) {
        coefs2[1,i] = coefs2[1,i] + logn;
      } else {
        coefs2[1,i] = coefs2[1,i] + lgeometric_series(coefs2[2,i] * tau, n);
      }
    }
  }

  //print("AFTER: coefs1 = ", coefs1);
  //print("AFTER: coefs2 = ", coefs2);

  lsystem1 = evolve_lsystem(3, Dt, coefs1, coefs1_map);
  lsystem2 = evolve_lsystem(2, Dt, coefs2, coefs2_map);

  //print("lsystem1 = ", lsystem1);
  //print("lsystem2 = ", lsystem2);

  // final system is the difference of the two solutions
  for(t in 1:num_elements(Dt)) {
    lsystem1[t,2] = log_diff_exp(lsystem1[t,2], lsystem2[t,1]);
    lsystem1[t,3] = log_diff_exp(lsystem1[t,3], lsystem2[t,2]);
  }

  return(lsystem1);
}

// same as above, but no depot cmt
matrix pk_1cmt_metabolite(vector lref, vector Dt, real lk1, real lk12, real lk20, real tau, int n) {
  matrix[2,4] coefs1;
  int coefs1_map[4];
  matrix[2,2] coefs2;
  int coefs2_map[2];
  // positive terms
  matrix[num_elements(Dt),2] lsystem1;
  // negative terms (due to Bateman function)
  matrix[num_elements(Dt),1] lsystem2;
  real ldeltaSq;
  real nk1;
  real nk20;

  ldeltaSq = 2*log_diff_exp_abs(lk1, lk20);

  nk1  = -exp(lk1);
  nk20 = -exp(lk20);

  // setup coefficient matrix and coef index vectors
  coefs1_map[1] = 1;
  coefs1[1,1] = lref[1];
  coefs1[2,1] = nk1;

  coefs1_map[2] = 2;
  coefs1[1,2] = lref[1] + lk12 - ldeltaSq + lk1;
  coefs1[2,2] = nk20;
  coefs1_map[3] = 2;
  coefs1[1,3] = lref[1] + lk12 - ldeltaSq + lk20;
  coefs1[2,3] = nk1;

  coefs1_map[4] = 2;
  coefs1[1,4] = lref[2];
  coefs1[2,4] = nk20;

  // for the negative terms we only use a two cmts; hence 2 is
  // relabeled to 1, and 3 to 2
  coefs2_map[1] = 1;
  coefs2[1,1] = lref[1] + lk12 - ldeltaSq + lk1;
  coefs2[2,1] = nk1;
  coefs2_map[2] = 1;
  coefs2[1,2] = lref[1] + lk12 - ldeltaSq + lk20;
  coefs2[2,2] = nk20;

  // in case the initial state is dosed in a regular pattern, we can
  // take advantage of the geometric series here by modifing the
  // coefficients
  if(n>1) {
    for(i in 1:4) {
      coefs1[1,i] = coefs1[1,i] + lgeometric_series(coefs1[2,i] * tau, n);
    }
    for(i in 1:2) {
      coefs2[1,i] = coefs2[1,i] + lgeometric_series(coefs2[2,i] * tau, n);
    }
  }

  //print("AFTER: coefs1 = ", coefs1);
  //print("AFTER: coefs2 = ", coefs2);

  lsystem1 = evolve_lsystem(2, Dt, coefs1, coefs1_map);
  lsystem2 = evolve_lsystem(1, Dt, coefs2, coefs2_map);

  //print("lsystem1 = ", lsystem1);
  //print("lsystem2 = ", lsystem2);

  // final system is the difference of the two solutions
  for(t in 1:num_elements(Dt)) {
    lsystem1[t,2] = log_diff_exp(lsystem1[t,2], lsystem2[t,1]);
  }

  return(lsystem1);
}

// calculates PK metrics, for now only SS concentration (per unit of dose administered)
vector pk_1cmt_metabolite_metrics(real tau, real lk1, real lk12, real lk20) {
  vector[3] metrics;
  real k1;
  real k20;

  k1  = exp(lk1);
  k20 = exp(lk20);

  // SS in main cmt
  metrics[1] = - log1m_exp(-k1*tau) - k1*tau;

  // SS in metabolite cmt
  metrics[2] = lk12
    - log_diff_exp_abs(lk1, lk20)
    + log_diff_exp_abs(-k20*tau  - log1m_exp(-k20*tau), -k1 * tau - log1m_exp(-k1*tau));

  // Rac in main cmt
  metrics[3] = - log1m_exp(-k1*tau);

  return(metrics);
}

// Oral 2cmt analytical solutions: Needs alignment with other
// functions!

/** ANALYTICAL SOLUTION:
 * Calculates the 2-cmt model for one patient given as input nonmem
 * type data and as input parameters the logarith of the micro rate
 * and micro constants
 *
 * returns the log-concentration of the central compartement
 *
 * No ADDL dosing supported right now; needs update
 **/
matrix pk_oral_2cmt(vector state0, vector Dt,
                    real lka, real lalphaR, real lbetaR, real lA, real lB) {
  real lstateRefOral; // ref state for the 2-cmt with oral cmt (only the oral cmt)
  real lstateRef[2];  // ref state for the 2-cmt without oral cmt
  int N;
  real alphaR;
  real betaR;
  real ka;
  real lAt;
  real lBt;
  matrix[num_elements(Dt),3] lstate;
  real lk12;
  real lk21;
  real lD;
  real lad2;
  real lbd2;
  real ltemp;

  N = num_elements(Dt);

  ka = exp(lka);
  alphaR = exp(lalphaR);
  betaR  = exp(lbetaR);

  // Bateman coefficients
  lAt = lA + lka - log_diff_exp_abs(lka, lalphaR);
  lBt = lB + lka - log_diff_exp_abs(lka, lbetaR );

  // needed constant for the unobserved peripheral cmt C
  lD = log_diff_exp(lalphaR, lbetaR);   // discriminant which is always positive
  ltemp = log_sum_exp(lB + lD, lbetaR);
  lk12 = log_diff_exp(log_sum_exp(lalphaR, lbetaR), log_sum_exp(2*ltemp, lalphaR + lbetaR) - ltemp );
  lk21 = log_diff_exp(lalphaR, lA + lD);

  lad2 = 2 * log_diff_exp_abs(lalphaR, lka);
  lbd2 = 2 * log_diff_exp_abs(lbetaR , lka);

  // by convention time starts just at the first observation
  lstateRefOral = state0[1];
  lstateRef[1]  = state0[2];
  lstateRef[2]  = state0[3];
  for(i in 1:N) {
    lstate[i,1] = lstateRefOral - ka * Dt[i];
    // solution for the concentration which is in the central and
    // peripheral cmt
    lstate[i,2] = lstateRef[1] + log_sum_exp(lA - alphaR * Dt[i], lB - betaR * Dt[i]);
    lstate[i,3] = lstateRef[2] + log_sum_exp(lB - alphaR * Dt[i], lA - betaR * Dt[i]);

    // other changes in the state can only meaningful be calculated
    // if Dt[i] is large enough to allow for diffusion to occur
    if(Dt[i] > 0.) {
      lstate[i,2] = log_sum_exp(lstate[i,2], lstateRef[2] + lD - lk12 + lA + lB + log_diff_exp(- betaR * Dt[i], - alphaR * Dt[i]) );
      lstate[i,3] = log_sum_exp(lstate[i,3], lstateRef[1] + lk12 - lD + log_diff_exp(-betaR * Dt[i], -alphaR * Dt[i]));


      // add in the part which stems from oral cmt which results in
      // the superposition of Bateman functions in the main cmt
      lstate[i,2] = log_sum_exp(lstate[i,2], lstateRefOral + log_sum_exp(lAt + log_diff_exp_abs( -alphaR * Dt[i], -ka * Dt[i]),
                                                                         lBt + log_diff_exp_abs( -betaR  * Dt[i], -ka * Dt[i]))
                                );
      // last, we add into the peripheral cmt the effect of the oral
      // cmt dosing
      //lstate[i,3] = log_sum_exp(lstate[i,3], lstateRefOral + lk12 + lka - lD - ka * Dt[i] + log( D*A2i*B2i + A2i * exp(- (alphaR-ka) * Dt[i]) - B2i * exp(-(betaR-ka)*Dt[i])   ) );
      // the huge expression below is a sign-sorted version of (brute force)
      // k12 * ka /[ (alpha-ka) * (beta-ka) ] * [ exp(-ka * t) - (ka-beta)/(alpha-beta) * exp(-alpha*t) + (ka-alpha)/(alpha-beta) * exp(-beta*t) ]
      lstate[i,3] = log_sum_exp(lstate[i,3], lstateRefOral + lk12 + lka - lD - lad2 - lbd2 +
                                log_diff_exp(log_sum_exp(log_sum_exp(lD + log_sum_exp(lalphaR + lbetaR, 2*lka) - ka * Dt[i], lalphaR + lbd2 - alphaR * Dt[i]), lka    + lad2 - betaR * Dt[i] ),
                                             log_sum_exp(log_sum_exp(lD + lka + log_sum_exp(lalphaR,   lbetaR) - ka * Dt[i], lka     + lbd2 - alphaR * Dt[i]), lbetaR + lad2 - betaR * Dt[i] )
                                             )
                                );
    }
  }

  return lstate;
}

  /*
   * transform from log of lka, CL, V1, Q, V2 to log of lka, alpha,
   * beta, A, B
   * NOTE: A & B are not scaled by any volume!
   */
  vector trans_oral2cmt_macro2micro(real lka, real lCL, real lV1, real lQ, real lV2) {
    // first project onto "k" parametrization and then onto micro
    // and macro constants
    vector[3] lk = [lCL - lV1,   // lk
                    lQ  - lV1,   // lk12
                    lQ  - lV2]'; // lk21
    real lkSum = log_sum_exp(lk);
    real lk1_times_k21 = lk[1] + lk[3];
    vector[5] mm;

    mm[1] = lka;

    // check that discriminat is for all patients real, i.e. that 
    // (k10 + k12 + k21)^2 - 4 k10 k21 > 0
    // otherwise the eigenvalues would be imaginary leading to oscillations
    if(2*lkSum < log(4.0) + lk1_times_k21) {
      reject("System discriminant must be real!");
    }

    // first and second rate constant are roots for pq-quadratic form
    // with p = -kSum, q = k21 * k = k[3] * k[1]
    
    // log of first rate constant alpha
    //mm[2] = log(0.5) + log_sum_exp(lkSum, 0.5 * log_diff_exp(2.0*lkSum, log(4.0) + lk[1] + lk[3]) );
    mm[2] = log(0.5) + log_sum_exp(lkSum, 0.5 * log_diff_exp(2.0*lkSum, log(4.0) + lk1_times_k21) );

    // log of second rate constant beta (Vieta)
    //mm[3] = lk[1] + lk[3] - mm[2];
    mm[3] = lk1_times_k21 - mm[2];

    // macro constants
    mm[4] = log_diff_exp_abs(lk[3], mm[2]) - log_diff_exp_abs(mm[2], mm[3]);
    mm[5] = log_diff_exp_abs(lk[3], mm[3]) - log_diff_exp_abs(mm[3], mm[2]);

    return(mm);
  }

// forward declare pk system functions
matrix pk_system(vector lref, vector Dt, vector theta, real[] x_r, int[] x_i);


matrix pk_system_addl(vector lref, vector Dt, int cmt, real lamt, real tau, int n, vector theta, real[] x_r, int[] x_i);

// model evaluation function taking dosing (and addl dosing) into
// account for a single patient. The function calculates always with
// respect to a reference state. The initial reference state is the
// initial time and initial_lstate. Upon the first dosing event past
// the initial time, the reference time and state is set to the
// respective dosing event.
matrix pk_model_fast(vector dose_lamt, int[] dose_cmt, vector dose_time, vector dose_tau, int[] dose_addl, int[] dose_next_obs,
                     vector init_lstate, real init_time,
                     vector obs_time, int[] obs_timeRank, int[] obs_dose_given,
                     vector theta,
                     vector lscale,
                     real[] x_r, int[] x_i)
{
  int D = num_elements(dose_lamt);
  int O = num_elements(obs_time);
  int d = 0;
  int o = 1;
  int active_addl = 0; // 0 = FALSE
  int init_ref = 1;
  int num_states = num_elements(init_lstate);
  vector[num_states] lstate_ref = init_lstate;
  matrix[O,num_states] lstate;

  /*
  print("dose_lamt = ", dose_lamt');
  print("dose_cmt = ", dose_cmt);
  print("dose_time = ", dose_time');
  print("dose_tau = ", dose_tau');
  print("dose_addl = ", dose_addl);
  print("dose_next_obs = ", dose_next_obs);
  print("obs_time = ", obs_time');
  print("obs_timeRank = ", obs_timeRank);
  print("obs_dose_given = ", obs_dose_given);
  */

  // skip all dosing records prior to init_time
  while(d < D && (dose_time[d+1] + dose_tau[d+1] * dose_addl[d+1]) < init_time) { d = d + 1; }
  // next, process all elements which are past active dosing
  //active_addl = dose_addl[d] > 0;
  while(o <= O) {
    // first update reference state to be just after the last
    // dosing
    //print("obs o = ", o, "; dose d = ", d);
    while(d != obs_timeRank[o]) {
      int nd = d + 1;
      //print("==> obs o = ", o, "; dose nd = ", nd);
      // the dose of the reference state is not yet the very
      // last dose given before this observation, add it
      if(init_ref) {
        lstate_ref = to_vector(pk_system(lstate_ref, rep_vector(dose_time[nd] - init_time, 1), theta,
                                          x_r, x_i)[1]);
      } else if(active_addl) {
        // in case of an active addl record, we have to super-impose
        // the developed reference state with the emerged dosing which
        // is handled by the _addl function
        lstate_ref = to_vector(pk_system_addl(lstate_ref, rep_vector(dose_time[nd] - dose_time[d], 1), dose_cmt[d], dose_lamt[d], dose_tau[d], dose_addl[d], theta,
                                               x_r, x_i)[1]);
      } else {
        lstate_ref = to_vector(pk_system(lstate_ref, rep_vector(dose_time[nd] - dose_time[d], 1), theta,
                                          x_r, x_i)[1]);
      }
      // the new reference time point is the dosing event; time is now
      // measured as time-after-dose

      // add in the dosing, but only if we have a simple dosing
      // event, i.e. no additional dosings
      active_addl = dose_addl[nd] > 0;
      if(!active_addl) {
        lstate_ref[dose_cmt[nd]] = log_sum_exp(lstate_ref[dose_cmt[nd]], dose_lamt[nd]);
      }
      d = nd;
      // at this point the initial time is not any more the reference
      // state
      init_ref = 0;
    }
    // ok, evolve from reference (last dose or initial) to current
    // observation...
    if(init_ref) {
      lstate[o] = pk_system(lstate_ref, segment(obs_time, o, 1) - init_time, theta,
                             x_r, x_i)[1];
      o = o + 1;
    } else if(active_addl) {
      // ...in case of addl dosing, the effect of the multiple
      // dosing has not yet been added
      // number of dosings given from the active addl records
      int ndose = obs_dose_given[o];
      int event_free = 0;
      // advance as far as we can by counting the number of
      // observations which have the same number of doses given, have
      // the same rank wrt to the time vector
      while((o + event_free) < O && obs_dose_given[o + event_free + 1] == ndose && obs_timeRank[o + event_free + 1] == obs_timeRank[o])
        event_free = event_free + 1;

      lstate[o:(o+event_free)] = pk_system_addl(lstate_ref, segment(obs_time, o, event_free + 1) - dose_time[d], dose_cmt[d], dose_lamt[d], dose_tau[d], ndose, theta,
                                                 x_r, x_i);
      o = o + event_free + 1;
    } else {
      // ... which is simple for non-addl dosing as dose is
      // already merged, evolve as long as no other dosing occurs
      int event_free = dose_next_obs[d];
      lstate[o:(o+event_free-1)] = pk_system(lstate_ref, segment(obs_time, o, event_free) - dose_time[d], theta,
                                              x_r, x_i);
      o = o + event_free;
    }
  }
  {
    row_vector[num_states] lscale_row = to_row_vector(lscale);
    for(i in 1:O)
      lstate[i] = lstate[i] - lscale_row;
  }
  return(lstate);
}

matrix pk_model(vector dose_lamt, int[] dose_cmt, vector dose_time, vector dose_tau, int[] dose_addl,
                vector init_lstate, real init_time,
                vector obs_time,
                vector theta,
                vector lscale,
                real[] x_r, int[] x_i) {
  int obs_timeRank[num_elements(obs_time)] = find_interval(obs_time, dose_time);
  check_addl_dosing(dose_time, dose_tau, dose_addl);
  check_order_increasing(dose_time);
  check_order_increasing(obs_time);
  if (init_time > dose_time[1])
    reject("Initial time must be at or before first dose!");
  return(pk_model_fast(dose_lamt, dose_cmt, dose_time, dose_tau, dose_addl, count_obs_event_free(obs_timeRank, size(dose_cmt)),
                       init_lstate, init_time,
                       obs_time, obs_timeRank,
                       count_dose_given(obs_time, dose_time, dose_tau, dose_addl),
                       theta,
                       lscale,
                       x_r, x_i));
}

matrix evaluate_model_fast(int[] dose_M, vector dose_lamt, int[] dose_cmt, vector dose_time, vector dose_tau, int[] dose_addl, int[] dose_next_obs,
                           matrix init_lstate, vector init_time,
                           int[] obs_M, vector obs_time, int[] obs_timeRank, int[] obs_dose_given,
                           matrix theta,
                           matrix lscale,
                           real[] x_r, int[] x_i) {
  matrix[num_elements(obs_time), cols(init_lstate)] lstate;
  int J = num_elements(dose_M);
  int d = 1;
  int o = 1;
  int S = cols(lscale);

  for(j in 1:J) {
    matrix[obs_M[j],S] lstate_j;
    int d_m = dose_M[j];
    int o_m = obs_M[j];
    //print("Processing patient ", j);
    lstate_j = pk_model_fast(segment(dose_lamt, d, d_m), segment(dose_cmt, d, d_m), segment(dose_time, d, d_m), segment(dose_tau, d, d_m), segment(dose_addl, d, d_m), segment(dose_next_obs, d, d_m)
                              ,to_vector(init_lstate[j]), init_time[j]
                              ,segment(obs_time, o, o_m), segment(obs_timeRank, o, o_m), segment(obs_dose_given, o, o_m)
                              ,to_vector(theta[j])
                              ,to_vector(lscale[j])
                              ,x_r, x_i);

    for(i in 1:o_m)
      lstate[i + o - 1] = lstate_j[i];

    d = d + d_m;
    o = o + o_m;
  }
  return(lstate);
}

// returns only at the observed cmts
vector evaluate_model_fast_cmt(int[] dose_M, vector dose_lamt, int[] dose_cmt, vector dose_time, vector dose_tau, int[] dose_addl, int[] dose_next_obs,
                               matrix init_lstate, vector init_time,
                               int[] obs_M, vector obs_time, int[] obs_timeRank, int[] obs_dose_given, int[] obs_cmt,
                               matrix theta,
                               matrix lscale,
                               real[] x_r, int[] x_i) {
  vector[num_elements(obs_time)] lstate;
  int J = num_elements(dose_M);
  int d = 1;
  int o = 1;
  int S = cols(lscale);

  for(j in 1:J) {
    matrix[obs_M[j],S] lstate_j;
    int d_m = dose_M[j];
    int o_m = obs_M[j];
    //print("Processing patient ", j);
    lstate_j = pk_model_fast(segment(dose_lamt, d, d_m), segment(dose_cmt, d, d_m), segment(dose_time, d, d_m), segment(dose_tau, d, d_m), segment(dose_addl, d, d_m), segment(dose_next_obs, d, d_m)
                              ,to_vector(init_lstate[j]), init_time[j]
                              ,segment(obs_time, o, o_m), segment(obs_timeRank, o, o_m), segment(obs_dose_given, o, o_m)
                              ,to_vector(theta[j])
                              ,to_vector(lscale[j])
                              ,x_r, x_i);

    for(i in 1:o_m)
      lstate[i + o - 1] = lstate_j[i, obs_cmt[i + o - 1] ];

    d = d + d_m;
    o = o + o_m;
  }
  return(lstate);
}

matrix evaluate_model(int[] dose_M, vector dose_lamt, int[] dose_cmt, vector dose_time, vector dose_tau, int[] dose_addl,
                      matrix init_lstate, vector init_time,
                      int[] obs_M, vector obs_time,
                      matrix theta,
                      matrix lscale,
                      real[] x_r, int[] x_i) {
  int obs_timeRank[num_elements(obs_time)] = find_interval_blocked(obs_M, obs_time, dose_M, dose_time);
  check_addl_dosing_blocked(dose_M, dose_time, dose_tau, dose_addl);
  check_order_increasing_blocked(dose_M, dose_time);
  check_order_increasing_blocked(obs_M, obs_time);
  return(evaluate_model_fast(dose_M, dose_lamt, dose_cmt, dose_time, dose_tau, dose_addl, count_obs_event_free_blocked(obs_M, obs_timeRank, dose_M),
                             init_lstate, init_time,
                             obs_M, obs_time, obs_timeRank,
                             count_dose_given_blocked(obs_M, obs_time, dose_M, dose_time, dose_tau, dose_addl),
                             theta,
                             lscale,
                             x_r, x_i));
}

matrix evaluate_model_nm(int[] id, vector time, int[] cmt, int[] evid, vector amt, vector tau, int[] addl, int[] mdv,
                         matrix init_lstate, vector init_time,
                         matrix theta, matrix lscale,
                         real[] x_r, int[] x_i) {
  int dose_ind[count_elem(evid, 1)] = which_elem(evid, 1);

  return(evaluate_model(rle_int(id[dose_ind]), log(amt[dose_ind]), cmt[dose_ind], time[dose_ind], tau[dose_ind], addl[dose_ind],
                        init_lstate, init_time,
                        rle_int(id), time,
                        theta,
                        lscale,
                        x_r, x_i));
}

// returns only the states for which the respective row was specified
// for
vector evaluate_model_nm_cmt(int[] id, vector time, int[] cmt, int[] evid, vector amt, vector tau, int[] addl, int[] mdv,
                             matrix init_lstate, vector init_time,
                             matrix theta, matrix lscale,
                             real[] x_r, int[] x_i) {
  int dose_ind[count_elem(evid, 1)] = which_elem(evid, 1);
  int N = num_elements(time);
  matrix[N,cols(init_lstate)] states;
  vector[N] obs_states;

  states = evaluate_model(rle_int(id[dose_ind]), log(amt[dose_ind]), cmt[dose_ind], time[dose_ind], tau[dose_ind], addl[dose_ind],
                          init_lstate, init_time,
                          rle_int(id), time,
                          theta,
                          lscale,
                          x_r, x_i);
  for(i in 1:N)
    obs_states[i] = states[i,cmt[i]];
  
  return(obs_states);
}

