library(rstan)
library(plyr)
library(RBesT)
library(metaSEM)
library(ggplot2)
library(data.table)

###-------------------------------------------------------###
### true values for the linear system in simulation study ###
###-------------------------------------------------------###
get_doses <- function() c(log(1500), log(1500) ) #log(1500), log(1500))
get_dose_times <- function() c(0, 24)
get_dose_cmts <- function() c(1, 1)
get_dose_taus <- function() c(0,24)

get_lka75 <- function() log(log(2)/2)
get_lv75 <- function() log(10)
get_lcl75 <- function() log(log(2)/12) + get_lv75() 
get_lomega_clearance <- function() log(log(1.5)/1.96) # allow 50% variation
get_lomega_volume <- function() log(log(1.5)/1.96) # allow 50% variation
get_alpha <- function() 0.75
get_beta <- function() 1.0
get_lsigma_y <- function()log(0.1)
get_mat_mean <- function() c(2.94, log(62.9/52.0))
get_mat_omega <- function() c(log(2)/1.96, log(2)/1.96)
get_mean <- function(ag)  (ag)*2.76+3.7 # mean for weight for paricular age
get_std <- function(ag) (ag)*0.5+0.5  # std for weight for particular age
mat <- function(age, h=NULL, lEC50=NULL, maturation=TRUE, doublehill=FALSE, EC50_1=0.7, h_1=6, EC50_2=3.5, h_2=7){ #h=2.94
  if(is.null(h)) {
    h=get_mat_mean()[1]
  }
  if(is.null(lEC50)){
    lEC50=get_mat_mean()[2]
  }
  if(maturation){
    if(doublehill){
      return( 0.5*inv_logit(h_1*(log(age)-log(EC50_1))) +0.5*inv_logit(h_2*(log(age)-log(EC50_2))) )
    } else {
      return(inv_logit(h*(log(age) - lEC50))) # defaults for sirolimus from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4999604/pdf/PSP4-5-411.pdf
    }
  } else {
    return(as.array(rep(1, length(age))))
  }
}

# means for 1-log of ka, 2-non-lin rel bio, 3-log of CL, 4-log of V, 5-theta_allo_peds, 6-lsigma_p, 7-lsigma_y, 8-lomega[1] (noise of clearance), 9-omega[2] (noise for volume)
get_simulated_prior_mean <- function() c(get_lka75(), -100, get_lcl75(), get_lv75(), 0.0, -1e6, get_lsigma_y(), get_lomega_clearance(), get_lomega_volume())
get_simulated_prior_sigma_diag <- function() c(log(2.0), log(1.05), log(1.05), log(1.15), log(8), log(1.2), log(1.05), log(1.1), log(1.1))/1.96
unit_covariance <- bdiagMat(list(c(1), rbind(c(1,-0.8,-0.7), c(-0.8,1,0.5), c(-0.7,0.5,1)), c(1), rbind(c(1,-0.8),c(-0.8,1)), rbind(c(1,0.5), c(0.5,1))))
get_simulated_prior_sigma <- function() get_simulated_prior_sigma_diag() %*% t(get_simulated_prior_sigma_diag()) * unit_covariance

get_known_for_simulation <- function() c(0,1,0,0,1,1,0,0,0)

get_at <- function() as.array(exp(seq(from=log(0.7), to=log(5.7), length.out = 100))) # as.array(seq(from=0.7, to=5.7, by=0.1))
get_at_lin <- function() as.array(seq(from=0.7, to=5.7, length.out = 100)) # as.array(seq(from=0.7, to=5.7, by=0.1))
lin_approx <- function(ref75, pow, w) ref75 + pow*(log(w) -log(75.0))

###-----------------------------------------------###
### Functions for generating data for simulations ###
###-----------------------------------------------###
sample_ages <- function(J, min=0.7, max=5.7) {
  # Samples J ages between min and max so that there is exactly one
  # sample in equally large bin of size (max-min)/J
  ag <- c()
  for (i in 1:J) { #divide to bins and sample one from each
    ag <- c(ag, runif(1, min=(i-1)/J, max=i/J) )
  }
  return((ag)*(max-min)+min) #scale it to between min max
}

gen_pediatric_data <- function(stan_pk_model, ag, wt, ke, use_maturation=TRUE, plotit=FALSE, num_doses=6, prob_measurement=1.0, num_observations_dense=1, double_hill=FALSE, random_effects=FALSE) {
  # Generates pediatric data in nonmem data format given the stan_pk_model,
  # ages of the patients, weights of the patients and individual log clearances
  # (without maturation) of patients. In addition, user can define if Hill function
  # maturation is used, if the results are plotted, number of doses per patient,
  # probability of measuring at each dose and how many measurements are made after
  # the first dose
  
  ## time unit is [h]
  
  J <- length(ag)
  ## these are num_dosed additional doses
  dose_addl <- c(0, num_doses-1)
  init_lstate <- rep(-25, 2)
  # Scale doses by weight (as we would do for adults)
  dose_lamts1 <- lin_approx(get_doses()[1], get_alpha(), wt)
  dose_lamts2 <- lin_approx(get_doses()[2], get_alpha(), wt)
  
  ## simulate per patient parameters of ke and V when assuming allometric scaling:
  ## simulate subject specific ke
  if(use_maturation) {
    ke <- ke + log(mat(ag,doublehill=double_hill))
  }
  
  lV_i <- lin_approx(get_lv75(), get_beta(), wt) + rnorm(J, 0, exp(get_lomega_volume()))
  
  ## assemble Stan model
  ## and export functions to R, used for simulation
  if(! is.na( stan_pk_model) ) {
    expose_stan_functions(stan_pk_model)
  }
  ## define design, i.e 6 observations on densely sampled days (first
  ## and last); in between we only measure the trough concentration
  ## (just before dosing)
  obs_time_dense <- c(seq(0.5,7, length=num_observations_dense))
  
  obs_time <- sort(c(obs_time_dense, seq(24, 24 * num_doses, by = 24) ))
  mean <- get_simulated_prior_mean()
  if(random_effects) {
    sigma <- get_simulated_prior_sigma()
    known <- get_known_for_simulation()
    nu <- mvrnorm(1, mean<-mean, Sigma<-sigma)
    nu[known] <- mean[known]
  } else {
    nu <- mean
  }
  phi <- nu[c(1,2,3,4,6,7)]
  phi[c(5,6)] <- exp(phi[c(5,6)])
  ###------------------------------------------###
  ### create Nonmem data sets for each patient ###
  ###------------------------------------------###
  nm <- data.frame()
  for(i in 1:J) {
    ## observation part for each patient
    time_t =c(obs_time[1], obs_time[-(1)][runif(length(obs_time)-1, min = 0, max = 1)<prob_measurement])
    nm_obs <- data.frame(time=time_t, cmt=2, evid=0, amt=0, tau=0, addl=0, mdv=0, lndv=0, id=i,age=ag[i],weight=wt[i],stud=0)
    ## dosing part for each patient
    nm_dose1 <- data.frame(time=get_dose_times()[1], cmt=1, evid=1, amt=exp(dose_lamts1[i]), tau=get_dose_taus()[1], addl=dose_addl[1], mdv=1, lndv=0, id=i,age=ag[i],weight=wt[i],stud=0)
    nm_dose2 <- data.frame(time=get_dose_times()[2], cmt=1, evid=1, amt=exp(dose_lamts2[i]), tau=get_dose_taus()[2], addl=dose_addl[2], mdv=1, lndv=0, id=i,age=ag[i],weight=wt[i],stud=0)
    
    nm_t <- arrange(rbind(nm_obs, nm_dose1, nm_dose2), time, evid, cmt)
    nm <- rbind(nm, nm_t)
  }
  
  
  ## sort by id, time, evid (first observation, then dosing), finally cmt
  nm <- arrange(nm, id, time, evid, cmt)
  ## simulate design
  vars <- list(N=nrow(nm), id=nm$id, lndv=nm$lndv, time=nm$time, amt=nm$amt, cmt=nm$cmt, mdv=nm$mdv, evid=nm$evid, addl=nm$addl, tau=nm$tau, phi=phi, Eta=lapply(1:J, FUN = function(i) c(ke[i], lV_i[i])))
  #print(vars)
  sim <- do.call(simulate_model_rng, vars)

  ## save simulated true mean conc in sdv
  nm$sdv[nm$evid==0] <- exp(sim)
  ## apply noise
  nm$lndv[nm$evid==0] <- sim
  nm$sdv[nm$cmt != 2] <- nm$lndv[nm$cmt != 2] <- 0
  
  return(nm)
}

get_stan_data <- function(nm, ag, wt, ke,deriv=T, max_age_deriv=3.7, max_age=5.7) {
  # Generates variables for the stan with prior distributions and
  # nonmem data and patient parameters like age and weight
  
  #Prior parameters for GP variance
  rho_a2 <- 5
  rho_b2 <- 1
  rho_p2 <- -1
  
  #Prior parameters for GP length scale
  rho_a <- 5
  rho_b <- 1
  rho_p <- -1
  if(deriv){
    rho_b <- 2
    rho_p <- -6    
  }
  
  # Virtual observations for discrepancy
  agef <- as.array(c(max_age))           #ages where we know the error exactly
  ageds <- as.array(exp(seq(from=log(0.7), to=log(max_age_deriv), length.out = 0)))
  agedf <- as.array(exp(seq(from=log(0.7), to=log(max_age_deriv), length.out = 0)))              #ages of fixed derivative observations
  if(deriv){
    ageds <- as.array(exp(seq(from=log(0.7), to=log(max_age_deriv), length.out = 8)))           #ages of derivative sign observations
    agedf <- as.array(c(5.7))                #ages of fixed derivative observations
  }
  
  Nf <- length(agef)            #number of fixed discrepancy values
  Nds <- length(ageds)             #number of derivative sign values
  Ndf <- length(agedf)              #number of fixed derivative values
  
  ef <- as.array(rep(0, Nf))                   #values of fixed error where we know it
  eds <- as.array(rep(1, Nds))                 #Values of derivative sign observations
  edf <- as.array(rep(0, Ndf))                 #Values of fixed derivative observations
  
  #Ages and weights for prediction of discrepancy
  aget <- get_at()
  weightt <- get_mean(aget)
  Nt <- length(aget)
  
  return(c(nm, list(N=nrow(nm),
                    alpha=get_alpha(),beta=get_beta(),delta_min=0.2,delta2=1e-4,
                    Nf=Nf,Nt=Nt,Nds=Nds,Ndf=Ndf,agef=agef,aget=aget,weightt=weightt,
                    ageds=ageds,agedf=agedf,ef=ef,eds=eds,edf=edf,rho_p=rho_p,
                    rho_b=rho_b,rho_a=rho_a,delta=1e-6,lin_min_mean=0.8,
                    lin_min_std=0.05, nu_gp=10,sigma_ef=0.05,sigma_edf=0.05,ag=ag,
                    wt=wt,mat_mean=get_mat_mean(),mat_omega=get_mat_omega(),
                    rho_p2=rho_p2,rho_b2=rho_b2,rho_a2=rho_a2,ke=ke,num_omega=2,
                    known=get_known_for_simulation(),
                    prior_mean=get_simulated_prior_mean(),
                    prior_Sigma=get_simulated_prior_sigma(), ref_weight_kid=20.0
                    )
          ))
}

gen_data_files <- function(J=25, num_measurements=7, seed=1234, maturation=TRUE, stan_file="stan/oral_1cmt_omni.stan", plotit=FALSE, double_hill=FALSE, deriv=TRUE) {
  # Generate and save data used to run the Stan model
  
  set.seed(seed)
  stan_pk_model = stan_generate_model(stan_file,enable_stan_next=T)
  # Weights and ages of patients
  ag <- seq(from=0.7, to=5.7, length.out=J)# sample_ages(J)# seq(from=0.7, to=5.7, length.out=J)# 
  wt <- get_mean(ag) 
  ke <- lin_approx(get_lcl75(), get_alpha(), wt) + rnorm(J, 0, exp(get_lomega_clearance()))
  nm <- gen_pediatric_data(stan_pk_model,ag,wt, ke, use_maturation=maturation,plotit=plotit, num_doses=num_measurements, prob_measurement=1.01, num_observations_dense=0, double_hill=double_hill)
  stan_data <- get_stan_data(nm, ag, wt, ke, deriv=deriv)
  
  ## save stan data
  dir.create("generated_files")
  stan_rdump(names(stan_data), "generated_files/oral_1cmt_allo_run.data.R", envir=list2env(stan_data))
  
  ## also serialize out Stan program with all includes applied
  
  cat("\n", file="generated_files/oral_1cmt_allo_run.stan", append=TRUE)
  ## finally save Nonmem data set
  write.csv(nm[1], file="generated_files/oral_1cmt_allo_run_nm.csv", quote=FALSE, row.names=FALSE)
}

stan_generate_model <- function(file, ..., postfix="_generated", enable_stan_next=FALSE) {
  # Generate the Stan model from the helper filess
  library(tools)
  generated_file <- gsub(".stan$", paste0(postfix, ".stan"), file)
  md5 <- md5sum(file)
  ## in case the file already exists, check if there is a MD5 hint
  if(file.exists(generated_file)) {
    stan_source_gen <- readLines(generated_file)
    md5_line <- grepl("^//MD5:", stan_source_gen)
    if(sum(md5_line) == 1) {
      md5_last <- gsub("^//MD5:", "", stan_source_gen[md5_line])
      if (md5_last == md5) {
        cat("Using existing file.\n")
        return(generated_file)
      }
    }
  }
  stan_model <- stanc_builder(file, ...)
  cat(stan_model$model_code, file=generated_file)
  if(enable_stan_next)
    system(paste0("sed -i_orig 's#\\/\\/stan_next:##g' ", generated_file))
  cat(paste0("\n\n//MD5:",md5, "\n\n"), file=generated_file, append=TRUE)
  generated_file
}

inits <- function() {
  #Gives inits for stan in the simulated data case
  theta <- log(c(log(2)/2, log(2)/12, 10))
  return(list(theta=rnorm(3, theta, 0.5),
       xi=matrix(theta[2:3], 2, J),
       omega=rlnorm(2, log(0.1), 1),
       sigma_y=rlnorm(1, log(0.1), 1)))
}

###-----------------------------###
### For visualizing the results ###
###-----------------------------###
limits <- function(samples) { #Samples are of size N times J, where N is the number of repetitions
  # Returns quantiles of the samples
  #print(apply(samples, 2, function(s) quantile(s, c(.022, .158, .5, 0.842, 0.978))))
  return(apply(samples, 2, function(s) quantile(na.omit(s), c(.022, .158, .5, 0.842, 0.978))))
}

plot_comparison <- function(ag, samples, trainag=NULL, trueag=NULL, true=NULL, pl=NULL, points=NULL, points_simple=NULL) {
  # Plots distribution of samples as a function of ag. If true is not null, then it is plotted as a solid line.
  # If points are given, their distributions are plotted. If points_simple are given, they are plotted.
  fnew = limits(samples)
  fnew <- as.data.frame(fnew)
  setattr(fnew, "row.names", c("1","2","CL", "4", "5"))
  eta <- data.frame('ag'=ag, 'id'=t(fnew))
  if(is.null(pl)) {
    pl = ggplot() 
  }
  if(!is.null(points)) {
    ety1 <- data.frame('ag'=trainag, 'elimination'=points)
    pl <- pl + geom_point(data=ety1, aes(x=ag,y=elimination), color="red", alpha=0.5)
  }
  if(!is.null(points_simple)) {
    ety <- data.frame('ag'=trainag, 'elimination'=points_simple)
    pl <- pl + geom_point(data=ety, aes(x=ag,y=elimination), color="blue", alpha=0.5)
  }
  if(is.null(true)) {
    pl = pl +
      geom_line(data=eta, aes(y=id.CL, x=ag), color="blue") + 
      geom_ribbon(data=eta, aes(ymin=id.1,ymax=id.5,x=ag), fill="blue", alpha=0.1) +
      geom_ribbon(data=eta, aes(ymin=id.2,ymax=id.4,x=ag), fill="blue", alpha=0.1)
    
  } else {
    discrepancy <- data.frame(ag=trueag, disc=true)
    pl = pl +
      geom_line(data=eta, aes(y=id.CL, x=ag), color="blue") + #s_legend 
      geom_line(data=discrepancy, aes(x=ag,y=disc), color="red") + #t_legend
      geom_ribbon(data=eta, aes(ymin=id.1,ymax=id.5,x=ag), fill="blue", alpha=0.1) +
      geom_ribbon(data=eta, aes(ymin=id.2,ymax=id.4,x=ag), fill="blue", alpha=0.1)
  }
  return(pl)
}

