library(ggplot2)
library(StanHeaders)
library(rstan)
library(reshape2)
theme_set(theme_bw())
source("utils.R", echo = FALSE)

##-------------------------------------------------------##
## Select data generation method and maturation function ##
##-------------------------------------------------------##
data = "hill" # 1="adult", 2="hill", 3="double-hill"
method = "parametric" # 1="linear", 2="parametric", 3="gp", 4="gp2"

##----------------##
## Generate files ##
##----------------##
if(data == 1 || data == "adult"){
  print("Using no maturation in data")
  data_title  = "No maturation in data"
  maturation = F
  data = "adult"
} else if(data == 2 || data == "hill") {
  print("Using hill ín data")
  data_title  = "Hill function in data"
  maturation = T
  double_hill = F
  data = "hill"
} else if(method == 3 || method == "double-hill") {
  print("Using double-hill in data")
  data_title  = "Double-hill function in data"
  maturation = T
  double_hill = T
  data = "double-hill"
}

gen_data_files(J=20, num_measurements=10, seed=12345, maturation=maturation, plotit=TRUE, double_hill=double_hill, deriv=F) 

##---------------------##
## Generate STAN model ##
##---------------------##
deriv=F
#Use one of the following lines:
if(method == 1 || method == "linear"){
  print("Using linear model for fitting")
  method_title = " and linear model"
  maturation_type=0
  method = "linear"
} else if(method == 2 || method == "parametric") {
  print("Using parametric model for fitting")
  method_title = " and parametric model"
  maturation_type=1  
  method = "parametric"
} else if(method == 3 || method == "gp") {
  print("Using naive gp model for fitting")
  method_title = " and naive GP model"
  maturation_type=2
  deriv=F
  method = "gp"
} else if(method == 4 || method == "gp2") {
  print("Using gp model for fitting")
  method_title = " and GP model"
  maturation_type=2
  deriv=T
  method = "gp2"
}

plot_title = paste(data_title, method_title, sep='')

stan_pk_model = stan_model(file=paste('./', stan_generate_model(file="stan/oral_1cmt_omni.stan",enable_stan_next=T), sep=''))

rstan_options(auto_write=TRUE)
options(mc.cores = 4)
set.seed(567583)

##-------------------------------------------##
## Load training data and run the Stan model ##
##-------------------------------------------##
source("generated_files/oral_1cmt_allo_run.data.R", echo=FALSE)

J <- max(id)
# The next line is SLOW (e.g ~1h)
fit <- sampling(stan_pk_model, chains=4, warmup=1500, iter=2000, init=inits, refresh=50, control = list(max_treedepth=10, adapt_delta = 0.9999)) #

##----------------------------------------##
## Extract MCMC samples from the Stan fit ##
##----------------------------------------##
samples <- extract(fit)


##-------------------------------------------##
## Visualize the detected maturation trend   ##
##-------------------------------------------##
trueag = seq(from=0.7, to=5.7, length.out=100)
pl = plot_comparison(aget, exp(samples$ft), trainag=ag, trueag=trueag, true=mat(trueag, maturation=maturation, doublehill=double_hill), points_simple= exp(apply(samples$y_cl,2,function(x) median(na.omit(x)))), points=exp(ke - lin_approx(get_lcl75(), get_alpha(), wt) + log(mat(ag, maturation=maturation,doublehill=double_hill))))
pl = pl  + xlab("Age [y]") + ylab("Maturation")   + #+ ggtitle("Real and estimated maturations")scale_y_continuous(breaks=c(-4, -3,-2,-1,0), labels=c("exp(-4)", "exp(-3)", "exp(-2)", "exp(-1)", "1"), limits=c(-3.9,0.5)) +
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1,1.25), labels=c(NULL,NULL,NULL,NULL,NULL,NULL)) +
  scale_x_continuous( limits=c(0.69,5.8),breaks=c(0.7,1.7,2.7,3.7,4.7,5.7), labels=c(0,1,2,3,4,5))+
  coord_trans(x = "log10", y = "identity", clip="on", limy=limitsy, limx=c(0.69,5.8)) +
  theme(legend.key = element_blank(), axis.ticks.y = element_blank(),
        legend.title = element_blank(), legend.position=c(0.82, 0.4), text = element_text(size=7), axis.text.x=element_text(size=6), axis.text.y=element_text(size=6), plot.title=element_text(size=7))  
pl
fr <- data.frame('time'=ageds, 'y'=vdy)
pl = pl + geom_point(data=fr, aes(y=y, x=time), shape=22, color="grey", fill="grey")
pl <- add_legend(pl,textsize=6, keysize = 0.035, ind=c(1,4,2,3,7), labels=c('True individual maturations','True maturation trend', 'Estimated individual maturations','Estimated maturation trend','Virtual derivative sign observations'), position=c(0.25,0.05), u_shapes=c(16, NA, 16, NA,22), u_linetypes=c('blank','solid','blank','solid','blank'))
pl <- pl + ggtitle(plot_title) +  theme(plot.title = element_text(hjust = 0.01, margin = margin(t = 10, b = -10)),  plot.margin=margin(l=0, r=0))
pl

##-----------------------------------------------------------------------##
## Save the environment for later use to avoid re-running the Stan model ##
##-----------------------------------------------------------------------##
remove(pl) # no need to save the generated plot
dir.create("backups")
save.image(paste("backups/", data, "_", method,"_simulated.RData",sep=''))

