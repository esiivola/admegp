// NM data block declaration
int<lower = 1> N; // number of lines of nm data set
int<lower=1, upper=N> id[N];
vector<lower=0>[N] time;
vector<lower=0>[N] amt;
int<lower=0> cmt[N];
int<lower=0, upper=1> mdv[N];
int<lower=0, upper=2> evid[N];
int<lower=0> addl[N];
vector<lower=0>[N] tau;
