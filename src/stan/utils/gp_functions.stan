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
    Sigma[1:N,1:N] = cov_exp_quad(x,alpha,length_scale) .* cov_lin(x,x,N,N,c,b,v);
    Sigma[N+1:N+Nd, 1:N] = dx1_cov_exp_quad(xd,x, Nd,N,alpha,length_scale) .* cov_lin(xd,x,Nd,N,c,b,v) + cov_exp_quad(xd,x,alpha,length_scale) .* dx1_cov_lin(xd,x,Nd,N,c,b,v);
    Sigma[1:N, N+1:N+Nd] = -1.0*dx1_cov_exp_quad(x,xd,N,Nd,alpha,length_scale) .* cov_lin(x,xd,N,Nd,c,b,v) + cov_exp_quad(x,xd,alpha,length_scale) .* dx2_cov_lin(x,xd,N,Nd,c,b,v);
    Sigma[N+1:N+Nd, N+1:N+Nd] = dx1dx2_cov_exp_quad(xd,xd,Nd,Nd,alpha, length_scale) .* cov_lin(xd,xd,Nd,Nd,c,b,v) + cov_exp_quad(xd,alpha,length_scale) .* dx1dx2_cov_lin(xd,xd,Nd,Nd,c,b,v) + dx1_cov_exp_quad(xd,xd,Nd,Nd,alpha,length_scale) .* dx2_cov_lin(xd,xd,Nd,Nd,c,b,v) -1.0*dx1_cov_exp_quad(xd,xd,Nd,Nd,alpha,length_scale) .* dx1_cov_lin(xd,xd,Nd,Nd,c,b,v);
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
  vector latent_f(real rho, real alpha, real delta, vector eta, real[] x, real[] xd, int n, int nd, real rho_fixed, real c, real b, real v) {
    // can we avoid the rho < 0 if statement here?
    matrix[n+nd, n+nd] K = cov_exp_quad_full(x, xd, n, nd, alpha, rho < 0 ? rho_fixed : rho, c,b,v);
    matrix[n+nd, n+nd] L_K;
    for (j in 1:n+nd)
    	K[j, j] += delta;
    L_K = cholesky_decompose(K);
    return L_K * eta;
  }
