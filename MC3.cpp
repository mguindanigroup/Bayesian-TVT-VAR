#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat sample_alpha1(int H, int N, int T, int P, arma::mat obs, arma::mat vectX_Khatri,
                        arma::mat gamma, arma::mat old_alpha1, arma::mat W1, double tau,
                        arma::vec phi, arma::vec sigma) {
  arma::mat D1(T-P,H);
  arma::mat y1_tilde;
  arma::field<arma::mat> Sigma1(H);
  arma::mat nu1(N,H);
  nu1.zeros();
  arma::mat new_alpha1(H,N);
  D1 = gamma.t()%vectX_Khatri;
  for (int h=0;h<H;h++) {
    y1_tilde = obs.rows(P,T-1);
    if (H>1) {
      if (h==0) {y1_tilde -= D1.cols(1,H-1)*old_alpha1.rows(1,H-1);}
      else if (h==(H-1)) {y1_tilde -= D1.cols(0,H-2)*new_alpha1.rows(0,H-2);}
      else {y1_tilde -= D1.cols(0,h-1)*new_alpha1.rows(0,h-1)+D1.cols(h+1,H-1)*old_alpha1.rows(h+1,H-1);}
    }
    Sigma1(h) = arma::inv_sympd(arma::diagmat(1/W1.row(h))/(tau*phi(h))+as_scalar(D1.col(h).t()*D1.col(h))*arma::diagmat(1/sigma));
    nu1.col(h) = Sigma1(h)*arma::diagmat(1/sigma)*y1_tilde.t()*D1.col(h);
    new_alpha1.row(h) = arma::mvnrnd(nu1.col(h),Sigma1(h)).t();
  }
  return(new_alpha1);
}

// [[Rcpp::export]]
arma::mat sample_alpha2(int H, int N, int T, int P, arma::mat obs, arma::mat gamma,
                        arma::mat alpha1, arma::mat old_alpha2, arma::mat alpha3,
                        arma::mat W2, double tau, arma::vec phi, arma::vec sigma) {
  arma::field<arma::mat> D2(T-P,H);
  arma::mat y2_tilde;
  arma::field<arma::mat> Sigma2(H);
  arma::mat nu2(N,H);
  nu2.zeros();
  arma::mat new_alpha2(H,N);
  arma::mat D_2;
  for (int h=0;h<H;h++) {
    Sigma2(h) = arma::mat(N,N).zeros();
    D_2 = alpha1.row(h).t()*alpha3.row(h);
    for (int t=P;t<T;t++) {
      D2(t-P,h) = arma::as_scalar(gamma(h,t-P))*D_2*arma::reverse(obs.rows(t-P,t-1));
      Sigma2(h) += D2(t-P,h).t()*arma::diagmat(1/sigma)*D2(t-P,h);
    }
    Sigma2(h) = arma::inv_sympd(arma::diagmat(1/W2.row(h))/(tau*phi(h))+Sigma2(h));
  }
  for (int h=0;h<H;h++) {
    y2_tilde = obs.rows(P,T-1);
    if (H>1) {
      if (h==0) {
        for (int t=P;t<T;t++) {
          for (int j=1;j<H;j++) {
            y2_tilde.row(t-P) -= old_alpha2.row(j)*D2(t-P,j).t();
          }
        }
      } else if (h==(H-1)) {
        for (int t=P;t<T;t++) {
          for (int j=0;j<H-1;j++) {
            y2_tilde.row(t-P) -= new_alpha2.row(j)*D2(t-P,j).t();
          }
        }
      } else {
        for (int t=P;t<T;t++) {
          for (int j=0;j<h;j++) {
            y2_tilde.row(t-P) -= new_alpha2.row(j)*D2(t-P,j).t();
          }
          for (int j=h+1;j<H;j++) {
            y2_tilde.row(t-P) -= old_alpha2.row(j)*D2(t-P,j).t();
          }
        }
      }
    }
    for (int t=P;t<T;t++) {
      nu2.col(h) += D2(t-P,h).t()*arma::diagmat(1/sigma)*y2_tilde.row(t-P).t();
    }
    nu2.col(h) = Sigma2(h)*nu2.col(h);
    new_alpha2.row(h) = arma::mvnrnd(nu2.col(h),Sigma2(h)).t();
  }
  return(new_alpha2);
}

// [[Rcpp::export]]
arma::mat sample_alpha3(int H, int N, int T, int P, arma::mat obs, arma::mat gamma, arma::mat alpha1, arma::mat alpha2,
                        arma::mat old_alpha3, arma::mat W3, double tau, arma::vec phi, arma::vec sigma) {
  arma::field<arma::mat> D3(T-P,H);
  arma::mat y3_tilde;
  arma::field<arma::mat> Sigma3(H);
  arma::mat nu3(P,H);
  nu3.zeros();
  arma::mat new_alpha3(H,P);
  for (int h=0;h<H;h++) {
    Sigma3(h) = arma::mat(P,P).zeros();
    for (int t=P;t<T;t++) {
      D3(t-P,h) = arma::as_scalar(gamma(h,t-P))*alpha1.row(h).t()*alpha2.row(h)*arma::reverse(obs.rows(t-P,t-1)).t();
      Sigma3(h) += D3(t-P,h).t()*arma::diagmat(1/sigma)*D3(t-P,h);
    }
    Sigma3(h) = arma::inv_sympd(arma::diagmat(1/W3.row(h))/(tau*phi(h))+Sigma3(h));
  }
  for (int h=0;h<H;h++) {
    y3_tilde = obs.rows(P,T-1);
    if (H>1) {
      if (h==0) {
        for (int t=P;t<T;t++) {
          for (int j=1;j<H;j++) {
            y3_tilde.row(t-P) -= old_alpha3.row(j)*D3(t-P,j).t();
          }
        }
      } else if (h==(H-1)) {
        for (int t=P;t<T;t++) {
          for (int j=0;j<H-1;j++) {
            y3_tilde.row(t-P) -= new_alpha3.row(j)*D3(t-P,j).t();
          }
        }
      } else {
        for (int t=P;t<T;t++) {
          for (int j=0;j<h;j++) {
            y3_tilde.row(t-P) -= new_alpha3.row(j)*D3(t-P,j).t();
          }
          for (int j=h+1;j<H;j++) {
            y3_tilde.row(t-P) -= old_alpha3.row(j)*D3(t-P,j).t();
          }
        }
      }
    }
    for (int t=P;t<T;t++) {
      nu3.col(h) += D3(t-P,h).t()*arma::diagmat(1/sigma)*y3_tilde.row(t-P).t();
    }
    nu3.col(h) = Sigma3(h)*nu3.col(h);
    new_alpha3.row(h) = arma::mvnrnd(nu3.col(h),Sigma3(h)).t();
  }
  return(new_alpha3);
}

double conditional_likelihood(int T, int P, arma::vec theta, arma::vec kappa, arma::vec z, int n) {
  double s;
  if (n == 0) {s = theta(n)+kappa(n)*z(n+1);}
  else if (n == T-P-1) {s = theta(n)+kappa(n-1)*z(n-1);}
  else {s = theta(n)+kappa(n)*z(n+1)+kappa(n-1)*z(n-1);}
  return(exp(z(n)*s)/(exp(s)+1));
}

// [[Rcpp::export]]
double psudo_likelihood(int T, int P, arma::vec theta, arma::vec kappa, arma::vec z) {
  double s=0;
  for (int t=0;t<(T-P);t++) {
    s -= log(conditional_likelihood(T,P,theta,kappa,z,t));
  }
  return(s);
}

// [[Rcpp::export]]
arma::vec NDARMA1(int T, int P, double p1, double p2) {
  IntegerVector choice_set(2);
  choice_set(0) = 0;
  choice_set(1) = 1;
  arma::vec prob_p1(2);
  prob_p1(0) = 1-p1;
  prob_p1(1) = p1;
  arma::vec prob_p2(2);
  prob_p2(0) = 1-p2;
  prob_p2(1) = p2;
  arma::vec NDARMA_samples(T-P);
  int a_t;
  NDARMA_samples(0) = RcppArmadillo::sample(choice_set, 1, true, prob_p2)(0);
  for (int ii=1;ii<(T-P);ii++) {
    a_t = RcppArmadillo::sample(choice_set, 1, true, prob_p1)(0);
    NDARMA_samples(ii) = a_t*NDARMA_samples(ii-1) + (1-a_t)*RcppArmadillo::sample(choice_set, 1, true, prob_p2)(0);
  }
  return(NDARMA_samples);
}

arma::vec NDARMA2(int T, int P, arma::vec p, arma::vec q) {
  IntegerVector choice_set(2);
  choice_set(0) = 0;
  choice_set(1) = 1;
  arma::vec prob_p(2);
  arma::vec prob_q(2);
  arma::vec NDARMA_samples(T-P);
  int a_t;
  prob_q(0) = 1-q(0);
  prob_q(1) = q(0);
  NDARMA_samples(0) = RcppArmadillo::sample(choice_set, 1, true, prob_q)(0);
  for (int ii=1;ii<(T-P);ii++) {
    prob_p(0) = 1-p(ii-1);
    prob_p(1) = p(ii-1);
    a_t = RcppArmadillo::sample(choice_set, 1, true, prob_p)(0);
    prob_q(0) = 1-q(ii);
    prob_q(1) = q(ii);
    NDARMA_samples(ii) = a_t*NDARMA_samples(ii-1) + (1-a_t)*RcppArmadillo::sample(choice_set, 1, true, prob_q)(0);
  }
  return(NDARMA_samples);
}

double fun_p1(double theta, double kappa) {
  return(exp(theta)*(exp(kappa)-1)/((exp(theta+kappa)+1)*(exp(theta)+1)));
}

double fun_p2(double theta, double kappa) {
  return(exp(theta)*(exp(theta+kappa)+1)/(exp(2*theta+kappa)+2*exp(theta)+1));
}

arma::mat fun_pq(arma::vec theta, arma::vec kappa) {
  int length = theta.n_elem;
  arma::mat pq(length, 2);
  double a = exp(theta(length-1));
  double b = exp(kappa(length-2));
  if (a > 1e+150) {
    pq(length-2,0) = 0;
    pq(length-1,1) = 1;
  } else {
    pq(length-2,0) = a*(b-1)/(a+1)/(a*b+1);
    pq(length-1,1) = 1-(a+1)/(a*a*b+2*a+1);
  }
  for (int ii=length-2;ii>=1;ii--) {
    if (a > 1e+150) {
      a = exp(theta(ii))*b;
    } else {
      a = exp(theta(ii))*(1+exp(log(pq(ii,0))-log(1-pq(ii,0))-log(a+1)+log(a*a*b+2*a+1)));
    }
    b = exp(kappa(ii-1));
    if (a > 1e+150) {
      pq(ii-1,0) = 0;
      pq(ii,1) = 1;
    } else {
      pq(ii-1,0) = a*(b-1)/(a+1)/(a*b+1);
      pq(ii,1) = 1-(a+1)/(a*a*b+2*a+1);
    }
  }
  if (a > 1e+150) {
    a = exp(theta(0))*b;
  } else {
    a = exp(theta(0))*(1+exp(log(pq(0,0))-log(1-pq(0,0))-log(a+1)+log(a*a*b+2*a+1)));
  }
  pq(0,1) = 1-1/(a+1);
  return(pq);
}

// [[Rcpp::export]]
arma::mat sample_gamma(int H, int N, int T, int P, arma::mat obs, arma::mat old_gamma,
                       arma::cube y_hat, arma::mat theta, arma::mat kappa,
                       arma::vec sigma) {
  arma::mat new_gamma(H, T-P);
  arma::mat y_tilde=obs.rows(P,T-1);
  arma::vec psi(T-P);
  arma::mat pq(T-P,2);
  arma::mat theta_kappa(T-P,2);
  for (int h=1;h<H;h++) {
    y_tilde -= arma::diagmat(old_gamma.row(h))*y_hat.slice(h);
  }
  for (int h=0;h<H;h++) {
    for (int t=0;t<T-P;t++) {
      psi(t) = arma::as_scalar(y_hat.slice(h).row(t)*arma::diagmat(1/sigma)*(y_tilde.row(t).t()-0.5*y_hat.slice(h).row(t).t()));
    }
    pq = fun_pq(psi+theta.row(h).t(), kappa.row(h).t());
    new_gamma.row(h) = NDARMA2(T, P, pq.col(0).rows(0,T-P-2), pq.col(1)).t();
    if (h<H-1) {
      y_tilde += arma::diagmat(old_gamma.row(h+1))*y_hat.slice(h+1)-arma::diagmat(new_gamma.row(h))*y_hat.slice(h);
    }
  }
  return(new_gamma);
}

// [[Rcpp::export]]
List sample_theta(int H, int N, int T, int P, arma::mat gamma, arma::vec old_theta, arma::vec old_kappa,
                  arma::mat mple, double theta_upper, double theta_lower, double kappa_upper,
                  Nullable<double> a_p1, Nullable<double> b_p1, Nullable<double> a_p2, Nullable<double> b_p2,
                  arma::mat old_auxiliary) {
  double sd = 1;
  double prop_theta;
  double prop_kappa;
  double gamma_sum;
  double gamma_diffsum;
  double old_auxiliary_sum;
  double old_auxiliary_diffsum;
  double prop_auxiliary_sum;
  double prop_auxiliary_diffsum;
  double old_theta_star;
  double prop_theta_star;
  double mple_theta_star;
  arma::vec prop_auxiliary;
  double MH_H;
  bool p1_bool = a_p1.isNull() or b_p1.isNull();
  bool p2_bool = a_p2.isNull() or b_p2.isNull();
  arma::vec new_theta(H);
  arma::vec new_kappa(H);
  arma::mat new_auxiliary(H, T-P);
  for (int h=0;h<H;h++) {
    prop_theta = arma::randn()*sd+old_theta(h);
    prop_kappa = arma::randn()*sd+old_kappa(h);
    if ((theta_upper > prop_theta) and (theta_lower < prop_theta) and (0 < prop_kappa) and (prop_kappa < kappa_upper)) {
      prop_auxiliary = NDARMA1(T, P, fun_p1(prop_theta, prop_kappa), fun_p2(prop_theta, prop_kappa));
      gamma_sum = sum(gamma.row(h).cols(1,T-P-2));
      old_auxiliary_sum = sum(old_auxiliary.row(h).cols(1,T-P-2));
      prop_auxiliary_sum = sum(prop_auxiliary.rows(1,T-P-2));
      gamma_diffsum = sum(gamma.row(h).cols(0,T-P-2) % gamma.row(h).cols(1,T-P-1));
      old_auxiliary_diffsum = sum(old_auxiliary.row(h).cols(0,T-P-2) % old_auxiliary.row(h).cols(1,T-P-1));
      prop_auxiliary_diffsum = sum(prop_auxiliary.rows(0,T-P-2) % prop_auxiliary.rows(1,T-P-1));
      old_theta_star = log(exp(old_theta(h))*(exp(old_theta(h))+1)/(exp(old_theta(h)+old_kappa(h))+1));
      prop_theta_star = log(exp(prop_theta)*(exp(prop_theta)+1)/(exp(prop_theta+prop_kappa)+1));
      mple_theta_star = log(exp(mple(h,0))*(exp(mple(h,0))+1)/(exp(mple(h,0)+mple(h,1))+1));
      MH_H = arma::as_scalar(exp(
        (prop_theta-old_theta(h))*(gamma.row(h).col(0)+gamma.row(h).col(T-P-1))+
          (prop_theta_star-old_theta_star)*gamma_sum+(prop_kappa-old_kappa(h))*gamma_diffsum+
          (mple(h,0)-prop_theta)*(prop_auxiliary.row(0)+prop_auxiliary.row(T-P-1))+
          (mple_theta_star-prop_theta_star)*prop_auxiliary_sum+(mple(h,1)-prop_kappa)*prop_auxiliary_diffsum+
          (old_theta(h)-mple(h,0))*(old_auxiliary.row(h).col(0)+old_auxiliary.row(h).col(T-P-1))+
          (old_theta_star-mple_theta_star)*old_auxiliary_sum+(old_kappa(h)-mple(h,1))*old_auxiliary_diffsum));
      if (!p1_bool) {
        double a_p1(a_p1);
        double b_p1(b_p1);
        double prop_p1;
        double old_p1;
        prop_p1 = fun_p1(prop_theta, prop_kappa);
        old_p1 = fun_p1(old_theta(h), old_kappa(h));
        MH_H += arma::as_scalar(exp(
          (a_p1-1)*(log(prop_p1)-log(old_p1))+(b_p1-1)*(log(1-prop_p1)-log(1-old_p1))));
      }
      if (!p2_bool) {
        double a_p2(a_p2);
        double b_p2(b_p2);
        double prop_p2;
        double old_p2;
        prop_p2 = fun_p2(prop_theta, prop_kappa);
        old_p2 = fun_p2(old_theta(h), old_kappa(h));
        MH_H += arma::as_scalar(exp(
          (a_p2-1)*(log(prop_p2)-log(old_p2))+(b_p2-1)*(log(1-prop_p2)-log(1-old_p2))));
      }
    } else {
      MH_H = 0;
    }
    if (arma::randu()<MH_H) {
      new_theta(h) = prop_theta;
      new_kappa(h) = prop_kappa;
      new_auxiliary.row(h) = prop_auxiliary.t();
    } else {
      new_theta(h) = old_theta(h);
      new_kappa(h) = old_kappa(h);
      new_auxiliary.row(h) = old_auxiliary.row(h);
    }
  }
  return(List::create(Named("theta") = new_theta, Named("kappa") = new_kappa, Named("auxiliary") = new_auxiliary));
}

// [[Rcpp::export]]
arma::cube dynamic_coefficients(int T, int N, int P, int H, arma::cube gamma, arma::cube alpha1, arma::cube alpha2, arma::cube alpha3) {
  int iter = gamma.n_rows;
  arma::field<arma::mat> coef(iter+1,H);
  arma::vec alpha1_vec;
  arma::vec alpha2_vec;
  arma::vec alpha3_vec;
  for (int it=0;it<iter;it++) {
    for (int h=0;h<H;h++) {
      alpha1_vec = alpha1.tube(it,h);
      alpha2_vec = alpha2.tube(it,h);
      alpha3_vec = alpha3.tube(it,h);
      coef(it,h) = arma::kron(alpha3_vec.t(), alpha1_vec*alpha2_vec.t());
    }
  }
  int length = gamma.n_slices;
  arma::cube A_time(length,N,N*P);
  A_time.fill(0);
  for (int t=0;t<length;t++) {
    for (int h=0;h<H;h++) {
      for (int it=0;it<iter;it++) {
        A_time.row(t) += gamma.slice(t)(it,h)*coef(it,h);
      }
    }
    A_time.row(t) = A_time.row(t)/iter;
  }
  return(A_time);
}

// [[Rcpp::export]]
List credible_interval(arma::cube alpha1, arma::cube alpha2, arma::cube alpha3, arma::cube gamma, double qt) {
  int iter = arma::size(alpha1)(0);
  int N = arma::size(alpha1)(2);
  int P = arma::size(alpha3)(2);
  int Time = arma::size(gamma)(2) + P;
  arma::vec quantile_vec = {qt/2, 1-qt/2};
  arma::cube interval_upper(Time-P, N, N*P);
  arma::cube interval_lower(Time-P, N, N*P);
  arma::vec tmp(iter);
  arma::vec quantile_tmp(2); 
  for (int tt=0;tt<Time-P;tt++) {
    for (int n1=0;n1<N;n1++) {
      for (int n2=0;n2<N;n2++) {
        for (int pp=0;pp<P;pp++) {
          for (int ii=0;ii<iter;ii++) {
            tmp(ii) = arma::accu(alpha1.slice(n1).row(ii) % alpha2.slice(n2).row(ii) % alpha3.slice(pp).row(ii));
          }
          quantile_tmp = arma::quantile(tmp, quantile_vec);
          interval_lower(tt,n1,pp*N+n2) = quantile_tmp(0);
          interval_upper(tt,n1,pp*N+n2) = quantile_tmp(1);
        }
      }
    }
  }
  return(List::create(Named("lower") = interval_lower, Named("upper") = interval_upper));
}
