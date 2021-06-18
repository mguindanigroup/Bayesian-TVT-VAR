#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat sample_alpha1(int H, int N, int T, int P, arma::mat obs, arma::mat vectX_Khatri, arma::mat gamma, arma::mat old_alpha1, arma::mat alpha2, 
                        arma::mat alpha3, arma::mat W1, double tau, arma::vec phi, arma::vec sigma) {
  arma::mat D1(T-P,H);
  arma::field<arma::mat> y1_tilde(H);
  arma::field<arma::mat> Sigma1(H);
  arma::mat nu1(N,H);
  nu1.zeros();
  arma::mat new_alpha1(H,N);
  D1 = gamma.t()%vectX_Khatri;
  for (int h=0;h<H;h++) {
    y1_tilde(h) = obs.rows(P,T-1);
    if (h==0) {y1_tilde(h) -= D1.cols(1,H-1)*old_alpha1.rows(1,H-1);}
    else if (h==(H-1)) {y1_tilde(h) -= D1.cols(0,H-2)*new_alpha1.rows(0,H-2);}
    else {y1_tilde(h) -= D1.cols(0,h-1)*new_alpha1.rows(0,h-1)+D1.cols(h+1,H-1)*old_alpha1.rows(h+1,H-1);}
    Sigma1(h) = arma::inv_sympd(arma::diagmat(1/W1.row(h))/(tau*phi(h))+as_scalar(D1.col(h).t()*D1.col(h))*arma::diagmat(1/sigma));
    nu1.col(h) = Sigma1(h)*arma::diagmat(1/sigma)*y1_tilde(h).t()*D1.col(h);
    new_alpha1.row(h) = arma::mvnrnd(nu1.col(h),Sigma1(h)).t();
  }
  return(new_alpha1);
}

// [[Rcpp::export]]
arma::mat sample_alpha2(int H, int N, int T, int P, arma::mat obs, arma::mat gamma, arma::mat alpha1, arma::mat old_alpha2, 
                        arma::mat alpha3, arma::mat W2, double tau, arma::vec phi, arma::vec sigma) {
  arma::field<arma::mat> D2(T-P,H);
  arma::field<arma::mat> y2_tilde(H);
  arma::field<arma::mat> Sigma2(H);
  arma::mat nu2(N,H);
  nu2.zeros();
  arma::mat new_alpha2(H,N);
  arma::mat D_2(N,H);
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
    y2_tilde(h) = obs.rows(P,T-1);
    if (h==0) {
      for (int t=P;t<T;t++) {
        for (int j=1;j<H;j++) {
          y2_tilde(h).row(t-P) -= old_alpha2.row(j)*D2(t-P,j).t(); 
        }
      }
    }
    else if (h==(H-1)) {
      for (int t=P;t<T;t++) {
        for (int j=0;j<H-1;j++) {
          y2_tilde(h).row(t-P) -= new_alpha2.row(j)*D2(t-P,j).t(); 
        }
      }
    }
    else {
      for (int t=P;t<T;t++) {
        for (int j=0;j<h;j++) {
          y2_tilde(h).row(t-P) -= new_alpha2.row(j)*D2(t-P,j).t();
        }
        for (int j=h+1;j<H;j++) {
          y2_tilde(h).row(t-P) -= old_alpha2.row(j)*D2(t-P,j).t();
        }
      }
    }
    for (int t=P;t<T;t++) {
      nu2.col(h) += D2(t-P,h).t()*arma::diagmat(1/sigma)*y2_tilde(h).row(t-P).t();
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
  arma::field<arma::mat> y3_tilde(H);
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
    y3_tilde(h) = obs.rows(P,T-1);
    if (h==0) {
      for (int t=P;t<T;t++) {
        for (int j=1;j<H;j++) {
          y3_tilde(h).row(t-P) -= old_alpha3.row(j)*D3(t-P,j).t(); 
        }
      }
    }
    else if (h==(H-1)) {
      for (int t=P;t<T;t++) {
        for (int j=0;j<H-1;j++) {
          y3_tilde(h).row(t-P) -= new_alpha3.row(j)*D3(t-P,j).t(); 
        }
      }
    }
    else {
      for (int t=P;t<T;t++) {
        for (int j=0;j<h;j++) {
          y3_tilde(h).row(t-P) -= new_alpha3.row(j)*D3(t-P,j).t();
        }
        for (int j=h+1;j<H;j++) {
          y3_tilde(h).row(t-P) -= old_alpha3.row(j)*D3(t-P,j).t();
        }
      }
    }
    for (int t=P;t<T;t++) {
      nu3.col(h) += D3(t-P,h).t()*arma::diagmat(1/sigma)*y3_tilde(h).row(t-P).t();
    }
    nu3.col(h) = Sigma3(h)*nu3.col(h);
    new_alpha3.row(h) = arma::mvnrnd(nu3.col(h),Sigma3(h)).t();
  }
  return(new_alpha3);
}

// [[Rcpp::export]]
double likelihood(int T, int P, arma::vec theta, arma::vec int_theta, arma::vec z) {
  double s=0;
  for (int t=0;t<(T-P);t++) {
    if (t<(T-P-1)) {
      s += theta(t)*z(t)+int_theta(t)*z(t)*z(t+1);
    } else {
      s += theta(t)*z(t);
    } 
  }
  return(exp(s));
}

double conditional_likelihood(int T, int P, arma::vec theta, arma::vec int_theta, arma::vec z, int n) {
  double s;
  if (n == 0) {s = theta(n)+int_theta(n)*z(n+1);}
  else if (n == T-P-1) {s = theta(n)+int_theta(n-1)*z(n-1);}
  else {s = theta(n)+int_theta(n)*z(n+1)+int_theta(n-1)*z(n-1);}
  return(exp(z(n)*s)/(exp(s)+1));
}

// [[Rcpp::export]]
double psudo_likelihood(int T, int P, arma::vec theta, arma::vec int_theta, arma::vec z) {
  double s=0;
  for (int t=0;t<(T-P);t++) {
    s -= log(conditional_likelihood(T,P,theta,int_theta,z,t));
  }
  return(s);
}

double full_conditional(int T, int P, arma::vec theta, arma::vec int_theta, arma::vec z, int n) {
  double s;
  if (n == 0) {s = theta(n)+int_theta(n)*z(n+1);}
  else if (n == T-P-1) {s = theta(n)+int_theta(n-1)*z(n-1);}
  else {s = theta(n)+int_theta(n)*z(n+1)+int_theta(n-1)*z(n-1);}
  return(exp(s)/(exp(s)+1));
}

// [[Rcpp::export]]
arma::vec exact_sampling(int T, int P, arma::vec theta, arma::vec int_theta) {
  int TT = 1;
  arma::vec upper(T-P);
  arma::vec lower(T-P);
  arma::Col<int> u_node(1048576);
  arma::vec u_prob(1048576);
  /* double u; */
  u_node(0) = arma::as_scalar(arma::randi(arma::distr_param(0,T-P-1)));
  u_prob(0) = arma::randu();
  while (TRUE) {
    upper.fill(1);
    lower.fill(0);
    for (int t=-TT;t<0;t++) {
      if (u_prob(-t-1) <= full_conditional(T, P, theta, int_theta, upper, u_node(-t-1))) upper(u_node(-t-1)) = 1; else upper(u_node(-t-1)) = 0; 
      if (u_prob(-t-1) <= full_conditional(T, P, theta, int_theta, lower, u_node(-t-1))) lower(u_node(-t-1)) = 1; else lower(u_node(-t-1)) = 0;
    }
    if (arma::approx_equal(upper, lower, "absdiff", 0.1)) break;
    u_node.rows(TT,2*TT-1) = arma::randi(TT, arma::distr_param(0,T-P-1));
    u_prob.rows(TT,2*TT-1) = arma::randu(TT);
    TT = TT*2;
  }
  return(upper);
}

// [[Rcpp::export]]
arma::mat sample_gamma(int H, int N, int T, int P, arma::mat obs, arma::mat old_gamma, arma::cube y_hat, arma::mat theta, arma::mat int_theta, arma::vec sigma) {
  arma::mat new_gamma(H, T-P);
  arma::mat y_tilde=obs.rows(P,T-1);
  arma::vec psi(T-P);
  for (int h=1;h<H;h++) {
    y_tilde -= arma::diagmat(old_gamma.row(h))*y_hat.slice(h);
  }
  for (int h=0;h<H;h++) {
    for (int t=0;t<T-P;t++) {
      psi(t) = arma::as_scalar(y_hat.slice(h).row(t)*arma::diagmat(1/sigma)*(y_tilde.row(t).t()-0.5*y_hat.slice(h).row(t).t()));
    }
    new_gamma.row(h) = exact_sampling(T, P, psi+theta.row(h).t(), int_theta.row(h).t()).t();
    if (h<H-1) {
      y_tilde += arma::diagmat(old_gamma.row(h+1))*y_hat.slice(h+1)-arma::diagmat(new_gamma.row(h))*y_hat.slice(h);
    }
  }
  return(new_gamma);
}

arma::vec theta_star(int T, int P, double theta, double int_theta) {
  arma::vec theta_output(T-P);
  theta_output.fill(log(exp(theta)*(exp(theta)+1)/(exp(theta+int_theta)+1)));
  theta_output(0) = theta;
  theta_output(T-P-1) = theta;
  return(theta_output);
}

// [[Rcpp::export]]
List sample_theta(int H, int N, int T, int P, arma::mat gamma, arma::vec old_theta, arma::vec old_int_theta, 
                  arma::mat mple, double theta_upper, double theta_lower, double int_theta_upper, arma::mat old_auxiliary) {
  double sd = 1;
  double prop_theta;
  double prop_int_theta;
  arma::vec prop_auxiliary;
  double MH_H;
  arma::vec new_theta(H);
  arma::vec new_int_theta(H);
  arma::mat new_auxiliary(H, T-P);
  for (int h=0;h<H;h++) {
    prop_theta = arma::randn()*sd+old_theta(h);
    prop_int_theta = arma::randn()*sd+old_int_theta(h);
    prop_auxiliary = exact_sampling(T, P, theta_star(T,P,prop_theta,prop_int_theta), arma::vec(T-P-1).fill(prop_int_theta));
    if ((theta_upper > prop_theta) and (theta_lower < prop_theta) and (0 < prop_int_theta) and (prop_int_theta < int_theta_upper)) {
      MH_H = likelihood(T, P, theta_star(T,P,prop_theta,prop_int_theta), arma::vec(T-P-1).fill(prop_int_theta), gamma.row(h).t())*
        likelihood(T, P, theta_star(T,P,mple(h,0),mple(h,1)), arma::vec(T-P-1).fill(mple(h,1)), prop_auxiliary)*
        likelihood(T, P, theta_star(T,P,old_theta(h),old_int_theta(h)), arma::vec(T-P-1).fill(old_int_theta(h)), old_auxiliary.row(h).t())/
          (likelihood(T, P, theta_star(T,P,old_theta(h),old_int_theta(h)), arma::vec(T-P-1).fill(old_int_theta(h)), gamma.row(h).t())*
            likelihood(T, P, theta_star(T,P,mple(h,0),mple(h,1)), arma::vec(T-P-1).fill(mple(h,1)), old_auxiliary.row(h).t())*
            likelihood(T, P, theta_star(T,P,prop_theta,prop_int_theta), arma::vec(T-P-1).fill(prop_int_theta), prop_auxiliary));
    } else {
      MH_H = 0;
    }
    if (arma::randu()<MH_H) {
      new_theta(h) = prop_theta;
      new_int_theta(h) = prop_int_theta;
      new_auxiliary.row(h) = prop_auxiliary.t();
    } else {
      new_theta(h) = old_theta(h);
      new_int_theta(h) = old_int_theta(h);
      new_auxiliary.row(h) = old_auxiliary.row(h);
    }
  }
  return(List::create(Named("theta") = new_theta, Named("int_theta") = new_int_theta, Named("auxiliary") = new_auxiliary));
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


