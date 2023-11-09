// Fdepth cpp 
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <iostream>
#include <string>
#include <stdio.h>
#include <time.h>

#define RCPP_ARMADILLO_FIX_Field

using namespace Rcpp;
using namespace arma;
using namespace std;

//[[Rcpp::export]]
vec testtt(arma::vec x, arma::vec p){
  vec res = quantile(x, p);
  return(res);
}

//[[Rcpp::export]]
List test(arma::vec v){
  // uvec stidx = arma::sort_index(arma::sort_index(v)) + 1;
  
  arma::vec a = linspace(0, 1, 100);
  arma::mat b(100, 1);
  b.col(0) = a;
  b.set_size(10, 10);
  List ll=List::create(find(v==min(v))+1, a, b);
  return(ll);
}

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string currentDateTime() {
  time_t     now = time(0);
  struct tm  tstruct;
  char       buf[80];
  tstruct = *localtime(&now);
  // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
  // for more information about date/time format
  strftime(buf, sizeof(buf), "%Y-%m-%d, %X", &tstruct);
  
  return buf;
}

//[[Rcpp::export]]
arma::vec f_nd(arma::vec& tGrid){
  int tt = tGrid.size();
  arma::vec res(tt, fill::zeros);
  for(int idx_t=0; idx_t<tt; idx_t++){
    double t = tGrid(idx_t);
    if(t>=-0.05 & t<0.05){
      res(idx_t) = 0.05*t;
    }else if(t>=0.05 & t<0.15){
      res(idx_t) = -0.05*(t-0.1);
    }else if(t>=0.15 & t<0.25){
      res(idx_t) = 0.05*(t-0.2);
    }else if(t>=0.25 & t<0.35){
      res(idx_t) = -0.05*(t-0.3);
    }else if(t>=0.35 & t<0.45){
      res(idx_t) = 0.05*(t-0.4);
    }else if(t>=0.45 & t<0.55){
      res(idx_t) = -0.05*(t-0.5);
    }else if(t>=0.55 & t<0.65){
      res(idx_t) = 0.05*(t-0.6);
    }else if(t>=0.65 & t<0.75){
      res(idx_t) = -0.05*(t-0.7);
    }else if(t>=0.75 & t<0.85){
      res(idx_t) = 0.05*(t-0.8);
    }else if(t>=0.85 & t<0.95){
      res(idx_t) = -0.05*(t-0.9);
    }else if(t>=0.95 & t<1.05){
      res(idx_t) = 0.05*(t-1);
    }else{
      res(idx_t) = 0;
    }
    
  }
  return(res);
}

//[[Rcpp::export]]
double f_nd_each(double t){
  double f;
  if(t>=-0.05 & t<0.05){
    f = 0.05*t;
  }else if(t>=0.05 & t<0.15){
    f = -0.05*(t-0.1);
  }else if(t>=0.15 & t<0.25){
    f = 0.05*(t-0.2);
  }else if(t>=0.25 & t<0.35){
    f = -0.05*(t-0.3);
  }else if(t>=0.35 & t<0.45){
    f = 0.05*(t-0.4);
  }else if(t>=0.45 & t<0.55){
    f = -0.05*(t-0.5);
  }else if(t>=0.55 & t<0.65){
    f = 0.05*(t-0.6);
  }else if(t>=0.65 & t<0.75){
    f = -0.05*(t-0.7);
  }else if(t>=0.75 & t<0.85){
    f = 0.05*(t-0.8);
  }else if(t>=0.85 & t<0.95){
    f = -0.05*(t-0.9);
  }else if(t>=0.95 & t<1.05){
    f = 0.05*(t-1);
  }else{
    f = 0;
  }
  return(f);
}

//[[Rcpp::export]]
List IQRout(arma::vec& x, double f){
  arma::vec pr = {0.25, 0.75};
  arma::vec q13 = quantile(x, pr);
  double iqr = as_scalar(diff(q13));
  double l = q13(0) - f*iqr;
  double u = q13(1) + f*iqr;
  double med = median(x);
  arma::vec stats = {l, q13(0), med, q13(1), u};
  uvec which = find(x<l || x>u);
  arma::vec out = x(which);
  
  List res = List::create(Named("stats") = stats,
                          Named("IQR") = iqr,
                          Named("out") = out,
                          Named("out.which") = which);
  
  return(res);
}

//[[Rcpp::export]]
List IQRout_vec(arma::vec& x, arma::vec& f_vec){
  arma::vec pr = {0.25, 0.75};
  arma::vec q13 = quantile(x, pr);
  double iqr = as_scalar(diff(q13));
  double med = median(x);
  
  int f_size = f_vec.size();
  arma::mat stat_mat(5, f_size, fill::zeros);
  field<uvec> which_field(f_size);
  field<arma::vec> out_field(f_size);
  
  for(int ii_f=0; ii_f<f_size; ++ii_f){
    double l = q13(0) - f_vec(ii_f)*iqr;
    double u = q13(1) + f_vec(ii_f)*iqr;
    
    stat_mat.col(ii_f) = {l, q13(0), med, q13(1), u};
    
    uvec which = find(x<l || x>u);
    which_field(ii_f) = which;
    out_field(ii_f) = x(which);
  }
  
  
  List res = List::create(Named("IQR") = iqr,
                          Named("stats") = stat_mat,
                          Named("out") = out_field,
                          Named("out.which") = which_field);
  
  return(res);
}

//[[Rcpp::export]]
NumericVector rank_rcpp(arma::vec& x, Rcpp::CharacterVector ties){
  // Obtaining namespace of Matrix package
  Environment pkg = Environment::namespace_env("base");
  
  // Picking up Matrix() function from Matrix package
  Function f = pkg["rank"];
  
  // NumericVector x_nv = NumericVector(x.begin(), x.end());
  
  return f(x, Named("ties.method", ties));
}




// 

// simulation of curves X
//[[Rcpp::export]]
arma::mat Xsim(int n, arma::vec& ev, arma::mat& ef, int xi_type=1, double nu_xi=1){
  int J = ev.size(); int tt = ef.n_cols;
  arma::mat X(n, tt, fill::zeros);
  if(xi_type==0){   // Gaussian processes
    for(int i=0; i<n; i++){
      rowvec xi_vec =  Rcpp::rnorm(J, 0, 1); 
      rowvec X_i(tt, fill::zeros);
      for(int j=0; j<J; j++){
        X_i += sqrt(ev(j))*ef.row(j)*xi_vec(j);
      }
      X.row(i) = X_i;
    }
  }else if(xi_type==1){   // normal
    for(int i=0; i<n; i++){
      rowvec xi_vec = R::rnorm(0, 1) * Rcpp::rnorm(J, 0, 1); 
      rowvec X_i(tt, fill::zeros);
      for(int j=0; j<J; j++){
        X_i += sqrt(ev(j))*ef.row(j)*xi_vec(j);
      }
      X.row(i) = X_i;
    }
  }else if(xi_type==2){   // t distribution
    double sd = sqrt(nu_xi / (nu_xi-2));
    for(int i=0; i<n; i++){
      rowvec xi_vec = (R::rt(nu_xi)/sd) * Rcpp::rnorm(J, 0, 1); 
      rowvec X_i(tt, fill::zeros);
      for(int j=0; j<J; j++){
        X_i += sqrt(ev(j))*ef.row(j)*xi_vec(j);
      }
      X.row(i) = X_i;
    }
  }else if(xi_type==3){   // unif distribution
    double sd = sqrt(3);
    for(int i=0; i<n; i++){
      rowvec xi_vec = R::runif(-sd,sd) * Rcpp::runif(J, -sd, sd); 
      rowvec X_i(tt, fill::zeros);
      for(int j=0; j<J; j++){
        X_i += sqrt(ev(j))*ef.row(j)*xi_vec(j);
      }
      X.row(i) = X_i;
    }
  }
  
  return(X);
}

//[[Rcpp::export]]
arma::mat Xout(int n0, arma::vec& tGrid, double eta, int type){
  int tt = tGrid.size();
  
  arma::vec U = (randu<arma::vec>(n0)+3) / 4;
  arma::vec W = 2*round(randu<arma::vec>(n0))-1;
  arma::vec V = 3*eta*U % W/4;
  
  arma::mat Xout(n0, tt, fill::zeros);
  
  if(type==1){                        // magnitude outlier
    arma::vec temp = 2 * eta * U % W;
    arma::vec vec1(tt, fill::ones);
    Xout = temp * trans(vec1);
  }else if(type==2){                  // jump outlier
    arma::vec Uloc = (8*randu<arma::vec>(n0)+1)/10;
    for(int i0=0; i0<n0; i0++){
      for(int idx_t=0; idx_t<tt; idx_t++){
        double t = tGrid(idx_t);
        if(t>Uloc(i0)){Xout(i0, idx_t) = V(i0);
          }else{Xout(i0, idx_t) = -V(i0);}
      }
    }
  }else if(type==3){                  // peak outlier
    arma::vec Uloc = 9*randu<arma::vec>(n0)/10;
    for(int i0=0; i0<n0; i0++){
      for(int idx_t=0; idx_t<tt; idx_t++){
        double t = tGrid(idx_t);
        if(t>Uloc(i0) & t<Uloc(i0)+0.1){Xout(i0, idx_t) = V(i0);
        }else{Xout(i0, idx_t) = -V(i0);}
      }
    }
  }else if(type==4){                  // wiggle outlier
    arma::vec Uloc = randu<arma::vec>(n0)/10;
    for(int i0=0; i0<n0; i0++){
      for(int idx_t=0; idx_t<tt; idx_t++){
        double t = tGrid(idx_t);
        Xout(i0, idx_t) = V(i0)*sin(10*datum::pi*(t - Uloc(i0)));
      }
    }
  }else if(type==5){                  // non-diff wiggle outlier
    arma::vec Uloc = randu<arma::vec>(n0)/20;
    for(int i0=0; i0<n0; i0++){
      for(int idx_t=0; idx_t<tt; idx_t++){
        double t = tGrid(idx_t);
        Xout(i0, idx_t) = 400*V(i0)*f_nd_each(t - Uloc(i0));
      }
    }
  }else if(type==6){                  // linear outlier
    Xout = V * trans(2*tGrid-1);
  }else if(type==7){                  // quadratic outlier
    Xout = V * trans(8*(tGrid-0.5) % (tGrid-0.5) - 1);
  }else if(type==8){                  // cubic outlier
    Xout = V * trans(12 * sqrt(3) * tGrid % (tGrid-0.5) % (tGrid-1.0));
  }else{
    stop("The input type should be one of integers from 1 to 8.");
  }
  return(Xout);
}

//[[Rcpp::export]]
cube Xout_all(int n0, arma::vec& tGrid, double eta){
  int tt = tGrid.size();
  
  arma::mat U = (randu<arma::mat>(n0,8)+3) / 4;
  arma::mat W = 2*round(randu<arma::mat>(n0,8))-1;
  arma::mat V = 3*eta*U % W/4;
  
  cube Xout(n0, tt, 8, fill::zeros);
  
  // magnitude outlier
  arma::vec temp = 2*4/3*V.col(0);
  arma::vec vec1(tt, fill::ones);
  Xout.slice(0) = temp * trans(vec1);
  
  // jump outlier
  arma::vec Uloc1 = (8*randu<arma::vec>(n0)+1)/10;
  for(int i0=0; i0<n0; i0++){
    for(int idx_t=0; idx_t<tt; idx_t++){
      double t = tGrid(idx_t);
      if(t>Uloc1(i0)){Xout(i0, idx_t,1) = V(i0,1);
      }else{Xout(i0, idx_t,1) = -V(i0,1);}
    }
  }
  
  // peak outlier
  arma::vec Uloc2 = 9*randu<arma::vec>(n0)/10;
  for(int i0=0; i0<n0; i0++){
    for(int idx_t=0; idx_t<tt; idx_t++){
      double t = tGrid(idx_t);
      if(t>Uloc2(i0) & t<Uloc2(i0)+0.1){Xout(i0, idx_t,2) = V(i0,2);
      }else{Xout(i0, idx_t,2) = -V(i0,2);}
    }
  }
  
  // wiggle outlier
  arma::vec Uloc3 = randu<arma::vec>(n0)/10;
  for(int i0=0; i0<n0; i0++){
    for(int idx_t=0; idx_t<tt; idx_t++){
      double t = tGrid(idx_t);
      Xout(i0, idx_t,3) = V(i0,3)*sin(10*datum::pi*(t - Uloc3(i0)));
    }
  }
  
  // non-diff wiggle outlier
  arma::vec Uloc4 = randu<arma::vec>(n0)/20;
  for(int i0=0; i0<n0; i0++){
    for(int idx_t=0; idx_t<tt; idx_t++){
      double t = tGrid(idx_t);
      Xout(i0, idx_t,4) = 400*V(i0,4)*f_nd_each(t - Uloc4(i0));
    }
  }
  
  // linear outlier
  Xout.slice(5) = V.col(5) * trans(2*tGrid-1);
  
  // quadratic outlier
  Xout.slice(6) = V.col(6) * trans(8*(tGrid-0.5) % (tGrid-0.5) - 1);
  
  // cubic outlier
  Xout.slice(7) = V.col(7) * trans(12 * sqrt(3) * tGrid % (tGrid-0.5) % (tGrid-1.0));
  
  return(Xout);
}

// FPCA
//[[Rcpp::export]]
List FPCA(arma::mat& X, arma::vec& tGrid, double fve=0.95){
  int tt = X.n_cols;
  double scal = range(tGrid) / tt;

  // FPCA

  arma::mat GaHat = cov(X); arma::vec gaHat; arma::mat phiHat;
  eig_sym(gaHat, phiHat, GaHat);
  gaHat = reverse(gaHat) * scal;
  phiHat = fliplr(phiHat) / sqrt(scal);

  int J = min(find(cumsum(gaHat) / sum(gaHat) >= fve))+1;
  
  List res = List::create(Named("ev") = gaHat,
                          Named("ef") = phiHat,
                          Named("J") = J);
  return(res);
}

//[[Rcpp::export]]
List RKHS(arma::vec& ev, arma::vec& J_vec, int M=1000){
  
  int J_size = J_vec.size();
  
  field<arma::mat> dir_field(J_size);
  arma::mat normRKHS_mat(M, J_size, fill::zeros);
  for(int ii_J=0; ii_J<J_size; ++ii_J){
    int J = J_vec(ii_J);
    arma::mat dir_mat(M, J, fill::zeros);
    for(int m=0; m<M; ++m){
      arma::vec v(J, fill::zeros);
      for(int j=0;j<J;++j){
        double z = randn<double>();
        v(j) = z * sqrt(ev(j));
      }
      
      double temp = sqrt(sum(v%v));
      v = v / temp;
      dir_mat.row(m) = trans(v);
      
      double normRKHS = sqrt(sum( v%v / ev.subvec(0, J-1)));
      normRKHS_mat(m, ii_J) = normRKHS;
    }
    dir_field(ii_J) = dir_mat;
  }
  
  List res = List::create(Named("normRKHS") = normRKHS_mat, 
                          Named("Dir") = dir_field);
  return(res);
}

//[[Rcpp::export]]
List sort_cpp(arma::vec& x){
  arma::vec s1 = sort(x, "descend");
  uvec s2 = sort_index(x, "descend");
  List res = List::create(s1, s2);
  return(res);
  }





//[[Rcpp::export]]
List RHD_inner(
    arma::mat& x, arma::mat& X, arma::vec& tGrid, arma::vec& ev, arma::mat& ef,
    arma::mat& ld_mat, arma::vec& J_vec, arma::vec& f_iqr_vec, arma::vec& M_vec, int M_type=0,
    Rcpp::CharacterVector ties="min",
    int out_type=1,
    int num_M=100, int stt=0){

  // M_vec: for each J in J_vec,
  // 1) The total number of drawn directions when M_type=0
  // 2) The number of accepted directions for the smallest lambda when M_type=1

  // Obtaining namespace of base package
  Environment pkg = Environment::namespace_env("base");
  // Picking up rank() function from base package
  Function rankk = pkg["rank"];

  int n = X.n_rows; int tt = X.n_cols;
  double scal = range(tGrid) / tt;
  int n0 = x.n_rows;

  int J_size = J_vec.size();
  int ld_size=ld_mat.n_cols;
  int f_iqr_size = f_iqr_vec.size();

  if(ld_mat.n_rows != J_size){
    stop("The number of columns of ld_mat should be equal to the length of J.");
  }
  if(ev.size() != ef.n_cols){
    stop("The number of eigenvalues should be equal to the number of eigenfunctions.");
  }
  if(J_size != M_vec.size()){
    stop(
      "J_vec and M_vec should have an equal size."
      );
  }

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  // =============================================
  // generate random directions and check their accpetance

  field<uvec> idx_acpt_field(J_size, ld_size);    // indices of accepted directions
  arma::mat acpt_rate_mat(J_size, ld_size);             // acceptance rates
  arma::vec M_total_vec(J_size);               // total number of drawn directions

  field<arma::mat> vCoef_field(J_size);        // generated directions
  field<arma::vec> normRKHS_field(J_size);     // RKHS norms of directions

  arma::mat M_used_mat(J_size, ld_size, fill::zeros);
  if(M_type==0){

    for(int ii_J=0;ii_J<J_size;++ii_J){

      int J = J_vec(ii_J); int J_idx = J-1;
      arma::vec ld_vec = trans(ld_mat.row(ii_J));
      arma::vec ld_vec_sort = sort(ld_vec);
      uvec ld_vec_sort_index = sort_index(ld_vec);
        
      int M = M_vec(ii_J);
      M_total_vec(ii_J) = M;
      
      // generates directions
      arma::vec normRKHS_vec(M, fill::zeros);
      arma::mat vCoef_mat(M, J, fill::zeros);
      
      for(int m=0; m<M; ++m){
        arma::vec v(J, fill::zeros);
        for(int j=0;j<J;++j){
          double z = randn<double>();
          v(j) = z * sqrt(ev(j));
        }
        
        double temp = sqrt(sum(v%v));
        v = v / temp;
        
        double normRKHS = sqrt(sum( v%v / ev.subvec(0, J_idx)));
        
        vCoef_mat.row(m) = trans(v);
        normRKHS_vec(m) = normRKHS;
      }
      vCoef_field(ii_J) = vCoef_mat;
      normRKHS_field(ii_J) = normRKHS_vec;
      
      for(int ii_ld=0;ii_ld<ld_size;++ii_ld){

        int ii_ld_sort = ld_vec_sort_index(ii_ld);
        double ld = ld_vec_sort(ii_ld);
        
        arma::vec idx_vec_zo(M, fill::zeros);
        
        for(int m=0; m<M; ++m){

          if(normRKHS_vec(m)<=ld){idx_vec_zo(m)=1;}
        }

        idx_acpt_field(ii_J, ii_ld_sort) = find(idx_vec_zo>0);
        acpt_rate_mat(ii_J, ii_ld_sort) = mean(idx_vec_zo);
        
        M_used_mat(ii_J, ii_ld_sort) = sum(idx_vec_zo);

      }
    }

  }else if(M_type==1){

    for(int ii_J=0;ii_J<J_size;++ii_J){

      int J = J_vec(ii_J); int J_idx = J-1;
      int M = M_vec(ii_J);
      arma::vec ld_vec = trans(ld_mat.row(ii_J));
      
      arma::vec ld_vec_sort = sort(ld_vec);   // if want, put "descend"   
      uvec ld_vec_sort_index = sort_index(ld_vec);
      
      // generate direcions for the largest lambda
      
      int ii_ld_extr = ld_vec_sort_index(0);
      double ld_extr = ld_vec_sort(0);
      
      // =====
      
      int M_draw=0;
      int M_draw_total=1;
      
      arma::vec normRKHS_vec(1, fill::zeros);
      arma::vec idx_vec_zo(1, fill::zeros);
      arma::mat vCoef_mat(1, J, fill::zeros);
      
      if(stt==1){
        Rcout << "Start to draw normal vectors " <<
          " when J = " << J << " and lambda = " << ld_extr <<
            " at " << currentDateTime() << "\n";
      }
      
      while(M_draw<M){
        
        M_draw_total = M_draw_total + 1;
        
        arma::vec v(J, fill::zeros);
        for(int j=0;j<J;++j){
          double z = randn<double>();
          v(j) = z * sqrt(ev(j));
        }
        
        double temp = sqrt(sum(v%v));
        v = v / temp;
        
        double normRKHS = sqrt(sum( v%v / ev.subvec(0, J_idx)));
        
        normRKHS_vec.resize(M_draw_total);
        normRKHS_vec(M_draw_total-1) = normRKHS;
        vCoef_mat.insert_rows(M_draw_total-1, 1);
        vCoef_mat.row(M_draw_total-1) = trans(v);
        
        idx_vec_zo.resize(M_draw_total);
        if(normRKHS <= ld_extr){
          idx_vec_zo(M_draw_total-1) = 1;
          M_draw = M_draw+1;
        }else{
          idx_vec_zo(M_draw_total-1) = 0;
        }
        
        if(stt == 1){
          if( (M_draw_total-1) % num_M == 0 ){
            Rcout << "Accepted M is " << M_draw+1 <<
              " among " << M_draw_total-1 << " total draws " <<
                " when J = " << J << " and lambda = " << ld_extr <<
                  " at " << currentDateTime() << "\n";
          }
        }
        
      }
      if(stt == 1){
        Rcout << "End to draw normal vectors with M_draw " <<
          M_draw << " and M_draw_total " << M_draw_total-1  <<
            " when J = " << J << " and lambda = " << ld_extr <<
              " at " << currentDateTime() << "\n";
      }
      
      M_draw_total = M_draw_total - 1;
      
      arma::mat vCoef_mat_record = vCoef_mat.rows(1, M_draw_total);
      arma::vec normRKHS_vec_record = normRKHS_vec.subvec(1, M_draw_total);
      
      arma::vec idx_vec_zo_record_extr = idx_vec_zo.subvec(1, M_draw_total);
      uvec idx_uvec_extr = find(idx_vec_zo_record_extr>0);
      
      double acpt_rate = M/(M_draw_total*1.0);
      
      idx_acpt_field(ii_J, ii_ld_extr) = idx_uvec_extr;
      acpt_rate_mat(ii_J, ii_ld_extr) = acpt_rate;
      M_used_mat(ii_J, ii_ld_extr) = sum(idx_vec_zo_record_extr);
      
      M_total_vec(ii_J) = M_draw_total;
      vCoef_field(ii_J) = vCoef_mat_record;
      normRKHS_field(ii_J) = normRKHS_vec_record;
      
      // ===
      // the rest of lambdas
      
      for(int ii_ld=1;ii_ld<ld_size;++ii_ld){
        
        int ii_ld_sort = ld_vec_sort_index(ii_ld);
        double ld = ld_vec_sort(ii_ld);
        
        arma::vec idx_vec_zo(M_draw_total, fill::zeros);
        
        for(int m=0; m<M_draw_total; ++m){
          
          if(normRKHS_vec_record(m)<=ld){idx_vec_zo(m)=1;}
        }
        
        idx_acpt_field(ii_J, ii_ld_sort) = find(idx_vec_zo>0);
        acpt_rate_mat(ii_J, ii_ld_sort) = mean(idx_vec_zo);
        
        M_used_mat(ii_J, ii_ld_sort) = sum(idx_vec_zo);
        
      }
      
    }

  }else{
    stop("M_type should be one of 0 and 1.");
  }

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  // =============================================
  // compute projections

  field<arma::mat> xCoef_field(J_size);        // (n0 times J) coefficients of x
  field<arma::mat> XCoef_field(J_size);        // (n times J) coefficients of X
  field<arma::mat> proj_x_field(J_size);       // (n0 times M_acpt) projections of x onto accepted directions
  field<arma::mat> proj_X_field(J_size);       // (n times M_acpt) projections of X onto accepted directions

  field<field<uvec>> out_all_field(J_size); // outliers of projections onto all drawn directions
  
  for(int ii_J=0;ii_J<J_size;++ii_J){
    
    int J = J_vec(ii_J); int J_idx = J-1;
    int M = M_vec(ii_J); int M_total = M_total_vec(ii_J);
    
    arma::mat vCoef_mat_record = vCoef_field(ii_J);
    arma::vec normRKHS_vec_record = normRKHS_field(ii_J);
    
    // coefficients and projections
    
    arma::mat temp_xCoef(n0, J, fill::zeros);
    arma::mat temp_proj_x(n0, M_total, fill::zeros);
    
    for(int i0=0; i0<n0; ++i0){
      arma::vec xtemp = trans(x.row(i0));
      arma::vec xCoef = trans(ef.cols(0, J_idx)) * (xtemp) * scal;
      arma::vec target_x = vCoef_mat_record * xCoef;
      
      temp_xCoef.row(i0) = trans(xCoef);
      temp_proj_x.row(i0) = trans(target_x);
    }
    
    xCoef_field(ii_J) = temp_xCoef;
    proj_x_field(ii_J) = temp_proj_x;
    
    arma::mat temp_XCoef(n, J, fill::zeros);
    arma::mat temp_proj_X(n, M_total, fill::zeros);
    
    for(int i=0; i<n; ++i){
      arma::vec Xtemp = trans(X.row(i));
      arma::vec XCoef = trans(ef.cols(0, J_idx)) * (Xtemp) * scal;
      arma::vec target_X = vCoef_mat_record * XCoef;
      
      temp_XCoef.row(i) = trans(XCoef);
      temp_proj_X.row(i) = trans(target_X);
    }
    
    XCoef_field(ii_J) = temp_XCoef;
    proj_X_field(ii_J) = temp_proj_X;
    
    // uni boxplots for all projections of x
    
    field<uvec> out_all_inner(M_total, f_iqr_size);
    for(int m=0; m<M_total; ++m){
      arma::vec proj_x = temp_proj_x.col(m);
      List resIQR = IQRout_vec(proj_x, f_iqr_vec);
      
      field<uvec> which_field = resIQR["out.which"];
      for(int ii_f = 0; ii_f<f_iqr_size; ++ii_f){
        out_all_inner(m, ii_f) = which_field(ii_f);
      }
    }
    out_all_field(ii_J) = out_all_inner;
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  // =============================================
  // compute depths and rankings and detect outliers
  
  
  field<arma::mat> depth_prob_field(J_size, ld_size);   // depth prob before min
  // indices of directions that give the depth minimum
  field<field<uvec>> idx_vCoef_min_field(J_size, ld_size);

  cube depth_cube(J_size, ld_size, n0, fill::zeros);
  cube rank_cube(J_size, ld_size, n0, fill::zeros);
  // indices of sample curves that achieve minimum depths
  field<uvec> idx_depth_min_field(J_size, ld_size); 
  // indices of outliers
  field<uvec> out_field(J_size, ld_size, f_iqr_size);
  
  for(int ii_J=0;ii_J<J_size;++ii_J){

    int J = J_vec(ii_J); int J_idx = J-1;
    int M = M_vec(ii_J); int M_total = M_total_vec(ii_J);
    
    arma::mat temp_proj_x = proj_x_field(ii_J);
    arma::mat temp_proj_X = proj_X_field(ii_J);
    
    field<uvec> out_all_inner = out_all_field(ii_J);
    
    for(int ii_ld=0;ii_ld<ld_size;++ii_ld){

      int M_used = M_used_mat(ii_J, ii_ld);
      double ld = ld_mat(ii_J, ii_ld);

      uvec idx_uvec = idx_acpt_field(ii_J, ii_ld);
      arma::mat temp_proj_x_sub = temp_proj_x.cols(idx_uvec);
      arma::mat temp_proj_X_sub = temp_proj_X.cols(idx_uvec);
      
      
      // depth

      arma::mat depth_prob_mat(n0, M_used, fill::zeros);
      arma::vec depth_vec(n0, fill::zeros);
      field<uvec> idx_vCoef_min_field_inner(n0);
      // field<uvec> idx_vCoef_min_field_inner_uvec(n0);
      for(int i0=0; i0<n0; ++i0){

        arma::vec target_x = trans(temp_proj_x_sub.row(i0));

        arma::mat idx_zero_one(n, idx_uvec.size(), fill::zeros);
        for(int i=0; i<n; ++i){
          arma::vec target_X = trans(temp_proj_X_sub.row(i));

          arma::vec target = target_X-target_x;

          for(int ll=0; ll<idx_uvec.size(); ++ll){
            if(target(ll)>=0){idx_zero_one(i, ll)=1;}else{idx_zero_one(i, ll)=0;}
          }

        }
        arma::vec idx_prob = trans(mean(idx_zero_one));

        depth_prob_mat.row(i0) = trans(idx_prob);

        double idx_prob_min = min(idx_prob);
        depth_vec(i0) = idx_prob_min;

        uvec temp_min_idx_uvec = find(idx_prob == idx_prob_min);
        idx_vCoef_min_field_inner(i0) = temp_min_idx_uvec;

      }

      depth_prob_field(ii_J, ii_ld) = depth_prob_mat;
      // depth_field(ii_J, ii_ld) = depth_vec;
      idx_vCoef_min_field(ii_J, ii_ld) = idx_vCoef_min_field_inner;

      // rank

      NumericVector rk_nv = rankk(depth_vec, Named("ties.method", ties));
      arma::vec rk_vec = as<vec>(rk_nv);
      // rank_field(ii_J, ii_ld) = rk_vec;
      for(int i0=0; i0<n0; ++i0){
        depth_cube(ii_J, ii_ld, i0) = depth_vec(i0);
        rank_cube(ii_J, ii_ld, i0) = rk_vec(i0);
      }
      
      
      uvec idx_depth_min_uvec = find(depth_vec == min(depth_vec));
      idx_depth_min_field(ii_J, ii_ld) = idx_depth_min_uvec+1;
      
      
      // =========================================
      // outlier detection
      
      // accepted directions
      field<uvec> out_all_inner_sub(idx_uvec.size(), f_iqr_size);
      for(int iii=0; iii<idx_uvec.size(); ++iii){
        uword iii0 = idx_uvec(iii);
        out_all_inner_sub.row(iii) = out_all_inner.row(iii0);
      }
      // field<uvec> out_all_inner_sub = out_all_inner.rows(idx_uvec);   
      
      
      
      for(int ii_f=0; ii_f<f_iqr_size; ++ii_f){
        
        uvec out_idx_stack_outer(1, fill::zeros);
        for(int ii_min=0; ii_min<idx_depth_min_uvec.size(); ++ii_min){
          
          // directions achieving min depth among accepted directions
          uword idx_min = idx_depth_min_uvec(ii_min);
          uvec temp_min_idx_vec = idx_vCoef_min_field_inner(idx_min);
          
          field<uvec> out_idx_inner(idx_uvec.size(), f_iqr_size);
          for(int iii=0; iii<temp_min_idx_vec.size(); ++iii){
            uword iii0 = temp_min_idx_vec(iii);
            out_idx_inner.row(iii) = out_all_inner_sub.row(iii0);
          }
          
          // field<uvec> out_idx_inner = out_all_inner_sub.rows(temp_min_idx_vec);
          
          // stacking outliers
          uvec out_idx_stack(1, fill::zeros);
          for(int m=0; m<out_idx_inner.n_rows; ++m){
            uvec out_idx_fine = out_idx_inner(m, ii_f);
            out_idx_fine = out_idx_fine + 1;
            out_idx_stack.insert_rows(0, out_idx_fine);
          }
          
          out_idx_stack_outer.insert_rows(0, out_idx_stack);
          
        }
        
        uvec out_unique = unique(out_idx_stack_outer);
        if(out_unique.size()==1){
          out_field(ii_J, ii_ld, ii_f) = {};
        }else{
          uvec out_final = out_unique.subvec(1,out_unique.size()-1);
          out_field(ii_J, ii_ld, ii_f) = out_final;
        }
        
         
      }
      
      
    }
  }

  List misc = List::create(Named("idx_acpt") = idx_acpt_field,
                           Named("acpt_rate") = acpt_rate_mat,
                           Named("M_total") = M_total_vec,
                           Named("M_used") = M_used_mat,
                           Named("vCoef") = vCoef_field,
                           Named("normRKHS") = normRKHS_field,
                           Named("xCoef") = xCoef_field,
                           Named("XCeof") = XCoef_field,
                           Named("proj_x") = proj_x_field,
                           Named("proj_X") = proj_X_field,
                           Named("depth_prob_each") = depth_prob_field,
                           Named("idx_vCoef_min") = idx_vCoef_min_field,
                           Named("idx_depth_min") = idx_depth_min_field,
                           Named("out_all") = out_all_field);

  List res = List::create(Named("depth") = depth_cube,
                          Named("rank") = rank_cube,
                          Named("out") = out_field,
                          Named("misc") = misc);
  return(res);
}









//[[Rcpp::export]]
List RHD_inner_prob(
    arma::mat& x, arma::mat& X, arma::vec& tGrid, arma::vec& ev, arma::mat& ef,
    arma::vec& p_vec, arma::vec& J_vec, arma::vec& f_iqr_vec, arma::vec& M_vec, 
    Rcpp::CharacterVector ties="min",
    int out_type=1,
    int num_M=100, int stt=0){
  
  // M_vec: for each J in J_vec,
  // 1) The total number of drawn directions when M_type=0
  // 2) The number of accepted directions for the smallest lambda when M_type=1
  
  // Obtaining namespace of base package
  Environment pkg = Environment::namespace_env("base");
  // Picking up rank() function from base package
  Function rankk = pkg["rank"];
  
  int n = X.n_rows; int tt = X.n_cols;
  double scal = range(tGrid) / tt;
  int n0 = x.n_rows;
  
  int J_size = J_vec.size();
  int ld_size = p_vec.size();
  int f_iqr_size = f_iqr_vec.size();
  
  // if(ld_mat.n_rows != J_size){
  //   stop("The number of columns of ld_mat should be equal to the length of J.");
  // }
  if(ev.size() != ef.n_cols){
    stop("The number of eigenvalues should be equal to the number of eigenfunctions.");
  }
  if(J_size != M_vec.size()){
    stop(
      "J_vec and M_vec should have an equal size."
    );
  }
  
  
  
  
  // =============================================
  // generate random directions and check their accpetance
  
  field<uvec> idx_acpt_field(J_size, ld_size);    // indices of accepted directions
  arma::mat acpt_rate_mat(J_size, ld_size);             // acceptance rates
  arma::vec M_total_vec(J_size);               // total number of drawn directions
  
  field<arma::mat> vCoef_field(J_size);        // generated directions
  field<arma::vec> normRKHS_field(J_size);     // RKHS norms of directions
  
  arma::mat M_used_mat(J_size, ld_size, fill::zeros);
  arma::mat ld_mat(J_size, ld_size);    // ld matrix from quantile arma::vec
  
  
  for(int ii_J=0;ii_J<J_size;++ii_J){
    
    int J = J_vec(ii_J); int J_idx = J-1;
    int M = M_vec(ii_J);
    M_total_vec(ii_J) = M;
    
    // generates directions
    arma::vec normRKHS_vec(M, fill::zeros);
    arma::mat vCoef_mat(M, J, fill::zeros);
    
    for(int m=0; m<M; ++m){
      arma::vec v(J, fill::zeros);
      for(int j=0;j<J;++j){
        double z = randn<double>();
        v(j) = z * sqrt(ev(j));
      }
      
      double temp = sqrt(sum(v%v));
      v = v / temp;
      
      double normRKHS = sqrt(sum( v%v / ev.subvec(0, J_idx)));
      
      vCoef_mat.row(m) = trans(v);
      normRKHS_vec(m) = normRKHS;
    }
    vCoef_field(ii_J) = vCoef_mat;
    normRKHS_field(ii_J) = normRKHS_vec;
    
    // regularization
    arma::vec ld_vec = quantile(normRKHS_vec, p_vec);
    ld_mat.row(ii_J) = trans(ld_vec);
    
    // arma::vec ld_vec = trans(ld_mat.row(ii_J));
    arma::vec ld_vec_sort = sort(ld_vec);
    uvec ld_vec_sort_index = sort_index(ld_vec);
    
    for(int ii_ld=0;ii_ld<ld_size;++ii_ld){
      
      int ii_ld_sort = ld_vec_sort_index(ii_ld);
      double ld = ld_vec_sort(ii_ld);
      
      arma::vec idx_vec_zo(M, fill::zeros);
      
      for(int m=0; m<M; ++m){
        
        if(normRKHS_vec(m)<=ld){idx_vec_zo(m)=1;}
      }
      
      idx_acpt_field(ii_J, ii_ld_sort) = find(idx_vec_zo>0);
      acpt_rate_mat(ii_J, ii_ld_sort) = mean(idx_vec_zo);
      
      M_used_mat(ii_J, ii_ld_sort) = sum(idx_vec_zo);
      
    }
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  // =============================================
  // compute projections
  
  field<arma::mat> xCoef_field(J_size);        // (n0 times J) coefficients of x
  field<arma::mat> XCoef_field(J_size);        // (n times J) coefficients of X
  field<arma::mat> proj_x_field(J_size);       // (n0 times M_acpt) projections of x onto accepted directions
  field<arma::mat> proj_X_field(J_size);       // (n times M_acpt) projections of X onto accepted directions
  
  field<field<uvec>> out_all_field(J_size); // outliers of projections onto all drawn directions
  
  for(int ii_J=0;ii_J<J_size;++ii_J){
    
    int J = J_vec(ii_J); int J_idx = J-1;
    int M = M_vec(ii_J); int M_total = M_total_vec(ii_J);
    
    arma::mat vCoef_mat_record = vCoef_field(ii_J);
    arma::vec normRKHS_vec_record = normRKHS_field(ii_J);
    
    // coefficients and projections
    
    arma::mat temp_xCoef(n0, J, fill::zeros);
    arma::mat temp_proj_x(n0, M_total, fill::zeros);
    
    for(int i0=0; i0<n0; ++i0){
      arma::vec xtemp = trans(x.row(i0));
      arma::vec xCoef = trans(ef.cols(0, J_idx)) * (xtemp) * scal;
      arma::vec target_x = vCoef_mat_record * xCoef;
      
      temp_xCoef.row(i0) = trans(xCoef);
      temp_proj_x.row(i0) = trans(target_x);
    }
    
    xCoef_field(ii_J) = temp_xCoef;
    proj_x_field(ii_J) = temp_proj_x;
    
    arma::mat temp_XCoef(n, J, fill::zeros);
    arma::mat temp_proj_X(n, M_total, fill::zeros);
    
    for(int i=0; i<n; ++i){
      arma::vec Xtemp = trans(X.row(i));
      arma::vec XCoef = trans(ef.cols(0, J_idx)) * (Xtemp) * scal;
      arma::vec target_X = vCoef_mat_record * XCoef;
      
      temp_XCoef.row(i) = trans(XCoef);
      temp_proj_X.row(i) = trans(target_X);
    }
    
    XCoef_field(ii_J) = temp_XCoef;
    proj_X_field(ii_J) = temp_proj_X;
    
    // uni boxplots for all projections of x
    
    field<uvec> out_all_inner(M_total, f_iqr_size);
    for(int m=0; m<M_total; ++m){
      arma::vec proj_x = temp_proj_x.col(m);
      List resIQR = IQRout_vec(proj_x, f_iqr_vec);
      
      field<uvec> which_field = resIQR["out.which"];
      for(int ii_f = 0; ii_f<f_iqr_size; ++ii_f){
        out_all_inner(m, ii_f) = which_field(ii_f);
      }
    }
    out_all_field(ii_J) = out_all_inner;
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  // =============================================
  // compute depths and rankings and detect outliers
  
  
  field<arma::mat> depth_prob_field(J_size, ld_size);   // depth prob before min
  // indices of directions that give the depth minimum
  field<field<uvec>> idx_vCoef_min_field(J_size, ld_size);
  
  cube depth_cube(J_size, ld_size, n0, fill::zeros);
  cube rank_cube(J_size, ld_size, n0, fill::zeros);
  // indices of sample curves that achieve minimum depths
  field<uvec> idx_depth_min_field(J_size, ld_size); 
  // indices of outliers
  field<uvec> out_field(J_size, ld_size, f_iqr_size);
  
  for(int ii_J=0;ii_J<J_size;++ii_J){
    
    int J = J_vec(ii_J); int J_idx = J-1;
    int M = M_vec(ii_J); int M_total = M_total_vec(ii_J);
    
    arma::mat temp_proj_x = proj_x_field(ii_J);
    arma::mat temp_proj_X = proj_X_field(ii_J);
    
    field<uvec> out_all_inner = out_all_field(ii_J);
    
    for(int ii_ld=0;ii_ld<ld_size;++ii_ld){
      
      int M_used = M_used_mat(ii_J, ii_ld);
      double ld = ld_mat(ii_J, ii_ld);
      
      uvec idx_uvec = idx_acpt_field(ii_J, ii_ld);
      arma::mat temp_proj_x_sub = temp_proj_x.cols(idx_uvec);
      arma::mat temp_proj_X_sub = temp_proj_X.cols(idx_uvec);
      
      
      // depth
      
      arma::mat depth_prob_mat(n0, M_used, fill::zeros);
      arma::vec depth_vec(n0, fill::zeros);
      field<uvec> idx_vCoef_min_field_inner(n0);
      // field<uvec> idx_vCoef_min_field_inner_uvec(n0);
      for(int i0=0; i0<n0; ++i0){
        
        arma::vec target_x = trans(temp_proj_x_sub.row(i0));
        
        arma::mat idx_zero_one(n, idx_uvec.size(), fill::zeros);
        for(int i=0; i<n; ++i){
          arma::vec target_X = trans(temp_proj_X_sub.row(i));
          
          arma::vec target = target_X-target_x;
          
          for(int ll=0; ll<idx_uvec.size(); ++ll){
            if(target(ll)>=0){idx_zero_one(i, ll)=1;}else{idx_zero_one(i, ll)=0;}
          }
          
        }
        arma::vec idx_prob = trans(mean(idx_zero_one));
        
        depth_prob_mat.row(i0) = trans(idx_prob);
        
        double idx_prob_min = min(idx_prob);
        depth_vec(i0) = idx_prob_min;
        
        uvec temp_min_idx_uvec = find(idx_prob == idx_prob_min);
        idx_vCoef_min_field_inner(i0) = temp_min_idx_uvec;
        
      }
      
      depth_prob_field(ii_J, ii_ld) = depth_prob_mat;
      // depth_field(ii_J, ii_ld) = depth_vec;
      idx_vCoef_min_field(ii_J, ii_ld) = idx_vCoef_min_field_inner;
      
      // rank
      
      NumericVector rk_nv = rankk(depth_vec, Named("ties.method", ties));
      arma::vec rk_vec = as<vec>(rk_nv);
      // rank_field(ii_J, ii_ld) = rk_vec;
      for(int i0=0; i0<n0; ++i0){
        depth_cube(ii_J, ii_ld, i0) = depth_vec(i0);
        rank_cube(ii_J, ii_ld, i0) = rk_vec(i0);
      }
      
      
      uvec idx_depth_min_uvec = find(depth_vec == min(depth_vec));
      idx_depth_min_field(ii_J, ii_ld) = idx_depth_min_uvec+1;
      
      
      // =========================================
      // outlier detection
      
      // accepted directions
      field<uvec> out_all_inner_sub(idx_uvec.size(), f_iqr_size);
      for(int iii=0; iii<idx_uvec.size(); ++iii){
        uword iii0 = idx_uvec(iii);
        out_all_inner_sub.row(iii) = out_all_inner.row(iii0);
      }
      // field<uvec> out_all_inner_sub = out_all_inner.rows(idx_uvec);   
      
      
      
      for(int ii_f=0; ii_f<f_iqr_size; ++ii_f){
        
        uvec out_idx_stack_outer(1, fill::zeros);
        for(int ii_min=0; ii_min<idx_depth_min_uvec.size(); ++ii_min){
          
          // directions achieving min depth among accepted directions
          uword idx_min = idx_depth_min_uvec(ii_min);
          uvec temp_min_idx_vec = idx_vCoef_min_field_inner(idx_min);
          
          field<uvec> out_idx_inner(idx_uvec.size(), f_iqr_size);
          for(int iii=0; iii<temp_min_idx_vec.size(); ++iii){
            uword iii0 = temp_min_idx_vec(iii);
            out_idx_inner.row(iii) = out_all_inner_sub.row(iii0);
          }
          
          // field<uvec> out_idx_inner = out_all_inner_sub.rows(temp_min_idx_vec);
          
          // stacking outliers
          uvec out_idx_stack(1, fill::zeros);
          for(int m=0; m<out_idx_inner.n_rows; ++m){
            uvec out_idx_fine = out_idx_inner(m, ii_f);
            out_idx_fine = out_idx_fine + 1;
            out_idx_stack.insert_rows(0, out_idx_fine);
          }
          
          out_idx_stack_outer.insert_rows(0, out_idx_stack);
          
        }
        
        uvec out_unique = unique(out_idx_stack_outer);
        if(out_unique.size()==1){
          out_field(ii_J, ii_ld, ii_f) = {};
        }else{
          uvec out_final = out_unique.subvec(1,out_unique.size()-1);
          out_field(ii_J, ii_ld, ii_f) = out_final;
        }
        
        
      }
      
      
    }
  }
  
  List misc = List::create(Named("idx_acpt") = idx_acpt_field,
                           Named("acpt_rate") = acpt_rate_mat,
                           Named("M_total") = M_total_vec,
                           Named("M_used") = M_used_mat,
                           Named("vCoef") = vCoef_field,
                           Named("normRKHS") = normRKHS_field,
                           Named("ld_mat") = ld_mat,
                           Named("xCoef") = xCoef_field,
                           Named("XCeof") = XCoef_field,
                           Named("proj_x") = proj_x_field,
                           Named("proj_X") = proj_X_field,
                           Named("depth_prob_each") = depth_prob_field,
                           Named("idx_vCoef_min") = idx_vCoef_min_field,
                           Named("idx_depth_min") = idx_depth_min_field,
                           Named("out_all") = out_all_field);
  
  List res = List::create(Named("depth") = depth_cube,
                          Named("rank") = rank_cube,
                          Named("out") = out_field,
                          Named("misc") = misc);
  return(res);
}





  
  
//
//
//
// // =========================================================
// simulations

//[[Rcpp::export]]
List RHDonly(
    arma::mat& x, arma::mat& X, arma::vec& tGrid, arma::vec& ev, arma::mat& ef,
    arma::mat& ld_mat, arma::vec& J_vec, arma::vec& M_vec, int M_type=0,
    Rcpp::CharacterVector ties="min",
    int num_M=100, int stt=0){

  // M_vec: for each J in J_vec,
  // 1) The total number of drawn directions when M_type=0
  // 2) The number of accepted directions for the smallest lambda when M_type=1
  
  // Obtaining namespace of base package
  Environment pkg = Environment::namespace_env("base");
  // Picking up rank() function from base package
  Function rankk = pkg["rank"];
  
  int n = X.n_rows; int tt = X.n_cols;
  double scal = range(tGrid) / tt;
  int n0 = x.n_rows;
  
  int J_size = J_vec.size();
  int ld_size=ld_mat.n_cols;
  // int f_iqr_size = f_iqr_vec.size();
  
  if(ld_mat.n_rows != J_size){
    stop("The number of columns of ld_mat should be equal to the length of J.");
  }
  if(ev.size() != ef.n_cols){
    stop("The number of eigenvalues should be equal to the number of eigenfunctions.");
  }
  if(J_size != M_vec.size()){
    stop(
      "J_vec and M_vec should have an equal size."
    );
  }
  
  // generate random directions and check their accpetance
  
  field<uvec> idx_acpt_field(J_size, ld_size);    // indices of accepted directions
  arma::mat acpt_rate_mat(J_size, ld_size);             // acceptance rates
  arma::vec M_total_vec(J_size);               // total number of drawn directions
  
  field<arma::mat> vCoef_field(J_size);        // generated directions
  field<arma::vec> normRKHS_field(J_size);     // RKHS norms of directions
  
  arma::mat M_used_mat(J_size, ld_size, fill::zeros);
  if(M_type==0){
    
    for(int ii_J=0;ii_J<J_size;++ii_J){
      
      int J = J_vec(ii_J); int J_idx = J-1;
      arma::vec ld_vec = trans(ld_mat.row(ii_J));
      arma::vec ld_vec_sort = sort(ld_vec);
      uvec ld_vec_sort_index = sort_index(ld_vec);
      
      int M = M_vec(ii_J);
      M_total_vec(ii_J) = M;
      
      // generates directions
      arma::vec normRKHS_vec(M, fill::zeros);
      arma::mat vCoef_mat(M, J, fill::zeros);
      
      for(int m=0; m<M; ++m){
        arma::vec v(J, fill::zeros);
        for(int j=0;j<J;++j){
          double z = randn<double>();
          v(j) = z * sqrt(ev(j));
        }
        
        double temp = sqrt(sum(v%v));
        v = v / temp;
        
        double normRKHS = sqrt(sum( v%v / ev.subvec(0, J_idx)));
        
        vCoef_mat.row(m) = trans(v);
        normRKHS_vec(m) = normRKHS;
      }
      vCoef_field(ii_J) = vCoef_mat;
      normRKHS_field(ii_J) = normRKHS_vec;
      
      for(int ii_ld=0;ii_ld<ld_size;++ii_ld){
        
        int ii_ld_sort = ld_vec_sort_index(ii_ld);
        double ld = ld_vec_sort(ii_ld);
        
        arma::vec idx_vec_zo(M, fill::zeros);
        
        for(int m=0; m<M; ++m){
          
          if(normRKHS_vec(m)<=ld){idx_vec_zo(m)=1;}
        }
        
        idx_acpt_field(ii_J, ii_ld_sort) = find(idx_vec_zo>0);
        acpt_rate_mat(ii_J, ii_ld_sort) = mean(idx_vec_zo);
        
        M_used_mat(ii_J, ii_ld_sort) = sum(idx_vec_zo);
        
      }
    }
    
  }else if(M_type==1){
    
    for(int ii_J=0;ii_J<J_size;++ii_J){
      
      int J = J_vec(ii_J); int J_idx = J-1;
      int M = M_vec(ii_J);
      arma::vec ld_vec = trans(ld_mat.row(ii_J));
      
      arma::vec ld_vec_sort = sort(ld_vec);   // if want, put "descend"   
      uvec ld_vec_sort_index = sort_index(ld_vec);
      
      // generate direcions for the largest lambda
      
      int ii_ld_extr = ld_vec_sort_index(0);
      double ld_extr = ld_vec_sort(0);
      
      // =====
      
      int M_draw=0;
      int M_draw_total=1;
      
      arma::vec normRKHS_vec(1, fill::zeros);
      arma::vec idx_vec_zo(1, fill::zeros);
      arma::mat vCoef_mat(1, J, fill::zeros);
      
      if(stt==1){
        Rcout << "Start to draw normal vectors " <<
          " when J = " << J << " and lambda = " << ld_extr <<
            " at " << currentDateTime() << "\n";
      }
      
      while(M_draw<M){
        
        M_draw_total = M_draw_total + 1;
        
        arma::vec v(J, fill::zeros);
        for(int j=0;j<J;++j){
          double z = randn<double>();
          v(j) = z * sqrt(ev(j));
        }
        
        double temp = sqrt(sum(v%v));
        v = v / temp;
        
        double normRKHS = sqrt(sum( v%v / ev.subvec(0, J_idx)));
        
        normRKHS_vec.resize(M_draw_total);
        normRKHS_vec(M_draw_total-1) = normRKHS;
        vCoef_mat.insert_rows(M_draw_total-1, 1);
        vCoef_mat.row(M_draw_total-1) = trans(v);
        
        idx_vec_zo.resize(M_draw_total);
        if(normRKHS <= ld_extr){
          idx_vec_zo(M_draw_total-1) = 1;
          M_draw = M_draw+1;
        }else{
          idx_vec_zo(M_draw_total-1) = 0;
        }
        
        if(stt == 1){
          if( (M_draw_total-1) % num_M == 0 ){
            Rcout << "Accepted M is " << M_draw+1 <<
              " among " << M_draw_total-1 << " total draws " <<
                " when J = " << J << " and lambda = " << ld_extr <<
                  " at " << currentDateTime() << "\n";
          }
        }
        
      }
      if(stt == 1){
        Rcout << "End to draw normal vectors with M_draw " <<
          M_draw << " and M_draw_total " << M_draw_total-1  <<
            " when J = " << J << " and lambda = " << ld_extr <<
              " at " << currentDateTime() << "\n";
      }
      
      M_draw_total = M_draw_total - 1;
      
      arma::mat vCoef_mat_record = vCoef_mat.rows(1, M_draw_total);
      arma::vec normRKHS_vec_record = normRKHS_vec.subvec(1, M_draw_total);
      
      arma::vec idx_vec_zo_record_extr = idx_vec_zo.subvec(1, M_draw_total);
      uvec idx_uvec_extr = find(idx_vec_zo_record_extr>0);
      
      double acpt_rate = M/(M_draw_total*1.0);
      
      idx_acpt_field(ii_J, ii_ld_extr) = idx_uvec_extr;
      acpt_rate_mat(ii_J, ii_ld_extr) = acpt_rate;
      M_used_mat(ii_J, ii_ld_extr) = sum(idx_vec_zo_record_extr);
      
      M_total_vec(ii_J) = M_draw_total;
      vCoef_field(ii_J) = vCoef_mat_record;
      normRKHS_field(ii_J) = normRKHS_vec_record;
      
      // ===
      // the rest of lambdas
      
      for(int ii_ld=1;ii_ld<ld_size;++ii_ld){
        
        int ii_ld_sort = ld_vec_sort_index(ii_ld);
        double ld = ld_vec_sort(ii_ld);
        
        arma::vec idx_vec_zo(M_draw_total, fill::zeros);
        
        for(int m=0; m<M_draw_total; ++m){
          
          if(normRKHS_vec_record(m)<=ld){idx_vec_zo(m)=1;}
        }
        
        idx_acpt_field(ii_J, ii_ld_sort) = find(idx_vec_zo>0);
        acpt_rate_mat(ii_J, ii_ld_sort) = mean(idx_vec_zo);
        
        M_used_mat(ii_J, ii_ld_sort) = sum(idx_vec_zo);
        
      }
      
    }
    
  }else{
    stop("M_type should be one of 0 and 1.");
  }
  
  // =======================================
  // compute projections
  
  field<arma::mat> xCoef_field(J_size);        // (n0 times J) coefficients of x
  field<arma::mat> XCoef_field(J_size);        // (n times J) coefficients of X
  field<arma::mat> proj_x_field(J_size);       // (n0 times M_acpt) projections of x onto accepted directions
  field<arma::mat> proj_X_field(J_size);       // (n times M_acpt) projections of X onto accepted directions
  
  for(int ii_J=0;ii_J<J_size;++ii_J){
    
    int J = J_vec(ii_J); int J_idx = J-1;
    int M = M_vec(ii_J); int M_total = M_total_vec(ii_J);
    
    arma::mat vCoef_mat_record = vCoef_field(ii_J);
    arma::vec normRKHS_vec_record = normRKHS_field(ii_J);
    
    // coefficients and projections
    
    arma::mat temp_xCoef(n0, J, fill::zeros);
    arma::mat temp_proj_x(n0, M_total, fill::zeros);
    
    for(int i0=0; i0<n0; ++i0){
      arma::vec xtemp = trans(x.row(i0));
      arma::vec xCoef = trans(ef.cols(0, J_idx)) * (xtemp) * scal;
      arma::vec target_x = vCoef_mat_record * xCoef;
      
      temp_xCoef.row(i0) = trans(xCoef);
      temp_proj_x.row(i0) = trans(target_x);
    }
    
    xCoef_field(ii_J) = temp_xCoef;
    proj_x_field(ii_J) = temp_proj_x;
    
    arma::mat temp_XCoef(n, J, fill::zeros);
    arma::mat temp_proj_X(n, M_total, fill::zeros);
    
    for(int i=0; i<n; ++i){
      arma::vec Xtemp = trans(X.row(i));
      arma::vec XCoef = trans(ef.cols(0, J_idx)) * (Xtemp) * scal;
      arma::vec target_X = vCoef_mat_record * XCoef;
      
      temp_XCoef.row(i) = trans(XCoef);
      temp_proj_X.row(i) = trans(target_X);
    }
    
    XCoef_field(ii_J) = temp_XCoef;
    proj_X_field(ii_J) = temp_proj_X;
    
  }
  
  
  
  
  
  // =======================================
  // compute depths and rankings and detect outliers
  
  
  field<arma::mat> depth_prob_field(J_size, ld_size);   // depth prob before min
  field<field<uvec>> idx_vCoef_min_field(J_size, ld_size);// indices of directions that achieve the depth minimum

  
  cube depth_cube(J_size, ld_size, n0, fill::zeros);
  cube rank_cube(J_size, ld_size, n0, fill::zeros);
  // field<uvec> out_field(J_size, ld_size, f_iqr_size);
  
  
  for(int ii_J=0;ii_J<J_size;++ii_J){
    
    int J = J_vec(ii_J); int J_idx = J-1;
    int M = M_vec(ii_J); int M_total = M_total_vec(ii_J);
    
    arma::mat temp_proj_x = proj_x_field(ii_J);
    arma::mat temp_proj_X = proj_X_field(ii_J);
    
    // field<uvec> out_all_inner = out_all_field(ii_J);
    
    for(int ii_ld=0;ii_ld<ld_size;++ii_ld){
      
      int M_used = M_used_mat(ii_J, ii_ld);
      double ld = ld_mat(ii_J, ii_ld);
      
      uvec idx_uvec = idx_acpt_field(ii_J, ii_ld);
      arma::mat temp_proj_x_sub = temp_proj_x.cols(idx_uvec);
      arma::mat temp_proj_X_sub = temp_proj_X.cols(idx_uvec);
      
      // depth
      
      arma::mat depth_prob_mat(n0, M_used, fill::zeros);
      arma::vec depth_vec(n0, fill::zeros);
      field<uvec> idx_vCoef_min_field_inner(n0);
      // field<uvec> idx_vCoef_min_field_inner_uvec(n0);
      for(int i0=0; i0<n0; ++i0){
        
        arma::vec target_x = trans(temp_proj_x_sub.row(i0));
        
        arma::mat idx_zero_one(n, idx_uvec.size(), fill::zeros);
        for(int i=0; i<n; ++i){
          arma::vec target_X = trans(temp_proj_X_sub.row(i));
          
          arma::vec target = target_X-target_x;
          
          for(int ll=0; ll<idx_uvec.size(); ++ll){
            if(target(ll)>=0){idx_zero_one(i, ll)=1;}else{idx_zero_one(i, ll)=0;}
          }
          
        }
        arma::vec idx_prob = trans(mean(idx_zero_one));
        
        depth_prob_mat.row(i0) = trans(idx_prob);
        
        double idx_prob_min = min(idx_prob);
        depth_vec(i0) = idx_prob_min;
        
        uvec temp_min_idx_uvec = find(idx_prob == idx_prob_min);
        idx_vCoef_min_field_inner(i0) = temp_min_idx_uvec;
        
      }
      
      depth_prob_field(ii_J, ii_ld) = depth_prob_mat;
      // depth_field(ii_J, ii_ld) = depth_vec;
      idx_vCoef_min_field(ii_J, ii_ld) = idx_vCoef_min_field_inner;
      
      // rank
      
      NumericVector rk_nv = rankk(depth_vec, Named("ties.method", ties));
      arma::vec rk_vec = as<arma::vec>(rk_nv);
      // rank_field(ii_J, ii_ld) = rk_vec;
      for(int i0=0; i0<n0; ++i0){
        depth_cube(ii_J, ii_ld, i0) = depth_vec(i0);
        rank_cube(ii_J, ii_ld, i0) = rk_vec(i0);
      }
      
      
    }
  }
  
  List misc = List::create(Named("idx_acpt") = idx_acpt_field,
                           Named("acpt_rate") = acpt_rate_mat,
                           Named("M_total") = M_total_vec,
                           Named("M_used") = M_used_mat,
                           Named("vCoef") = vCoef_field,
                           Named("normRKHS") = normRKHS_field,
                           Named("xCoef") = xCoef_field,
                           Named("XCeof") = XCoef_field,
                           Named("proj_x") = proj_x_field,
                           Named("proj_X") = proj_X_field,
                           Named("depth_prob_each") = depth_prob_field,
                           Named("idx_vCoef_min") = idx_vCoef_min_field);
  
  // List out = List::create(Named("idx_out") = idx_out_field,
  //                         Named("x_out") = x_out_field);
  
  List res = List::create(Named("depth") = depth_cube,
                          Named("rank") = rank_cube,
                          Named("misc") = misc);
  return(res);
}










//[[Rcpp::export]]
List RHDonly_prob(
    arma::mat& x, arma::mat& X, arma::vec& tGrid, arma::vec& ev, arma::mat& ef,
    arma::vec& p_vec, arma::vec& J_vec, arma::vec& M_vec, 
    Rcpp::CharacterVector ties="min",
    int num_M=100, int stt=0){
  
  // M_vec: for each J in J_vec,
  // 1) The total number of drawn directions when M_type=0
  // 2) The number of accepted directions for the smallest lambda when M_type=1
  
  // Obtaining namespace of base package
  Environment pkg = Environment::namespace_env("base");
  // Picking up rank() function from base package
  Function rankk = pkg["rank"];
  
  int n = X.n_rows; int tt = X.n_cols;
  double scal = range(tGrid) / tt;
  int n0 = x.n_rows;
  
  int J_size = J_vec.size();
  int ld_size = p_vec.size();
  // int f_iqr_size = f_iqr_vec.size();
  
  // if(ld_mat.n_rows != J_size){
  //   stop("The number of columns of ld_mat should be equal to the length of J.");
  // }
  if(ev.size() != ef.n_cols){
    stop("The number of eigenvalues should be equal to the number of eigenfunctions.");
  }
  if(J_size != M_vec.size()){
    stop(
      "J_vec and M_vec should have an equal size."
    );
  }
  
  
  
  
  // =============================================
  // generate random directions and check their accpetance
  
  field<uvec> idx_acpt_field(J_size, ld_size);    // indices of accepted directions
  arma::mat acpt_rate_mat(J_size, ld_size);             // acceptance rates
  arma::vec M_total_vec(J_size);               // total number of drawn directions
  
  field<arma::mat> vCoef_field(J_size);        // generated directions
  field<arma::vec> normRKHS_field(J_size);     // RKHS norms of directions
  
  arma::mat M_used_mat(J_size, ld_size, fill::zeros);
  arma::mat ld_mat(J_size, ld_size, fill::zeros);    // ld matrix from quantile arma::vec
  
  
  for(int ii_J=0;ii_J<J_size;++ii_J){
    
    int J = J_vec(ii_J); int J_idx = J-1;
    int M = M_vec(ii_J);
    M_total_vec(ii_J) = M;
    
    // generates directions
    arma::vec normRKHS_vec(M, fill::zeros);
    arma::mat vCoef_mat(M, J, fill::zeros);
    
    for(int m=0; m<M; ++m){
      arma::vec v(J, fill::zeros);
      for(int j=0;j<J;++j){
        double z = randn<double>();
        v(j) = z * sqrt(ev(j));
      }
      
      double temp = sqrt(sum(v%v));
      v = v / temp;
      
      double normRKHS = sqrt(sum( v%v / ev.subvec(0, J_idx)));
      
      vCoef_mat.row(m) = trans(v);
      normRKHS_vec(m) = normRKHS;
    }
    vCoef_field(ii_J) = vCoef_mat;
    normRKHS_field(ii_J) = normRKHS_vec;
    
    // regularization
    arma::vec ld_vec = quantile(normRKHS_vec, p_vec);
    ld_mat.row(ii_J) = trans(ld_vec);

    // arma::vec ld_vec = trans(ld_mat.row(ii_J));
    arma::vec ld_vec_sort = sort(ld_vec);
    uvec ld_vec_sort_index = sort_index(ld_vec);

    for(int ii_ld=0;ii_ld<ld_size;++ii_ld){

      int ii_ld_sort = ld_vec_sort_index(ii_ld);
      double ld = ld_vec_sort(ii_ld);

      arma::vec idx_vec_zo(M, fill::zeros);

      for(int m=0; m<M; ++m){

        if(normRKHS_vec(m)<=ld){idx_vec_zo(m)=1;}
      }

      idx_acpt_field(ii_J, ii_ld_sort) = find(idx_vec_zo>0);
      acpt_rate_mat(ii_J, ii_ld_sort) = mean(idx_vec_zo);

      M_used_mat(ii_J, ii_ld_sort) = sum(idx_vec_zo);

    }

  }
  
  
  
  
  
  
  
  
  
  

  // =======================================
  // compute projections

  field<arma::mat> xCoef_field(J_size);        // (n0 times J) coefficients of x
  field<arma::mat> XCoef_field(J_size);        // (n times J) coefficients of X
  field<arma::mat> proj_x_field(J_size);       // (n0 times M_acpt) projections of x onto accepted directions
  field<arma::mat> proj_X_field(J_size);       // (n times M_acpt) projections of X onto accepted directions

  for(int ii_J=0;ii_J<J_size;++ii_J){

    int J = J_vec(ii_J); int J_idx = J-1;
    int M = M_vec(ii_J); int M_total = M_total_vec(ii_J);

    arma::mat vCoef_mat_record = vCoef_field(ii_J);
    arma::vec normRKHS_vec_record = normRKHS_field(ii_J);

    // coefficients and projections

    arma::mat temp_xCoef(n0, J, fill::zeros);
    arma::mat temp_proj_x(n0, M_total, fill::zeros);

    for(int i0=0; i0<n0; ++i0){
      arma::vec xtemp = trans(x.row(i0));
      arma::vec xCoef = trans(ef.cols(0, J_idx)) * (xtemp) * scal;
      arma::vec target_x = vCoef_mat_record * xCoef;

      temp_xCoef.row(i0) = trans(xCoef);
      temp_proj_x.row(i0) = trans(target_x);
    }

    xCoef_field(ii_J) = temp_xCoef;
    proj_x_field(ii_J) = temp_proj_x;

    arma::mat temp_XCoef(n, J, fill::zeros);
    arma::mat temp_proj_X(n, M_total, fill::zeros);

    for(int i=0; i<n; ++i){
      arma::vec Xtemp = trans(X.row(i));
      arma::vec XCoef = trans(ef.cols(0, J_idx)) * (Xtemp) * scal;
      arma::vec target_X = vCoef_mat_record * XCoef;

      temp_XCoef.row(i) = trans(XCoef);
      temp_proj_X.row(i) = trans(target_X);
    }

    XCoef_field(ii_J) = temp_XCoef;
    proj_X_field(ii_J) = temp_proj_X;

  }





  // =======================================
  // compute depths and rankings and detect outliers


  field<arma::mat> depth_prob_field(J_size, ld_size);   // depth prob before min
  field<field<uvec>> idx_vCoef_min_field(J_size, ld_size);// indices of directions that achieve the depth minimum


  cube depth_cube(J_size, ld_size, n0, fill::zeros);
  cube rank_cube(J_size, ld_size, n0, fill::zeros);
  // field<uvec> out_field(J_size, ld_size, f_iqr_size);


  for(int ii_J=0;ii_J<J_size;++ii_J){

    int J = J_vec(ii_J); int J_idx = J-1;
    int M = M_vec(ii_J); int M_total = M_total_vec(ii_J);

    arma::mat temp_proj_x = proj_x_field(ii_J);
    arma::mat temp_proj_X = proj_X_field(ii_J);

    // field<uvec> out_all_inner = out_all_field(ii_J);

    for(int ii_ld=0;ii_ld<ld_size;++ii_ld){

      int M_used = M_used_mat(ii_J, ii_ld);
      double ld = ld_mat(ii_J, ii_ld);

      uvec idx_uvec = idx_acpt_field(ii_J, ii_ld);
      arma::mat temp_proj_x_sub = temp_proj_x.cols(idx_uvec);
      arma::mat temp_proj_X_sub = temp_proj_X.cols(idx_uvec);

      // depth

      arma::mat depth_prob_mat(n0, M_used, fill::zeros);
      arma::vec depth_vec(n0, fill::zeros);
      field<uvec> idx_vCoef_min_field_inner(n0);
      // field<uvec> idx_vCoef_min_field_inner_uvec(n0);
      for(int i0=0; i0<n0; ++i0){

        arma::vec target_x = trans(temp_proj_x_sub.row(i0));

        arma::mat idx_zero_one(n, idx_uvec.size(), fill::zeros);
        for(int i=0; i<n; ++i){
          arma::vec target_X = trans(temp_proj_X_sub.row(i));

          arma::vec target = target_X-target_x;

          for(int ll=0; ll<idx_uvec.size(); ++ll){
            if(target(ll)>=0){idx_zero_one(i, ll)=1;}else{idx_zero_one(i, ll)=0;}
          }

        }
        arma::vec idx_prob = trans(mean(idx_zero_one));

        depth_prob_mat.row(i0) = trans(idx_prob);

        double idx_prob_min = min(idx_prob);
        depth_vec(i0) = idx_prob_min;

        uvec temp_min_idx_uvec = find(idx_prob == idx_prob_min);
        idx_vCoef_min_field_inner(i0) = temp_min_idx_uvec;

      }

      depth_prob_field(ii_J, ii_ld) = depth_prob_mat;
      // depth_field(ii_J, ii_ld) = depth_vec;
      idx_vCoef_min_field(ii_J, ii_ld) = idx_vCoef_min_field_inner;

      // rank

      NumericVector rk_nv = rankk(depth_vec, Named("ties.method", ties));
      arma::vec rk_vec = as<arma::vec>(rk_nv);
      // rank_field(ii_J, ii_ld) = rk_vec;
      for(int i0=0; i0<n0; ++i0){
        depth_cube(ii_J, ii_ld, i0) = depth_vec(i0);
        rank_cube(ii_J, ii_ld, i0) = rk_vec(i0);
      }


    }
  }

  List misc = List::create(Named("idx_acpt") = idx_acpt_field,
                           Named("acpt_rate") = acpt_rate_mat,
                           Named("M_total") = M_total_vec,
                           Named("M_used") = M_used_mat,
                           Named("vCoef") = vCoef_field,
                           Named("normRKHS") = normRKHS_field,
                           Named("ld_mat") = ld_mat,
                           Named("xCoef") = xCoef_field,
                           Named("XCeof") = XCoef_field,
                           Named("proj_x") = proj_x_field,
                           Named("proj_X") = proj_X_field,
                           Named("depth_prob_each") = depth_prob_field,
                           Named("idx_vCoef_min") = idx_vCoef_min_field);

  // List out = List::create(Named("idx_out") = idx_out_field,
  //                         Named("x_out") = x_out_field);

  List res = List::create(Named("depth") = depth_cube,
                          Named("rank") = rank_cube,
                          Named("misc") = misc);
  return(res);
}
















// =============================================================
// sim1: non-degeneracy for Gaussian case

//[[Rcpp::export]]
List MC1nd(int nX, int nY, arma::vec& tGrid, arma::vec& ev, arma::mat& ef,
           arma::mat& ld_mat, arma::vec& J_vec,
           arma::vec& M_vec, int M_type=0,
           Rcpp::CharacterVector ties="min"){
  // int tt = 100;
  // arma::vec tGrid = linspace(0, 1, tt);
  int tt = tGrid.size();
  arma::mat X=Xsim(nX, ev, ef, 0);
  arma::mat Y=Xsim(nY, ev, ef, 0);

  List resFPCA = FPCA(X, tGrid);
  arma::vec gaHat = resFPCA["ev"];
  arma::mat phiHat = resFPCA["ef"];

  List resRHD = RHDonly(Y, X, tGrid, gaHat, phiHat,
                        ld_mat, J_vec, M_vec, M_type, ties);

  return(resRHD);
}

//[[Rcpp::export]]
List MC1nd_prob(int nX, int nY, arma::vec& tGrid, arma::vec& ev, arma::mat& ef,
            arma::vec& prob_vec, arma::vec& J_vec,
            arma::vec& M_vec, int M_type=0,
            Rcpp::CharacterVector ties="min"){
  // int tt = 100;
  // arma::vec tGrid = linspace(0, 1, tt);
  int tt = tGrid.size();
  arma::mat X=Xsim(nX, ev, ef, 0);
  arma::mat Y=Xsim(nY, ev, ef, 0);

  List resFPCA = FPCA(X, tGrid);
  arma::vec gaHat = resFPCA["ev"];
  arma::mat phiHat = resFPCA["ef"];

  List resRKHS = RKHS(gaHat, J_vec);
  arma::mat normRKHS = resRKHS["normRKHS"];
  arma::mat ld_mat = trans(quantile(normRKHS, prob_vec));

  for(int ii_pr=0; ii_pr<prob_vec.size(); ++ii_pr){
    if(prob_vec(ii_pr)==1){
      uvec idx_uvec = find(gaHat>0);
      arma::vec gaHatPlus = gaHat(idx_uvec);
      double ldMax = 1/sqrt(gaHatPlus.min());
      arma::vec temp(J_vec.size(), fill::ones);
      ld_mat.col(ii_pr) = temp*ldMax;
    }
  }

  List resRHD = RHDonly(Y, X, tGrid, gaHat, phiHat,
                        ld_mat, J_vec, M_vec, M_type, ties);

  return(resRHD);
}

// sim2: depth and ranking values for one outlier

//[[Rcpp::export]]
List MC2one(int n, int n_out, arma::vec& tGrid, arma::vec& ev, arma::mat& ef,
            arma::mat& ld_mat, arma::vec& J_vec,
            arma::vec& M_vec, int M_type=0,
            Rcpp::CharacterVector ties="min"){
  int n_in = n-n_out;

  arma::mat Xin = Xsim(n_in, ev, ef, 3, 1);

  double eta = 3*sqrt(sum(ev));
  cube Xout = Xout_all(n_out, tGrid, eta);
  int out_size = Xout.n_slices / n_out;
  field<cube> depth_field(out_size);
  field<cube> rank_field(out_size);
  for(int l=0; l<out_size; ++l){
    arma::mat X = join_cols(Xin, Xout.slice(l));

    List resFPCA = FPCA(X, tGrid);
    arma::vec gaHat = resFPCA["ev"];
    arma::mat phiHat = resFPCA["ef"];

    List resRKHS = RKHS(gaHat, J_vec);
    arma::mat normRKHS = resRKHS["normRKHS"];
    List resRHD = RHDonly(X, X, tGrid, gaHat, phiHat,
                          ld_mat, J_vec, M_vec, M_type, ties);
    cube depth_cube = resRHD["depth"];
    cube rank_cube = resRHD["rank"];

    depth_field(l) = depth_cube;
    rank_field(l) = rank_cube;
  }

  List res = List::create(Named("depth") = depth_field,
                          Named("rank") = rank_field);

  return(res);
}

//[[Rcpp::export]]
List MC2one_prob(int n, int n_out, arma::vec& tGrid, arma::vec& ev, arma::mat& ef,
         arma::vec& prob_vec, arma::vec& J_vec,
         arma::vec& M_vec, int M_type=0,
         Rcpp::CharacterVector ties="min"){
  int n_in = n-n_out;

  arma::mat Xin = Xsim(n_in, ev, ef, 3, 1);

  double eta = 3*sqrt(sum(ev));
  cube Xout = Xout_all(n_out, tGrid, eta);
  int out_size = Xout.n_slices / n_out;
  field<cube> depth_field(out_size);
  field<cube> rank_field(out_size);
  for(int l=0; l<out_size; ++l){
    arma::mat X = join_cols(Xin, Xout.slice(l));

    List resFPCA = FPCA(X, tGrid);
    arma::vec gaHat = resFPCA["ev"];
    arma::mat phiHat = resFPCA["ef"];

    List resRKHS = RKHS(gaHat, J_vec);
    arma::mat normRKHS = resRKHS["normRKHS"];
    arma::mat ld_mat = trans(quantile(normRKHS, prob_vec));

    for(int ii_pr=0; ii_pr<prob_vec.size(); ++ii_pr){
      if(prob_vec(ii_pr)==1){
        uvec idx_uvec = find(gaHat>0);
        arma::vec gaHatPlus = gaHat(idx_uvec);
        double ldMax = 1/sqrt(gaHatPlus.min());
        arma::vec temp(J_vec.size(), fill::ones);
        ld_mat.col(ii_pr) = temp*ldMax;
      }
    }

    List resRHD = RHDonly(X, X, tGrid, gaHat, phiHat,
                          ld_mat, J_vec, M_vec, M_type, ties);
    cube depth_cube = resRHD["depth"];
    cube rank_cube = resRHD["rank"];

    depth_field(l) = depth_cube;
    rank_field(l) = rank_cube;
  }

  List res = List::create(Named("depth") = depth_field,
                          Named("rank") = rank_field);

  return(res);
}

// sim3: the ratios of correctly/fasely detected outliers

// sim based choice of f


//[[Rcpp::export]]
cube Fadj_inner(int nSim, arma::mat& covMat, arma::vec& tGrid,
                arma::mat& ld_mat, arma::vec& J_vec, arma::vec& f_iqr_vec,
                arma::vec& M_vec, int M_type=0){

  int tt = tGrid.size();
  int J_size = J_vec.size();
  int ld_size=ld_mat.n_cols;
  int f_iqr_size = f_iqr_vec.size();

  // arma::mat GaHatCor = cor(X);
  arma::vec mu0(tt, fill::zeros);
  arma::mat XsimGau = trans(mvnrnd(mu0, covMat, nSim));

  List resFPCA = FPCA(XsimGau, tGrid);
  arma::vec gaHat = resFPCA["ev"];
  arma::mat phiHat = resFPCA["ef"];

  List resRHD = RHD_inner(XsimGau, XsimGau, tGrid, gaHat, phiHat,
                          ld_mat, J_vec, f_iqr_vec, M_vec, M_type);
  
  field<uvec> out_field = resRHD["out"];
  List misc_List = resRHD["misc"];
  field<uvec> idx_depth_min_field = misc_List["idx_depth_min"];
  
  out_field.set_size(J_size, ld_size, f_iqr_size);
  idx_depth_min_field.set_size(J_size, ld_size);
  
  cube res_cube(J_size, ld_size, f_iqr_size, fill::zeros);
  for(int ii_J = 0; ii_J<J_size; ++ii_J){
    for(int ii_ld = 0; ii_ld<ld_size; ++ii_ld){
      
      uvec idx_min_depth = idx_depth_min_field(ii_J, ii_ld);
      
      for(int ii_f = 0; ii_f<f_iqr_size; ++ii_f){
        
        uvec idx_out = out_field(ii_J, ii_ld, ii_f);
        
        uvec idx_uvec = intersect(idx_out, idx_min_depth);
        
        if(idx_uvec.size()==0){
          res_cube(ii_J, ii_ld, ii_f)=1;
        }
        
        // Rcout << "ii_J=" << ii_J << ", ii_ld=" << ii_ld <<
        //   ", ii_f=" << ii_f << "\n";
      }
    }
  }
  
  return(res_cube);

}


//[[Rcpp::export]]
cube Fadj(int M_adj, int nSim, arma::mat& covMat, arma::vec& tGrid,
          arma::mat& ld_mat, arma::vec& J_vec, arma::vec& f_iqr_vec,
          arma::vec& M_vec, int M_type=0, int num_M=100, int stt=0){

  int tt = tGrid.size();
  int J_size = J_vec.size();
  int ld_size=ld_mat.n_cols;
  int f_iqr_size = f_iqr_vec.size();

  cube res_cube(J_size, ld_size, f_iqr_size, fill::zeros);
  if(stt==1){
    Rcout << "m=" << 0 <<
      " at " << currentDateTime() << "\n";
    }
  for(int m=0; m<M_adj; ++m){
    cube res_innder = Fadj_inner(nSim, covMat, tGrid, 
                                 ld_mat, J_vec, f_iqr_vec,
                                 M_vec, M_type);
    res_cube += res_innder;

    if(stt==1){
      if( (m+1) % num_M == 0 ){
        Rcout << "m=" << m+1 <<
          " at " << currentDateTime() << "\n";
      }
    }

  }
  res_cube = res_cube/M_adj;

  return(res_cube);
}




//[[Rcpp::export]]
cube Fadj_inner_prob(int nSim, arma::mat& covMat, arma::vec& tGrid,
                arma::vec& prob_vec, arma::vec& J_vec, arma::vec& f_iqr_vec,
                arma::vec& M_vec){

  int tt = tGrid.size();
  int J_size = J_vec.size();
  int ld_size = prob_vec.size();
  int f_iqr_size = f_iqr_vec.size();

  // arma::mat GaHatCor = cor(X);
  arma::vec mu0(tt, fill::zeros);
  arma::mat XsimGau = trans(mvnrnd(mu0, covMat, nSim));

  List resFPCA = FPCA(XsimGau, tGrid);
  arma::vec gaHat = resFPCA["ev"];
  arma::mat phiHat = resFPCA["ef"];

  // List resRKHS = RKHS(gaHat, J_vec);
  // arma::mat normRKHS = resRKHS["normRKHS"];
  // arma::mat ld_mat = trans(quantile(normRKHS, prob_vec));
  // 
  // for(int ii_pr=0; ii_pr<prob_vec.size(); ++ii_pr){
  //   if(prob_vec(ii_pr)==1){
  //     uvec idx_uvec = find(gaHat>0);
  //     arma::vec gaHatPlus = gaHat(idx_uvec);
  //     double ldMax = 1/sqrt(gaHatPlus.min());
  //     arma::vec temp(J_vec.size(), fill::ones);
  //     ld_mat.col(ii_pr) = temp*ldMax;
  //   }
  // }

  List resRHD = RHD_inner_prob(
    XsimGau, XsimGau, tGrid, gaHat, phiHat,
    prob_vec, J_vec, f_iqr_vec, M_vec
    );
  
  field<uvec> out_field = resRHD["out"];
  List misc_List = resRHD["misc"];
  field<uvec> idx_depth_min_field = misc_List["idx_depth_min"];
  
  out_field.set_size(J_size, ld_size, f_iqr_size);
  idx_depth_min_field.set_size(J_size, ld_size);
  
  cube res_cube(J_size, ld_size, f_iqr_size, fill::zeros);
  for(int ii_J = 0; ii_J<J_size; ++ii_J){
    for(int ii_ld = 0; ii_ld<ld_size; ++ii_ld){
      
      uvec idx_min_depth = idx_depth_min_field(ii_J, ii_ld);
      
      for(int ii_f = 0; ii_f<f_iqr_size; ++ii_f){
        
        uvec idx_out = out_field(ii_J, ii_ld, ii_f);
        
        uvec idx_uvec = intersect(idx_out, idx_min_depth);
        
        if(idx_uvec.size()==0){
          res_cube(ii_J, ii_ld, ii_f)=1;
        }
        
        // Rcout << "ii_J=" << ii_J << ", ii_ld=" << ii_ld <<
        //   ", ii_f=" << ii_f << "\n";
      }
    }
  }
  
  return(res_cube);

}


//[[Rcpp::export]]
cube Fadj_prob(int M_adj, int nSim, arma::mat& covMat, arma::vec& tGrid,
          arma::vec& prob_vec, arma::vec& J_vec, arma::vec& f_iqr_vec,
          arma::vec& M_vec, int num_M=100){

  int tt = tGrid.size();
  int J_size = J_vec.size();
  int ld_size = prob_vec.size();
  int f_iqr_size = f_iqr_vec.size();

  cube res_cube(J_size, ld_size, f_iqr_size, fill::zeros);
  Rcout << "m=" << 0 <<
    " at " << currentDateTime() << "\n";
  for(int m=0; m<M_adj; ++m){
    cube res_innder = Fadj_inner_prob(
      nSim, covMat, tGrid, 
      prob_vec, J_vec, f_iqr_vec,
      M_vec
    );
    res_cube += res_innder;

    if( (m+1) % num_M == 0 ){
      Rcout << "m=" << m+1 <<
        " at " << currentDateTime() << "\n";
    }

  }
  res_cube = res_cube/M_adj;

  return(res_cube);
}







