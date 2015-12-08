// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]

arma::mat alphpow(double x, arma::mat mat){
  // construct matrix of powers of x
  arma::mat alphpowmat = arma::exp(log(x) * mat);
  return(alphpowmat);
}


// [[Rcpp::export]]

List cmat(Rcpp::NumericVector ctimes, double alpha, Rcpp::String corrmod, 
                   Rcpp::String diffmeth, double h){
    
   // input data
   arma::vec times = Rcpp::as<arma::vec>(ctimes);
    
   // sparse matrices 
   unsigned int ntimes = times.n_elem;
   arma::sp_mat icmat(ntimes,ntimes), gicmat(ntimes,ntimes), ggicmat(ntimes,ntimes); 

   if(ntimes > 1){

   // two or more time points
   double lalpha, ualpha, lphi, uphi;
   lalpha = alpha; ualpha = alpha; lphi = alpha; uphi = alpha;

    if(corrmod == "ar1"){
     
     // construct icmat
     arma::sp_mat mlag1(ntimes,ntimes), mlag2(ntimes,ntimes);
     arma::mat diff1(1, ntimes), diff2(1, ntimes), lag1, lag2, a1, a2;
     arma::mat la1, ua1, la2, ua2;
     mlag1.diag(0) = -1 * arma::ones(ntimes);
     mlag1.diag(-1) = arma::ones(ntimes - 1);
     diff1 = times.t() * mlag1;
     lag1 = diff1.cols(0, ntimes - 2);
     a1 = pow(1 - alphpow(alpha, 2 * lag1), -1);
     if(diffmeth == "numeric"){
      lphi = log(alpha) - log(1 - alpha) - h;
      uphi = log(alpha) - log(1 - alpha) + h;
      lalpha = exp(lphi) / (1 + exp(lphi));
      ualpha = exp(uphi) / (1 + exp(uphi));
      la1 = pow(1 - alphpow(lalpha, 2 * lag1), -1);
      ua1 = pow(1 - alphpow(ualpha, 2 * lag1), -1);
     }
     if(ntimes > 2){
      mlag2.diag(0) = -1 * arma::ones(ntimes);
      mlag2.diag(-2) = arma::ones(ntimes - 2);
      diff2 = times.t() * mlag2;
      lag2 = diff2.cols(0, ntimes - 3);
      a2 = pow(1 - alphpow(alpha, 2 * lag2), -1);
      if(diffmeth == "numeric"){
       la2 = pow(1 - alphpow(lalpha, 2 * lag2), -1);
       ua2 = pow(1 - alphpow(ualpha, 2 * lag2), -1);
      }
     } else {
      a2.reset();
      if(diffmeth == "numeric"){
       la2.reset();
       ua2.reset();
      }
     }
     arma::mat diagnum = arma::join_rows(arma::join_rows(arma::ones(1), pow(a2, -1)), arma::ones(1));
     arma::mat diagdenom1 = arma::trans(arma::repmat(pow(a1.t(), -1), 1, 2));
     diagdenom1.reshape(1, 2 *(ntimes - 1));
     arma::mat diagdenom = arma::join_rows(arma::join_rows(arma::ones(1), diagdenom1), arma::ones(1));
     diagdenom.reshape(2, ntimes);
     icmat.diag(0) = diagnum / prod(diagdenom, 0);
     icmat.diag(1) = - a1 % alphpow(alpha, lag1);
     icmat.diag(-1) = - a1 % alphpow(alpha, lag1);


     if(diffmeth == "analytic"){

      // construct gicmat
      arma::mat gicmatodiag = - pow(a1, 2) % lag1 % alphpow(alpha, lag1 - 1) % (1 + alphpow(alpha, 2 * lag1));
      arma::mat gicmatediag = 2 * pow(a1, 2) % lag1 % alphpow(alpha, 2 * lag1 - 1);
      arma::mat gicmatcdiag;
      if(ntimes > 2){
       gicmatcdiag = gicmatediag.cols(0, ntimes - 3) % (a1.cols(1, ntimes - 2) / a2)
                       - 2 * a1.cols(0, ntimes - 3) % pow(a1.cols(1, ntimes - 2), 2)
                       % alphpow(alpha, 2 * lag1.cols(1, ntimes - 2) + 2 * lag2 - 1)
                       % (lag2 % (alphpow(alpha, - 2 * lag1.cols(1, ntimes - 2)) - 1)
                       - lag1.cols(1, ntimes - 2) % (alphpow(alpha, - 2 * lag2) - 1));
      } else {
       gicmatcdiag.reset();
      }
      gicmat.diag(0) = arma::join_rows(arma::join_rows(gicmatediag.cols(0, 0), gicmatcdiag), 
                                                              gicmatediag.cols(ntimes - 2,ntimes - 2));
      gicmat.diag(1) = gicmatodiag;
      gicmat.diag(-1) = gicmatodiag;

      // construct ggicmat
      arma::mat ggicmatodiag = - pow(a1, 3) % lag1 % (alphpow(alpha, lag1 - 2)) % ((lag1 + 1) % alphpow(alpha, 4 * lag1)
                                       + 6 * lag1 % alphpow(alpha, 2 * lag1) + (lag1 - 1));
      arma::mat ggicmatediag = 2 * pow(a1, 3) % lag1 % (alphpow(alpha, 2 * lag1 - 2))
                                             % ((2 * lag1 - 1) + (2 * lag1 + 1) % alphpow(alpha, 2 * lag1));
      arma::mat ggicmatcdiag;
      if(ntimes > 2){
       ggicmatcdiag = ggicmatediag.cols(0, ntimes - 3) % (a1.cols(1, ntimes - 2) / a2)
                       - 4 * gicmatediag.cols(0, ntimes - 3) % pow(a1.cols(1, ntimes - 2), 2) 
                       % alphpow(alpha, 2 * lag1.cols(1, ntimes - 2) + 2 * lag2 - 1)
                       % (lag2 % (alphpow(alpha, - 2 * lag1.cols(1, ntimes - 2))- 1) 
                       - lag1.cols(1, ntimes - 2) % (alphpow(alpha, -2 * lag2) - 1))
                       - 2 * a1.cols(0, ntimes - 3) % pow(a1.cols(1, ntimes - 2), 3) 
                       % alphpow(alpha, 2 * lag1.cols(1, ntimes - 2) + 2 * lag2 - 2)
                       % (lag2 % (2 * lag2 - 1) % (1 + (4 * lag1.cols(1, ntimes - 2) / (2 * lag2 - 1) - 1) 
                       % alphpow(alpha, 2 * lag1.cols(1, ntimes - 2))) % (alphpow(alpha, - 2 * lag1.cols(1, ntimes - 2)) - 1) 
                       - lag1.cols(1, ntimes - 2) % (2 * lag1.cols(1, ntimes - 2) - 1)
                       % (1 + (4 * lag1.cols(1, ntimes - 2) / (2 * lag1.cols(1, ntimes - 2) - 1) - 1) 
                       % alphpow(alpha, 2 * lag1.cols(1, ntimes - 2))) % (alphpow(alpha, -2 * lag2) - 1));

      } else {
       ggicmatcdiag.reset();
      }
      ggicmat.diag(0) = arma::join_rows(arma::join_rows(ggicmatediag.cols(0, 0), ggicmatcdiag), 
                                                              ggicmatediag.cols(ntimes - 2, ntimes - 2));
      ggicmat.diag(1) = ggicmatodiag;
      ggicmat.diag(-1) = ggicmatodiag;

      } else if(diffmeth == "numeric"){

      // upper and lower icmat
      arma::mat ldiagnum = arma::join_rows(arma::join_rows(arma::ones(1), pow(la2, -1)), arma::ones(1));
      arma::mat ldiagdenom1 = arma::trans(arma::repmat(pow(la1.t(), -1), 1, 2));
      ldiagdenom1.reshape(1, 2 *(ntimes - 1));
      arma::mat ldiagdenom = arma::join_rows(arma::join_rows(arma::ones(1), ldiagdenom1), arma::ones(1));
      ldiagdenom.reshape(2, ntimes);
      arma::mat udiagnum = arma::join_rows(arma::join_rows(arma::ones(1), pow(ua2, -1)), arma::ones(1));
      arma::mat udiagdenom1 = arma::trans(arma::repmat(pow(ua1.t(), -1), 1, 2));
      udiagdenom1.reshape(1, 2 *(ntimes - 1));
      arma::mat udiagdenom = arma::join_rows(arma::join_rows(arma::ones(1), udiagdenom1), arma::ones(1));
      udiagdenom.reshape(2, ntimes);
      gicmat.diag(0) = ldiagnum / prod(ldiagdenom, 0);
      gicmat.diag(1) = - la1 % alphpow(lalpha, lag1);
      gicmat.diag(-1) = - la1 % alphpow(lalpha, lag1); 
      ggicmat.diag(0) = udiagnum / prod(udiagdenom, 0);
      ggicmat.diag(1) = - ua1 % alphpow(ualpha, lag1);
      ggicmat.diag(-1) = - ua1 % alphpow(ualpha, lag1);

      }

     } else if(corrmod == "uniform") {
      
     // construct icmat
     double lphi = log(alpha) - log(1 - alpha) - h;
     double uphi = log(alpha) - log(1 - alpha) + h;
     double lalpha = exp(lphi) / (1 + exp(lphi));
     double ualpha = exp(uphi) / (1 + exp(uphi));
     arma::sp_mat oneone(1,1);
     oneone(0,0) = 1;
     arma::sp_mat onemat = arma::repmat(arma::repmat(oneone, ntimes, 1), 1, ntimes);
     arma::sp_mat cmat = alpha * onemat;
     cmat.diag(0) = - (1 + (ntimes - 2) * alpha) * arma::ones(ntimes);
     icmat = pow((ntimes - 1) * (alpha - 1) * (alpha + pow(ntimes - 1, -1)), -1) * cmat;

     if(diffmeth == "analytic"){

      // construct gicmat
      arma::sp_mat gcmat = - (1 + (ntimes -1) * pow(alpha, 2)) * onemat;
      gcmat.diag(0) = (alpha * (ntimes - 1) * (2 + (ntimes - 2) * alpha)) * arma::ones(ntimes);
      gicmat = pow((ntimes - 1) * (alpha - 1) * (alpha + pow(ntimes - 1, -1)), -2) * gcmat;

      // construct ggicmat
      arma::sp_mat ggcmat = (2 * alpha * (ntimes - 1) * ((ntimes - 1) * pow(alpha, 2) + 3) - 2 * (ntimes - 2)) * onemat;
      ggcmat.diag(0) = -(2 * (ntimes - 1) *(pow(alpha, 3) * (ntimes - 2) * 
                         (ntimes - 1) + 3 * pow(alpha, 2) * (ntimes - 1) + 1)) * arma::ones(ntimes);
      ggicmat = pow((ntimes - 1) * (alpha - 1) * (alpha + pow(ntimes - 1, -1)), -3) * ggcmat;

     } else if(diffmeth == "numeric"){

     // upper and lower icmat
     arma::sp_mat lcmat = (lalpha) * onemat;
     arma::sp_mat ucmat = (ualpha) * onemat;
     lcmat.diag(0) = - (1 + (ntimes - 2) * (lalpha)) * arma::ones(ntimes);
     ucmat.diag(0) = - (1 + (ntimes - 2) * (ualpha)) * arma::ones(ntimes);
     gicmat = pow((ntimes - 1) * ((lalpha) - 1) * ((lalpha) + pow(ntimes - 1, -1)), -1) * lcmat;
     ggicmat = pow((ntimes - 1) * ((ualpha) - 1) * ((ualpha) + pow(ntimes - 1, -1)), -1) * ucmat;

     }

    }

   } else {

    // one time point only

     // construct icmat, gicmat and ggicmat
     icmat(0,0) = 1;
     gicmat(0,0) = 1;
     ggicmat(0,0) = 1;

   }

  if(diffmeth == "analytic"){

  // output data
  return Rcpp::List::create(Rcpp::Named("icmat")=icmat,
                            Rcpp::Named("gicmat")=gicmat,
                            Rcpp::Named("ggicmat")=ggicmat);
  } else {

  // output data
  return Rcpp::List::create(Rcpp::Named("icmat")=icmat,
                            Rcpp::Named("licmat")=gicmat,
                            Rcpp::Named("uicmat")=ggicmat);

  }

}

