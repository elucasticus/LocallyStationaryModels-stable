// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppEigen.h>
#include "crippadecarlo.hpp"
using namespace cd;
using namespace LBFGSpp;
using namespace std::chrono;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Rcpp::List fullmodel(const Eigen::VectorXd &y, const Eigen::MatrixXd &d, const double& epsilon, const unsigned int& h) {
    auto start = high_resolution_clock::now();
    
    
    //d = stampante::caricamatrice("d.csv");
   //y = stampante::caricavettore("y.csv");
   
    vectorind positions({14,17,37,87,107,117,23,34});
    
    matrixptr dd = std::make_shared<matrix>(d);
    vectorptr yy = std::make_shared<vector>(y);
    

    crippadecarlo CD(dd, yy, epsilon, h, positions);
    vector newpos = dd->row(0);
    double delta = CD.get_delta();
    double epsilon_ = CD.get_epsilon();
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
   
    
    
    return Rcpp::List::create(Rcpp::Named("positions")=positions,
                              Rcpp::Named("values")=y,
                              Rcpp::Named("ypredicted")=CD.predict_y(newpos),
                              Rcpp::Named("predictedmean")=CD.predict_mean(newpos),
                              Rcpp::Named("delta")=delta,
                              Rcpp::Named("epsilon")=epsilon_);
    
    
    
 
  
    
}

