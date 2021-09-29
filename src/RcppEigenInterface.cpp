// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppEigen.h>
#include "crippadecarlo.hpp"
using namespace cd;
using namespace LBFGSpp;
using namespace std::chrono;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Rcpp::List fullmodel(const Eigen::VectorXd &y, const Eigen::MatrixXd &d, const double& epsilon, const unsigned int& n_angles, const unsigned int& n_intervals) {
    auto start = high_resolution_clock::now();
    
    
    //d = stampante::caricamatrice("d.csv");
   //y = stampante::caricavettore("y.csv");
   

    matrixptr dd = std::make_shared<matrix>(d);
    vectorptr yy = std::make_shared<vector>(y);
    

    crippadecarlo CD(dd, yy, dd ,epsilon, n_angles,n_intervals);
    vector newpos = dd->row(0);
    double delta = CD.get_delta();
    double epsilon_ = CD.get_epsilon();
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
   
    
    
    return Rcpp::List::create(Rcpp::Named("values")=y,
                              Rcpp::Named("ypredicted")=CD.predict_y(newpos),
                              Rcpp::Named("predictedmean")=CD.predict_mean(newpos),
                              Rcpp::Named("delta")=delta,
                              Rcpp::Named("epsilon")=epsilon_);
    
    
    
 
  
    
}

