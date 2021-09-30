// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppEigen.h>
#include "crippadecarlo.hpp"
using namespace cd;
using namespace LBFGSpp;
using namespace std::chrono;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Rcpp::List fullmodel(const Eigen::VectorXd &y, const Eigen::MatrixXd &d, const Eigen::MatrixXd &anchorpoints, const double& epsilon, const unsigned int& n_angles, const unsigned int& n_intervals) {
    auto start = high_resolution_clock::now();
    
    
    //d = stampante::caricamatrice("d.csv");
   //y = stampante::caricavettore("y.csv");
   

    matrixptr dd = std::make_shared<matrix>(d);
    vectorptr yy = std::make_shared<vector>(y);
    matrixptr anchorpointsptr = std::make_shared<matrix>(anchorpoints);
    
    

    crippadecarlo CD(dd, yy, anchorpointsptr ,epsilon, n_angles,n_intervals);
    vector newpos = dd->row(0);
    double delta = CD.get_delta();
    double epsilon_ = CD.get_epsilon();
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
   
    
    
    return Rcpp::List::create(Rcpp::Named("anchorpoints")=anchorpoints,
                              Rcpp::Named("values")=y,
                              Rcpp::Named("kernel")=*(CD.get_kernel()),
                              Rcpp::Named("grid")=*(CD.get_grid()),
                              Rcpp::Named("empiricvariogram")=*(CD.get_empiricvariogram()),
                              Rcpp::Named("solutions")=*(CD.get_solutions()),
                              Rcpp::Named("ypredicted")=CD.predict_y(newpos),
                              Rcpp::Named("predictedmean")=CD.predict_mean(newpos),
                              Rcpp::Named("delta")=delta,
                              Rcpp::Named("epsilon")=epsilon_);
    
    
    
 
  
    
}

