// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppEigen.h>
#include "crippadecarlo.hpp"
#include "ancora.hpp"
using namespace cd;
using namespace LBFGSpp;
using namespace std::chrono;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Rcpp::List fullmodel(const Eigen::VectorXd &y, const Eigen::MatrixXd &d, const Eigen::MatrixXd &anchorpoints, const Eigen::VectorXd &parameters, const double& epsilon, const unsigned int& n_angles, 
    const unsigned int& n_intervals, const std::string &kernel_id, const std::string &variogram_id) {
    auto start = high_resolution_clock::now(); 

    matrixptr dd = std::make_shared<matrix>(d);
    vectorptr yy = std::make_shared<vector>(y);
    matrixptr anchorpointsptr = std::make_shared<matrix>(anchorpoints);
    
    crippadecarlo CD(dd, yy, anchorpointsptr , parameters,epsilon, n_angles,n_intervals, kernel_id, variogram_id);
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
                              Rcpp::Named("ypredicted")=CD.predict_ys(anchorpoints),
                              Rcpp::Named("predictedmean")=CD.predict_means(anchorpoints),
                              Rcpp::Named("delta")=delta,
                              Rcpp::Named("epsilon")=epsilon_);   
}


// [[Rcpp::export]]
Rcpp::List predikt(const Eigen::VectorXd &y, const Eigen::MatrixXd &d, const Eigen::MatrixXd &anchorpoints, const double& epsilon, const double &delta, const Eigen::MatrixXd &solutions,
    const Eigen::MatrixXd &positions) {

    matrixptr dd = std::make_shared<matrix>(d);
    vectorptr yy = std::make_shared<vector>(y);
    matrixptr solutionsptr = std::make_shared<matrix>(solutions);
    matrixptr anchorpointsptr = std::make_shared<matrix>(anchorpoints);

    crippadecarlo CD(dd, yy, anchorpointsptr, epsilon, delta, solutionsptr, "esponenziale");

    return Rcpp::List::create(Rcpp::Named("ypredicted")=CD.predict_ys(positions),
                              Rcpp::Named("predictedmean")=CD.predict_means(positions));    
}

// [[Rcpp::export]]

Rcpp::List find_anchorpoints(const Eigen::MatrixXd &d, const unsigned int& n_cubotti) {
    auto start = high_resolution_clock::now();
    
    
    //d = stampante::caricamatrice("d.csv");
   //y = stampante::caricavettore("y.csv");
   

    matrixptr dd = std::make_shared<matrix>(d);

    ancora a(dd, n_cubotti);

    cd::matrix anchorpos = a.find_anchorpoints();
   
    
    
    return Rcpp::List::create(Rcpp::Named("anchorpoints")=anchorpos,
                            Rcpp::Named("center_x")=a.get_center().first,
                            Rcpp::Named("center_y")=a.get_center().second,
                            Rcpp::Named("width")=a.get_dimensionecubotti().first,
                            Rcpp::Named("height")=a.get_dimensionecubotti().second); 
}


// [[Rcpp::export]]
Rcpp::List rawmodel(const Eigen::VectorXd &y, const Eigen::MatrixXd &d, const Eigen::MatrixXd &anchorpoints, const Eigen::VectorXd &parameters, const double& epsilon, const unsigned int& n_angles, 
    const unsigned int& n_intervals, const std::string &kernel_id, const std::string &variogram_id) {
    auto start = high_resolution_clock::now(); 

    matrixptr dd = std::make_shared<matrix>(d);
    vectorptr yy = std::make_shared<vector>(y);
    matrixptr anchorpointsptr = std::make_shared<matrix>(anchorpoints);
    
    crippadecarlo CD(dd, yy, anchorpointsptr , parameters,epsilon, n_angles,n_intervals, kernel_id, variogram_id);
    double delta = CD.get_delta();
    double epsilon_ = CD.get_epsilon();
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
  
    return Rcpp::List::create(Rcpp::Named("anchorpoints")=anchorpoints,
                              Rcpp::Named("kernel")=*(CD.get_kernel()),
                              Rcpp::Named("grid")=*(CD.get_grid()),
                              Rcpp::Named("empiricvariogram")=*(CD.get_empiricvariogram()),
                              Rcpp::Named("solutions")=*(CD.get_solutions()),
                              Rcpp::Named("delta")=delta,
                              Rcpp::Named("epsilon")=epsilon_);   
}


// [[Rcpp::export]]
Rcpp::List fullmodelCV(const Eigen::VectorXd &y, const Eigen::MatrixXd &d, const Eigen::MatrixXd &anchorpoints, const Eigen::VectorXd &parameters, const double& epsilonmin, const double& epsilonmax, const unsigned int& nepsilons
                      , const unsigned int& n_angles, 
                     const unsigned int& n_intervals, const std::string &kernel_id, const std::string &variogram_id) {
  auto start = high_resolution_clock::now(); 
  
  matrixptr dd = std::make_shared<matrix>(d);
  vectorptr yy = std::make_shared<vector>(y);
  matrixptr anchorpointsptr = std::make_shared<matrix>(anchorpoints);
  
  crippadecarlo CD(dd, yy, anchorpointsptr , parameters, epsilonmin, epsilonmax, nepsilons, n_angles,n_intervals, kernel_id, variogram_id);
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
                            Rcpp::Named("ypredicted")=CD.predict_ys(anchorpoints),
                            Rcpp::Named("predictedmean")=CD.predict_means(anchorpoints),
                            Rcpp::Named("delta")=delta,
                            Rcpp::Named("epsilon")=epsilon_);   
}


// [[Rcpp::export]]
Rcpp::List smoothing(const Eigen::MatrixXd solutions, const Eigen::MatrixXd &anchorpoints, const double &delta, const Eigen::MatrixXd &positions)
{
    matrixptr solutionsptr = std::make_shared<matrix>(solutions);
    matrixptr anchorpointsptr = std::make_shared<matrix>(anchorpoints);
    
    smt smt_(solutionsptr, anchorpointsptr, delta);
    
    Eigen::MatrixXd result(positions.rows(), solutions.cols());
    for (size_t i=0; i<positions.rows(); ++i)
        result.row(i)=smt_.smooth_vector(positions.row(i));

    return Rcpp::List::create(Rcpp::Named("parameters")=result);
}