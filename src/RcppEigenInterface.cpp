// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

/// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
/// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>
/// Under MIT license

#include <RcppEigen.h>
#include "cvinterface.hpp"
#include "ancora.hpp"
using namespace cd;
using namespace LBFGSpp;
using namespace std::chrono;

// [[Rcpp::depends(RcppEigen)]]

/**
 * \brief               finds the anchorpoints given the position of the points in the initial dataset
 * \param d             a matrix with the coordinates of the points in the original dataset
 * \param n_cubotti     the number of cells per row and coloumn in the grid of the anchorpoints
*/
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

/**
 * \brief                   calculate the empiric variogram in each anchor points
 * \param y                 a vector with the values of Y for each point in the dataset d
 * \param d                 a matrix with the coordinates of the points in the original dataset
 * \param anchorpoints      a matrix with the coordinates of each anchorpoints
 * \param epsilon           the value of the parameter epsilon regulating the kernel
 * \param n_angles          the number of the angles for the grid
 * \param n_intervals       the number of intervals for the grid
 * \param kernel_id         the type of kernel to be used
*/
// [[Rcpp::export]]
Rcpp::List variogramlsm(const Eigen::VectorXd &y, const Eigen::MatrixXd &d, const Eigen::MatrixXd &anchorpoints, const double& epsilon, const unsigned int& n_angles, 
    const unsigned int& n_intervals, const std::string &kernel_id) {
    
    auto start = high_resolution_clock::now();
  
    matrixptr dd = std::make_shared<matrix>(d);
    vectorptr yy = std::make_shared<vector>(y);
    matrixptr anchorpointsptr = std::make_shared<matrix>(anchorpoints);

    samplevar samplevar_(kernel_id, n_angles, n_intervals, epsilon);
    samplevar_.build_samplevar(dd, anchorpointsptr, yy);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    
    Rcpp::Rcout << "task successfully completed in " << duration.count() << "ms" << std::endl;

    return Rcpp::List::create(Rcpp::Named("kernel")=*(samplevar_.get_kernel()),
                              Rcpp::Named("grid")=*(samplevar_.get_grid()),
                              Rcpp::Named("mean.x")=*(samplevar_.get_x()),
                              Rcpp::Named("mean.y")=*(samplevar_.get_y()),
                              Rcpp::Named("squaredweigths")=*(samplevar_.get_squaredweights()),
                              Rcpp::Named("empiricvariogram")=*(samplevar_.get_variogram()),
                              Rcpp::Named("anchorpoints")=anchorpoints,
                              Rcpp::Named("epsilon")=epsilon
                              );
}

/**
 * \brief                       for each anchorpoints solves a problem of nonlinear optimization and returns the results
 * \param anchorpoints          a matrix with the coordinates of each anchorpoints
 * \param empiricvariogram      the empiric variogram returned by the previour function
 * \param squaredweights        squared weigths returned by the previous function
 * \param x                     mean.x returned by the previous function
 * \param y                     mean.y returned by the previous function
 * \param variogram_id          the variogram to be used
 * \param parameters            the starting position to be given to the optimizer
 * \param epsilon               the value of epsilon regulating the kernel
*/
//[[Rcpp::export]]
Rcpp::List findsolutionslsm(const Eigen::MatrixXd &anchorpoints, const Eigen::MatrixXd &empiricvariogram, const Eigen::MatrixXd &squaredweights, const Eigen::VectorXd &x, const Eigen::VectorXd &y, std::string &variogram_id,
    const std::string &kernel_id, const Eigen::VectorXd &parameters, const double &epsilon) {
    
    auto start = high_resolution_clock::now();
  
    matrixptr empiricvariogramptr = std::make_shared<matrix>(empiricvariogram);
    matrixptr squaredweightsptr = std::make_shared<matrix>(squaredweights);
    vectorptr xptr = std::make_shared<vector>(x);
    vectorptr yptr = std::make_shared<vector>(y);
    matrixptr anchorpointsptr = std::make_shared<matrix>(anchorpoints);
    
    opt opt_(empiricvariogramptr, squaredweightsptr, xptr,  yptr, variogram_id, parameters);
    opt_.findallsolutions();

    smt smt_(opt_.get_solutions(), anchorpointsptr, epsilon/10, epsilon/2, kernel_id);

    double delta_ottimale = smt_.get_optimal_delta();

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    
    Rcpp::Rcout << "task successfully completed in " << duration.count() << "ms" << std::endl;

    return Rcpp::List::create(Rcpp::Named("solutions")=*(opt_.get_solutions()),
                              Rcpp::Named("delta")=delta_ottimale,
                              Rcpp::Named("epsilon")=epsilon,
                              Rcpp::Named("anchorpoints")=anchorpoints
                              );    
}

/**
 * \brief                   predict the mean value and the punctual value of Yù
 * \param y                 a vector with the values of Y for each point in the dataset d
 * \param d                 a matrix with the coordinates of the points in the original dataset
 * \param anchorpoints      a matrix with the coordinates of each anchorpoints
 * \param epsilon           epsilon regulating the kernel
 * \param delta             delta regulating the smoothing
 * \param solutions         the solution of the nonlinear optimization problem returned by the previous function
 * \param positions         the position in which to perform the kriging
 * \param variogram_id      the variogram to be used
 * \param kernel_id         the kernel to be used inside the smoother
*/
// [[Rcpp::export]]
Rcpp::List predikt(const Eigen::VectorXd &y, const Eigen::MatrixXd &d, const Eigen::MatrixXd &anchorpoints, const double& epsilon, const double &delta, const Eigen::MatrixXd &solutions,
    const Eigen::MatrixXd &positions, const std::string &variogram_id, const std::string &kernel_id) {

    auto start = high_resolution_clock::now();
  
    matrixptr dd = std::make_shared<matrix>(d);
    vectorptr yy = std::make_shared<vector>(y);
    matrixptr solutionsptr = std::make_shared<matrix>(solutions);
    matrixptr anchorpointsptr = std::make_shared<matrix>(anchorpoints);

    smt smt_(solutionsptr, anchorpointsptr, delta, kernel_id);
    predictor predictor_(variogram_id, yy, smt_, epsilon, dd);
    
    vector predicted_ys(predictor_.predict_y<cd::matrix, cd::vector>(positions));
    vector predicted_means(predictor_.predict_mean<cd::matrix, cd::vector>(positions));
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    
    Rcpp::Rcout << predicted_ys.size() << " pairs of values predicted in " << duration.count() << " ms" << std::endl;

    return Rcpp::List::create(Rcpp::Named("ypredicted")=predicted_ys,
                              Rcpp::Named("predictedmean")=predicted_means);    
}

/**
 * \brief                find the value of the parameters regulating the variogram
 * \param solutions      the solution of the nonlinear optimization problem returned by the previous function
 * \param anchorpoints   the coordinates of the anchorpoints in which the optimization problem has been solved
 * \param delta          the value of delta regulating the smoothing
 * \param positions      where to smooth the parameters
 * \param kernel_id      the kernel to be used inside the smoother
*/
// [[Rcpp::export]]
Rcpp::List smoothing(const Eigen::MatrixXd solutions, const Eigen::MatrixXd &anchorpoints, const double &delta, const Eigen::MatrixXd &positions, const std::string &kernel_id)
{
    matrixptr solutionsptr = std::make_shared<matrix>(solutions);
    matrixptr anchorpointsptr = std::make_shared<matrix>(anchorpoints);
    
    smt smt_(solutionsptr, anchorpointsptr, delta, kernel_id);
    
    Eigen::MatrixXd result(positions.rows(), solutions.cols());
    for (size_t i=0; i<positions.rows(); ++i)
        result.row(i)=smt_.smooth_vector(positions.row(i));

    return Rcpp::List::create(Rcpp::Named("parameters")=result);
}

// DA SISTEMARE E/O COMMENTARE
// [[Rcpp::export]]
Rcpp::List fullmodelCV(const Eigen::VectorXd &y, const Eigen::MatrixXd &d, const Eigen::MatrixXd &anchorpoints, const Eigen::VectorXd &parameters, const double& epsilonmin, const double& epsilonmax, const unsigned int& nepsilons
                      , const unsigned int& n_angles, 
                     const unsigned int& n_intervals, const std::string &kernel_id, const std::string &variogram_id) {
  Rcpp::Rcout << "WARNING: the execution of this function may require a lot of time, proceed only at your own risk" << std::endl;
  auto start = high_resolution_clock::now(); 
  
  matrixptr dd = std::make_shared<matrix>(d);
  vectorptr yy = std::make_shared<vector>(y);
  matrixptr anchorpointsptr = std::make_shared<matrix>(anchorpoints);
  
  cvinterface CD(dd, yy, anchorpointsptr , parameters, epsilonmin, epsilonmax, nepsilons, n_angles,n_intervals, kernel_id, variogram_id);
  double delta = CD.get_delta();
  double epsilon_ = CD.get_epsilon();
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<milliseconds>(stop - start);
  
  Rcpp::Rcout << "task successfully completed in " << duration.count() << "ms" << std::endl;
  
  return Rcpp::List::create(Rcpp::Named("anchorpoints")=anchorpoints,
                            Rcpp::Named("values")=y,
                            Rcpp::Named("kernel")=*(CD.get_kernel()),
                            Rcpp::Named("grid")=*(CD.get_grid()),
                            Rcpp::Named("empiricvariogram")=*(CD.get_empiricvariogram()),
                            Rcpp::Named("solutions")=*(CD.get_solutions()),
                            Rcpp::Named("ypredicted")=CD.predict_y<cd::matrix, cd::vector>(anchorpoints),
                            Rcpp::Named("predictedmean")=CD.predict_mean<cd::matrix, cd::vector>(anchorpoints),
                            Rcpp::Named("delta")=delta,
                            Rcpp::Named("epsilon")=epsilon_);   
}
