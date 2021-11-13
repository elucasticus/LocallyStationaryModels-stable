/// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
/// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>
/// Under MIT license

#include "variogramfit.hpp"
#include <cmath>
#include <cfloat>
#include <iostream>

using namespace cd;
using namespace LBFGSpp;

cd::scalar funzionedaottimizzare::operator() (const cd::vector &params)
{
    variogramfunction &gammaiso = *(gammaisoptr);
    vector w = squaredweights->row(x0);
    vector truegamma(empiricvariogram->rows());

    for (unsigned int h = 0; h < truegamma.size(); ++h)
    {
        truegamma[h] = gammaiso(params, mean_x->operator[](h), mean_y->operator[](h));
    }
    vector empiricgamma = empiricvariogram->col(x0);
    return w.dot((truegamma - empiricgamma).cwiseProduct(truegamma - empiricgamma));
}

cd::scalar funzionedaottimizzare::operator() (const cd::vector &params, vector &grad)
{
    variogramfunction &gammaiso = *(gammaisoptr);
    vector w = squaredweights->row(x0);
    vector truegamma(empiricvariogram->rows());

    for (unsigned int h = 0; h < truegamma.size(); ++h)
    {
        truegamma[h] = gammaiso(params, mean_x->operator[](h), mean_y->operator[](h));
    }
    
    for (unsigned int i=0; i<params.size(); ++i)
    {
        vector paramsdeltaplus(params);
        vector paramsdeltaminus(params);

        double increment = 10e-8*params[i];

        paramsdeltaplus[i] += increment;
        paramsdeltaminus[i] -= increment;

        grad[i] = (funzionedaottimizzare::operator()(paramsdeltaplus) - funzionedaottimizzare::operator()(paramsdeltaminus))/(2*increment);
    }
    
    vector empiricgamma = empiricvariogram->col(x0);
    return w.dot((truegamma - empiricgamma).cwiseProduct(truegamma - empiricgamma));
}

funzionedaottimizzare::funzionedaottimizzare(const cd::matrixptr empiricvariogram_, const cd::matrixptr squaredweights_, const cd::vectorptr mean_x_, const cd::vectorptr mean_y_, unsigned int x0_, 
    const std::string &id): empiricvariogram(empiricvariogram_), squaredweights(squaredweights_), mean_x(mean_x_), mean_y(mean_y_), x0(x0_), gammaisoptr(make_variogramiso(id)) {};

opt::opt(const cd::matrixptr empiricvariogram_, const cd::matrixptr squaredweights_, const cd::vectorptr mean_x_, const cd::vectorptr mean_y_, const std::string &id_, 
    const cd::vector &initialparameters_, const cd::vector &lowerbound_, const cd::vector &upperbound_): 
    empiricvariogram(empiricvariogram_), squaredweights(squaredweights_), mean_x(mean_x_), mean_y(mean_y_), id(id_), initialparameters(initialparameters_), lowerbound(lowerbound_),
    upperbound(upperbound_)
{
    solutions = std::make_shared<matrix>(matrix::Zero(empiricvariogram->cols(),initialparameters.size()));
};

vector opt::findonesolution(const unsigned int pos) const
{
    funzionedaottimizzare fun(empiricvariogram, squaredweights, mean_x,  mean_y, pos, id);

    // Set up parameters
    LBFGSBParam<double> param;
    param.epsilon = 1e-6;
    param.max_iterations = 1000000;

    // Create solver and function object
    LBFGSBSolver<double> solver(param);

    // Bounds
    Eigen::VectorXd lb(lowerbound);
    Eigen::VectorXd ub(upperbound);
    
    cd::vector x(initialparameters);
    // x will be overwritten to be the best point found
    double fx;

    try
    {
        /*int niter = */solver.minimize(fun, x, fx, lb, ub);
    } catch (std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        x = initialparameters;
    }
    return x;
}

void opt::findallsolutions()
{
    #pragma omp parallel for
    for (unsigned int i=0; i<empiricvariogram->cols(); ++i)
    {
        vector sol = findonesolution(i);
        for (unsigned int j = 0; j < initialparameters.size(); ++j)
        {
            solutions->operator()(i, j) = sol[j];
        }
    }
}

cd::matrixptr opt::get_solutions() const {return solutions;}
