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
        truegamma[h] = gammaiso(params, x->operator[](h), y->operator[](h));
    }
    vector empiricgamma = empiricvariogram->col(x0);
    return w.transpose() * (truegamma - empiricgamma).cwiseProduct(truegamma - empiricgamma);
}

cd::scalar funzionedaottimizzare::operator() (const cd::vector &params, vector &grad)
{
    variogramfunction &gammaiso = *(gammaisoptr);
    vector w = squaredweights->row(x0);
    vector truegamma(empiricvariogram->rows());

    for (unsigned int h = 0; h < truegamma.size(); ++h)
    {
        truegamma[h] = gammaiso(params, x->operator[](h), y->operator[](h));
    }
    
    for (unsigned int i=0; i<params.size(); ++i)
    {
        vector paramsdeltaplus(params);
        vector paramsdeltaminus(params);

        double increment = 10e-6*params[i];

        paramsdeltaplus[i] += increment;
        paramsdeltaminus[i] -= increment;

        grad[i] = (funzionedaottimizzare::operator()(paramsdeltaplus) - funzionedaottimizzare::operator()(paramsdeltaminus))/(2*increment);
    }
    
    vector empiricgamma = empiricvariogram->col(x0);
    return w.transpose() * (truegamma - empiricgamma).cwiseProduct(truegamma - empiricgamma);
}

funzionedaottimizzare::funzionedaottimizzare(const cd::matrixptr empiricvariogram_, const cd::matrixptr squaredweights_, const cd::vectorptr x_, const cd::vectorptr y_, unsigned int x0_, 
    const std::string &id): empiricvariogram(empiricvariogram_), squaredweights(squaredweights_), x(x_), y(y_), x0(x0_), gammaisoptr(make_variogramiso(id)) {};

opt::opt(const cd::matrixptr empiricvariogram_, const cd::matrixptr squaredweights_, const cd::vectorptr x_, const cd::vectorptr y_, const std::string &id_, const cd::vector &initialparameters_): 
    empiricvariogram(empiricvariogram_), squaredweights(squaredweights_), x(x_), y(y_), id(id_), initialparameters(initialparameters_) 
{
    solutions = std::make_shared<matrix>(matrix::Zero(empiricvariogram->cols(),initialparameters.size()));
};

vector opt::findonesolution(const unsigned int pos) const
{
    funzionedaottimizzare fun(empiricvariogram, squaredweights, x,  y, pos, id);

    // Set up parameters
    LBFGSBParam<double> param;
    param.epsilon = 1e-6;
    param.max_iterations = 1000000;

    // Create solver and function object
    LBFGSBSolver<double> solver(param);

    // Bounds
    Eigen::VectorXd lb(cd::vector::Zero(initialparameters.size()));
    
    for(size_t e=0; e<lb.size(); e++)
        lb(e)=1e-08;
   
   
    auto inf = std::numeric_limits<double>::infinity();
    Eigen::VectorXd ub(cd::vector::Zero(initialparameters.size()));
    for(size_t e=0; e<ub.size(); e++)
        ub(e)=inf;
    ub(2) = M_PI_2;
    
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
    #pragma omp parallel for //RIMUOVERE IL COMMENTO FA CRASHARE R
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
