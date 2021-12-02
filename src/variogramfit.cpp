/// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
/// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>
/// Under MIT license

#include "variogramfit.hpp"
#include <cmath>
#include <cfloat>
#include <iostream>

namespace LocallyStationaryModels
{
using namespace cd;
using namespace LBFGSpp;

cd::scalar FunzioneDaOttimizzare::operator() (const cd::vector &params)
{
    VariogramFunction &gammaiso = *(m_gammaisoptr);
    vector w = m_squaredweights->row(m_x0);
    vector truegamma(m_empiricvariogram->rows());

    for (unsigned int h = 0; h < truegamma.size(); ++h)
    {
        truegamma[h] = gammaiso(params, m_mean_x->operator[](h), m_mean_y->operator[](h));
    }
    vector empiricgamma = m_empiricvariogram->col(m_x0);
    return w.dot((truegamma - empiricgamma).cwiseProduct(truegamma - empiricgamma));
}

cd::scalar FunzioneDaOttimizzare::operator() (const cd::vector &params, vector &grad)
{
    VariogramFunction &gammaiso = *(m_gammaisoptr);
    vector w = m_squaredweights->row(m_x0);
    vector truegamma(m_empiricvariogram->rows());

    for (unsigned int h = 0; h < truegamma.size(); ++h)
    {
        truegamma[h] = gammaiso(params, m_mean_x->operator[](h), m_mean_y->operator[](h));
    }
    
    double c = 10e-8;
    // we update the gradient of the function
    // partial derivative are calcutated with central differences method
    // the step for the numerical estimation of the gradient is chosen proportionally to the parameter with respect to which we are
    // calculating the derivative
    for (unsigned int i=0; i<params.size(); ++i)
    {
        vector paramsdeltaplus(params);
        vector paramsdeltaminus(params);

        double increment = c*params[i];

        paramsdeltaplus[i] += increment;
        paramsdeltaminus[i] -= increment;

        grad[i] = (FunzioneDaOttimizzare::operator()(paramsdeltaplus) - FunzioneDaOttimizzare::operator()(paramsdeltaminus))/(2*increment);
    }
    
    vector empiricgamma = m_empiricvariogram->col(m_x0);
    return w.dot((truegamma - empiricgamma).cwiseProduct(truegamma - empiricgamma));
}

FunzioneDaOttimizzare::FunzioneDaOttimizzare(const cd::matrixptr empiricvariogram, const cd::matrixptr squaredweights, const cd::vectorptr mean_x, const cd::vectorptr mean_y, unsigned int x0, 
    const std::string &id): m_empiricvariogram(empiricvariogram), m_squaredweights(squaredweights), m_mean_x(mean_x), m_mean_y(mean_y), m_x0(x0), m_gammaisoptr(make_variogramiso(id)) {};

Opt::Opt(const cd::matrixptr empiricvariogram, const cd::matrixptr squaredweights, const cd::vectorptr mean_x, const cd::vectorptr mean_y, const std::string &id, 
    const cd::vector &initialparameters, const cd::vector &lowerbound, const cd::vector &upperbound): 
    m_empiricvariogram(empiricvariogram), m_squaredweights(squaredweights), m_mean_x(mean_x), m_mean_y(mean_y), m_id(id), m_initialparameters(initialparameters), m_lowerbound(lowerbound),
    m_upperbound(upperbound)
{
    m_solutions = std::make_shared<matrix>(matrix::Zero(m_empiricvariogram->cols(),m_initialparameters.size()));
};

vector Opt::findonesolution(const unsigned int pos) const
{
    FunzioneDaOttimizzare fun(m_empiricvariogram, m_squaredweights, m_mean_x,  m_mean_y, pos, m_id);

    // Set up parameters
    LBFGSBParam<double> param;
    param.epsilon = 1e-6;
    param.max_iterations = 1000000;

    // Create solver and function object
    LBFGSBSolver<double> solver(param);

    // Bounds
    Eigen::VectorXd lb(m_lowerbound);
    Eigen::VectorXd ub(m_upperbound);
    
    cd::vector x(m_initialparameters);
    // x will be overwritten to be the best point found
    double fx;

    try
    {
        /*int niter = */solver.minimize(fun, x, fx, lb, ub);
    } catch (std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        x = m_initialparameters;
    }
    return x;
}

void Opt::findallsolutions()
{
    #pragma omp parallel for
    for (unsigned int i=0; i<m_empiricvariogram->cols(); ++i)
    {
        vector sol = findonesolution(i);
        for (unsigned int j = 0; j < m_initialparameters.size(); ++j)
        {
            m_solutions->operator()(i, j) = sol[j];
        }
    }
}

cd::matrixptr Opt::get_solutions() const {return m_solutions;}
}; // namespace LocallyStationaryModels
