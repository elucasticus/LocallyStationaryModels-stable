/// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
/// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>
/// Under MIT license

#ifndef LOCALLY_STATIONARY_MODELS_GRADIENT
#define LOCALLY_STATIONARY_MODELS_GRADIENT

#include "traits.hpp"
#include "LBFGS/LBFGSB.h"
#include "variogramfunctions.hpp"
#include <vector>

/**
 * \brief functor to pass to the optimizer that contains the wls to be minimized
*/
struct funzionedaottimizzare
{
    const cd::matrixptr empiricvariogram;
    const cd::matrixptr squaredweights;
    const cd::vectorptr mean_x;
    const cd::vectorptr mean_y;
    unsigned int x0;
    std::unique_ptr<variogramfunction> gammaisoptr;

    /**
     * \brief                       constructor
     * \param empiricvariogram_     a shared pointer to the empiric variogram
     * \param squaredweights_       a shared pointer to the squared weights
     * \param mean_x_                    a shared pointer to the vector of the abscissas of the centers
     * \param mean_y_                    a shared pointer to the vector of the ordinates of the centers
     * \param x0_                   the index of the position x0
     * \param id                    the name of the variogram of your choice
    */
    funzionedaottimizzare(const cd::matrixptr empiricvariogram_, const cd::matrixptr squaredweights_, const cd::vectorptr mean_x_, const cd::vectorptr mean_y_, unsigned int x0_, 
    const std::string &id);

    /**
     * \param params    a vector containing the previous value of the parameters of the function (lambda1, lambda2, phi, sigma, etc.)
     * \param grad      a vector containing the previous value of the gradient which will be updated
    */
    cd::scalar operator() (const cd::vector &params, cd::vector &grad);
    cd::scalar operator() (const cd::vector &params);
};


/**
 * \brief a class to estimate the value of the parameters of the variogram in each point by optimizing the correspondent funzionedaottimizzare relying on the library LBFGSpp
*/
class opt
{
private:
    cd::matrixptr empiricvariogram;
    cd::matrixptr squaredweights;
    cd::vectorptr mean_x;
    cd::vectorptr mean_y;
    std::string id;
    cd::vector initialparameters;
    cd::vector lowerbound;
    cd::vector upperbound;
    cd::matrixptr solutions = nullptr;

    /**
     * \brief       find the optimal solution for the point in position pos
     * \param pos   the index of the position in which find the optimal solution
    */
    cd::vector findonesolution(const unsigned int pos) const;

public:
    /**
     * \brief                       constructor
     * \param empiricvariogram_     a shared pointer to the empiric variogram
     * \param squaredweights_       a shared pointer to the squared weights
     * \param mean_x_                    a shared pointer to the vector of the abscissas of the centers
     * \param mean_y_                    a shared pointer to the vector of the ordinates of the centers
     * \param id                    the name of the variogram of your choice
     * \param initialparameters_    the initial value of the parameters required from the optimizer to start the search for a minimum
     * \param lowerbound_           the lower bounds for the parameters in the nonlinear optimization problem
     * \param upperbound_           the upper bounds for the parameters in the nonlinear optimization problem
    */
    opt(const cd::matrixptr empiricvariogram_, const cd::matrixptr squaredweights_, const cd::vectorptr mean_x_, const cd::vectorptr mean_y_,
    const std::string &id_, const cd::vector &initialparameters_, const cd::vector &lowerbound_, const cd::vector &upperbound_);

    /**
     * \brief   find the optimal solution in all the position
    */
    void findallsolutions();

    /**
     * \brief   return the solutions found by solving the problem of nonlinear optimization
    */
    cd::matrixptr get_solutions() const;
};

#endif //LOCALLY_STATIONARY_MODELS_GRADIENT
