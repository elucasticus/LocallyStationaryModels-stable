/// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
/// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>
/// Under MIT license

#ifndef LOCALLY_STATIONARY_MODELS_GRADIENT
#define LOCALLY_STATIONARY_MODELS_GRADIENT

#include "LBFGS/LBFGSB.h"
#include "traits.hpp"
#include "variogramfunctions.hpp"

namespace LocallyStationaryModels
{
/**
 * \brief functor to pass to the optimizer that contains the wls to be minimized
*/
struct FunzioneDaOttimizzare
{
    const cd::matrixptr m_empiricvariogram; ///< sample variogram matrix
    const cd::matrixptr m_squaredweights; ///< matrix of the squared weights
    const cd::vectorptr m_mean_x; ///< vector with the x of each cell of the grid (mean of the x of all the pairs inside)
    const cd::vectorptr m_mean_y; ///< vector with the y of each cell of the grid (mean of the y of all the pairs inside)
    size_t m_x0; ///< index of the position where to evaluate gammaisoptr
    std::shared_ptr<VariogramFunction> m_gammaisoptr; ///< pointer to the variogram function

    /**
     * \brief                       constructor
     * \param empiricvariogram      a shared pointer to the empiric variogram
     * \param squaredweights        a shared pointer to the squared weights
     * \param mean_x                a shared pointer to the vector of the abscissas of the centers
     * \param mean_y                a shared pointer to the vector of the ordinates of the centers
     * \param x0                    the index of the position x0
     * \param id                    the name of the variogram of your choice
    */
    FunzioneDaOttimizzare(const cd::matrixptr &empiricvariogram, const cd::matrixptr &squaredweights, const cd::vectorptr &mean_x, 
    const cd::vectorptr &mean_y, const size_t &x0, const std::string &id);

    /**
     * \param params    a vector containing the previous value of the parameters of the function (lambda1, lambda2, phi, sigma, etc.)
     * \param grad      a vector containing the previous value of the gradient which is updated at each iteration
    */
    double operator() (const cd::vector &params, cd::vector &grad);
    double operator() (const cd::vector &params);
}; // struct FunzioneDaOttimizzare


/**
 * \brief a class to estimate the value of the parameters of the variogram in each point by optimizing the correspondent funzionedaottimizzare relying on the library LBFGSpp
*/
class Opt
{
private:
    cd::matrixptr m_empiricvariogram; ///< sample variogram matrix
    cd::matrixptr m_squaredweights; ///< matrix with the squared weights
    cd::vectorptr m_mean_x; ///< vector with the x of each cell of the grid (mean of the x of all the pairs inside)
    cd::vectorptr m_mean_y; ///< vector with the y of each cell of the grid (mean of the y of all the pairs inside)
    std::string m_id; ///< name of the chosen variogram
    cd::vector m_initialparameters; ///< initial parameters for the optimizer
    cd::vector m_lowerbound; ///< lower bounds for the optimizer
    cd::vector m_upperbound; ///< upper bounds for the optimizer
    cd::matrixptr m_solutions = nullptr; ///< matrix with the solution in all the anchor points

    /**
     * \brief       find the optimal solution for the point in position pos
     * \param pos   the index of the position in which find the optimal solution
    */
    cd::vector findonesolution(const size_t &pos) const;

public:
    /**
     * \brief                     constructor
     * \param empiricvariogram    a shared pointer to the empiric variogram
     * \param squaredweights      a shared pointer to the squared weights
     * \param mean_x              a shared pointer to the vector of the abscissas of the centers
     * \param mean_y              a shared pointer to the vector of the ordinates of the centers
     * \param id                  the name of the variogram of your choice
     * \param initialparameters   the initial value of the parameters required from the optimizer to start the search for a minimum
     * \param lowerbound          the lower bounds for the parameters in the nonlinear optimization problem
     * \param upperbound          the upper bounds for the parameters in the nonlinear optimization problem
    */
    Opt(const cd::matrixptr &empiricvariogram, const cd::matrixptr &squaredweights, const cd::vectorptr &mean_x, 
    const cd::vectorptr &mean_y, const std::string &id, const cd::vector &initialparameters, const cd::vector &lowerbound, 
    const cd::vector &upperbound);

    /**
     * \brief   find the optimal solution in all the position
    */
    void findallsolutions();

    /**
     * \return  the solutions found by solving the problem of nonlinear optimization
    */
    cd::matrixptr get_solutions() const;
}; // class Opt
} // namespace LocallyStationaryModels

#endif //LOCALLY_STATIONARY_MODELS_GRADIENT
