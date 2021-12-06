/// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
/// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>
/// Under MIT license

#ifndef LOCALLY_STATIONARY_MODELS_SAMPLEVAR
#define LOCALLY_STATIONARY_MODELS_SAMPLEVAR

#include "kernel.hpp"
#include "grid.hpp"
#include "traits.hpp"

namespace LocallyStationaryModels
{
/**
 * \brief a class to build and store the empiric variogram in all the anchor points
*/
class SampleVar
{
private:
    cd::matrixptr m_variogram = nullptr; /// sample variogram matrix
    cd::matrixptr m_denominators = nullptr; /// a matrix with the denominators necessary to compute the squared weights
    cd::matrixptr m_squaredweights = nullptr; /// matrix with the squared weights
    Kernel m_kernel; /// kernel
    Grid m_grid; /// grid
    unsigned int m_n_angles; /// number of angles of the grid
    unsigned int m_n_intervals; /// number of intervals per angle of the grid

    /**
     * \brief a "helper" function which built the squared weights for the wls problem needed by the optimizer
    */
    void build_squaredweights();
    
public:
    /**
	 * \brief               constructor
	 * \param kernel_id     the name of the function you want to use for the kernel
     * \param n_angles      the number of angles to be passed to the grid
     * \param n_intervals   the number of inervals to be passed to the grid
     * \param epsilon       the bandwidth parameter regulating the kernel
	*/
    SampleVar(const std::string &kernel_id, const unsigned int &n_angles, const unsigned int &n_intervals, const cd::scalar &epsilon);
    
    /**
     * \brief a default constructor for the class which calls the default constructors for both the kernel and the grid
    */
    SampleVar();

    /**
     * \brief                build the matrix of the empiric variogram
     * \param data           a shared pointer to the matrix of the coordinates of the original dataset
     * \param anchorpoints   a shared pointer to the matrix of the coordinates of the anchor poitns 
     * \param z              a shared pointer to the vector of the value of Z
    */
    void build_samplevar(const cd::matrixptr &data, const cd::matrixptr &anchorpoints, const cd::vectorptr &z);

    /**
     * \brief   return a shared pointer to the sample variogram
    */
    const cd::matrixptr get_variogram() const;
    /**
     * \brief   return a shared pointer to the matrix of the denominators
    */
    const cd::matrixptr get_denominators() const;
    /**
     * \brief   return a shared pointers to the squaredweigths required to evaluate the function to be optimized
    */
    const cd::matrixptr get_squaredweights() const;
    /**
     * \brief   return m_grid.mean_x
    */
    const cd::vectorptr get_x() const;
    /**
     * \brief   return m_grid.mean_y
    */
    const cd::vectorptr get_y() const;
    /**
     * \brief   return m_kernel.k
    */
    const cd::matrixptr get_kernel() const;
    /**
     * \brief   return m_grid.g
    */
    const cd::matrixIptr get_grid() const;
    /**
     * \brief   return m_grid.normh
    */
    const cd::vectorptr get_normh() const;
}; // class SampleVar
} // namespace LocallyStationaryModels

#endif // LOCALLY_STATIONARY_MODELS_SAMPLEVAR
