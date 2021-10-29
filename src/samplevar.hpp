/// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
/// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>
/// Under MIT license

#ifndef LOCALLY_STATIONARY_MODELS_SAMPLEVAR
#define LOCALLY_STATIONARY_MODELS_SAMPLEVAR

#include "traits.hpp"
#include "kernel.hpp"
#include "grid.hpp"

/**
 * \brief a class to build and store the empiric variogram build from data
*/
class samplevar
{
private:
    cd::matrixptr variogram = nullptr;
    cd::matrixptr denominators = nullptr;
    cd::matrixptr squaredweights = nullptr;
    kernel kernel_;
    grid grid_;
    unsigned int n_angles;
    unsigned int n_intervals;

    /**
     * \brief a "helper" function which built the squared weights needed by the optimizer
    */
    void build_squaredweights();
    
public:
    /**
	 * \brief               constructor
	 * \param kernel_id     the name of the function you want to use for the kernel
     * \param n_angles      the number of angles to be passed to the grid
     * \param n_intervals   the number of inervals to be passed to the grid
     * \param epsilon       the parameter regulating the kernel
	*/
    samplevar(const std::string &kernel_id, const unsigned int &n_angles_, const unsigned int &n_intervals_, const cd::scalar &epsilon);
    
    /**
     * \brief a default constructor for the class which calls the default constructors for both the kernel and the grid
    */
    samplevar();

    /**
     * \brief                   builds the matrix of the empiric variogram
     * \param dptr              a shared pointer to the matrix of the coordinates of the original dataset
     * \param anchorpointsptr   a shared pointer to the matrix of the coordinates of the anchorpoitns 
     * \param yptr              a shared pointer to the vector of the value of Y
    */
    void build_samplevar(const cd::matrixptr &dptr, const cd::matrixptr &anchorpointsptr, const cd::vectorptr &yptr);

    /**
     * \brief   returns a shared pointer to the sample variogram
    */
    const cd::matrixptr get_variogram() const;
    /**
     * \brief   returns a shared pointer to the matrix of the denominators
    */
    const cd::matrixptr get_denominators() const;
    /**
     * \brief   returns a shared pointers to the squaredweigths required to evaluate the function to be optimized
    */
    const cd::matrixptr get_squaredweights() const;
    /**
     * \brief   returns grid_.x
    */
    const cd::vectorptr get_x() const;
    /**
     * \brief   returns grid_.y
    */
    const cd::vectorptr get_y() const;
    /**
     * \brief   returns kernel_.k
    */
    const cd::matrixptr get_kernel() const;
    /**
     * \brief   returns grid_.g
    */
    const cd::matrixIptr get_grid() const;
    /**
     * \brief   returns grid_.normh
    */
    const cd::vectorptr get_normh() const;
};

#endif // LOCALLY_STATIONARY_MODELS_SAMPLEVAR