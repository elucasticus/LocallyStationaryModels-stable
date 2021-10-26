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

    const cd::matrixptr get_variogram() const;
    const cd::matrixptr get_denominators() const;
    const cd::matrixptr get_squaredweights() const;
    const cd::vectorptr get_x() const;
    const cd::vectorptr get_y() const;
    const cd::matrixptr get_kernel() const;
    const cd::matrixIptr get_grid() const;
    const cd::vectorptr get_normh() const;
};

#endif // LOCALLY_STATIONARY_MODELS_SAMPLEVAR