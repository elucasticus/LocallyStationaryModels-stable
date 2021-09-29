#ifndef SAMPLEVAR
#define SAMPLEVAR

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
    grid<2> grid_;

    unsigned int n_angles;
    unsigned int n_intervals;

    /**
     * \brief a "helper" function which built the squared weights needed by the optimizer
    */
    void build_squaredweights();
    
public:
    /**
	 * \brief constructor
	 * \param kernel_id     the name of the function you want to use for the kernel
     * \param grid_id       the name of the partition you want for the grid
     * \param h             a parameter proportional to the number of the cells of the grid
     * \param epsilon       the parameter regulating the kernel
	*/
    samplevar(const std::string &kernel_id, const unsigned int &n_angles_, const unsigned int &n_intervals_, const cd::scalar &epsilon);
    
    /**
     * \brief a default constructor for the class which calls the default constructors for both the kernel and the grid
    */
    samplevar();

    /**
     * \brief builds the matrix of the empiric variogram
     * \param dptr  a shared pointer to the matrix of the coordinates
     * \param yptr  a shared pointer to the vector of the value of Y
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

#endif // SAMPLEVAR