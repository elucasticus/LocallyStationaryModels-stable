/// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
/// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>
/// Under MIT license

#ifndef LOCALLY_STATIONARY_MODELS_KERNEL
#define LOCALLY_STATIONARY_MODELS_KERNEL

#include "traits.hpp"
#include <string>

/**
 * \return  e^(-norm(x-y)^2/epsilon^2)
*/
cd::scalar gaussian(const cd::vector &x, const cd::vector &y, const cd::scalar &epsilon);

/**
 * \brief       allow to select between pre-built kernel functions
 * \param id	name of the function of choice
*/
cd::kernelfunction make_kernel(const std::string &id);

/**
 * \brief   class to compute the kernel matrix
*/
class kernel
{
private:
    cd::scalar epsilon;
    cd::kernelfunction f;
    cd::matrixptr k = std::make_shared<cd::matrix>(0,0);

public:
    /**
	 * \brief               constructor
	 * \param id 	        name of the kernel function
     * \param epsilon_      value of epsilon
	*/
	kernel(const std::string &id, const cd::scalar &epsilon_);

    /** 
     * \brief   default constuctor with f gaussian and epsilon equal to -4.162009e-07
    */
    kernel();

    /**
	 * \return  f(x ,y) where f is the kernel function
	*/
    cd::scalar operator()(const cd::vector &x, const cd::vector &y) const;

    /**
     * \brief                   build the "star" version of the kernel needed for the empiric variogram
     * \param coordinates       a shared pointer to the matrix with the coordinates of the original dataset
     * \param anchorpoints      a shared pointer to the matrix with the coordinates of the anchorpoints
    */
    void build_kernel(const cd::matrixptr &coordinates, const cd::matrixptr &anchorpoints);

    /**
     * \brief               these functions build the "standard" version of the kernel needed for smoothing
     * \param coordinates   a shared pointer to the matrix with all coordinates
     * \param epsilon_      if needed, replace the old epsilon with a new value
    */
    void build_simple_kernel(const cd::matrixptr &coordinates);
    void build_simple_kernel(const cd::matrixptr &coordinates, const cd::scalar &epsilon_);

    /**
     * \brief   return a shared pointer to the matrix pointed by k
    */
    const cd::matrixptr get_kernel() const;
};

#endif // LOCALLY_STATIONARY_MODELS_KERNEL