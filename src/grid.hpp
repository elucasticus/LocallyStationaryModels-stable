/// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
/// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>
/// Under MIT license

#ifndef LOCALLY_STATIONARY_MODELS_GRID
#define LOCALLY_STATIONARY_MODELS_GRID

#include "traits.hpp"

/**
 * \brief               this function builds a 2D-grid using a "a fette di pizza" (slices-of-pizza like) algorithm to partition the domain
 * \param d             a shared pointer to the matrix of the coordinates
 * \param n_angles      number of slices of the pizza
 * \param n_intervals   number of the pieces for each slice of the pizza
 * \param epsilon       the same epsilon regulating the kernel
*/
cd::matrixIptr Pizza(const cd::matrixptr &d, const unsigned int &n_angles, const unsigned int &n_intervals, const double &epsilon);

/**
 * \brief       allow to select between the preferred method to build the grid
 * \param id	name of the function of choice
*/
cd::gridfunction make_grid(const std::string &id);

/**
 * \brief   class to build the grid
*/
class grid
{
private:
    cd::gridfunction f; /// grid function
    cd::matrixIptr g = std::make_shared<cd::matrixI>(0,0); /// grid matrix
    cd::vectorptr normh = nullptr; /// vector with the norm of each cell of the grid (mean of the norm of all the pairs inside)
    cd::vectorptr mean_x = nullptr; /// vector with the x of each cell of the grid (mean of the x of all the pairs inside)
    cd::vectorptr mean_y = nullptr; /// vector with the y of each cell of the grid (mean of the y of all the pairs inside)
    double epsilon; /// bandwidth parameter inside

    /**
     * \brief           a "helper" function to build the vector containing the position of the centers of the cells of the grid. 
     * Each pair of coordinates is assigned to a position of the grid.
     * \param data      a shared pointer to the matrix of the coordinates
    */
    void build_normh(const cd::matrixptr &data);

public:
    /**
	 * \brief           constructor
	 * \param id 	    name of the grid function
     * \param epsilon   the same epsilon regulating the kernel
	*/
	grid(const std::string &id, const double epsilon_);

    /**
     * \brief   use the Pizza style by default
    */
    grid();

    /**
     * \brief               build the grid
     * \param d             a shared pointer to the matrix of the coordinates
     * \param n_angles      number of slices of the pizza
     * \param n_intervals   number of the pieces for each slice of the pizza
    */
    void build_grid(const cd::matrixptr &d, const unsigned int &n_angles, const unsigned int &n_intervals);

    /**
     * \brief   return a shared pointer to the grid
    */
    const cd::matrixIptr get_grid() const;
    /**
     * \brief   return a shared pointer to normh
    */
    const cd::vectorptr get_normh() const;
    /**
     * \brief   return a pointer to the vector containing the xs of the centers of the cells of the grid
    */
    const cd::vectorptr get_x() const;
    /**
     * \brief   return a pointer to the vector containing the ys of the centers of the cells of the grid
    */
    const cd::vectorptr get_y() const;
};

#endif // LOCALLY_STATIONARY_MODELS_GRID
