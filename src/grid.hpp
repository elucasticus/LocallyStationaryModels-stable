#ifndef GRID
#define GRID

#include "traits.hpp"

/**
 * \brief these function build a 2D-grid using a "a fette di pizza" (slices-of-pizza like) algorithm to partition the domain
 * \param d     a shared pointer to the matrix of the coordinates
 * \param h     a integer number proportional to the number of the slices
 * it is highly suggeted to use this partition instead of the "a-cubotti" one
*/
cd::matrixIptr Pizza(const cd::matrixptr &d, const unsigned int &n_angles, const unsigned int &n_intervals, const double &epsilon);


/**
 * \brief allow to select between the preferred method to build the grid
 * \param id	name of the function of choice
*/
template<unsigned int T>
cd::gridfunction make_grid(const std::string &id);
template<>
cd::gridfunction make_grid<2>(const std::string &id);


/**
 * \brief a class containing all the info related to the grid
*/
template<unsigned int n>
class grid
{};


/**
 * \brief the two-dimensional version of the grid class
*/
template<>
class grid<2>
{
private:
    cd::gridfunction f;
    cd::matrixIptr g = std::make_shared<cd::matrixI>(0,0);
    cd::vectorptr normh = std::make_shared<cd::vector>(0);
    cd::vectorptr x = std::make_shared<cd::vector>(0);
    cd::vectorptr y = std::make_shared<cd::vector>(0);
    double epsilon;

    /**
     * \brief a "helper" function to build the vector containing the position of the centers of the cells of the grid
     * \param data      a shared pointer to the matrix of the coordinates
    */
    void build_normh(const cd::matrixptr &data);

public:
    /**
	 * \brief constructor
	 * \param id 	name of the grid function
	*/
	grid(const std::string &id, const double epsilon_);

    /**
     * \brief use the Pizza style by default
    */
    grid();

    /**
     * \brief build the grid
     * \param d     a shared pointer to the matrix of the coordinates
     * \param h     a integer number proportional to the number of cells of the grid
    */
    void build_grid(const cd::matrixptr &d, const unsigned int &n_angles, const unsigned int &n_intervals);

    const cd::matrixIptr get_grid() const;
    const cd::vectorptr get_normh() const;
    const cd::vectorptr get_x() const;
    const cd::vectorptr get_y() const;
};






#endif // GRID