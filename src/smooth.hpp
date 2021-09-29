#ifndef SMOOTH
#define SMOOTH

#include "kernel.hpp"
#include "traits.hpp"
#include <vector>


/**
 * \brief a class to smooth the solutions of the optimization and obtain cooler results
*/
class smt
{
private:
    cd::matrixptr solutions = nullptr;
    cd::matrixptr anchorpos = nullptr;

    kernel kernel_;

    double optimal_delta = 0;

    /**
     * \brief smooth a single parameter for a point in position pos
     * \param pos   a vector of coordinates or the index of the position of the point where to find the smoothed value of the parameter
     * \param n     the index of the parameter to obtain
    */
    double smooth_value(const unsigned int &pos, const unsigned int &n) const;
    double smooth_value(const cd::vector &pos, const unsigned int &n) const;

public:
    /**
     * \brief constructor
     * \param solutions_    a shared pointer to the solutions of the optimization
     * \param anchorpos_    a vector containing the indeces of the anchor position obtained by clustering
     * \param d             a shared pointer to the matrix of the coordinates
     * \param min_delta     the minimum exponent for the cross-validation of delta
     * \param max_delta     the maximum exponent for the cross-validation of delta
    */
    smt(const cd::matrixptr solutions_, const cd::matrixptr &anchorpos_, const cd::scalar &min_delta, const cd::scalar &max_delta);
    /**
     * \brief constructor
     * \param solutions_    a shared pointer to the solutions of the optimization
     * \param anchorpos_    a vector containing the indeces of the anchor position obtained by clustering
     * \param d             a shared pointer to the matrix of the coordinates
     * \param delta         a user-chosen value for delta
    */
    smt(const cd::matrixptr solutions_, const cd::matrixptr &anchorpos_, const double delta);
    /**
     * \brief calls the default constructor for kernel_
    */
    smt();

    /**
     * \brief smooth all the parameters for a point in position pos
     * \param pos   a vector of coordinates or the index of the position of the point where to find the smoothed value of the parameters
    */
    cd::vector smooth_vector(const unsigned int &pos) const;
    cd::vector smooth_vector(const cd::vector &pos) const;

    const cd::matrixptr get_solutions() const;

    double get_optimal_delta() const;

    const cd::matrixptr get_anchorpos() const;
};


#endif //SMOOTH