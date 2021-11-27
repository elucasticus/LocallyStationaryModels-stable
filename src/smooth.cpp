/// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
/// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>
/// Under MIT license

#include "smooth.hpp"
#include  <algorithm>
#include <iostream>

using namespace cd;

double smt::smooth_value(const unsigned int &pos, const unsigned int &n) const
{
    const matrix &K = *(kernel_.get_kernel());

    double numerator = 0;
    double denominator = 0;

    for(size_t i=0; i<anchorpos->rows(); ++i)
    {
        numerator += K(pos, i) * solutions->operator()(i, n);
        denominator += K(pos, i);
    }
    if (denominator < std::numeric_limits<double>::min())
        return 0;
    return numerator/denominator;
}

double smt::smooth_value(const cd::vector &pos, const unsigned int &n) const
{
    double numerator = 0;
    double denominator = 0;

    for(size_t i=0; i<anchorpos->rows(); ++i)
    {
        numerator += kernel_(pos, anchorpos->row(i)) * solutions->operator()(i, n);
        denominator += kernel_(pos, anchorpos->row(i));
    }
    if (denominator < std::numeric_limits<double>::min())
        return 0;
    return numerator/denominator;
}

smt::smt(const cd::matrixptr solutions_, const matrixptr &anchorpos_, const cd::scalar &min_delta, const cd::scalar &max_delta, const std::string &kernel_id): 
    anchorpos(anchorpos_), solutions(solutions_), kernel_(kernel_id, min_delta)
{
    double min_error= std::numeric_limits<double>::infinity();
    optimal_delta = (max_delta-min_delta)/2;
    const unsigned int n_deltas = 1000;
    // find the optimal value of delta via cross-validation
    for (size_t i=0; i<=n_deltas; i++)
    {
        double delta = min_delta + i*(max_delta-min_delta)/n_deltas;
        // build a new kernel with bandwidth parameter equal to delta
        kernel_.build_simple_kernel(anchorpos_, delta);
        const matrix &Kk = *(kernel_.get_kernel());
        
        double error = 0;
        // find the value of the error function for the current value of delta
        #pragma omp parallel for reduction(+:error)
        for(size_t j=0; j<anchorpos->rows(); ++j)
        {
            cd::vector Kkrow = Kk.row(j);
            double predicted_value = smooth_value(j, 3);
            double real_value = solutions->operator()(j, 3);
            double weightk2= (1-Kk(j,j)/Kkrow.sum())*(1-Kk(j,j)/Kkrow.sum());

            error += (real_value - predicted_value)*(real_value - predicted_value)/weightk2;
        }
        if (error < min_error)
        {
            optimal_delta = delta;
            min_error = error;
        }
    }
    // build the final kernel with the optimal value of delta
    kernel_.build_simple_kernel(anchorpos, optimal_delta);
}

smt::smt(const cd::matrixptr solutions_, const matrixptr &anchorpos_, const double delta, const std::string &kernel_id): anchorpos(anchorpos_),  solutions(solutions_),
    kernel_(kernel_id, delta), optimal_delta(delta)
{
    kernel_.build_simple_kernel(anchorpos);
}

const cd::matrixptr smt::get_solutions() const
{
    return solutions;
}

double smt::get_optimal_delta() const
{
    return optimal_delta;
}

const cd::matrixptr smt::get_anchorpos() const 
{
    return anchorpos;
}

smt::smt(): kernel_() {};
