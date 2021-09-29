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
        if (i!=pos)
        {
            numerator += K(pos, i) * solutions->operator()(i, n);
            denominator += K(pos, i);
        }
    }
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
    return numerator/denominator;
}


// SISTEMA CROSS-VALIDATION DI SIGMA!!!!!!!!!!!!
smt::smt(const cd::matrixptr solutions_, const matrixptr &anchorpos_, const cd::scalar &min_delta, const cd::scalar &max_delta): anchorpos(anchorpos_), solutions(solutions_), kernel_()
{
    double min_error;
    const unsigned int n_deltas = 100;

    for (int i=0; i<n_deltas; i++)
    {
        double delta = min_delta + i*(max_delta-min_delta)/n_deltas;

        kernel_.build_simple_kernel(anchorpos_, delta);

        double error = 0;

        for(size_t j=0; j<anchorpos->rows(); ++j)
        {
            double predicted_value = smooth_value(j, 3);
            double real_value = solutions->operator()(j, 3);

            error += (real_value - predicted_value)*(real_value - predicted_value);
        }
        if (i==0 || error < min_error)
        {
            optimal_delta = delta;
            min_error = error;
        }
    }

    kernel_.build_simple_kernel(anchorpos, optimal_delta);
}



smt::smt(const cd::matrixptr solutions_, const matrixptr &anchorpos_, const double delta): anchorpos(anchorpos_),  solutions(solutions_),
    kernel_("gaussian", delta), optimal_delta(delta)
{
    kernel_.build_simple_kernel(anchorpos);
}


cd::vector smt::smooth_vector(const unsigned int &pos) const
{
    vector result(solutions->cols());
    for (unsigned int i=0; i<solutions->cols(); ++i)
        result(i) = smooth_value(pos, i);

    return result;
}

cd::vector smt::smooth_vector(const vector &pos) const
{
    vector result(solutions->cols());
    for (unsigned int i=0; i<solutions->cols(); ++i)
        result(i) = smooth_value(pos, i);

    return result;
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