#include "smooth.hpp"
#include  <algorithm>
#include <iostream>

using namespace cd;

double smt::smooth_value(const unsigned int &pos, const unsigned int &n) const
{
    const matrix &K = *(kernel_.get_kernel());

    double numerator = 0;
    double denominator = 0;

    for(size_t i=0; i<anchorpos.size(); ++i)
    {
        if (anchorpos[i]!=pos)
        {
            numerator += K(pos, anchorpos[i]) * solutions->operator()(anchorpos[i], n);
            denominator += K(pos, anchorpos[i]);
        }
    }
    return numerator/denominator;
}

double smt::smooth_value(const cd::vector &pos, const unsigned int &n) const
{
    double numerator = 0;
    double denominator = 0;

    for(size_t i=0; i<anchorpos.size(); ++i)
    {
        numerator += kernel_(pos, data->row(anchorpos[i])) * solutions->operator()(anchorpos[i], n);
        denominator += kernel_(pos, data->row(anchorpos[i]));
    }
    return numerator/denominator;
}



smt::smt(const cd::matrixptr solutions_, const vectorind &anchorpos_, const matrixptr &d, const int min_delta, const int max_delta): anchorpos(anchorpos_), solutions(solutions_), data(d), kernel_()
{
    double min_error;

    for (int i=min_delta; i<max_delta; ++i)
    {
        double delta = -exp(i);

        kernel_.build_simple_kernel(d, delta);

        double error = 0;

        for(size_t j=0; j<anchorpos.size(); ++j)
        {
            double predicted_value = smooth_value(anchorpos[j], 3);
            double real_value = solutions->operator()(anchorpos[j], 3);

            error += (real_value - predicted_value)*(real_value - predicted_value);
        }
        if (i==min_delta || error < min_error)
        {
            optimal_delta = delta;
            min_error = error;
        }
    }

    kernel_.build_simple_kernel(d, optimal_delta);
}



smt::smt(const cd::matrixptr solutions_, const vectorind &anchorpos_, const matrixptr &d, const double delta): anchorpos(anchorpos_),  solutions(solutions_),
    data(d), kernel_("gaussian", delta), optimal_delta(delta)
{
    kernel_.build_simple_kernel(d);
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


void smt::smooth_solutions()
{
    for (unsigned int i=0; i<solutions->rows(); ++i)
    {
        if (std::find(anchorpos.begin(), anchorpos.end(), i) == anchorpos.end())
        {
            vector result = smooth_vector(i);
            solutions->row(i) = result;
        }
    }
}

const cd::matrixptr smt::get_solutions() const
{
    return solutions;
}

double smt::get_optimal_delta() const
{
    return optimal_delta;
}

const cd::matrixptr smt::get_d() const 
{
    return data;
}

smt::smt(): kernel_() {};