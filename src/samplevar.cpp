#include "samplevar.hpp"
#include <iostream>

using namespace cd;

void samplevar::build_samplevar(const cd::matrixptr &dptr, const cd::vectorptr &yptr)
{
    grid_.build_grid(dptr, h_);

    kernel_.build_kernel(dptr);

    
    const matrix &d = *(dptr);
    const vector &y = *(yptr);
    const matrixIptr g = grid_.get_grid();
    const matrixptr K = kernel_.get_kernel();

    if (d.rows() != y.size())
        throw std::length_error("void samplevar::build_samplevar(const cd::matrixptr &dptr, const cd::vectorptr &yptr): d and y must have the same number of rows");
    else
    {
        size_t n = g->rows();

        size_t hh = g->maxCoeff();

        variogram = std::make_shared<matrix>(matrix::Zero(hh+1, n));
        denominators = std::make_shared<matrix>(matrix::Zero(hh+1, n));


        matrix Y(n, n);

        #pragma omp for
        for (size_t i = 0; i < n; ++i)
        {
            for (size_t j = 1; j < n; ++j)
            {
                Y(i, j) = (y(i) - y(j)) * (y(i) - y(j));
                Y(j, i) = Y(i, j);
            }
        }

        #pragma omp parallel
        {
            int k = 0;
            /// for every location in d
            #pragma omp for
            for (size_t l=0; l < n; ++l)
            {
                Eigen::VectorXi counters = Eigen::VectorXi::Zero(hh+1);
                /// for every couple of locations in d
                for (size_t i = 0; i < n-1; ++i)
                {
                    for (size_t j = i+1; j < n; ++j)
                    {
                        /// if the vector between i and j belongs to the cell k
                        k = g->operator()(i, j);
                        if (k >= 0)
                        {
                            variogram->operator()(k, l) += K->operator()(l, i) * K->operator()(l, j) * Y(i, j);
                            denominators->operator()(k, l) += K->operator()(l, i) * K->operator()(l, j);
                            counters[k]++;
                        }   
                    }
                }
                for (size_t u = 0; u < hh+1; ++u)
                {
                    if (counters[u] != 0)
                        variogram->operator()(u, l) /= (2*denominators->operator()(u, l));
                }
            }
        }
        build_squaredweights();
    }
}


void samplevar::build_squaredweights()
{
    const matrixIptr g = grid_.get_grid(); 
    const vectorptr normh = grid_.get_normh();

    size_t n = g->rows();
    size_t htot = normh->rows();

    squaredweights = std::make_shared<matrix>(matrix::Zero(n, htot));

    #pragma omp parallel for
    for (size_t k = 0; k < n ; ++k)
    {
        for (size_t h = 0 ; h < htot; ++h)
        {
            if (normh->operator()(h) != 0)
                squaredweights->operator()(k,h) = denominators->operator()(h,k)/normh->operator()(h);
        }
    } 
}

samplevar::samplevar(const std::string kernel_id, const std::string grid_id, const unsigned int h, const scalar &epsilon): kernel_(kernel_id, epsilon), grid_(grid_id, epsilon), h_(h) {};

samplevar::samplevar(): kernel_(), grid_() {};

const matrixptr samplevar::get_variogram() const
{
    return variogram;
}

const matrixptr samplevar::get_denominators() const
{
    return denominators;
}

const matrixptr samplevar::get_squaredweights() const
{
    return squaredweights;
}

const vectorptr samplevar::get_x() const
{
    return grid_.get_x();
}

const vectorptr samplevar::get_y() const
{
    return grid_.get_y();
}

const matrixptr samplevar::get_kernel() const
{
    return kernel_.get_kernel();
}

const matrixIptr samplevar::get_grid() const
{
    return grid_.get_grid();
}

const vectorptr samplevar::get_normh() const
{
    return grid_.get_normh();
}