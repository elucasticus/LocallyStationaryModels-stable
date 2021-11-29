/// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
/// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>
/// Under MIT license

#include "samplevar.hpp"
#include <iostream>

namespace LocallyStationaryModels
{
using namespace cd;

void samplevar::build_samplevar(const cd::matrixptr &dptr, const cd::matrixptr &anchorpointsptr, const cd::vectorptr &zptr)
{
    grid_.build_grid(dptr, n_angles, n_intervals);

    kernel_.build_kernel(dptr, anchorpointsptr);

    // d is the matrix with the coordinates of the initial points
    const matrix &d = *(dptr);
    // z is the vector with z(d)
    const vector &z = *(zptr);
    // a is the matrix with the coordinates of the anchor points
    const matrix &a = *(anchorpointsptr);
    const matrixIptr g = grid_.get_grid();
    const matrix &K = *(kernel_.get_kernel());

    size_t n = g->rows();

    size_t max_index = g->maxCoeff();

    size_t N = a.rows();

    variogram = std::make_shared<matrix>(matrix::Zero(max_index+1, N));
    denominators = std::make_shared<matrix>(matrix::Zero(max_index+1, N));

    // Z is the matrix that contains the squared norm of the difference between each possible pair z_i and z_j
    matrix Z(n, n);

    #pragma omp parallel for
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 1; j < n; ++j)
        {
            Z(i, j) = (z(i) - z(j)) * (z(i) - z(j));
            Z(j, i) = Z(i, j);
        }
    }

    #pragma omp parallel
    {
        int k = 0;
        // for every location in d
        #pragma omp for
        for (size_t l=0; l < N; ++l)
        {
            Eigen::VectorXi counters = Eigen::VectorXi::Zero(max_index+1);
            // for every couple of locations in d
            for (size_t i = 0; i < n-1; ++i)
            {
                for (size_t j = i+1; j < n; ++j)
                {
                    // if the vector between i and j belongs to the cell k
                    k = g->operator()(i, j);
                    if (k >= 0)
                    {
                        scalar prodotto = K(l, i) * K(l, j);
                        variogram->operator()(k, l) += prodotto * Z(i, j);
                        denominators->operator()(k, l) += prodotto;
                        counters[k]++;
                    }   
                }
            }
            for (size_t u = 0; u < max_index+1; ++u)
            {
                if (counters[u] != 0)
                    variogram->operator()(u, l) /= (2*denominators->operator()(u, l));
            }
        }
    }
    build_squaredweights();
}

void samplevar::build_squaredweights()
{
    const matrixIptr g = grid_.get_grid(); 
    const vectorptr normh = grid_.get_normh();

    size_t N = denominators->cols();
    size_t htot = normh->rows();

    squaredweights = std::make_shared<matrix>(matrix::Zero(N, htot));

    #pragma omp parallel for
    for (size_t k = 0; k < N ; ++k)
    {
        for (size_t h = 0 ; h < htot; ++h)
        {
            if (normh->operator()(h) != 0)
                squaredweights->operator()(k,h) = denominators->operator()(h,k)/normh->operator()(h);
        }
    } 
}

samplevar::samplevar(const std::string &kernel_id, const unsigned int &n_angles_, const unsigned int &n_intervals_, const scalar &epsilon): kernel_(kernel_id, epsilon), grid_("Pizza", epsilon), n_angles(n_angles_), n_intervals(n_intervals_) {};

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
}; // namespace LocallyStationaryModels
