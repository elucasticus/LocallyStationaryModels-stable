/// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
/// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>
/// Under MIT license

#include "grid.hpp"
#include <iostream>
#include <unordered_map>

namespace LocallyStationaryModels
{
using namespace cd;

matrixIptr Pizza(const matrixptr &d, const unsigned int &n_angles, const unsigned int &n_intervals, const double &epsilon)
{
    double pi = 4*std::atan(1.);
    // create a square matrix of dimension d->rows()^2 and fill it with -1
    matrixIptr grid(std::make_shared<matrixI>(matrixI::Constant(d->rows(), d->rows(), -1)));
    // b is the radius of locally stationary neighbourhood as function of bandwidth parameter epsilon
    double b = 2*epsilon;
    double cell_length = b / n_intervals;
    double cell_angle = pi / (n_angles);
    
    // for every couple of points i and j in d compute the position of the vector (j - i) in the grid and fill grid(i, j) accordingly
    // since grid is symmetric we only need to fill the upper triangular part of the matrix 
    #pragma omp parallel for
    for (unsigned int i = 0; i < d->rows() - 1; ++i)
    {
        for (unsigned int j = i+1; j < d->rows(); ++j)
        {
            scalar deltax =  d->operator()(j, 0) - d->operator()(i, 0);
            scalar deltay =  d->operator()(j, 1) - d->operator()(i, 1);
            scalar radius =  std::sqrt( deltax*deltax + deltay*deltay );

            if (radius >= b)
                grid->operator()(i, j) = -1;
            else if (deltax!=0)
                grid->operator()(i, j) = floor( radius / cell_length ) + n_intervals *  floor( (pi/2 + std::atan( deltay / deltax )) / cell_angle );
            else 
                grid->operator()(i, j) = floor( radius / cell_length );
        }
    }
    return grid;
}

gridfunction make_grid(const std::string &id)
{
    return Pizza;
}

void grid::build_grid(const matrixptr &d, const unsigned int &n_angles, const unsigned int &n_intervals)
{
    g = f(d, n_angles, n_intervals, epsilon);
    build_normh(d);
}

grid::grid(const std::string &id, const double epsilon_): f(make_grid(id)), epsilon(epsilon_){};

grid::grid(): grid("Pizza", 1.) {};

const matrixIptr grid::get_grid() const {return g;}

const vectorptr grid::get_normh() const {return normh;}

const vectorptr grid::get_x() const {return mean_x;}

const vectorptr grid::get_y() const {return mean_y;}

void grid::build_normh(const matrixptr &data)
{
    const matrix &d = *(data);
    // n is the number of rows of grid which is equal to the number of points in d
    size_t n = g->rows();
    // max_index is the maximum index assigned to any pair of points in the grid
    size_t max_index = g->maxCoeff()+1;

    //mean_x will be the vector with the x of each cell of the grid (mean of the x of all the pairs inside)
    mean_x = std::make_shared<vector>(vector::Zero(max_index));
    // mean_y will be the vector with the y of each cell of the grid (mean of the y of all the pairs inside)
    mean_y = std::make_shared<vector>(vector::Zero(max_index));
    // normh will be the vector with the norm of each cell of the grid (mean of the norm of all the pairs inside)
    normh = std::make_shared<vector>(vector::Zero(max_index));
    // nn[k] will count how many times we will update the k-th element of mean_x, mean_y and normh
    // which corresponds to the number of vectors which fell in the k-th cell of the grid
    Eigen::VectorXi nn = Eigen::VectorXi::Zero(max_index);
    // k is a variable used to index elements according to their position in the grid
    // it has to be an integer since we have initialized each component of the grid with -1
    int k = 0;

    // for every couple of index i and j update mean_x, mean_y and normh in poistion k if g->operator(i, j) == k
    for (size_t i = 0; i < n-1 ; ++i)
    {
        for (size_t j = i+1 ; j < n; ++j)
        {
            k = g->operator()(i, j);

            if (k >= 0)
            {
                normh->operator()(k) += sqrt((d(j, 0) - d(i, 0))*(d(j, 0) - d(i, 0)) + (d(j,1)-d(i,1))*(d(j,1)-d(i,1)) );
                // because of the way we have constructed the grid we need to add the absolute value of pairs in the first and third quadrant
                // and subtract in the second and fourth ones
                if ((d(j, 0) - d(i, 0))*(d(j,1) - d(i,1)) < 0)
                    mean_x->operator()(k) -= std::abs(d(j, 0) - d(i, 0));
                else
                    mean_x->operator()(k) += std::abs(d(j, 0) - d(i, 0));  
                mean_y->operator()(k) += std::abs(d(j,1) - d(i,1));
                nn[k]++;
            }
        }
    }
    // now divide element by element normh, mean_x and mean_y by nn to get the sample mean
    for (size_t u = 0; u < max_index; ++u) 
    {
        if (nn[u] != 0)
        {
            normh->operator[](u) /= nn[u];
            mean_x->operator[](u) /= nn[u];
            mean_y->operator[](u) /= nn[u];
        }
    }
}
}; // namespace LocallyStationaryModels
