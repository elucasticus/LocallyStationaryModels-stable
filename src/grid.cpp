#include "grid.hpp"
#include <iostream>
#include <unordered_map>

using namespace cd;

///-------------------------------------------------
/// GRID FUNCTIONS
///-------------------------------------------------

matrixIptr Pizza(const matrixptr &d, const unsigned int &n_angles, const unsigned int &n_intervals, const double &epsilon)
{
    if (n_angles == 0 || n_intervals == 0)
        throw std::length_error("matrixptr squares(const matrixptr &d, const unsigned int &h): h cannot be 0");
    else if (d->cols() != 2)
        throw std::domain_error("matrixptr squares(const matrixptr &d, const unsigned int &h): d.cols() must be equal to 2");
    else if (d->rows() <= 0)
        throw std::length_error("matrixptr squares(const matrixptr &d, const unsigned int &h): d.rows() must be greater than 0");
    else
    {
        double pi = 4*std::atan(1.);
      
        matrixIptr grid(std::make_shared<matrixI>(matrixI::Constant(d->rows(), d->rows(), -1)));
        double b = 2*epsilon;
        double cell_length = b / n_intervals;
        double cell_angle = pi / (n_angles);
        

        #pragma omp parallel for
        for (unsigned int i = 0; i < d->rows() - 1; ++i)
        {
            for (unsigned int j = i+1; j < d->rows(); ++j)
            {
                scalar deltax =  d->operator()(j, 0) - d->operator()(i, 0);
                scalar deltay =  d->operator()(j, 1) - d->operator()(i, 1);
                scalar radius =  std::sqrt( deltax*deltax + deltay*deltay );
                //scalar angle  =  std::atan( deltay / deltax );

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
}


template<>
gridfunction make_grid<2>(const std::string &id)
{
    return Pizza;
}

///--------------------------------------------------------------------------

///-------------------------------------------------
/// GRID CLASS
///-------------------------------------------------

void grid<2>::build_grid(const matrixptr &d, const unsigned int &n_angles, const unsigned int &n_intervals)
{
    g = f(d, n_angles, n_intervals, epsilon);
    build_normh(d);
}

grid<2>::grid(const std::string &id, const double epsilon_): f(make_grid<2>(id)), epsilon(epsilon_){};

grid<2>::grid(): grid("Pizza",-1.162009e-07) {};

const matrixIptr grid<2>::get_grid() const {return g;}

const vectorptr grid<2>::get_normh() const {return normh;}

const vectorptr grid<2>::get_x() const {return x;}

const vectorptr grid<2>::get_y() const {return y;}

void grid<2>::build_normh(const matrixptr &data)
{
    const matrix &d = *(data);
    size_t n = g->rows();
    size_t hh = g->maxCoeff()+1;

    x->resize(hh);
    y->resize(hh);
    normh->resize(hh);
  
    vector nn = vector::Zero(hh);
    int k = 0;

    for (size_t i = 0; i < n-1 ; ++i)
    {
        for (size_t j = i+1 ; j < n; ++j)
        {
            k = g->operator()(i, j);

            if (k >= 0)
            {
                normh->operator()(k) += sqrt((d(j, 0) - d(i, 0))*(d(j, 0) - d(i, 0)) + (d(j,1)-d(i,1))*(d(j,1)-d(i,1)) );
                if ((d(j, 0) - d(i, 0))*(d(j,1) - d(i,1)) < 0)
                    x->operator()(k) -= std::abs(d(j, 0) - d(i, 0));
                else
                    x->operator()(k) += std::abs(d(j, 0) - d(i, 0));  
                y->operator()(k) += std::abs(d(j,1) - d(i,1));
                nn[k]++;
            }
        }
    }
    #pragma omp parallel for
    for (size_t u = 0; u < hh; ++u) 
    {
        if (nn[u] != 0)
        {
            normh->operator[](u) /= nn[u];
            x->operator[](u) /= nn[u];
            y->operator[](u) /= nn[u];
        }
    }
}