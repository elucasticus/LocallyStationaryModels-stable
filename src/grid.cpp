#include "grid.hpp"
#include <iostream>
#include <unordered_map>

using namespace cd;

///-------------------------------------------------
/// GRID FUNCTIONS
///-------------------------------------------------

matrixIptr cubottiA(const matrixptr &d, const unsigned int &h, const double &epsilon)
{
    if (h == 0)
        throw std::length_error("matrixptr squares(const matrixptr &d, const unsigned int &h): h cannot be 0");
    else if (d->cols() != 2)
        throw std::domain_error("matrixptr squares(const matrixptr &d, const unsigned int &h): d.cols() must be equal to 2");
    else if (d->rows() <= 0)
        throw std::length_error("matrixptr squares(const matrixptr &d, const unsigned int &h): d.rows() must be greater than 0");
    else
    {
        matrixIptr grid(std::make_shared<matrixI>(matrixI::Constant(d->rows(), d->rows(), -1)));

        scalar grid_length = (d->col(0).maxCoeff() - d->col(0).minCoeff()) * 1.01;
        scalar grid_height = (d->col(1).maxCoeff() - d->col(1).minCoeff()) * 1.01;

        scalar cell_length = grid_length / h;
        scalar cell_height = grid_height / h;

        #pragma omp parallel for
        for (unsigned int i = 0; i < d->rows() - 1; ++i)
        {
            for (unsigned int j = i+1; j < d->rows(); ++j)
            {
                if (d->operator()(j, 0) < d->operator()(i, 0))
                    grid->operator()(i, j) = ceil( -(d->operator()(j, 0) - d->operator()(i, 0)) / cell_length ) + h * (h + floor( -(d->operator()(j, 1) - d->operator()(i, 1)) / cell_height) );
                else
                    grid->operator()(i, j) = ceil( (d->operator()(j, 0) - d->operator()(i, 0)) / cell_length ) + h * (h + floor( (d->operator()(j, 1) - d->operator()(i, 1)) / cell_height) );
            }
        }
        return grid;
    }
}


matrixIptr cubottiB(const matrixptr &d, const unsigned int &h, const double &epsilon)
{
    if (h == 0)
        throw std::length_error("matrixptr squares(const matrixptr &d, const unsigned int &h): h cannot be 0");
    else if (d->cols() != 2)
        throw std::domain_error("matrixptr squares(const matrixptr &d, const unsigned int &h): d.cols() must be equal to 2");
    else if (d->rows() <= 0)
        throw std::length_error("matrixptr squares(const matrixptr &d, const unsigned int &h): d.rows() must be greater than 0");
    else
    {
        matrixIptr grid(std::make_shared<matrixI>(matrixI::Constant(d->rows(), d->rows(), -1)));

        scalar grid_length = (d->col(0).maxCoeff() - d->col(0).minCoeff()) * 1.01;
        scalar grid_height = (d->col(1).maxCoeff() - d->col(1).minCoeff()) * 1.01;

        scalar cell_length = grid_length / h;
        scalar cell_height = grid_height / h;

        std::unordered_map<unsigned int, unsigned int> map_;
        unsigned int counter = 0;

        #pragma omp parallel for ordered
        for (unsigned int i = 0; i < d->rows() - 1; ++i)
        {
            for (unsigned int j = i+1; j < d->rows(); ++j)
            {
                unsigned int oldposition = 0;
                if (d->operator()(j, 0) < d->operator()(i, 0))
                    oldposition = ceil( -(d->operator()(j, 0) - d->operator()(i, 0)) / cell_length ) + h * (h + floor( -(d->operator()(j, 1) - d->operator()(i, 1)) / cell_height) );
                else
                    oldposition = ceil( (d->operator()(j, 0) - d->operator()(i, 0)) / cell_length ) + h * (h + floor( (d->operator()(j, 1) - d->operator()(i, 1)) / cell_height) );
                #pragma omp ordered
                {
                    auto iterator = map_.find(oldposition);
                    if (iterator == map_.end())
                    {
                        map_.insert(std::make_pair(oldposition, counter));
                        grid->operator()(i, j) = counter;
                        counter++;
                    }
                    else
                        grid->operator()(i, j) = map_.at(oldposition);
                }
            }
        }
        return grid;
    }
}

matrixIptr Pizza(const matrixptr &d, const unsigned int &h, const double &epsilon)
{
    if (h == 0)
        throw std::length_error("matrixptr squares(const matrixptr &d, const unsigned int &h): h cannot be 0");
    else if (d->cols() != 2)
        throw std::domain_error("matrixptr squares(const matrixptr &d, const unsigned int &h): d.cols() must be equal to 2");
    else if (d->rows() <= 0)
        throw std::length_error("matrixptr squares(const matrixptr &d, const unsigned int &h): d.rows() must be greater than 0");
    else
    {
        matrixIptr grid(std::make_shared<matrixI>(matrixI::Constant(d->rows(), d->rows(), -1)));
        scalar b = 2*epsilon;
        scalar cell_length = b / h;
        scalar cell_angle = 3.14159265358979323846 / (h);

        //#pragma omp parallel for collapse(2)
        for (unsigned int i = 0; i < d->rows() - 1; ++i)
        {
            for (unsigned int j = i+1; j < d->rows(); ++j)
            {
                scalar deltax =  d->operator()(j, 0) - d->operator()(i, 0);
                scalar deltay =  d->operator()(j, 1) - d->operator()(i, 1);
                scalar radius =  std::sqrt( deltax*deltax + deltay*deltay );
                scalar angle  =  std::atan( deltay / deltax );

                if (radius > b)
                    grid->operator()(i, j) = -1;
                else
                    grid->operator()(i, j) = floor( (3.14159265358979323846/2 + angle) / cell_angle ) + h * ( floor( (radius) / cell_length) );
            }
        }
        return grid;
    }
}


template<>
gridfunction make_grid<2>(const std::string &id)
{
    if (id == "cubottiA")
        return cubottiA;
    else if (id == "cubottiB")
        return cubottiB;
    else if (id == "Pizza")
        return Pizza;
    else
        return Pizza;
}

///--------------------------------------------------------------------------

///-------------------------------------------------
/// GRID CLASS
///-------------------------------------------------

void grid<2>::build_grid(const matrixptr &d, const unsigned int &h)
{
    g = f(d, h, sqrt(-0.5/epsilon));
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