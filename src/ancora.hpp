/// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
/// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>
/// Under MIT license

#ifndef LOCALLY_STATIONARY_MODELS_ANCORA
#define LOCALLY_STATIONARY_MODELS_ANCORA

#include "traits.hpp"
#include <vector>
#include <algorithm>

namespace LocallyStationaryModels
{
/**
 * \brief   a simple class to find the anchor points given the data
*/
class ancora
{
private:
    cd::matrixptr data; /// matrix to generate the anchor points
    double n_cubotti; /// number of tiles per row and column of the grid
    double larghezza = 0; /// total width of the grid
    double altezza = 0; /// total height of the grid
    double larghezza_cubo = 0; /// width of each tile
    double altezza_cubo = 0; /// height of each tile
    double center_x = 0; /// x of the origin of the grid
    double center_y = 0; /// y of the origin of the gird

    /**
     * \brief   return the index of the position in the grid of each of the points of the dataset "data"
    */
    Eigen::VectorXi cubotti() 
    {
        unsigned int n = data->rows();

        center_x = (data->col(0)).minCoeff()*0.999999;
        center_y = (data->col(1)).minCoeff()*0.999999;

        larghezza = (data->col(0)).maxCoeff()*1.000001 - center_x;
        altezza = (data->col(1)).maxCoeff()*1.000001 - center_y;
        larghezza_cubo = larghezza/n_cubotti;
        altezza_cubo = altezza/n_cubotti;

        // fill a vector with the position of each point
        Eigen::VectorXi result(n);
        for (unsigned int i=0; i<n; ++i)
        {
            cd::vector coordinates = data->row(i);
            result(i) = ceil((coordinates(0)-center_x)/larghezza_cubo) + n_cubotti*floor((coordinates(1)-center_y)/altezza_cubo);
        }
        return result;
    }

public:
    /**
     * \brief           constructor
     * \param data_     shared pointer to the matrix with the coordinates of the dataset points
     * \param h_        the number of squares per row and column of the grid
    */
    ancora(const cd::matrixptr &data_, const double h_): data(data_), n_cubotti(h_){};

    /**
     * \brief   this function returns the coordinates of the anchor points in a way such that every anchor point has at least one point of the domain in its neighbourhood
    */
    const cd::matrix find_anchorpoints()
    {
        unsigned int n = data->rows();
        Eigen::VectorXi results = cubotti();

        // build a new vector without duplicates
        std::vector<unsigned int> positions;
        for (unsigned int i=0; i<n; ++i)
        {
            unsigned int pos = results(i);
            if (std::find(positions.begin(), positions.end(), pos) == positions.end())
                positions.push_back(pos);
        }

        // fill a new matrix with the coordinates of each anchorpoins
        cd::matrix anchorpos(positions.size(), data->cols());
        for (unsigned int i=0; i<anchorpos.rows(); ++i)
        {
            unsigned int I = positions[i];
            anchorpos(i,0) = center_x + (I - floor((I*0.999999)/n_cubotti)*n_cubotti)*larghezza_cubo - larghezza_cubo/2;
            anchorpos(i,1) = center_y + ceil((I*0.999999)/n_cubotti)*altezza_cubo - altezza_cubo/2;
        }
        return anchorpos;
    }

    /**
     * \brief   return the coordinates of the origin of the grid
    */
    std::pair<double, double> get_center() const{return std::make_pair(center_x, center_y);}
    /**
     * \brief   return the dimensions (height and width) of each cell of the grid
    */
    std::pair<double, double> get_dimensionecubotti() const{return std::make_pair(larghezza_cubo, altezza_cubo);}
};
}; // namespace LocallyStationaryModels

#endif //LOCALLY_STATIONARY_MODELS_ANCORA