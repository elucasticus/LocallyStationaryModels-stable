/// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
/// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>
/// Under MIT license

#ifndef LOCALLY_STATIONARY_MODELS_ANCORA
#define LOCALLY_STATIONARY_MODELS_ANCORA

#include "traits.hpp"

namespace LocallyStationaryModels
{
/**
 * \brief   a simple class to find the anchor points given the data
*/
class Ancora
{
private:
    cd::matrixptr m_data; /// matrix to generate the anchor points
    double m_n_cubotti; /// number of tiles per row and column of the grid
    double m_larghezza = 0; /// total width of the grid
    double m_altezza = 0; /// total height of the grid
    double m_larghezza_cubo = 0; /// width of each tile
    double m_altezza_cubo = 0; /// height of each tile
    double m_center_x = 0; /// x of the origin of the grid
    double m_center_y = 0; /// y of the origin of the gird

    /**
     * \brief   return the index of the position in the grid of each of the points of the dataset "m_data"
    */
    Eigen::VectorXi cubotti() 
    {
        unsigned int n = m_data->rows();

        m_center_x = (m_data->col(0)).minCoeff()*(1 - Tolerances::anchor_tollerance);
        m_center_y = (m_data->col(1)).minCoeff()*(1 - Tolerances::anchor_tollerance);

        m_larghezza = (m_data->col(0)).maxCoeff()*(1 + Tolerances::anchor_tollerance) - m_center_x;
        m_altezza = (m_data->col(1)).maxCoeff()*(1 + Tolerances::anchor_tollerance) - m_center_y;
        m_larghezza_cubo = m_larghezza/m_n_cubotti;
        m_altezza_cubo = m_altezza/m_n_cubotti;

        // fill a vector with the position of each point
        Eigen::VectorXi result(n);
        for (unsigned int i=0; i<n; ++i)
        {
            cd::vector coordinates = m_data->row(i);
            result(i) = ceil((coordinates(0)-m_center_x)/m_larghezza_cubo) + m_n_cubotti*floor((coordinates(1)-m_center_y)/m_altezza_cubo);
        }
        return result;
    }

public:
    /**
     * \brief          constructor
     * \param data     shared pointer to the matrix with the coordinates of the dataset points
     * \param h        the number of squares per row and column of the grid
    */
    Ancora(const cd::matrixptr &data, const double h): m_data(data), m_n_cubotti(h){};

    /**
     * \brief   this function returns the coordinates of the anchor points in a way such that every anchor point has at least one point of the domain in its neighbourhood
    */
    const cd::matrix find_anchorpoints()
    {
        unsigned int n = m_data->rows();
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
        cd::matrix anchorpos(positions.size(), m_data->cols());
        for (unsigned int i=0; i<anchorpos.rows(); ++i)
        {
            unsigned int I = positions[i];
            anchorpos(i,0) = m_center_x + (I - floor((I*0.999999)/m_n_cubotti)*m_n_cubotti)*m_larghezza_cubo - m_larghezza_cubo/2;
            anchorpos(i,1) = m_center_y + ceil((I*0.999999)/m_n_cubotti)*m_altezza_cubo - m_altezza_cubo/2;
        }
        return anchorpos;
    }

    /**
     * \brief   return the coordinates of the origin of the grid
    */
    std::pair<double, double> get_center() const{return std::make_pair(m_center_x, m_center_y);}
    /**
     * \brief   return the dimensions (height and width) of each cell of the grid
    */
    std::pair<double, double> get_dimensionecubotti() const{return std::make_pair(m_larghezza_cubo, m_altezza_cubo);}
}; // class Ancora
}; // namespace LocallyStationaryModels

#endif //LOCALLY_STATIONARY_MODELS_ANCORA