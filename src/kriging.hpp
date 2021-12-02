/// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
/// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>
/// Under MIT license

#ifndef LOCALLY_STATIONARY_MODELS_KRIGING
#define LOCALLY_STATIONARY_MODELS_KRIGING

#include "traits.hpp"
#include "smooth.hpp"
#include "variogramfit.hpp"

namespace LocallyStationaryModels
{
/**
 * \brief   class to perform kriging on the data
*/
class Predictor
{
private:
    std::shared_ptr<VariogramFunction> m_gammaisoptr; /// variogram function
    cd::vectorptr m_z = nullptr; /// z(d)
    Smt m_smt; /// smoother
    double m_b; /// cutoff-radius of locally stationary neighbourhood
    cd::vectorptr m_means = nullptr; /// vector with the mean predicted in each anchor point
    cd::matrixptr m_d = nullptr; /// dataset with the initial points

    /**
     * \brief       build a vector with the index of the points in the neighbourhood of radius b of the point in position pos
     * \param pos   a vector of coordinates or the index of the position of the center of the neighbourhood
    */
    cd::vectorind build_neighbourhood(const cd::vector &pos) const;
    cd::vectorind build_neighbourhood(const unsigned int &pos) const;

    /**
     * \brief                   build the vector eta necessary to perform kriging on the mean of Y in a point
     * \param params            the params obtained by smoothing in the center of the neighbourhood
     * \param neighbourhood     a "neighbourhood" vector build with the previous functions
    */
    cd::vector build_eta(cd::vector &params, cd::vectorind &neighbourhood) const;

    /**
     * \brief                   build the vector eta necessary to perform kriging on Y in a point
     * \param params            the params obtained by smoothing in the center of the neighbourhood
     * \param neighbourhood     a "neighbourhood" vector build with the previous functions
    */
    std::pair<cd::vector, double> build_etakriging(const cd::vector &params,const cd::vector &pos) const;
    
public:
    /**
     * \brief           constructor
     * \param id        name of the variogram function associated with the problem
     * \param z         the vector with the value of the function Y in the known points
     * \param mysmt     the one used to previously smooth the variogram
     * \param b         the radius of the neighbourhood of the point where to perform kriging
     * \param d         a shared pointer to the matrix with the coordinates of the original dataset
    */
    Predictor(const std::string &id, const cd::vectorptr &z, const Smt &mysmt, const double b, const cd::matrixptr &d);
    /**
     * \brief   gammaiso set by default to exponential
    */
    Predictor();

    /**
     * \brief   predict the mean
    */
    template<typename Input, typename Output>    
    Output predict_mean(const Input &pos) const;

    /**
     * \brief   predict Z
    */
    template<typename Input, typename Output>
    Output predict_z(const Input &pos) const;
}; // class Predictor
}; // namespace LocallyStationaryModels

#endif //LOCALLY_STATIONARY_MODELS_KRIGING
