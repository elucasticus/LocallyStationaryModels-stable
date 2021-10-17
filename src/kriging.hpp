#ifndef KRIGING
#define KRIGING

#include "traits.hpp"
#include "smooth.hpp"
#include "variogramfit.hpp"

/**
 * \brief class to perform kriging on the data
*/
class predictor
{
private:
    cd::variogramfunction gammaiso;
    cd::vectorptr y = nullptr;
    smt smt_;
    double b;
    cd::vectorptr means = nullptr;
    cd::matrixptr d = nullptr;

    /**
     * \brief build a vector with the index of the points in the neighbourhood of radius b of the point in position pos
     * \param pos   a vector of coordinates or the index of the position of the center of the neighbourhood
    */
    cd::vectorind build_neighbourhood(const cd::vector &pos) const;
    cd::vectorind build_neighbourhood(const unsigned int &pos) const;

    /**
     * \brief the vector eta necessary to perform kriging on the data on the neighbourhood of a point
     * \param params     the params obtained by smoothing in the center of the neighbourhood
     * \param neighbourhood     a "neighbourhood" vector build with the previous functions
    */
    cd::vector build_eta(cd::vector &params, cd::vectorind &neighbourhood) const;

public:
    /**
     * \brief constructor
     * \param id    name of the variogram function associated with the problem
     * \param y_    the vector with the value of the function Y in the known points
     * \param smt_  the one used to previously smooth the variogram
     * \param b_    the radius of the neighbourhood of the point where to perform kriging
    */
    predictor(const std::string &id, const cd::vectorptr &y_, const smt &smt__, const double b_, const cd::matrixptr &d_);
    /**
     * \brief gammaiso set by default to exponential
    */
    predictor();

    /**
     * \brief predict the mean value of Y in position pos
     * \param pos   a vector of coordinates or the index of the position
    */
    double predict_mean(const cd::vector &pos) const;
    double predict_mean(const unsigned int &pos) const;

    cd::vector predict_means(const cd::matrix &pos) const;

    /**
     * \brief predict the mean value of Y in position pos
     * \param pos   a vector of coordinates or the index of the position
    */
    double predict_y(const cd::vector &pos) const;

    cd::vector predict_ys(const cd::matrix &pos) const;
};








#endif //KRIGIGN