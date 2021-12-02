/// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
/// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>
/// Under MIT license

#ifndef LOCALLY_STATIONARY_MODES_VARIOGRAMFUNCTIONS
#define LOCALLY_STATIONARY_MODES_VARIOGRAMFUNCTIONS

#include "traits.hpp"

namespace LocallyStationaryModels
{
class VariogramFunction
{
protected:
    /**
     * \brief   convert the isotropic variogram in the equivalent anisotropic one calculating the norm of the spatial lag rotated and
     * expanded according to the eigenvalues and eigenvector of the anisotropy matrix
    */
    cd::scalar compute_anisotropic_h(const cd::scalar &lambda1, const cd::scalar &lambda2, const cd::scalar &phi, const cd::scalar &x, const cd::scalar &y);
public:
    VariogramFunction() = default;
    /**
     * \brief   return f(params, x, y)
    */
    virtual cd::scalar operator()(const cd::vector &params, const cd::scalar &x, const cd::scalar &y) = 0;
}; // class VariogramFunction

class Exponential: public VariogramFunction
{
public:
    Exponential() = default;
    /**
     * \brief           return sigma * sigma * (1 - exp(-h))
     * \param params    a vector with lambda1, lambda2, phi and sigma in this exact order
    */
    cd::scalar operator()(const cd::vector &params, const cd::scalar &x, const cd::scalar &y) override;
}; // class Exponential

class Matern: public VariogramFunction
{
public:
    Matern() = default;
    /**
     * \brief           return sigma * sigma *(1 - std::pow(std::sqrt(2*nu)*h, nu)*std::cyl_bessel_k(nu, std::sqrt(2*nu)*h)/(std::tgamma(nu)*std::pow(2,nu-1)))
     * \param params    a vector with lambda1, lambda2, phi, sigma and nu in this exact order
    */
    cd::scalar operator()(const cd::vector &params, const cd::scalar &x, const cd::scalar &y) override;
}; // class Matern

class MaternNuFixed: public VariogramFunction
{
private:
    double m_NU = 0.5;
public:
    MaternNuFixed(const double &NU): m_NU(NU) {};
    /**
     * \brief           return sigma * sigma *(1 - std::pow(std::sqrt(2*nu)*h, nu)*std::cyl_bessel_k(nu, std::sqrt(2*nu)*h)/(std::tgamma(nu)*std::pow(2,nu-1)))
     * \param params    a vector with lambda1, lambda2, phi and sigma in this exact order
    */
    cd::scalar operator()(const cd::vector &params, const cd::scalar &x, const cd::scalar &y) override;
}; // class MaternNuFixed

class Gaussian: public VariogramFunction
{
public:
    Gaussian() = default;
    /**
     * \brief           return sigma * sigma * (1 - exp(-h*h))
     * \param params    a vector with lambda1, lambda2, phi and sigma in this exact order
    */
    cd::scalar operator()(const cd::vector &params, const cd::scalar &x, const cd::scalar &y) override;
}; // class Gaussian

/**
 * \brief       allow to select between different functions for the variogram
 * \param id    the name of chosen variogram
*/
std::shared_ptr<VariogramFunction> make_variogramiso(const std::string &id);
}; // namespace LocallyStationaryModels

#endif //LOCALLY_STATIONARY_MODES_VARIOGRAM_FUNCTIONS