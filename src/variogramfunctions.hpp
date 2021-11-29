/// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
/// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>
/// Under MIT license

#ifndef LOCALLY_STATIONARY_MODES_VARIOGRAMFUNCTIONS
#define LOCALLY_STATIONARY_MODES_VARIOGRAMFUNCTIONS

#include "traits.hpp"

namespace LocallyStationaryModels
{
class variogramfunction
{
protected:
    /**
     * \brief   convert the isotropic variogram in the equivalent anisotropic one calculating the norm of the spatial lag rotated and
     * expanded according to the eigenvalues and eigenvector of the anisotropy matrix
    */
    cd::scalar compute_anisotropic_h(const cd::scalar &lambda1, const cd::scalar &lambda2, const cd::scalar &phi, const cd::scalar &x, const cd::scalar &y);
public:
    variogramfunction() = default;
    /**
     * \brief   return f(params, x, y)
    */
    virtual cd::scalar operator()(const cd::vector &params, const cd::scalar &x, const cd::scalar &y) = 0;
};

class exponential: public variogramfunction
{
public:
    exponential() = default;
    /**
     * \brief           return sigma * sigma * (1 - exp(-h))
     * \param params    a vector with lambda1, lambda2, phi and sigma in this exact order
    */
    cd::scalar operator()(const cd::vector &params, const cd::scalar &x, const cd::scalar &y) override;
};

class matern: public variogramfunction
{
public:
    matern() = default;
    /**
     * \brief           return sigma * sigma *(1 - std::pow(std::sqrt(2*nu)*h, nu)*std::cyl_bessel_k(nu, std::sqrt(2*nu)*h)/(std::tgamma(nu)*std::pow(2,nu-1)))
     * \param params    a vector with lambda1, lambda2, phi, sigma and nu in this exact order
    */
    cd::scalar operator()(const cd::vector &params, const cd::scalar &x, const cd::scalar &y) override;
};

class maternNuFixed: public variogramfunction
{
private:
    double NU = 0.5;
public:
    maternNuFixed(const double &NU_): NU(NU_) {};
    /**
     * \brief           return sigma * sigma *(1 - std::pow(std::sqrt(2*nu)*h, nu)*std::cyl_bessel_k(nu, std::sqrt(2*nu)*h)/(std::tgamma(nu)*std::pow(2,nu-1)))
     * \param params    a vector with lambda1, lambda2, phi and sigma in this exact order
    */
    cd::scalar operator()(const cd::vector &params, const cd::scalar &x, const cd::scalar &y) override;
};

class gaussian: public variogramfunction
{
public:
    gaussian() = default;
    /**
     * \brief           return sigma * sigma * (1 - exp(-h*h))
     * \param params    a vector with lambda1, lambda2, phi and sigma in this exact order
    */
    cd::scalar operator()(const cd::vector &params, const cd::scalar &x, const cd::scalar &y) override;
};

/**
 * \brief       allow to select between different functions for the variogram
 * \param id    the name of chosen variogram
*/
std::shared_ptr<variogramfunction> make_variogramiso(const std::string &id);
}; // namespace LocallyStationaryModels

#endif //LOCALLY_STATIONARY_MODES_VARIOGRAM_FUNCTIONS