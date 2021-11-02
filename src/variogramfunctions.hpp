/// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
/// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>
/// Under MIT license

#ifndef LOCALLY_STATIONARY_MODES_VARIOGRAMFUNCTIONS
#define LOCALLY_STATIONARY_MODES_VARIOGRAMFUNCTIONS

#include "traits.hpp"

class variogramfunction
{
protected:
    /**
     * \brief   convert the isotropic variogram in the equivalent anisotropic one calculating the norm of the spatial lag rotated and
     *          expanded according to the eigenvalues and eigenvector of the anisotropy matrix
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
     * \brief   return sigma * sigma * (1 - exp(-h))
    */
    cd::scalar operator()(const cd::vector &params, const cd::scalar &x, const cd::scalar &y) override;
};

class matern: public variogramfunction
{
public:
    matern() = default;
    /**
     * \brief   return sigma * sigma *(1 - std::pow(std::sqrt(2*nu)*h, nu)*std::cyl_bessel_k(nu, std::sqrt(2*nu)*h)/(std::tgamma(nu)*std::pow(2,nu-1)))
    */
    cd::scalar operator()(const cd::vector &params, const cd::scalar &x, const cd::scalar &y) override;
};

class gaussian: public variogramfunction
{
public:
    gaussian() = default;
    /**
     * \brief   return sigma * sigma * (1 - exp(-h*h))
    */
    cd::scalar operator()(const cd::vector &params, const cd::scalar &x, const cd::scalar &y) override;
};

/**
 * \brief       allow to select between different functions for the variogram
 * \param id    the name of your favourite variogram
*/
std::unique_ptr<variogramfunction> make_variogramiso(const std::string &id);

#endif //LOCALLY_STATIONARY_MODES_VARIOGRAM_FUNCTIONS