/// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
/// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>
/// Under MIT license

#include "variogramfunctions.hpp"

namespace LocallyStationaryModels
{
using namespace cd;

cd::scalar VariogramFunction::compute_anisotropic_h(const cd::scalar &lambda1, const cd::scalar &lambda2, const cd::scalar &phi, const cd::scalar &x, const cd::scalar &y)
{
    scalar xx = x * x;
    scalar yy = y * y;
    scalar xy = x * y;
    
    return sqrt((lambda2 * lambda2 * xx * cos(phi) * cos(phi) + lambda1 * lambda1 * yy * cos(phi) * cos(phi)
                     + lambda1 * lambda1 * xx * sin(phi) * sin(phi) + lambda2 * lambda2 * yy * sin(phi) * sin(phi)
                     + lambda1 * lambda1 * xy * sin(2 * phi) - lambda2 * lambda2 * xy * sin(2 * phi)) / (lambda1 * lambda1 * lambda2 * lambda2));
}

cd::scalar Exponential::operator()(const cd::vector &params, const cd::scalar &x, const cd::scalar &y)
{
    double lambda1 = params[0];
    double lambda2 = params[1];
    double phi = params[2];
    double sigma = params[3];

    scalar h = compute_anisotropic_h(lambda1, lambda2, phi, x, y);
    return sigma * sigma * (1 - exp(-h));
}

cd::scalar Matern::operator()(const cd::vector &params, const cd::scalar &x, const cd::scalar &y)
{
    double lambda1 = params[0];
    double lambda2 = params[1];
    double phi = params[2];
    double sigma = params[3];
    double nu = params[4];

    if (std::abs(x) < Tolerances::min_norm && std::abs(y) < Tolerances::min_norm)
    {
        return Tolerances::infinity;
    }
    
    scalar h = compute_anisotropic_h(lambda1, lambda2, phi, x, y);
    return sigma * sigma *(1 - std::pow(std::sqrt(2*nu)*h, nu)*std::cyl_bessel_k(nu, std::sqrt(2*nu)*h)/(std::tgamma(nu)*std::pow(2,nu-1)));
}

cd::scalar MaternNuFixed::operator()(const cd::vector &params, const cd::scalar &x, const cd::scalar &y)
{
    double lambda1 = params[0];
    double lambda2 = params[1];
    double phi = params[2];
    double sigma = params[3];
    double nu = m_NU;

    if (std::abs(x) < Tolerances::min_norm && std::abs(y) < Tolerances::min_norm)
    {
        return Tolerances::infinity;
    }
    
    scalar h = compute_anisotropic_h(lambda1, lambda2, phi, x, y);
    return sigma * sigma *(1 - std::pow(std::sqrt(2*nu)*h, nu)*std::cyl_bessel_k(nu, std::sqrt(2*nu)*h)/(std::tgamma(nu)*std::pow(2,nu-1)));
}

cd::scalar Gaussian::operator()(const cd::vector &params, const cd::scalar &x, const cd::scalar &y)
{
    double lambda1 = params[0];
    double lambda2 = params[1];
    double phi = params[2];
    double sigma = params[3];
    
    scalar h = compute_anisotropic_h(lambda1, lambda2, phi, x, y);
    return sigma * sigma * (1 - exp(-h*h));
}

std::shared_ptr<VariogramFunction> make_variogramiso(const std::string &id)
{
    if(id == "exponential" || id == "esponenziale")
    {
        return std::make_shared<Exponential>();
    }
    if(id == "matern" || id == "Matern")
    {
        return std::make_shared<Matern>();
    }
    if(id == "gaussian" || id == "Gaussian")
    {
        return std::make_shared<Gaussian>();
    }
    // using the following method we can set directly from R passing a string a constant value for nu
    if(id.substr(0, 13) == "maternNuFixed")
    {
        try
        {
            double NU = std::stod(id.substr(14));
            return std::make_shared<MaternNuFixed>(NU);
        }
        catch (std::exception &e)
        {
            return std::make_shared<Exponential>();
        }
    }
    return std::make_shared<Exponential>();
}
} // namespace LocallyStationaryModels
