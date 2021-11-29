/// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
/// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>
/// Under MIT license

#include "kriging.hpp"

namespace LocallyStationaryModels
{
using namespace cd;

vectorind predictor::build_neighbourhood(const cd::vector &pos) const
{
    vectorind n;
    for (unsigned int i=0; i< d->rows(); ++i)
    {
        vector datapos = d->row(i);
        // if datapos is in a neighbourhood of radius b
        if ((pos - datapos).norm() < b)
            n.push_back(i);
    }
    return n;
}

vectorind predictor::build_neighbourhood(const unsigned int &pos) const
{
    vectorind n;
    const vector &pospos = d->row(pos);
    for (unsigned int i=0; i< d->rows(); ++i)
    {
        const vector &posi =  d->row(i);
        // if pos is in a neighbourhood of radius b
        if ((pospos - posi).norm() < b)
            n.push_back(i);
    }
    return n;
}

cd::vector predictor::build_eta(cd::vector &params, vectorind &neighbourhood) const
{
    unsigned int n = neighbourhood.size();
    matrix gamma(n,n);
    variogramfunction &gammaiso = *(gammaisoptr);
    // compute gamma
    for (unsigned int i=0; i<n; ++i)
    {
        for (unsigned int j=0; j<n; ++j)
        {
            const vector &posi = d->row(neighbourhood[i]);
            const vector &posj = d->row(neighbourhood[j]);
            cd::vector s = posi - posj;
            gamma(i, j) = gammaiso(params, s[0], s[1]);
        }
    }

    vector ones = vector::Ones(n);
    // if gamma is not invertible, return ones/n
    if (std::abs(gamma.determinant()) < 1e-12)
        return ones/n;

    // compute eta
    vector gammaones(n);
    gammaones = gamma.fullPivHouseholderQr().solve(ones);
    double denominator = ones.dot(gammaones);
    vector eta = (gammaones) / denominator;
    return eta;
}

std::pair<cd::vector, double> predictor::build_etakriging(const cd::vector &params,const cd::vector &pos) const
{
    unsigned int n = d->rows();
    
    vector etakriging(n);
    vector C0(n);
    matrix correlationmatrix(n,n);
    double sigma2 = params[3]*params[3];
    variogramfunction &gammaiso = *(gammaisoptr);
    // compute the corralation matrix and C0
    for (unsigned int i=0; i<n; ++i)
    {
        const vector &posi = d->row(i);
        for (unsigned int j=0; j<n; ++j)
        {
            const vector &posj = d->row(j);
            cd::vector s = posi - posj;
            correlationmatrix(i, j) = sigma2-gammaiso(params, s[0], s[1]);
        }
        cd::vector s0 = posi - pos;
        C0(i) = sigma2-gammaiso(params, s0[0], s0[1]);
    }
    // compute etakriging
    etakriging =  correlationmatrix.colPivHouseholderQr().solve(C0);
    // compute the variance
    double krigingvariance = params(3)*params(3) - C0.transpose()*etakriging;
    return std::make_pair(etakriging, krigingvariance);
}

template<>
double predictor::predict_mean<cd::vector, double>(const cd::vector &pos) const
{
    // find the value of the parameters in pos
    cd::vector params = smt_.smooth_vector(pos);
    // find the anchorpoints in its neighbourhood
    vectorind neighbourhood = build_neighbourhood(pos);
    unsigned int n = neighbourhood.size();
    // build eta
    vector eta(build_eta(params, neighbourhood));

    double result = 0;
    // compute the mean of z in pos
    for (unsigned int i=0; i<n; ++i)
        result += eta(i) * z->operator()(neighbourhood[i]);
    
    return result;
}

template<>
double predictor::predict_mean<unsigned int, double>(const unsigned int &pos) const
{
    // find the value of the parameters relative to the anchorpoint in row pos
    cd::vector params = smt_.smooth_vector(d->row(pos));
    // find the anchropoints in its neighbourhood
    vectorind neighbourhood = build_neighbourhood(pos);
    unsigned int n = neighbourhood.size();
    // build eta
    vector eta(build_eta(params, neighbourhood));

    double result = 0;
    // compute the mean of z in position pos
    for (unsigned int i=0; i<n; ++i)
        result += eta(i) * z->operator()(neighbourhood[i]);
    
    return result;
}

template<>
cd::vector predictor::predict_mean<cd::matrix, cd::vector>(const cd::matrix &pos) const
{
    vector result(pos.rows());
    #pragma omp parallel for
    for (size_t i=0; i<pos.rows(); ++i)
        result(i) = predict_mean<cd::vector, double>(pos.row(i));
    return result;
}

template<>
std::pair<double,double> predictor::predict_z<cd::vector, std::pair<double,double>>(const cd::vector &pos) const
{
    unsigned int n=d->rows();
    // predict the mean of z in pos
    double m0 = predict_mean<cd::vector, double>(pos);
    double result = m0;
    // find the value of the parameters in pos
    cd::vector params = smt_.smooth_vector(pos);
    // build etakriging and calculate the variance
    std::pair<vector, double> fulletakriging(build_etakriging(params, pos));
    vector &etakriging = fulletakriging.first;
    // predict the value of z(pos)
    for (unsigned int i=0; i<n; ++i)
        result += etakriging(i)*(z->operator()(i)-means->operator()(i)); 
    // return z(pos) and the kriging variance
    return std::make_pair(result, fulletakriging.second);
}

template<>
cd::matrix predictor::predict_z<cd::matrix, cd::matrix>(const cd::matrix &pos) const
{
    matrix result(pos.rows(), 2);
    #pragma omp parallel for
    for (size_t i=0; i<pos.rows(); ++i)
    {
        std::pair<double, double> prediction = predict_z<cd::vector, std::pair<double,double>>(pos.row(i));
        result(i, 0) = prediction.first;
        result(i, 1) = prediction.second;
    }
    return result;
}

predictor::predictor(const std::string &id, const cd::vectorptr &z_, const smt &mysmt, const double b_, const cd::matrixptr &d_): gammaisoptr(make_variogramiso(id)), z(z_), smt_(mysmt), b(b_), d(d_) 
{
    means = std::make_shared<vector>(z_->size());
    // build a vector with the prediction of the mean of z in every anchorpoint to speed up the next computations
    #pragma omp parallel for
    for (unsigned int i=0; i<means->size(); ++i)
        means->operator()(i) = predict_mean<unsigned int, double>(i);
};

predictor::predictor(): gammaisoptr(make_variogramiso("esponenziale")) {}
}; // namespace LocallyStationaryModels
