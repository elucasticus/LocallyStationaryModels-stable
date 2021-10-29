/// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
/// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>
/// Under MIT license

#include "kriging.hpp"

using namespace cd;

vectorind predictor::build_neighbourhood(const cd::vector &pos) const
{
    vectorind n;
    for (unsigned int i=0; i< d->rows(); ++i)
    {
        vector datapos = d->row(i);
        /// if datapos is in a neighbourhood of radius b
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
        /// if pos is in a neighbourhood of radius b
        if ((pospos - posi).norm() < b)
            n.push_back(i);
    }
    return n;
}


cd::vector predictor::build_eta(cd::vector &params, vectorind &neighbourhood) const
{
    unsigned int n = neighbourhood.size();
    matrix gamma(n,n);

    #pragma omp parallel for
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

    if (std::abs(gamma.determinant()) < 1e-12)
        return ones/n;
/*
    vector gammaones = gamma.colPivHouseholderQr().solve(ones); // QUESTO COMANDO FA CRASHARE R

    double denominator = ones.dot(gammaones);

    vector eta = (gammaones) / denominator;
*/

    gamma = gamma.inverse();

    double denominator = ones.transpose() * gamma * ones;

    vector eta = (gamma * ones) / denominator;

    return eta;
}

cd::vector predictor::build_etakriging(const cd::vector &params,const cd::vector &pos) const
{
    unsigned int n = d->rows();
    
    vector etakriging(n);
    vector C0(n);
    matrix correlationmatrix(n,n);
    double sigma2 = params[3]*params[3];
    
    #pragma omp parallel for
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
    
    etakriging =  correlationmatrix.colPivHouseholderQr().solve(C0);
    //etakriging = correlationmatrix.inverse()*C0;
    
    //double krigingvariance = params(3) - C0.transpose()*etakriging;
    return etakriging;
}

template<>
double predictor::predict_mean<cd::vector, double>(const cd::vector &pos) const
{
    cd::vector params = smt_.smooth_vector(pos);

    vectorind neighbourhood = build_neighbourhood(pos);
    unsigned int n = neighbourhood.size();

    vector eta(build_eta(params, neighbourhood));

    double result = 0;

    #pragma omp parallel for reduction(+:result)
    for (unsigned int i=0; i<n; ++i)
        result += eta(i) * y->operator()(neighbourhood[i]);
    
    return result;
}

template<>
double predictor::predict_mean<unsigned int, double>(const unsigned int &pos) const
{
    cd::vector params = smt_.smooth_vector(d->row(pos));

    vectorind neighbourhood = build_neighbourhood(pos);
    unsigned int n = neighbourhood.size();

    vector eta(build_eta(params, neighbourhood));

    double result = 0;

    #pragma omp parallel for reduction(+:result)
    for (unsigned int i=0; i<n; ++i)
        result += eta(i) * y->operator()(neighbourhood[i]);
    
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
double predictor::predict_y<cd::vector, double>(const cd::vector &pos) const
{
    unsigned int n=d->rows();
    double m0 = predict_mean<cd::vector, double>(pos);
    double result = m0;

    cd::vector params = smt_.smooth_vector(pos);

    vector etakriging(build_etakriging(params, pos));

    #pragma omp parallel for reduction(+:result)
    for (unsigned int i=0; i<n; ++i)
        result += etakriging(i)*(y->operator()(i)-means->operator()(i)); //sistemare predict mean 
    
    return result;
}

template<>
cd::vector predictor::predict_y<cd::matrix, cd::vector>(const cd::matrix &pos) const
{
    vector result(pos.rows());
    #pragma omp parallel for
    for (size_t i=0; i<pos.rows(); ++i)
        result(i) = predict_y<cd::vector, double>(pos.row(i));
    return result;
}

predictor::predictor(const std::string &id, const cd::vectorptr &y_, const smt &mysmt, const double b_, const cd::matrixptr &d_): gammaiso(make_variogramiso(id)), y(y_), smt_(mysmt), b(b_), d(d_) 
{
    means = std::make_shared<vector>(y_->size());
    #pragma omp parallel for
    for (unsigned int i=0; i<means->size(); ++i)
        means->operator()(i) = predict_mean<unsigned int, double>(i);
};

predictor::predictor(): gammaiso(make_variogramiso("esponenziale")) {}