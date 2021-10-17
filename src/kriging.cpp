#include "kriging.hpp"

using namespace cd;

vectorind predictor::build_neighbourhood(const cd::vector &pos) const
{
    vectorind n;

    #pragma omp parallel for
    for (unsigned int i=0; i< d->rows(); ++i)
    {
        vector datapos = d->row(i);
        #pragma omp critical
        if ((pos - datapos).norm() < b)
            n.push_back(i);
    }

    return n;
}


vectorind predictor::build_neighbourhood(const unsigned int &pos) const
{
    vectorind n;

    #pragma omp parallel for
    for (unsigned int i=0; i< d->rows(); ++i)
    {
        #pragma omp critical
        if ((d->row(pos) - d->row(i)).norm() < b && pos != i)
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
            cd::vector s = d->row(i) - d->row(j);
            gamma(i, j) = gammaiso(params, s[0], s[1]);
        }
    }

    gamma = gamma.inverse();

    vector ones = vector::Ones(n);

    double denominator = ones.transpose() * gamma * ones;

    vector eta = (gamma * ones) / denominator;

    return eta;
}


double predictor::predict_mean(const cd::vector &pos) const
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


double predictor::predict_mean(const unsigned int &pos) const
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



double predictor::predict_y(const cd::vector &pos) const
{
    double m0 = predict_mean(pos);
    double result = m0;

    vectorind neighbourhood = build_neighbourhood(pos);
    unsigned int n = neighbourhood.size();

    cd::vector params = smt_.smooth_vector(pos);

    vector eta(build_eta(params, neighbourhood));

    #pragma omp parallel for reduction(+:result)
    for (unsigned int i=0; i<n; ++i)
        result += eta(i) * (y->operator()(neighbourhood[i]) - means->operator()(neighbourhood[i]));
    
    return result;
}


predictor::predictor(const std::string &id, const cd::vectorptr &y_, const smt &smt__, const double b_, const cd::matrixptr &d_): gammaiso(make_variogramiso(id)), y(y_), smt_(smt__), b(b_), d(d_) 
{
    means = std::make_shared<vector>(y_->size());
    #pragma omp parallel for
    for (unsigned int i=0; i<means->size(); ++i)
        means->operator()(i) = predict_mean(i);
};

predictor::predictor(): gammaiso(make_variogramiso("esponenziale")) {}

cd::vector predictor::predict_means(const cd::matrix &pos) const
{
    vector result(pos.rows());
    #pragma omp parallel for
    for (size_t i=0; i<pos.rows(); ++i)
        result(i) = predict_mean(pos.row(i));
    return result;
}

cd::vector predictor::predict_ys(const cd::matrix &pos) const
{
    vector result(pos.rows());
    #pragma omp parallel for
    for (size_t i=0; i<pos.rows(); ++i)
        result(i) = predict_y(pos.row(i));
    return result;
}