/// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
/// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>
/// Under MIT license

#include "kriging.hpp"

namespace LocallyStationaryModels {
using namespace cd;

vectorind Predictor::build_neighbourhood(const cd::vector& pos) const
{
    vectorind n;
    for (size_t i = 0; i < m_data->rows(); ++i) {
        vector datapos = m_data->row(i);
        // if datapos is in a neighbourhood of radius m_b
        if ((pos - datapos).norm() < m_b) {
            n.push_back(i);
        }
    }
    return n;
}

vectorind Predictor::build_neighbourhood(const size_t& pos) const
{
    vectorind n;
    const vector& pospos = m_data->row(pos);
    for (size_t i = 0; i < m_data->rows(); ++i) {
        const vector& posi = m_data->row(i);
        // if pos is in a neighbourhood of radius m_b
        if ((pospos - posi).norm() < m_b) {
            n.push_back(i);
        }
    }
    return n;
}

cd::vector Predictor::build_eta(cd::vector& params, vectorind& neighbourhood) const
{
    size_t n = neighbourhood.size();
    matrix gamma(n, n);
    VariogramFunction& gammaiso = *(m_gammaisoptr);
    // compute gamma
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            const vector& posi = m_data->row(neighbourhood[i]);
            const vector& posj = m_data->row(neighbourhood[j]);
            cd::vector s = posi - posj;
            gamma(i, j) = gammaiso(params, s[0], s[1]);
        }
    }

    vector ones = vector::Ones(n);
    // if gamma is not invertible, return ones/n
    if (std::abs(gamma.determinant()) < Tolerances::min_determinant) {
        return ones / n;
    }

    // compute eta
    vector gammaones(n);
    gammaones = gamma.fullPivHouseholderQr().solve(ones);
    double denominator = ones.dot(gammaones);
    vector eta = (gammaones) / denominator;
    return eta;
}

std::pair<cd::vector, double> Predictor::build_etakriging(const cd::vector& params, const cd::vector& pos) const
{
    size_t n = m_data->rows();

    vector etakriging(n);
    vector C0(n);
    matrix correlationmatrix(n, n);
    double sigma2 = params[3] * params[3];
    VariogramFunction& gammaiso = *(m_gammaisoptr);
    // compute the corralation matrix and C0
    for (size_t i = 0; i < n; ++i) {
        const vector& posi = m_data->row(i);
        for (size_t j = 0; j < n; ++j) {
            const vector& posj = m_data->row(j);
            cd::vector s = posi - posj;
            correlationmatrix(i, j) = sigma2 - gammaiso(params, s[0], s[1]);
        }
        cd::vector s0 = posi - pos;
        C0(i) = sigma2 - gammaiso(params, s0[0], s0[1]);
    }
    // compute etakriging
    etakriging = correlationmatrix.colPivHouseholderQr().solve(C0);
    // compute the variance
    double krigingvariance = params(3) * params(3) - C0.transpose() * etakriging;
    return std::make_pair(etakriging, krigingvariance);
}

template <> double Predictor::predict_mean<cd::vector, double>(const cd::vector& pos) const
{
    // find the value of the parameters in pos
    cd::vector params = m_smt.smooth_vector(pos);
    // find the anchorpoints in its neighbourhood
    vectorind neighbourhood = build_neighbourhood(pos);
    size_t n = neighbourhood.size();
    // build eta
    vector eta(build_eta(params, neighbourhood));

    double result = 0;
    // compute the mean of z in pos
    for (size_t i = 0; i < n; ++i) {
        result += eta(i) * m_z->operator()(neighbourhood[i]);
    }

    return result;
}

template <> double Predictor::predict_mean<size_t, double>(const size_t& pos) const
{
    // find the value of the parameters relative to the anchorpoint in row pos
    cd::vector params = m_smt.smooth_vector(m_data->row(pos));
    // find the anchropoints in its neighbourhood
    vectorind neighbourhood = build_neighbourhood(pos);
    size_t n = neighbourhood.size();
    // build eta
    vector eta(build_eta(params, neighbourhood));

    double result = 0;
    // compute the mean of z in position pos
    for (size_t i = 0; i < n; ++i) {
        result += eta(i) * m_z->operator()(neighbourhood[i]);
    }

    return result;
}

template <> cd::vector Predictor::predict_mean<cd::matrix, cd::vector>(const cd::matrix& pos) const
{
    vector result(pos.rows());
    #pragma omp parallel for
    for (size_t i = 0; i < pos.rows(); ++i) {
        result(i) = predict_mean<cd::vector, double>(pos.row(i));
    }
    return result;
}

template <>
std::pair<double, double> Predictor::predict_z<cd::vector, std::pair<double, double>>(const cd::vector& pos) const
{
    size_t n = m_data->rows();
    // predict the mean of z in pos
    double m0 = predict_mean<cd::vector, double>(pos);
    double result = m0;
    // find the value of the parameters in pos
    cd::vector params = m_smt.smooth_vector(pos);
    // build etakriging and calculate the variance
    std::pair<vector, double> fulletakriging(build_etakriging(params, pos));
    vector& etakriging = fulletakriging.first;
    // predict the value of z(pos)
    for (size_t i = 0; i < n; ++i) {
        result += etakriging(i) * (m_z->operator()(i) - m_means->operator()(i));
    }
    // return z(pos) and the kriging variance
    return std::make_pair(result, fulletakriging.second);
}

template <> cd::matrix Predictor::predict_z<cd::matrix, cd::matrix>(const cd::matrix& pos) const
{
    matrix result(pos.rows(), 2);
    #pragma omp parallel for
    for (size_t i = 0; i < pos.rows(); ++i) {
        std::pair<double, double> prediction = predict_z<cd::vector, std::pair<double, double>>(pos.row(i));
        result(i, 0) = prediction.first;
        result(i, 1) = prediction.second;
    }
    return result;
}

Predictor::Predictor(
    const std::string& id, const cd::vectorptr& z, const Smt& mysmt, const double& b, const cd::matrixptr& data)
    : m_gammaisoptr(make_variogramiso(id))
    , m_z(z)
    , m_smt(mysmt)
    , m_b(b)
    , m_data(data)
{
    m_means = std::make_shared<vector>(z->size());
    // build a vector with the prediction of the mean of z in every anchorpoint to speed up the next computations
    #pragma omp parallel for
    for (size_t i = 0; i < m_means->size(); ++i) {
        m_means->operator()(i) = predict_mean<size_t, double>(i);
    }
};

Predictor::Predictor()
    : m_gammaisoptr(make_variogramiso("esponenziale"))
{
}
} // namespace LocallyStationaryModels
