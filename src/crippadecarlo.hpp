/// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
/// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>
/// Under MIT license

#ifndef LOCALLY_STATIONARY_MODELS_CRIPPADECARLO
#define LOCALLY_STATIONARY_MODELS_CRIPPADECARLO

#include "samplevar.hpp"
#include "variogramfit.hpp"
#include <iostream>
#include <chrono>
#include "smooth.hpp"
#include "kriging.hpp"

class crippadecarlo
{
private:
    cd::matrixptr d;
    cd::vectorptr y;
    cd::matrixptr anchorpoints;
    double epsilon_ottimale;
    double delta_ottimale;

    predictor predictor_;

    cd::matrixptr solutions = nullptr;
    cd::matrixptr empvar = nullptr;
    cd::matrixptr kernelmatrix = nullptr;
    cd::matrixIptr gridptr = nullptr;

public:
    crippadecarlo(const cd::matrixptr &d_, const cd::vectorptr &y_, cd::matrixptr anchorpoints_, const cd::vector &parameters, const double epsilon, const unsigned int n_angles, const unsigned int n_intervals,
        const std::string &kernel_id, const std::string &variogram_id);
    crippadecarlo(const cd::matrixptr &d_, const cd::vectorptr &y_, cd::matrixptr anchorpoints_, const cd::vector &parameters, const double min_epsilon, const double max_epsilon, const unsigned int& nepsilons,  const unsigned int n_angles, 
        const unsigned int n_intervals, const std::string &kernel_id, const std::string &variogram_id);
    crippadecarlo(const cd::matrixptr &d_, const cd::vectorptr &y_, cd::matrixptr anchorpoints_, const double epsilon, const double delta, const cd::matrixptr &solutions_, const std::string &variogram_id);
    
    template<class Input, class Output>
    Output predict_mean(const Input &pos) const{return predictor_.predict_mean<Input, Output>(pos);};

    template<class Input, class Output>
    Output predict_y(const Input &pos) const{return predictor_.predict_y<Input, Output>(pos);};

    double get_epsilon() const;
    double get_delta() const;
    const cd::matrixptr get_solutions() const;
    const cd::matrixptr get_empiricvariogram() const;
    const cd::matrixptr get_kernel() const;
    const cd::matrixIptr get_grid() const;
};








#endif //LOCALLY_STATIONARY_MODELS_CRIPPADECARLO
