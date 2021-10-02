#ifndef CRIPPADECARLO
#define CRIPPADECARLO

#include "samplevar.hpp"
#include "variogramfit.hpp"
#include <iostream>
#include "LBFGS.h"
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

    xatu xatu_;

    cd::matrixptr solutions = nullptr;
    cd::matrixptr empvar = nullptr;
    cd::matrixptr kernelmatrix = nullptr;
    cd::matrixIptr gridptr = nullptr;

public:
    crippadecarlo(const cd::matrixptr &d_, const cd::vectorptr &y_, cd::matrixptr anchorpoints_, const double epsilon, const unsigned int n_angles, const unsigned int n_intervals);
    crippadecarlo(const cd::matrixptr &d_, const cd::vectorptr &y_, cd::matrixptr anchorpoints_, const double min_epsilon, const double max_epsilon, const unsigned int n_angles, const unsigned int n_intervals);
    crippadecarlo(const cd::matrixptr &d_, const cd::vectorptr &y_, cd::matrixptr anchorpoints_, const double epsilon, const double delta, const cd::matrixptr &solutions_);
    
    double predict_mean(const cd::vector &pos) const;
    double predict_y(const cd::vector &pos) const;

    double get_epsilon() const;
    double get_delta() const;
    const cd::matrixptr get_solutions() const;
    const cd::matrixptr get_empiricvariogram() const;
    const cd::matrixptr get_kernel() const;
    const cd::matrixIptr get_grid() const;
};








#endif //CRIPPADECARLO
