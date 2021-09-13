#ifndef CRIPPADECARLO
#define CRIPPADECARLO

#include "samplevar.hpp"
#include "gradient.hpp"
#include "help.hpp"
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
    double epsilon_ottimale;
    double delta_ottimale;

    xatu xatu_;

public:
    crippadecarlo(const cd::matrixptr &d_, const cd::vectorptr &y_, const double epsilon, const unsigned int h_, const cd::vectorind &indici);
    crippadecarlo(const cd::matrixptr &d_, const cd::vectorptr &y_, const double min_epsilon, const double max_epsilon, const unsigned int h_, const cd::vectorind &indici);
    
    double predict_mean(const cd::vector &pos) const;
    double predict_y(const cd::vector &pos) const;

    double get_epsilon() const;
    double get_delta() const;
};








#endif //CRIPPADECARLO
