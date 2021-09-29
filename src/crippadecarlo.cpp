#include "crippadecarlo.hpp"

using namespace cd;
using namespace LBFGSpp;
using namespace std::chrono;

crippadecarlo::crippadecarlo(const cd::matrixptr &d_, const cd::vectorptr &y_, cd::matrixptr anchorpoints_, const double epsilon, const unsigned int n_angles, const unsigned int n_intervals): d(d_), y(y_), anchorpoints(anchorpoints_), epsilon_ottimale(epsilon) 
{
    samplevar samplevar_("gaussian", n_angles, n_intervals, epsilon_ottimale);
    samplevar_.build_samplevar(d, anchorpoints, y);

    opt opt_(samplevar_.get_variogram(), samplevar_.get_squaredweights(), samplevar_.get_x(),  samplevar_.get_y(), "esponenziale");
    opt_.findallsolutions();

    smt smt_(opt_.get_solutions(), d, exp(-20), exp(-18));

    delta_ottimale = smt_.get_optimal_delta();

    xatu_ = xatu("esponenziale", y, smt_, epsilon, d);
}




crippadecarlo::crippadecarlo(const cd::matrixptr &d_, const cd::vectorptr &y_, cd::matrixptr anchorpoints_, const double min_epsilon, const double max_epsilon, const unsigned int n_angles, const unsigned int n_intervals): d(d_), y(y_), anchorpoints(anchorpoints_)
{
    double min_error;

    for(int k=min_epsilon; k<max_epsilon; ++k)
    {
        double epsilon = -exp(k);

        std::cout << "provo con epsilon = " << epsilon << std::endl;

        double delta;
        double error = 0;
        for (unsigned int i=0; i<d->rows(); ++i)
        {
            matrixptr di = std::make_shared<matrix>(d->rows()-1, d->cols());
            vectorptr yi = std::make_shared<vector>(d->rows()-1);

            unsigned int counter=0;
            for (unsigned int j=0; j<d->rows(); ++j)
            {
                if (j!=i)
                {
                    di->row(counter) = d->row(j);
                    yi->operator()(counter) = y->operator()(j);
                    counter++;
                }
            }

            crippadecarlo CDi(di, yi, anchorpoints, epsilon, n_angles, n_intervals);
            double prediction = CDi.predict_y(d->row(i));
            double real = y->operator()(i);
            error += (prediction - real) * (prediction - real);
            delta = CDi.delta_ottimale;
        }
        if (k == min_epsilon || error < min_error)
        {
            epsilon_ottimale = epsilon;
            min_error = error;
            delta_ottimale = delta;
        }
    }

    samplevar samplevar_("gaussian", n_angles, n_intervals, epsilon_ottimale);
    samplevar_.build_samplevar(d, anchorpoints, y);

    opt opt_(samplevar_.get_variogram(), samplevar_.get_squaredweights(), samplevar_.get_x(),  samplevar_.get_y(), "esponenziale");
    opt_.findallsolutions();

    smt smt_(opt_.get_solutions(), d, delta_ottimale);

    xatu_ = xatu("esponenziale", y, smt_, epsilon_ottimale, d);
}





double crippadecarlo::predict_mean(const cd::vector &pos) const
{
    return xatu_.predict_mean(pos);
}

double crippadecarlo::predict_y(const cd::vector &pos) const
{
    return xatu_.predict_y(pos);
}

double crippadecarlo::get_epsilon() const
{
    return epsilon_ottimale;
}

double crippadecarlo::get_delta() const
{
    return delta_ottimale;
}