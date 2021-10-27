/// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
/// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>
/// Under MIT license

#include "cvinterface.hpp"

using namespace cd;
using namespace LBFGSpp;
using namespace std::chrono;

cvinterface::cvinterface(const cd::matrixptr &d_, const cd::vectorptr &y_, cd::matrixptr anchorpoints_, const cd::vector &parameters, const double epsilon, const unsigned int n_angles, 
    const unsigned int n_intervals, const std::string &kernel_id, const std::string &variogram_id): d(d_), y(y_), anchorpoints(anchorpoints_), epsilon_ottimale(epsilon) 
{
    samplevar samplevar_(kernel_id, n_angles, n_intervals, epsilon_ottimale);
    samplevar_.build_samplevar(d, anchorpoints, y);

    kernelmatrix = samplevar_.get_kernel();
    gridptr = samplevar_.get_grid();

    empvar = samplevar_.get_variogram();

    opt opt_(samplevar_.get_variogram(), samplevar_.get_squaredweights(), samplevar_.get_x(),  samplevar_.get_y(), variogram_id, parameters);
    opt_.findallsolutions();

    solutions = opt_.get_solutions();

    smt smt_(opt_.get_solutions(), anchorpoints, epsilon/10, epsilon/2);

    delta_ottimale = smt_.get_optimal_delta();

    predictor_ = predictor(variogram_id, y, smt_, epsilon, d);
}

cvinterface::cvinterface(const cd::matrixptr &d_, const cd::vectorptr &y_, cd::matrixptr anchorpoints_, const double epsilon, const double delta, const cd::matrixptr &solutions_, 
    const std::string &variogram_id):
    d(d_), y(y_), anchorpoints(anchorpoints_), epsilon_ottimale(epsilon), delta_ottimale(delta)
{
    smt smt_(solutions_, anchorpoints, delta_ottimale);
    predictor_ = predictor(variogram_id, y, smt_, epsilon, d);
}



cvinterface::cvinterface(const cd::matrixptr &d_, const cd::vectorptr &y_, cd::matrixptr anchorpoints_, const cd::vector &parameters, const double min_epsilon, const double max_epsilon, const unsigned int& nepsilons, 
const unsigned int n_angles, const unsigned int n_intervals, const std::string &kernel_id, const std::string &variogram_id): d(d_), y(y_), anchorpoints(anchorpoints_)
{
    double min_error;

    for(double k=min_epsilon; k<=max_epsilon; k+=(max_epsilon-min_epsilon)/(nepsilons-1))
    {
        double epsilon = k;

        std::cout << "provo con epsilon = " << epsilon << std::endl;

        double error = 0;
        #pragma omp parallel for reduction(+:error)
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

            cvinterface CDi(di, yi, anchorpoints, parameters, epsilon, n_angles, n_intervals, kernel_id, variogram_id);
            double prediction = CDi.predict_y<cd::vector, double>(d->row(i));
            double real = y->operator()(i);
            error += (prediction - real) * (prediction - real);
        }
        if (k == min_epsilon || error < min_error)
        {
            epsilon_ottimale = epsilon;
            min_error = error;
        }
    }

    samplevar samplevar_(kernel_id, n_angles, n_intervals, epsilon_ottimale);
    samplevar_.build_samplevar(d, anchorpoints, y);
    
    kernelmatrix = samplevar_.get_kernel();
    gridptr = samplevar_.get_grid();
    
    empvar = samplevar_.get_variogram();
    
    opt opt_(samplevar_.get_variogram(), samplevar_.get_squaredweights(), samplevar_.get_x(),  samplevar_.get_y(), variogram_id, parameters);
    opt_.findallsolutions();
    
    solutions = opt_.get_solutions();
    
    smt smt_(opt_.get_solutions(), anchorpoints, epsilon_ottimale/10, 10*epsilon_ottimale);
    
    delta_ottimale = smt_.get_optimal_delta();
    
    predictor_ = predictor(variogram_id, y, smt_, epsilon_ottimale, d);
}



/*
template<class Input, class Output>
Output cvinterface::predict_mean(const Input &pos) const
{
    return predictor_.predict_mean<Input, Output>(pos);
}


template<class Input, class Output>
Output cvinterface::predict_y(const Input &pos) const
{
    return predictor_.predict_y<Input, Output>(pos);
}
*/

double cvinterface::get_epsilon() const
{
    return epsilon_ottimale;
}

double cvinterface::get_delta() const
{
    return delta_ottimale;
}

const cd::matrixptr cvinterface::get_solutions() const
{
    return solutions;
}

const cd::matrixptr cvinterface::get_empiricvariogram() const
{
    return empvar;
}

const cd::matrixptr cvinterface::get_kernel() const
{
    return kernelmatrix;
}

const cd::matrixIptr cvinterface::get_grid() const
{
    return gridptr;
}