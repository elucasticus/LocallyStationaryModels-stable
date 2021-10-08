#include "crippadecarlo.hpp"

using namespace cd;
using namespace LBFGSpp;
using namespace std::chrono;

crippadecarlo::crippadecarlo(const cd::matrixptr &d_, const cd::vectorptr &y_, cd::matrixptr anchorpoints_, const cd::vector &parameters, const double epsilon, const unsigned int n_angles, 
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

    smt smt_(opt_.get_solutions(), anchorpoints, epsilon/10, 10*epsilon);

    delta_ottimale = smt_.get_optimal_delta();

    xatu_ = xatu(variogram_id, y, smt_, epsilon, d);
}

crippadecarlo::crippadecarlo(const cd::matrixptr &d_, const cd::vectorptr &y_, cd::matrixptr anchorpoints_, const double epsilon, const double delta, const cd::matrixptr &solutions_, 
    const std::string &variogram_id):
    d(d_), y(y_), anchorpoints(anchorpoints_), epsilon_ottimale(epsilon), delta_ottimale(delta)
{
    smt smt_(solutions_, anchorpoints, delta_ottimale);
    xatu_ = xatu(variogram_id, y, smt_, epsilon, d);
}



crippadecarlo::crippadecarlo(const cd::matrixptr &d_, const cd::vectorptr &y_, cd::matrixptr anchorpoints_, const cd::vector &parameters, const double min_epsilon, const double max_epsilon, const unsigned int& nepsilons, 
const unsigned int n_angles, const unsigned int n_intervals, const std::string &kernel_id, const std::string &variogram_id): d(d_), y(y_), anchorpoints(anchorpoints_)
{
    double min_error;

    for(int k=min_epsilon; k<=max_epsilon; k+=(max_epsilon-min_epsilon)/(nepsilons-1))
    {
        double epsilon = k;

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

            crippadecarlo CDi(di, yi, anchorpoints, parameters, epsilon, n_angles, n_intervals, kernel_id, variogram_id);
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
    
    xatu_ = xatu(variogram_id, y, smt_, epsilon_ottimale, d);
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

const cd::matrixptr crippadecarlo::get_solutions() const
{
    return solutions;
}

const cd::matrixptr crippadecarlo::get_empiricvariogram() const
{
    return empvar;
}

const cd::matrixptr crippadecarlo::get_kernel() const
{
    return kernelmatrix;
}

const cd::matrixIptr crippadecarlo::get_grid() const
{
    return gridptr;
}


cd::vector crippadecarlo::predict_means(const cd::matrix &pos) const
{
    return xatu_.predict_means(pos);
}

cd::vector crippadecarlo::predict_ys(const cd::matrix &pos) const
{
    return xatu_.predict_ys(pos);
}