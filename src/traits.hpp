/// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
/// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>
/// Under MIT license

#ifndef LOCALLY_STATIONARY_MODELS_TRAITS
#define LOCALLY_STATIONARY_MODELS_TRAITS

#include <functional>
#include <memory>
#include <omp.h>
#include <vector>

#include "Eigen/Dense"

namespace cd
{
    /// defining basic types
    using scalar = double;
    using vector = Eigen::VectorXd;
    using matrix = Eigen::MatrixXd;
    using matrixI = Eigen::MatrixXi;
    using vectorptr = std::shared_ptr<vector>;
    using matrixptr = std::shared_ptr<matrix>;
    using matrixIptr = std::shared_ptr<matrixI>;
    using vectorind = std::vector<unsigned int>;
    using vofvs = std::vector<std::vector<unsigned int>>;
    using vofvsptr = std::shared_ptr<vofvs>;

    /// defining function types
    using kernelfunction = std::function<scalar(const vector&, const vector&, const scalar&)>;
    using gridfunction = std::function<matrixIptr(const matrixptr&, const unsigned int&, const unsigned int&, const double &)>;
    using variogramfunction = std::function<scalar(const vector&, const scalar&, const scalar&)>;
     
}; // namespace cd

#endif //LOCALLY_STATIONARY_MODELS_TRAITS