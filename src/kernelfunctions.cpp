#include "kernelfunctions.hpp"

namespace LocallyStationaryModels {
namespace kf {
    using namespace cd;

    double gaussian(const vector& x, const vector& y, const double& epsilon)
    {
        return std::exp(-(x - y).squaredNorm() / (2 * epsilon * epsilon));
    }

    double identity(const vector& x, const vector& y, const double& epsilon)
    {
        if ((x - y).squaredNorm() > epsilon * epsilon) {
            return 0;
        } else {
            return 1;
        }
    }

    kernelfunction make_kernel(const std::string& id) { return gaussian; }
} // namespace kf
} // namespace LocallyStationaryModels
