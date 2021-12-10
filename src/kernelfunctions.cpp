#include "kernelfunctions.hpp"

namespace LocallyStationaryModels
{
namespace kf
{
    using namespace cd;

    double gaussian(const vector &x, const vector &y, const double &epsilon)
    {
        return std::exp(-(x - y).squaredNorm()/(2*epsilon*epsilon));
    }

    kernelfunction make_kernel(const std::string &id)
    {
        return gaussian;
    }
} // namespace kf
} // namespace LocallyStationaryModels
