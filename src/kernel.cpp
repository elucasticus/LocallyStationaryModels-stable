/// Copyright (C) Luca Crippa <luca7.crippa@mail.polimi.it>
/// Copyright (C) Giacomo De Carlo <giacomo.decarlo@mail.polimi.it>
/// Under MIT license

#include "kernel.hpp"

namespace LocallyStationaryModels
{
using namespace cd;

scalar gaussian(const vector &x, const vector &y, const scalar &epsilon)
{
	return std::exp(-(x - y).squaredNorm()/(2*epsilon*epsilon));
}

kernelfunction make_kernel(const std::string &id)
{
 	return gaussian;
}

kernel::kernel(const std::string &id, const scalar &epsilon_): epsilon(epsilon_), f(make_kernel(id)) {};

kernel::kernel(): kernel("Gaussian", 1.) {};

scalar kernel::operator()(const vector &x, const vector &y) const
{
	return f(x, y, epsilon);
}

void kernel::build_kernel(const matrixptr &d, const matrixptr &anchorpoints)
{
	size_t n = d->rows();

	size_t N = anchorpoints->rows();

	k->resize(N, n);

	// fill each component of k with the value of the kernel function evaluated between the i-th anchor point and the j-th initial point
	#pragma omp parallel for
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < n; ++j)
		{
			k->operator()(i, j) = this->operator()(anchorpoints->row(i), d->row(j));
		}
	}

	// create a vector with the sum on each row of k
	vector sums(N);
	#pragma omp parallel for
	for (size_t i=0; i < N; ++i)
	{
		sums(i) = (k->row(i)).sum();
	}

	// divide each element of k by the sum of the elements of its row to obtained the normalized version of the kernel matrix K*
	#pragma omp parallel for
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < n; ++j)
		{
			k->operator()(i, j) /= sums(i);
		}
	}
}

void kernel::build_simple_kernel(const matrixptr &d)
{
	size_t n = d->rows();
	k->resize(n, n);
	// fill each component of k with the kernel function evaluated between the i-th and the j-th point of d
	#pragma omp parallel for
	for (size_t i = 0; i < n; ++i)
	{
		for (size_t j = i; j < n; ++j)
		{
			k->operator()(i, j) = this->operator()(d->row(i), d->row(j));
			k->operator()(j, i) = k->operator()(i, j);
		}
	}
}

void kernel::build_simple_kernel(const matrixptr &d, const scalar &epsilon_)
{
	epsilon = epsilon_;
	build_simple_kernel(d);
}

const matrixptr kernel::get_kernel() const {return k;}
}; // namespace LocallyStationaryModels
