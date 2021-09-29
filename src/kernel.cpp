#include "kernel.hpp"

using namespace cd;

///-------------------------------------------------
/// KERNEL FUNCTIONS
///-------------------------------------------------

scalar gaussian(const vector &x, const vector &y, const scalar &epsilon)
{
	return std::exp((x - y).squaredNorm()/(2*epsilon*epsilon));
}

kernelfunction make_kernel(const std::string &id)
{
 	return gaussian;
}

///--------------------------------------------------------------------------

///-------------------------------------------------
/// KERNEL CLASS
///-------------------------------------------------

kernel::kernel(const std::string &id, const scalar &epsilon_): epsilon(epsilon_), f(make_kernel(id)) {};

kernel::kernel(): kernel("Gaussian", 1.) {};

scalar kernel::operator()(const vector &x, const vector &y) const
{
	if (x.size() != y.size())
		throw std::out_of_range("scalar kernel::operator()(const vector &x, const vector &y) const: vectors dimensions do not coincide");
	else
		return f(x, y, epsilon);
}


void kernel::build_kernel(const matrixptr &d, const matrixptr &anchorpoints)
{
	size_t n = d->rows();

	size_t N = anchorpoints->rows();

	k->resize(N, n);

	#pragma omp parallel for
	for (size_t i = 0; i < N; ++i)
	{
		for (size_t j = 0; j < n; ++j)
		{
			k->operator()(i, j) = this->operator()(anchorpoints->row(i), d->row(j));
		}
	}

	/// create a vector with the sum on each row of k
	vector sums(N);

	#pragma omp parallel for
	for (size_t i=0; i < N; ++i)
	{
		sums(i) = (k->row(i)).sum();
	}

	/// fill k with the new values
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

	#pragma omp parallel for
	for (size_t i = 0; i < n; ++i)
	{
		for (size_t j = i+1; j < n; ++j)
		{
			k->operator()(i, j) = this->operator()(d->row(i), d->row(j));
			k->operator()(j, i) = k->operator()(i, j);
		}
	}
}


void kernel::build_simple_kernel(const matrixptr &d, const scalar &epsilon_)
{
	epsilon = epsilon_;
	size_t n = d->rows();
	k->resize(n, n);

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

const matrixptr kernel::get_kernel() const {return k;}