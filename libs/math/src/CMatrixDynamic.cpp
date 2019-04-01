/* +------------------------------------------------------------------------+
|                     Mobile Robot Programming Toolkit (MRPT)            |
|                          https://www.mrpt.org/                         |
|                                                                        |
| Copyright (c) 2005-2019, Individual contributors, see AUTHORS file     |
| See: https://www.mrpt.org/Authors - All rights reserved.               |
| Released under BSD License. See: https://www.mrpt.org/License          |
+------------------------------------------------------------------------+ */

#include "math-precomp.h"  // Precompiled headers

#include <mrpt/math/CMatrixDynamic.h>
#include <Eigen/Dense>

using namespace mrpt::math;

template <typename T>
CMatrixDynamic<T> CMatrixDynamic<T>::inverse() const
{
	ASSERT_(cols() == rows());
	const auto N = cols();
	CMatrixDynamic<T> inv(N, N);
	inv.asEigen() = this->asEigen().inverse().eval();
	return inv;
}

template <typename T>
CMatrixDynamic<T> CMatrixDynamic<T>::inverseLLt() const
{
	ASSERT_(cols() == rows());
	const auto N = cols();

	using MatX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
	const auto I = MatX::Identity(N, N);
	CMatrixDynamic<T> inv(N, N);
	inv.asEigen() = this->asEigen().llt().solve(I).eval();
	return inv;
}

// Template instantiation:
template class mrpt::math::CMatrixDynamic<float>;
template class mrpt::math::CMatrixDynamic<double>;
#ifdef HAVE_LONG_DOUBLE
template class mrpt::math::CMatrixDynamic<long double>;
#endif
