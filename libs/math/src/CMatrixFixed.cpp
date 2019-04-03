/* +------------------------------------------------------------------------+
|                     Mobile Robot Programming Toolkit (MRPT)            |
|                          https://www.mrpt.org/                         |
|                                                                        |
| Copyright (c) 2005-2019, Individual contributors, see AUTHORS file     |
| See: https://www.mrpt.org/Authors - All rights reserved.               |
| Released under BSD License. See: https://www.mrpt.org/License          |
+------------------------------------------------------------------------+ */

#include "math-precomp.h"  // Precompiled headers

#include <mrpt/math/CMatrixFixed.h>
#include <Eigen/Dense>
#include "MatrixVectorBase_impl.h"

using namespace mrpt::math;

template <typename T, std::size_t ROWS, std::size_t COLS>
template <typename T2>
CMatrixFixed<T2, ROWS, COLS> CMatrixFixed<T, ROWS, COLS>::cast() const
{
	CMatrixFixed<T, ROWS, COLS> r(rows(), cols());
	r.asEigen() = asEigen().template cast<T2>();
	return r;
}

// Template instantiations:
#define DO_MATFIXED_INSTANTIATION_NM(T_, N_, M_)         \
	template class mrpt::math::CMatrixFixed<T_, N_, M_>; \
	template class mrpt::math::MatrixVectorBase<         \
		T_, mrpt::math::CMatrixFixed<T_, N_, M_>>;

#define DO_MATFIXED_INSTANTIATION(T_)      \
	DO_MATFIXED_INSTANTIATION_NM(T_, 2, 2) \
	DO_MATFIXED_INSTANTIATION_NM(T_, 3, 3) \
	DO_MATFIXED_INSTANTIATION_NM(T_, 4, 4) \
	DO_MATFIXED_INSTANTIATION_NM(T_, 6, 6)

DO_MATFIXED_INSTANTIATION(float);
DO_MATFIXED_INSTANTIATION(double);
