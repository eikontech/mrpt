/* +------------------------------------------------------------------------+
|                     Mobile Robot Programming Toolkit (MRPT)            |
|                          https://www.mrpt.org/                         |
|                                                                        |
| Copyright (c) 2005-2019, Individual contributors, see AUTHORS file     |
| See: https://www.mrpt.org/Authors - All rights reserved.               |
| Released under BSD License. See: https://www.mrpt.org/License          |
+------------------------------------------------------------------------+ */

#include "math-precomp.h"  // Precompiled headers

#include <mrpt/math/CVectorFixed.h>
#include <Eigen/Dense>
#include "MatrixVectorBase_impl.h"

// Explicit instantiation of "MatrixVectorBase_impl.h" methods:

#define EXPL_INST_VECTOR_FIX(T_)                       \
	template class mrpt::math::CMatrixFixed<T_, 2, 1>; \
	template class mrpt::math::CMatrixFixed<T_, 3, 1>; \
	template class mrpt::math::CMatrixFixed<T_, 4, 1>; \
	template class mrpt::math::CMatrixFixed<T_, 5, 1>; \
	template class mrpt::math::CMatrixFixed<T_, 6, 1>

EXPL_INST_VECTOR_FIX(float);
EXPL_INST_VECTOR_FIX(double);
