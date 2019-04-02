/* +------------------------------------------------------------------------+
|                     Mobile Robot Programming Toolkit (MRPT)            |
|                          https://www.mrpt.org/                         |
|                                                                        |
| Copyright (c) 2005-2019, Individual contributors, see AUTHORS file     |
| See: https://www.mrpt.org/Authors - All rights reserved.               |
| Released under BSD License. See: https://www.mrpt.org/License          |
+------------------------------------------------------------------------+ */

#include "math-precomp.h"  // Precompiled headers

#include <mrpt/math/CVectorDynamic.h>
#include <Eigen/Dense>
#include "MatrixVectorBase_impl.h"

// Explicit instantiation of "MatrixVectorBase_impl.h" methods:

template class mrpt::math::CVectorDynamic<float>;
template class mrpt::math::CVectorDynamic<double>;
