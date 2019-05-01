/* +------------------------------------------------------------------------+
   |                     Mobile Robot Programming Toolkit (MRPT)            |
   |                          https://www.mrpt.org/                         |
   |                                                                        |
   | Copyright (c) 2005-2019, Individual contributors, see AUTHORS file     |
   | See: https://www.mrpt.org/Authors - All rights reserved.               |
   | Released under BSD License. See: https://www.mrpt.org/License          |
   +------------------------------------------------------------------------+ */
#pragma once

#include <mrpt/math/CMatrixFixed.h>
#include <Eigen/Dense>

namespace mrpt::math
{

template <typename T, std::size_t ROWS, std::size_t COLS>
CMatrixFixed<float, ROWS, COLS> CMatrixFixed<T, ROWS, COLS>::cast_float() const
{
	CMatrixFixed<float, ROWS, COLS> r(rows(), cols());
	r.asEigen() = asEigen().template cast<float>();
	return r;
}

template <typename T, std::size_t ROWS, std::size_t COLS>
CMatrixFixed<double, ROWS, COLS> CMatrixFixed<T, ROWS, COLS>::cast_double()
const
{
	CMatrixFixed<double, ROWS, COLS> r(rows(), cols());
	r.asEigen() = asEigen().template cast<double>();
	return r;
}


}  // namespace mrpt::math
