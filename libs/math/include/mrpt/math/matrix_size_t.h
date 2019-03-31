/* +------------------------------------------------------------------------+
   |                     Mobile Robot Programming Toolkit (MRPT)            |
   |                          https://www.mrpt.org/                         |
   |                                                                        |
   | Copyright (c) 2005-2019, Individual contributors, see AUTHORS file     |
   | See: https://www.mrpt.org/Authors - All rights reserved.               |
   | Released under BSD License. See: https://www.mrpt.org/License          |
   +------------------------------------------------------------------------+ */
#pragma once

#include <array>

namespace mrpt::math
{
/** Auxiliary class used in CMatrixDynamic:size(), CMatrixDynamic::resize(),
 * CMatrixFixed::size(), CMatrixFixed::resize(), to mimic the
 * behavior of STL-containers.
 * \ingroup mrpt_math_grp
 */
struct matrix_size_t : public std::array<size_t, 2>
{
	/** Cast to size_t as the overall number of matrix/vector elements */
	operator std::size_t() const { return at(0) * at(1); }
};

}  // namespace mrpt::math
