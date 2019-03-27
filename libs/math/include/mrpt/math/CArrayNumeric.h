/* +------------------------------------------------------------------------+
   |                     Mobile Robot Programming Toolkit (MRPT)            |
   |                          https://www.mrpt.org/                         |
   |                                                                        |
   | Copyright (c) 2005-2019, Individual contributors, see AUTHORS file     |
   | See: https://www.mrpt.org/Authors - All rights reserved.               |
   | Released under BSD License. See: https://www.mrpt.org/License          |
   +------------------------------------------------------------------------+ */
#pragma once

#include <mrpt/math/CMatrixFixedNumeric.h>
#include <mrpt/typemeta/TTypeName.h>
#include <mrpt/typemeta/num_to_string.h>

namespace mrpt
{
namespace math
{
/** CArrayNumeric is an array for numeric types supporting several mathematical
 * operations (actually, just a wrapper on Eigen::Matrix<T,N,1>)
 * \sa CArrayFloat, CArrayDouble, CArray
 */
template <typename T, std::size_t N>
using CArrayNumeric = CMatrixFixedNumeric<T, N, 1>;

/** Specialization of CArrayNumeric for float numbers \sa CArrayNumeric  */
template <std::size_t N>
using CArrayFloat = CArrayNumeric<float, N>;

/** Specialization of CArrayNumeric for double numbers \sa CArrayNumeric  */
template <std::size_t N>
using CArrayDouble = CArrayNumeric<double, N>;

}  // namespace math

namespace typemeta
{
// Extensions to mrpt::typemeta::TTypeName for matrices:
template <typename T, size_t N>
struct TTypeName<mrpt::math::CArrayNumeric<T, N>>
{
	constexpr static auto get()
	{
		return literal("CArrayNumeric<") + TTypeName<T>::get() + literal(",") +
			   literal(num_to_string<N>::value) + literal(">");
	}
};
template <size_t N>
struct TTypeName<mrpt::math::CArrayDouble<N>>
{
	constexpr static auto get()
	{
		return literal("CArrayDouble<") + literal(num_to_string<N>::value) +
			   literal(">");
	}
};
template <size_t N>
struct TTypeName<mrpt::math::CArrayFloat<N>>
{
	constexpr static auto get()
	{
		return literal("CArrayFloat<") + literal(num_to_string<N>::value) +
			   literal(">");
	}
};
}  // namespace typemeta
}  // namespace mrpt
