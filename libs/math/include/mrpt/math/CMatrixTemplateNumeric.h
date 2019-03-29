/* +------------------------------------------------------------------------+
   |                     Mobile Robot Programming Toolkit (MRPT)            |
   |                          https://www.mrpt.org/                         |
   |                                                                        |
   | Copyright (c) 2005-2019, Individual contributors, see AUTHORS file     |
   | See: https://www.mrpt.org/Authors - All rights reserved.               |
   | Released under BSD License. See: https://www.mrpt.org/License          |
   +------------------------------------------------------------------------+ */
#pragma once

#include <mrpt/core/exceptions.h>
#include <mrpt/math/CMatrixTemplate.h>
#include <mrpt/math/point_poses2vectors.h>  // MRPT_MATRIX_CONSTRUCTORS_FROM_POSES()
#include <mrpt/typemeta/TTypeName.h>
#include <Eigen/Dense>

namespace mrpt
{
namespace math
{
/**  A matrix of dynamic size.
 *   This class is a container for 2D matrix data, stored in RowMajor order.
 *
 * \note For a complete introduction to Matrices and vectors in MRPT, see:
 * https://www.mrpt.org/Matrices_vectors_arrays_and_Linear_Algebra_MRPT_and_Eigen_classes
 * \ingroup mrpt_math_grp
 */
template <typename T>
using CMatrixTemplateNumeric = CMatrixTemplate<T>;

/** Declares a matrix of float numbers (non serializable).
 *  For a serializable version, use math::CMatrix
 *  \sa CMatrixDouble, CMatrix, CMatrixD
 */
using CMatrixFloat = CMatrixTemplateNumeric<float>;

/** Declares a matrix of double numbers (non serializable).
 *  For a serializable version, use math::CMatrixD
 *  \sa CMatrixFloat, CMatrix, CMatrixD
 */
using CMatrixDouble = CMatrixTemplateNumeric<double>;

/** Declares a matrix of unsigned ints (non serializable).
 *  \sa CMatrixDouble, CMatrixFloat
 */
using CMatrixUInt = CMatrixTemplateNumeric<unsigned int>;

#ifdef HAVE_LONG_DOUBLE
/** Declares a matrix of "long doubles" (non serializable), or of "doubles" if
 * the compiler does not support "long double".
 *  \sa CMatrixDouble, CMatrixFloat
 */
using CMatrixLongDouble = CMatrixTemplateNumeric<long double>;
#else
/** Declares a matrix of "long doubles" (non serializable), or of "doubles" if
 * the compiler does not support "long double".
 *  \sa CMatrixDouble, CMatrixFloat
 */
using CMatrixLongDouble = CMatrixTemplateNumeric<double>;
#endif

}  // namespace math

namespace typemeta
{
// Extensions to mrpt::typemeta::TTypeName for matrices:
template <typename T>
struct TTypeName<mrpt::math::CMatrixTemplateNumeric<T>>
{
	static auto get()
	{
		return literal("CMatrixTemplateNumeric<") + TTypeName<T>::get() +
			   literal(">");
	}
};
}  // namespace typemeta

}  // namespace mrpt
