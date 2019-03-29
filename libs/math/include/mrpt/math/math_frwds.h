/* +------------------------------------------------------------------------+
   |                     Mobile Robot Programming Toolkit (MRPT)            |
   |                          https://www.mrpt.org/                         |
   |                                                                        |
   | Copyright (c) 2005-2019, Individual contributors, see AUTHORS file     |
   | See: https://www.mrpt.org/Authors - All rights reserved.               |
   | Released under BSD License. See: https://www.mrpt.org/License          |
   +------------------------------------------------------------------------+ */
#pragma once

#include <mrpt/config.h>
#include <string>
#include <type_traits>

/*! \file math_frwds.h
 * Forward declarations of all mrpt::math classes related to vectors, arrays
 * and matrices.
 * Many of the function implementations are in ops_matrices.h, others in
 * ops_containers.h
 */

// Minimum Eigen forward declarations for use in MRPT templates:
namespace Eigen
{
template <typename PlainObjectType, int MapOptions, typename StrideType>
class Map;
template <int Value>
class InnerStride;
template <
	typename _Scalar, int _Rows, int _Cols, int _Options,
	int _MaxRows /*= _Rows*/, int _MaxCols /* = _Cols*/>
class Matrix;
// For reference: _Options =
// /*AutoAlign*/ 0 | ((_Rows == 1 && _Cols != 1) ? /*Eigen::RowMajor*/ 1
//  : (_Cols == 1 && _Rows != 1) ? /*Eigen::ColMajor*/ 0 : /*default: Col*/ 0)
// ==> That is: _Options=1 for RowMajor as it's the default in MRPT matrices.

template <typename Derived>
class MatrixBase;
}  // namespace Eigen

namespace mrpt::system
{
std::string MRPT_getVersion();
}
namespace mrpt::math
{
struct TPoseOrPoint;
// class CMatrix;  // mrpt-binary-serializable matrix
// class CMatrixD;  // mrpt-binary-serializable matrix

/** For usage in one of the constructors of CMatrixFixedNumeric or
   CMatrixTemplate (and derived classes), if it's not required
	 to fill it with zeros at the constructor to save time. */
enum TConstructorFlags_Matrices
{
	UNINITIALIZED_MATRIX = 0
};

// ---------------- Forward declarations: Classes ----------------
template <class T>
class CMatrixTemplate;
template <class T>
class CMatrixTemplateObjects;
template <class T>
class CQuaternion;

/** ContainerType<T>::element_t exposes the value of any STL or Eigen container.
 *  Default specialization works for STL and MRPT containers, there is another
 * one for Eigen in <mrpt/math/eigen_frwds.h> */
template <typename CONTAINER>
struct ContainerType
{
	using element_t = typename CONTAINER::value_type;
};

#define MRPT_MATRIX_CONSTRUCTORS_FROM_POSES(_CLASS_)                          \
	template <                                                                \
		class TPOSE, typename = std::enable_if_t<                             \
						 std::is_base_of_v<mrpt::math::TPoseOrPoint, TPOSE>>> \
	explicit inline _CLASS_(const TPOSE& p)                                   \
	{                                                                         \
		mrpt::math::containerFromPoseOrPoint(*this, p);                       \
	}                                                                         \
	template <class CPOSE, int = CPOSE::is_3D_val>                            \
	explicit inline _CLASS_(const CPOSE& p)                                   \
	{                                                                         \
		mrpt::math::containerFromPoseOrPoint(*this, p);                       \
	}

template <class CONTAINER1, class CONTAINER2>
void cumsum(const CONTAINER1& in_data, CONTAINER2& out_cumsum);

template <class CONTAINER>
inline typename CONTAINER::Scalar norm(const CONTAINER& v);
template <class CONTAINER>
inline typename CONTAINER::Scalar norm_inf(const CONTAINER& v);

template <class MAT_A, class SKEW_3VECTOR, class MAT_OUT>
void multiply_A_skew3(const MAT_A& A, const SKEW_3VECTOR& v, MAT_OUT& out);
template <class SKEW_3VECTOR, class MAT_A, class MAT_OUT>
void multiply_skew3_A(const SKEW_3VECTOR& v, const MAT_A& A, MAT_OUT& out);

/** Conversion of poses (TPose2D,TPoint2D,...,
 * mrpt::poses::CPoint2D,CPose3D,...) to MRPT containers (vector/matrix) */
template <class CONTAINER, class POINT_OR_POSE>
CONTAINER& containerFromPoseOrPoint(CONTAINER& C, const POINT_OR_POSE& p);

}  // namespace mrpt::math
