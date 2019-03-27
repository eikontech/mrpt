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

// Frwd decls for MRPT matrices:
namespace Eigen
{
template <typename PlainObjectType, int MapOptions, typename StrideType>
class Map;
template <int Value>
class InnerStride;
}  // namespace Eigen

namespace mrpt::system
{
std::string MRPT_getVersion();
}
namespace mrpt::math
{
struct TPoseOrPoint;
class CMatrix;  // mrpt-binary-serializable matrix
class CMatrixD;  // mrpt-binary-serializable matrix

namespace detail
{
/** Internal resize which compiles to nothing on fixed-size matrices. */
template <typename MAT, int TypeSizeAtCompileTime>
struct TAuxResizer
{
	static inline void internal_resize(MAT&, size_t, size_t) {}
	static inline void internal_resize(MAT&, size_t) {}
};
template <typename MAT>
struct TAuxResizer<MAT, -1>
{
	static inline void internal_resize(MAT& obj, size_t row, size_t col)
	{
		obj.derived().conservativeResize(row, col);
	}
	static inline void internal_resize(MAT& obj, size_t nsize)
	{
		obj.derived().conservativeResize(nsize);
	}
};
}  // namespace detail

/*! Selection of the number format in CMatrixTemplate::saveToTextFile
 */
enum TMatrixTextFileFormat
{
	/** engineering format '%e' */
	MATRIX_FORMAT_ENG = 0,
	/** fixed floating point '%f' */
	MATRIX_FORMAT_FIXED = 1,
	/** intergers '%i' */
	MATRIX_FORMAT_INT = 2
};

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
