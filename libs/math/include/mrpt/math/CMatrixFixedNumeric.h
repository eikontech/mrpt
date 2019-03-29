/* +------------------------------------------------------------------------+
   |                     Mobile Robot Programming Toolkit (MRPT)            |
   |                          https://www.mrpt.org/                         |
   |                                                                        |
   | Copyright (c) 2005-2019, Individual contributors, see AUTHORS file     |
   | See: https://www.mrpt.org/Authors - All rights reserved.               |
   | Released under BSD License. See: https://www.mrpt.org/License          |
   +------------------------------------------------------------------------+ */
#pragma once

#include <mrpt/core/alignment_req.h>
#include <mrpt/core/exceptions.h>
#include <mrpt/math/math_frwds.h>  // Forward declarations
#include <mrpt/math/point_poses2vectors.h>  // MRPT_MATRIX_CONSTRUCTORS_FROM_POSES()
#include <mrpt/typemeta/TTypeName.h>
#include <mrpt/typemeta/num_to_string.h>
#include <array>
#include <cstddef>  // std::size_t

namespace mrpt::math
{
/** A compile-time fixed-size numeric matrix container.
 * It uses a RowMajor element memory layout.
 *
 * \sa CMatrixTemplateNumeric (for dynamic-size matrices)
 * \note For a complete introduction to Matrices and vectors in MRPT, see:
 * https://www.mrpt.org/Matrices_vectors_arrays_and_Linear_Algebra_MRPT_and_Eigen_classes
 * \ingroup mrpt_math_grp
 */
template <typename T, std::size_t ROWS, std::size_t COLS>
class CMatrixFixedNumeric
{
   public:
	/** @name Constructors
	 *  @{ */

	/** Default constructor, initializes all elements to zero */
	inline CMatrixFixedNumeric() { fill(0); }

	/** Constructor which leaves the matrix uninitialized.
	 *  Example of usage: CMatrixFixedNumeric<double,3,2>
	 * M(mrpt::math::UNINITIALIZED_MATRIX);
	 */
	inline CMatrixFixedNumeric(TConstructorFlags_Matrices) {}

	MRPT_MATRIX_CONSTRUCTORS_FROM_POSES(CMatrixTemplateNumeric)

	/** @} */

	/** @name Matrix element access & modifiers
	 *  @{ */

	/** Get as an Eigen-compatible Eigen::Map object  */
	template <
		typename EIGEN_MATRIX = Eigen::Matrix<T, ROWS, COLS>,
		typename EIGEN_MAP = Eigen::Map<
			EIGEN_MATRIX, MRPT_MAX_ALIGN_BYTES, Eigen::InnerStride<1>>>
	EIGEN_MAP asEigen()
	{
		return EIGEN_MAP(&m_data[0], ROWS, COLS);
	}
	/** \overload (const version) */
	template <
		typename EIGEN_MATRIX = Eigen::Matrix<
			T, ROWS, COLS, (ROWS != 1 && COLS == 1) ? 0 /*rowMajor*/ : 1, ROWS,
			COLS>,
		typename EIGEN_MAP = Eigen::Map<
			const EIGEN_MATRIX, MRPT_MAX_ALIGN_BYTES, Eigen::InnerStride<1>>>
	EIGEN_MAP asEigen() const
	{
		return EIGEN_MAP(&m_data[0], ROWS, COLS);
	}

	/** Access (row,col), without out-of-bounds check (except in Debug builds)
	 */
	inline T& operator()(int row, int col)
	{
		ASSERTDEB_(row < ROWS);
		ASSERTDEB_(col < COLS);
		return m_data[row * COLS + col];
	}
	inline const T& operator()(int row, int col) const
	{
		ASSERTDEB_(row < ROWS);
		ASSERTDEB_(col < COLS);
		return m_data[row * COLS + col];
	}

	/** Access the [i-th] element (for 1xN or Nx1 matrices) */
	inline T& operator[](int i)
	{
		ASSERT_(ROWS == 1 || COLS == 1);
		ASSERTDEB_(i < ROWS * COLS);
		return m_data[i];
	}
	inline const T& operator[](int i) const
	{
		ASSERT_(ROWS == 1 || COLS == 1);
		ASSERTDEB_(i < ROWS * COLS);
		return m_data[i];
	}

	void fill(const T& value) { m_data.fill(value); }
	void setZero() { m_data.fill(0); }
	void setIdentity()
	{
		for (std::size_t r = 0; r < ROWS; r++)
			for (std::size_t c = 0; c < COLS; c++)
				(*this)(r, c) = (r == c) ? 1 : 0;
	}

	/** @} */
	template <typename VECTOR>
	void loadFromArray(const VECTOR& vals)
	{
		constexpr auto LEN = std::size(vals);
		static_assert(LEN == ROWS * COLS, "Array of incorrect size.");
		for (size_t r = 0, i = 0; r < ROWS; r++)
			for (size_t c = 0; c < COLS; c++) m_data[r * COLS + c] = vals[i];
	}

   private:
	/** RowMajor matrix data */
	alignas(MRPT_MAX_ALIGN_BYTES) std::array<T, ROWS * COLS> m_data;
};
}  // namespace mrpt::math

namespace mrpt::typemeta
{
// Extensions to mrpt::typemeta::TTypeName for matrices:
template <typename T, std::size_t N, std::size_t M>
struct TTypeName<mrpt::math::CMatrixFixedNumeric<T, N, M>>
{
	constexpr static auto get()
	{
		return literal("CMatrixFixedNumeric<") + TTypeName<T>::get() +
			   literal(",") + literal(num_to_string<N>::value) + literal(",") +
			   literal(num_to_string<M>::value) + literal(">");
	}
};
}  // namespace mrpt::typemeta
