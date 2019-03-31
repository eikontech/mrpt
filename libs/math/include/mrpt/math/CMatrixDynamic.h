/* +------------------------------------------------------------------------+
   |                     Mobile Robot Programming Toolkit (MRPT)            |
   |                          https://www.mrpt.org/                         |
   |                                                                        |
   | Copyright (c) 2005-2019, Individual contributors, see AUTHORS file     |
   | See: https://www.mrpt.org/Authors - All rights reserved.               |
   | Released under BSD License. See: https://www.mrpt.org/License          |
   +------------------------------------------------------------------------+ */
#pragma once

#include <mrpt/core/aligned_std_basicstring.h>
#include <mrpt/core/exceptions.h>  // ASSERT_()
#include <mrpt/core/format.h>
#include <mrpt/math/math_frwds.h>  // forward declarations
#include <mrpt/math/matrix_size_t.h>
#include <mrpt/typemeta/TTypeName.h>
#include <algorithm>  // swap()
#include <array>
#include <cstring>  // memset()
#include <type_traits>

namespace mrpt::math
{
/**  This template class provides the basic functionality for a general 2D
 *any-size, resizable container of numerical or non-numerical elements.
 * NOTES:
 * - This class is not serializable since it is a template. For using
 *serialization, see mrpt::math::CMatrixD.
 * - First row or column index is "0".
 * - This class includes range checks with ASSERT_() if compiling with "_DEBUG"
 *or "MRPT_ALWAYS_CHECKS_DEBUG_MATRICES=1".
 * - Use asEigen() to get an `Eigen::Map<>` object and to access full Algebra
 *functionality.
 *
 * \sa CMatrixTemplateNumeric
 * \ingroup mrpt_math_grp
 */
template <class T>
class CMatrixDynamic
{
   public:
	/** @name Matrix type definitions
	 * @{ */
	/** The type of the matrix elements */
	using value_type = T;
	using Scalar = T;
	using Index = int;
	using reference = T&;
	using const_reference = const T&;
	using size_type = int;
	using difference_type = std::ptrdiff_t;
	constexpr static int RowsAtCompileTime = -1;
	constexpr static int ColsAtCompileTime = -1;
	constexpr static int SizeAtCompileTime = -1;
	constexpr static int is_mrpt_type = 1;
	/** @} */

   private:
	/** RowMajor matrix data */
	mrpt::aligned_std_basicstring<T> m_data;
	size_t m_Rows{0}, m_Cols{0};

	/** Internal use only: It reallocs the memory for the 2D matrix, maintaining
	 * the previous contents if posible.
	 */
	void realloc(size_t row, size_t col, bool newElementsToZero = false)
	{
		if (row == m_Rows && col == m_Cols) return;
		const auto old_rows = m_Rows, old_cols = m_Cols;
		m_Rows = row;
		m_Cols = col;

		// New buffer:
		decltype(m_data) newData;
		newData.resize(m_Rows * m_Cols);
		// Copy old content:
		for (size_t r = 0; r < old_rows; r++)
		{
			if constexpr (std::is_trivial_v<T>)
				::memcpy(
					&newData[r * m_Cols], &m_data[r * old_cols],
					sizeof(T) * old_cols);
			else
				for (size_t c = 0; c < old_cols; c++)
					newData[r * m_Cols + c] = m_data[r * old_cols + c];
		}
		// New rows to zero?
		if (m_Rows > old_rows)
		{
			if constexpr (std::is_trivial_v<T>)
				::memset(
					&newData[old_rows * m_Cols], 0,
					sizeof(T) * (m_Rows - old_rows));
			else
				for (size_t r = old_rows; r < m_Rows; r++)
					for (size_t c = 0; c < m_Cols; c++)
						newData[r * m_Cols + c] = T();
		}
		// New cols to zero?
		if (m_Cols > old_cols)
		{
			for (size_t r = 0; r < old_rows; r++)
				if constexpr (std::is_trivial_v<T>)
					::memset(
						&newData[r * m_Cols + old_cols], 0,
						sizeof(T) * (m_Cols - old_cols));
				else
					for (size_t c = old_cols; c < m_Cols; c++)
						newData[r * m_Cols + c] = T();
		}
		// Swap:
		m_data.swap(newData);
	}

   public:
	/*! Fill all the elements with a given value (Note: named "fillAll" since
	 * "fill" will be used by child classes) */
	void fill(const T& val) { std::fill(m_data.begin(), m_data.end(), val); }

	inline void setZero() { fill(0); }
	inline void setZero(size_t nrows, size_t ncols)
	{
		realloc(nrows, ncols);
		fill(0);
	}

	/** Swap with another matrix very efficiently (just swaps a pointer and two
	 * integer values). */
	inline void swap(CMatrixDynamic<T>& o)
	{
		std::swap(m_data, o.m_data);
		std::swap(m_Rows, o.m_Rows);
		std::swap(m_Cols, o.m_Cols);
	}

	/** Constructors */
	CMatrixDynamic(const CMatrixDynamic& m) { (*this) = m; }

	CMatrixDynamic(size_t row = 1, size_t col = 1) { realloc(row, col); }

	/** Copy (casting from if needed) from another matrix  */
	template <typename U>
	explicit CMatrixDynamic(const CMatrixDynamic<U>& m)
	{
		(*this) = m;
	}

	/** Convert from Eigen matrix */
	template <class Derived>
	explicit CMatrixDynamic(const Eigen::MatrixBase<Derived>& m)
	{
		*this = m;
	}

	/** Copy constructor & crop from another matrix
	 */
	CMatrixDynamic(
		const CMatrixDynamic& m, const size_t cropRowCount,
		const size_t cropColCount)
	{
		ASSERT_(m.m_Rows >= cropRowCount);
		ASSERT_(m.m_Cols >= cropColCount);
		realloc(cropRowCount, cropColCount);
		for (size_t i = 0; i < m_Rows; i++)
			for (size_t j = 0; j < m_Cols; j++) (*this)(i, j) = m(i, j);
	}

	/** Constructor from a given size and a C array. The array length must match
	 *cols x row.
	 * \code
	 *  const double numbers[] = {
	 *    1,2,3,
	 *    4,5,6 };
	 *  CMatrixDouble   M(3,2, numbers);
	 * \endcode
	 */
	template <typename V, size_t N>
	CMatrixDynamic(size_t row, size_t col, V (&theArray)[N])
	{
		static_assert(N != 0, "Empty array!");
		realloc(row, col);
		if (m_Rows * m_Cols != N)
			THROW_EXCEPTION(format(
				"Mismatch between matrix size %lu x %lu and array of "
				"length %lu",
				static_cast<long unsigned>(m_Rows),
				static_cast<long unsigned>(m_Cols),
				static_cast<long unsigned>(N)));
		size_t idx = 0;
		for (size_t i = 0; i < m_Rows; i++)
			for (size_t j = 0; j < m_Cols; j++)
				(*this)(i, j) = static_cast<T>(theArray[idx++]);
	}

	/** Constructor from a given size and a STL container (std::vector,
	 * std::list,...) with the initial values. The vector length must match cols
	 * x row.
	 */
	template <typename V>
	CMatrixDynamic(size_t row, size_t col, const V& theVector)
	{
		const size_t N = theVector.size();
		realloc(row, col);
		if (m_Rows * m_Cols != N)
			THROW_EXCEPTION(format(
				"Mismatch between matrix size %lu x %lu and array of "
				"length %lu",
				static_cast<long unsigned>(m_Rows),
				static_cast<long unsigned>(m_Cols),
				static_cast<long unsigned>(N)));
		typename V::const_iterator it = theVector.begin();
		for (size_t i = 0; i < m_Rows; i++)
			for (size_t j = 0; j < m_Cols; j++)
				(*this)(i, j) = static_cast<T>(*(it++));
	}

	virtual ~CMatrixDynamic();

	template <class MAT>
	void setFromMatrixLike(const MAT& m)
	{
		MRPT_START
		setSize(m.rows(), m.cols());
		for (Index r = 0; r < rows(); r++)
			for (Index c = 0; c < cols(); c++) (*this)(r, c) = m(r, c);
		MRPT_END
	}

	/** Assignment operator from another matrix (possibly of a different type)
	 */
	template <typename U>
	CMatrixDynamic& operator=(const CMatrixDynamic<U>& m)
	{
		MRPT_START
		setFromMatrixLike(m);
		return *this;
		MRPT_END
	}

	/** Assignment from an Eigen matrix */
	template <class Derived>
	CMatrixDynamic& operator=(const Eigen::MatrixBase<Derived>& m)
	{
		MRPT_START
		setFromMatrixLike(m);
		return *this;
		MRPT_END
	}
	/** Assignment from a fixed matrix */
	template <typename U, std::size_t ROWS, std::size_t COLS>
	CMatrixDynamic& operator=(const CMatrixFixed<U, ROWS, COLS>& m)
	{
		MRPT_START
		setFromMatrixLike(m);
		return *this;
		MRPT_END
	}

	/** Assignment operator for initializing from a C array (The matrix must be
	 *set to the correct size before invoking this asignament)
	 * \code
	 *	 CMatrixDouble   M(3,2);
	 *  const double numbers[] = {
	 *    1,2,3,
	 *    4,5,6 };
	 *  M = numbers;
	 * \endcode
	 *  Refer also to the constructor with initialization data
	 *CMatrixDynamic::CMatrixDynamic
	 */
	template <typename V, size_t N>
	CMatrixDynamic& operator=(V (&theArray)[N])
	{
		static_assert(N != 0, "Empty array!");
		if (m_Rows * m_Cols != N)
		{
			THROW_EXCEPTION(format(
				"Mismatch between matrix size %lu x %lu and array of "
				"length %lu",
				m_Rows, m_Cols, N));
		}
		size_t idx = 0;
		for (size_t i = 0; i < m_Rows; i++)
			for (size_t j = 0; j < m_Cols; j++)
				(*this)(i, j) = static_cast<T>(theArray[idx++]);
		return *this;
	}

	/** Number of rows in the matrix \sa rows() */
	inline size_type rows() const { return m_Rows; }

	/** Number of columns in the matrix \sa rows() */
	inline size_type cols() const { return m_Cols; }

	/** Get a 2-vector with [NROWS NCOLS] (as in MATLAB command size(x)) */
	inline matrix_size_t size() const
	{
		matrix_size_t dims;
		dims[0] = m_Rows;
		dims[1] = m_Cols;
		return dims;
	}

	/** Changes the size of matrix, maintaining the previous contents. */
	void setSize(size_t row, size_t col, bool zeroNewElements = false)
	{
		realloc(row, col, zeroNewElements);
	}

	/** Resize the matrix */
	inline void resize(const matrix_size_t& siz, bool zeroNewElements = false)
	{
		setSize(siz[0], siz[1], zeroNewElements);
	}

	// These ones are to make template code compatible with Eigen & mrpt:
	CMatrixDynamic& derived() { return *this; }
	const CMatrixDynamic& derived() const { return *this; }
	void conservativeResize(size_t row, size_t col) { setSize(row, col); }

	/** Subscript operator to get/set individual elements
	 */
	inline T& operator()(size_t row, size_t col)
	{
#if defined(_DEBUG) || (MRPT_ALWAYS_CHECKS_DEBUG_MATRICES)
		if (row >= m_Rows || col >= m_Cols)
			THROW_EXCEPTION(format(
				"Indexes (%lu,%lu) out of range. Matrix is %lux%lu",
				static_cast<unsigned long>(row),
				static_cast<unsigned long>(col),
				static_cast<unsigned long>(m_Rows),
				static_cast<unsigned long>(m_Cols)));
#endif
		return m_data[row * m_Cols + col];
	}

	/** Subscript operator to get individual elements
	 */
	inline const T& operator()(size_t row, size_t col) const
	{
#if defined(_DEBUG) || (MRPT_ALWAYS_CHECKS_DEBUG_MATRICES)
		if (row >= m_Rows || col >= m_Cols)
			THROW_EXCEPTION(format(
				"Indexes (%lu,%lu) out of range. Matrix is %lux%lu",
				static_cast<unsigned long>(row),
				static_cast<unsigned long>(col),
				static_cast<unsigned long>(m_Rows),
				static_cast<unsigned long>(m_Cols)));
#endif
		return m_data[row * m_Cols + col];
	}

	/** Subscript operator to get/set an individual element from a row or column
	 * matrix.
	 * \exception std::exception If the object is not a column or row matrix.
	 */
	inline T& operator[](size_t ith)
	{
#if defined(_DEBUG) || (MRPT_ALWAYS_CHECKS_DEBUG_MATRICES)
		ASSERT_(m_Rows == 1 || m_Cols == 1);
#endif
		if (m_Rows == 1)
		{
// A row matrix:
#if defined(_DEBUG) || (MRPT_ALWAYS_CHECKS_DEBUG_MATRICES)
			if (ith >= m_Cols)
				THROW_EXCEPTION_FMT(
					"Index %u out of range!", static_cast<unsigned>(ith));
#endif
		}
		else
		{
// A columns matrix:
#if defined(_DEBUG) || (MRPT_ALWAYS_CHECKS_DEBUG_MATRICES)
			if (ith >= m_Rows)
				THROW_EXCEPTION_FMT(
					"Index %u out of range!", static_cast<unsigned>(ith));
#endif
		}
		return m_data[ith];
	}

	/** Subscript operator to get/set an individual element from a row or column
	 * matrix.
	 * \exception std::exception If the object is not a column or row matrix.
	 */
	inline const T& operator[](size_t ith) const
	{
#if defined(_DEBUG) || (MRPT_ALWAYS_CHECKS_DEBUG_MATRICES)
		ASSERT_(m_Rows == 1 || m_Cols == 1);
#endif
		if (m_Rows == 1)
		{
// A row matrix:
#if defined(_DEBUG) || (MRPT_ALWAYS_CHECKS_DEBUG_MATRICES)
			if (ith >= m_Cols)
				THROW_EXCEPTION_FMT(
					"Index %u out of range!", static_cast<unsigned>(ith));
#endif
		}
		else
		{
// A columns matrix:
#if defined(_DEBUG) || (MRPT_ALWAYS_CHECKS_DEBUG_MATRICES)
			if (ith >= m_Rows)
				THROW_EXCEPTION_FMT(
					"Index %u out of range!", static_cast<unsigned>(ith));
#endif
		}
		return m_data[ith];
	}

	/** Fast but unsafe method to write a value in the matrix
	 */
	inline void set_unsafe(size_t row, size_t col, const T& v)
	{
#ifdef _DEBUG
		if (row >= m_Rows || col >= m_Cols)
			THROW_EXCEPTION(format(
				"Indexes (%lu,%lu) out of range. Matrix is %lux%lu",
				static_cast<unsigned long>(row),
				static_cast<unsigned long>(col),
				static_cast<unsigned long>(m_Rows),
				static_cast<unsigned long>(m_Cols)));
#endif
		m_data[row][col] = v;
	}

	/** Fast but unsafe method to read a value from the matrix
	 */
	inline const T& get_unsafe(size_t row, size_t col) const
	{
#ifdef _DEBUG
		if (row >= m_Rows || col >= m_Cols)
			THROW_EXCEPTION(format(
				"Indexes (%lu,%lu) out of range. Matrix is %lux%lu",
				static_cast<unsigned long>(row),
				static_cast<unsigned long>(col),
				static_cast<unsigned long>(m_Rows),
				static_cast<unsigned long>(m_Cols)));
#endif
		return m_data[row][col];
	}

	/** Fast but unsafe method to get a reference from the matrix
	 */
	inline T& get_unsafe(size_t row, size_t col)
	{
#ifdef _DEBUG
		if (row >= m_Rows || col >= m_Cols)
			THROW_EXCEPTION(format(
				"Indexes (%lu,%lu) out of range. Matrix is %lux%lu",
				static_cast<unsigned long>(row),
				static_cast<unsigned long>(col),
				static_cast<unsigned long>(m_Rows),
				static_cast<unsigned long>(m_Cols)));
#endif
		return m_data[row][col];
	}

	/** Fast but unsafe method to obtain a pointer to a given row of the matrix
	 * (Use only in time critical applications)
	 */
	inline T* get_unsafe_row(size_t row)
	{
#ifdef _DEBUG
		if (row >= m_Rows)
			THROW_EXCEPTION(format(
				"Row index %lu out of range. Matrix is %lux%lu",
				static_cast<unsigned long>(row),
				static_cast<unsigned long>(m_Rows),
				static_cast<unsigned long>(m_Cols)));
#endif
		return m_data[row];
	}

	/** Fast but unsafe method to obtain a pointer to a given row of the matrix
	 * (Use only in critical applications)
	 */
	inline const T* get_unsafe_row(size_t row) const { return m_data[row]; }

	/** Appends a new row to the MxN matrix from a 1xN vector.
	 *  The lenght of the vector must match the width of the matrix, unless
	 * it's empty: in that case the matrix is resized to 1xN.
	 *  \code
	 *    CMatrixDouble  M(0,0);
	 *    CVectorDouble  v(7),w(7);
	 *    // ...
	 *    M.appendRow(v);
	 *    M.appendRow(w);
	 *  \endcode
	 * \exception std::exception On incorrect vector length.
	 * \sa extractRow
	 * \sa appendCol
	 */
	template <typename VECTOR>
	void appendRow(const VECTOR& in)
	{
		if (m_Cols == 0 || m_Rows == 0)
			ASSERT_(!in.empty());
		else
			ASSERT_(in.size() == m_Cols);
		const auto row = m_Rows;
		realloc(row + 1, m_Cols = in.size());
		for (size_t i = 0; i < m_Cols; i++) m_data[row][i] = in[i];
	}

	template <typename VECTOR>
	void setRow(const Index row, const VECTOR& v)
	{
		ASSERT_EQUAL_(cols(), v.size());
		for (Index c = 0; c < cols(); c++) (*this)(row, c) = v[c];
	}

	template <typename VECTOR>
	void setCol(const Index col, const VECTOR& v)
	{
		ASSERT_EQUAL_(rows(), v.size());
		for (Index r = 0; r < rows(); r++) (*this)(r, col) = v[r];
	}

	/** Appends a new column to the matrix from a vector.
	 * The length of the vector must match the number of rows of the matrix,
	 * unless it is (0,0).
	 * \exception std::exception On size mismatch.
	 * \sa extractCol
	 * \sa appendRow
	 */
	template <typename VECTOR>
	void appendCol(const VECTOR& in)
	{
		size_t r = m_Rows, c = m_Cols;
		if (m_Cols == 0 || m_Rows == 0)
		{
			ASSERT_(!in.empty());
			r = in.size();
			c = 0;
		}
		else
			ASSERT_(in.size() == m_Rows);
		realloc(r, c + 1);
		for (size_t i = 0; i < m_Rows; i++) m_data[i][m_Cols - 1] = in[i];
	}

	/** Returns a vector containing the matrix's values.
	 */
	template <typename VECTOR>
	void getAsVector(VECTOR& out) const
	{
		out.clear();
		out.reserve(m_Rows * m_Cols);
		for (size_t i = 0; i < m_Rows; i++)
			out.insert(out.end(), &(m_data[i][0]), &(m_data[i][m_Cols]));
	}

	/** Get as an Eigen-compatible Eigen::Map object  */
	template <
		typename EIGEN_MATRIX = Eigen::Matrix<T, -1, -1, 1, -1, -1>,
		typename EIGEN_MAP = Eigen::Map<
			EIGEN_MATRIX, MRPT_MAX_ALIGN_BYTES, Eigen::InnerStride<1>>>
	EIGEN_MAP asEigen()
	{
		return EIGEN_MAP(&m_data[0], m_Rows, m_Cols);
	}
	/** \overload (const version) */
	template <
		typename EIGEN_MATRIX = Eigen::Matrix<T, -1, -1, 1, -1, -1>,
		typename EIGEN_MAP = Eigen::Map<
			const EIGEN_MATRIX, MRPT_MAX_ALIGN_BYTES, Eigen::InnerStride<1>>>
	EIGEN_MAP asEigen() const
	{
		return EIGEN_MAP(&m_data[0], m_Rows, m_Cols);
	}

};  // end of class CMatrixDynamic

/** Declares a matrix of booleans (non serializable).
 *  \sa CMatrixDouble, CMatrixFloat, CMatrixB */
using CMatrixBool = CMatrixDynamic<bool>;

/** Declares a matrix of float numbers (non serializable).
 *  For a serializable version, use math::CMatrixF
 *  \sa CMatrixDouble, CMatrixF, CMatrixD
 */
using CMatrixFloat = CMatrixDynamic<float>;

/** Declares a matrix of double numbers (non serializable).
 *  For a serializable version, use math::CMatrixD
 *  \sa CMatrixFloat, CMatrixF, CMatrixD
 */
using CMatrixDouble = CMatrixDynamic<double>;

/** Declares a matrix of unsigned ints (non serializable).
 *  \sa CMatrixDouble, CMatrixFloat
 */
using CMatrixUInt = CMatrixDynamic<unsigned int>;

#ifdef HAVE_LONG_DOUBLE
/** Declares a matrix of "long doubles" (non serializable), or of "doubles" if
 * the compiler does not support "long double".
 *  \sa CMatrixDouble, CMatrixFloat
 */
using CMatrixLongDouble = CMatrixDynamic<long double>;
#else
/** Declares a matrix of "long doubles" (non serializable), or of "doubles" if
 * the compiler does not support "long double".
 *  \sa CMatrixDouble, CMatrixFloat
 */
using CMatrixLongDouble = CMatrixDynamic<double>;
#endif
}  // namespace mrpt::math

namespace mrpt::typemeta
{
// Extensions to mrpt::typemeta::TTypeName for matrices:
template <typename T>
struct TTypeName<mrpt::math::CMatrixDynamic<T>>
{
	static auto get()
	{
		return literal("CMatrixDynamic<") + TTypeName<T>::get() + literal(">");
	}
};
}  // namespace mrpt::typemeta
