/* +------------------------------------------------------------------------+
   |                     Mobile Robot Programming Toolkit (MRPT)            |
   |                          https://www.mrpt.org/                         |
   |                                                                        |
   | Copyright (c) 2005-2019, Individual contributors, see AUTHORS file     |
   | See: https://www.mrpt.org/Authors - All rights reserved.               |
   | Released under BSD License. See: https://www.mrpt.org/License          |
   +------------------------------------------------------------------------+ */
#pragma once

#include <mrpt/core/aligned_allocator.h>
#include <mrpt/core/exceptions.h>  // ASSERT_()
#include <mrpt/core/format.h>
#include <mrpt/math/math_frwds.h>  // forward declarations
#include <mrpt/math/matrix_size_t.h>
#include <algorithm>  // swap()
#include <array>
#include <cstring>  // memset()
#include <type_traits>
#include <vector>

namespace mrpt::math
{
/**  This template class provides the basic functionality for a general 2D
 *any-size, resizable container of numerical or non-numerical elements.
 * NOTES:
 * - This class is not serializable since it is a template. For using
 *serialization, see mrpt::math::CMatrixNumeric
 * - First row or column index is "0".
 * - This class includes range checks with ASSERT_() if compiling with "_DEBUG"
 *or "MRPT_ALWAYS_CHECKS_DEBUG_MATRICES=1".
 *
 * \note Memory blocks for each row are 16-bytes aligned (since MRPT 0.7.0).
 * \note For a complete introduction to Matrices and vectors in MRPT, see:
 *http://www.mrpt.org/Matrices_vectors_arrays_and_Linear_Algebra_MRPT_and_Eigen_classes
 * \sa CMatrixTemplateNumeric
 * \ingroup mrpt_math_grp
 */
template <class T>
class CMatrixTemplate
{
   public:
	// type definitions
	/** The type of the matrix elements */
	using value_type = T;
	using reference = T&;
	using const_reference = const T&;
	using size_type = std::size_t;
	using difference_type = std::ptrdiff_t;

   protected:
	T** m_Val;
	size_t m_Rows{0}, m_Cols{0};

	/** Internal use only: It reallocs the memory for the 2D matrix, maintaining
	 * the previous contents if posible.
	 */
	void realloc(size_t row, size_t col, bool newElementsToZero = false)
	{
		MRPT_TODO("Refactor storage to simple T*, not T** -> Map<> compatib.");

		if (row != m_Rows || col != m_Cols || m_Val == nullptr)
		{
			size_t r;
			bool doZeroColumns = newElementsToZero && (col > m_Cols);

			// If we are reducing rows, free that memory:
			for (r = row; r < m_Rows; r++) mrpt::aligned_free(m_Val[r]);

			// Realloc the vector of pointers:
			if (!row)
			{
				mrpt::aligned_free(m_Val);
				m_Val = nullptr;
			}
			else
				m_Val = static_cast<T**>(
					mrpt::aligned_realloc(m_Val, sizeof(T*) * row, 16));

			// How many new rows/cols?
			size_t row_size = col * sizeof(T);

			// Alloc new ROW pointers & resize previously existing rows, as
			// required:
			for (r = 0; r < row; r++)
			{
				if (r < m_Rows)
				{
					// This was an existing row: Resize the memory:
					m_Val[r] = static_cast<T*>(
						mrpt::aligned_realloc(m_Val[r], row_size, 16));

					if (doZeroColumns)
					{
						// Fill with zeros:
						if constexpr (std::is_trivial_v<T>)
							::memset(
								&m_Val[r][m_Cols], 0,
								sizeof(T) * (col - m_Cols));
						else
							for (size_t k = m_Cols; k < col; k++)
								m_Val[r][k] = T();
					}
				}
				else
				{
					// This is a new row, alloc the memory for the first time:
					m_Val[r] =
						static_cast<T*>(mrpt::aligned_malloc(row_size, 16));
					if constexpr (std::is_trivial_v<T>)
						::memset(m_Val[r], 0, row_size);
					else
						for (size_t k = 0; k < col; k++) m_Val[r][k] = T();
				}
			}
			// Done!
			m_Rows = row;
			m_Cols = col;
		}
	}

   public:
	/*! Fill all the elements with a given value (Note: named "fillAll" since
	 * "fill" will be used by child classes) */
	void fill(const T& val)
	{
		for (size_t r = 0; r < m_Rows; r++)
			for (size_t c = 0; c < m_Cols; c++) m_Val[r][c] = val;
	}

	/** Swap with another matrix very efficiently (just swaps a pointer and two
	 * integer values). */
	inline void swap(CMatrixTemplate<T>& o)
	{
		std::swap(m_Val, o.m_Val);
		std::swap(m_Rows, o.m_Rows);
		std::swap(m_Cols, o.m_Cols);
	}

	/** Constructors */
	CMatrixTemplate(const CMatrixTemplate& m)
		: m_Val(nullptr), m_Rows(0), m_Cols(0)
	{
		(*this) = m;
	}

	CMatrixTemplate(size_t row = 1, size_t col = 1) : m_Val(nullptr)
	{
		realloc(row, col);
	}

	/** Copy (casting from if needed) from another matrix  */
	template <typename U>
	explicit CMatrixTemplate(const CMatrixTemplate<U>& m)
		: m_Val(nullptr), m_Rows(0), m_Cols(0)
	{
		(*this) = m;
	}

	/** Copy constructor & crop from another matrix
	 */
	CMatrixTemplate(
		const CMatrixTemplate& m, const size_t cropRowCount,
		const size_t cropColCount)
		: m_Val(nullptr), m_Rows(0), m_Cols(0)
	{
		ASSERT_(m.m_Rows >= cropRowCount);
		ASSERT_(m.m_Cols >= cropColCount);
		realloc(cropRowCount, cropColCount);
		for (size_t i = 0; i < m_Rows; i++)
			for (size_t j = 0; j < m_Cols; j++) m_Val[i][j] = m.m_Val[i][j];
	}

	/** Constructor from a given size and a C array. The array length must match
	 *cols x row.
	 * \code
	 *  const double numbers[] = {
	 *    1,2,3,
	 *    4,5,6 };
	 *	 CMatrixDouble   M(3,2, numbers);
	 * \endcode
	 */
	template <typename V, size_t N>
	CMatrixTemplate(size_t row, size_t col, V (&theArray)[N])
		: m_Val(nullptr), m_Rows(0), m_Cols(0)
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
				m_Val[i][j] = static_cast<T>(theArray[idx++]);
	}

	/** Constructor from a given size and a STL container (std::vector,
	 * std::list,...) with the initial values. The vector length must match cols
	 * x row.
	 */
	template <typename V>
	CMatrixTemplate(size_t row, size_t col, const V& theVector)
		: m_Val(nullptr), m_Rows(0), m_Cols(0)
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
				m_Val[i][j] = static_cast<T>(*(it++));
	}

	/** Destructor */
	virtual ~CMatrixTemplate() { realloc(0, 0); }

	/** Assignment operator from another matrix (possibly of a different type)
	 */
	template <typename U>
	CMatrixTemplate& operator=(const CMatrixTemplate<U>& m)
	{
		realloc(m.rows(), m.cols());
		for (size_t i = 0; i < m_Rows; i++)
			for (size_t j = 0; j < m_Cols; j++)
				(*this)(i, j) = static_cast<T>(m(i, j));
		return *this;
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
	 *CMatrixTemplate::CMatrixTemplate
	 */
	template <typename V, size_t N>
	CMatrixTemplate& operator=(V (&theArray)[N])
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
				m_Val[i][j] = static_cast<T>(theArray[idx++]);
		return *this;
	}

	/** Number of rows in the matrix
	 * \sa rows(), getColCount, nr, nc
	 */
	inline size_t rows() const { return m_Rows; }
	/** Number of columns in the matrix
	 * \sa rows(), getColCount, nr, nc
	 */
	inline size_t cols() const { return m_Cols; }
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

	/** This method just checks has no effects in this class, but raises an
	 * exception if the expected size does not match */
	inline void resize(const matrix_size_t& siz, bool zeroNewElements = false)
	{
		setSize(siz[0], siz[1], zeroNewElements);
	}

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
		return m_Val[row][col];
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
		return m_Val[row][col];
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
			return m_Val[0][ith];
		}
		else
		{
// A columns matrix:
#if defined(_DEBUG) || (MRPT_ALWAYS_CHECKS_DEBUG_MATRICES)
			if (ith >= m_Rows)
				THROW_EXCEPTION_FMT(
					"Index %u out of range!", static_cast<unsigned>(ith));
#endif
			return m_Val[ith][0];
		}
	}

	/** Subscript operator to get/set an individual element from a row or column
	 * matrix.
	 * \exception std::exception If the object is not a column or row matrix.
	 */
	inline T operator[](size_t ith) const
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
			return m_Val[0][ith];
		}
		else
		{
// A columns matrix:
#if defined(_DEBUG) || (MRPT_ALWAYS_CHECKS_DEBUG_MATRICES)
			if (ith >= m_Rows)
				THROW_EXCEPTION_FMT(
					"Index %u out of range!", static_cast<unsigned>(ith));
#endif
			return m_Val[ith][0];
		}
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
		m_Val[row][col] = v;
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
		return m_Val[row][col];
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
		return m_Val[row][col];
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
		return m_Val[row];
	}

	/** Fast but unsafe method to obtain a pointer to a given row of the matrix
	 * (Use only in critical applications)
	 */
	inline const T* get_unsafe_row(size_t row) const { return m_Val[row]; }

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
	void appendRow(const std::vector<T>& in)
	{
		if (m_Cols == 0 || m_Rows == 0)
			ASSERT_(!in.empty());
		else
			ASSERT_(in.size() == m_Cols);

		const auto row = m_Rows;
		realloc(row + 1, m_Cols = in.size());
		for (size_t i = 0; i < m_Cols; i++) m_Val[row][i] = in[i];
	}

	/** Appends a new column to the matrix from a vector.
	 * The length of the vector must match the number of rows of the matrix,
	 * unless it is (0,0).
	 * \exception std::exception On size mismatch.
	 * \sa extractCol
	 * \sa appendRow
	 */
	void appendCol(const std::vector<T>& in)
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
		for (size_t i = 0; i < m_Rows; i++) m_Val[i][m_Cols - 1] = in[i];
	}

	/** Returns a vector containing the matrix's values.
	 */
	void getAsVector(std::vector<T>& out) const
	{
		out.clear();
		out.reserve(m_Rows * m_Cols);
		for (size_t i = 0; i < m_Rows; i++)
			out.insert(out.end(), &(m_Val[i][0]), &(m_Val[i][m_Cols]));
	}

	/** Get as an Eigen-compatible Eigen::Map object  */
	template <
		typename EIGEN_MATRIX,
		typename EIGEN_MAP = Eigen::Map<
			EIGEN_MATRIX, MRPT_MAX_ALIGN_BYTES, Eigen::InnerStride<1>>>
	EIGEN_MAP asEigen()
	{
		return EIGEN_MAP(m_Val, m_Rows, m_Cols);
	}
	/** \overload (const version) */
	template <
		typename EIGEN_MATRIX,
		typename EIGEN_MAP = Eigen::Map<
			const EIGEN_MATRIX, MRPT_MAX_ALIGN_BYTES, Eigen::InnerStride<1>>>
	EIGEN_MAP asEigen() const
	{
		return EIGEN_MAP(m_Val, m_Rows, m_Cols);
	}

};  // end of class CMatrixTemplate

/** Declares a matrix of booleans (non serializable).
 *  \sa CMatrixDouble, CMatrixFloat, CMatrixB */
using CMatrixBool = CMatrixTemplate<bool>;

}  // namespace mrpt::math
