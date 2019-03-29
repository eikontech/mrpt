/* +------------------------------------------------------------------------+
   |                     Mobile Robot Programming Toolkit (MRPT)            |
   |                          https://www.mrpt.org/                         |
   |                                                                        |
   | Copyright (c) 2005-2019, Individual contributors, see AUTHORS file     |
   | See: https://www.mrpt.org/Authors - All rights reserved.               |
   | Released under BSD License. See: https://www.mrpt.org/License          |
   +------------------------------------------------------------------------+ */
#pragma once

#include <mrpt/core/aligned_std_vector.h>
#include <mrpt/core/exceptions.h>  // ASSERT_()
#include <mrpt/math/matrix_size_t.h>
#include <mrpt/serialization/serialization_frwds.h>
#include <array>
#include <cstring>  // memset()
#include <type_traits>

namespace mrpt::math
{
/** Template for column vectors of dynamic size, compatible with Eigen.
 *
 * \note For a complete introduction to Matrices and vectors in MRPT, see:
 *http://www.mrpt.org/Matrices_vectors_arrays_and_Linear_Algebra_MRPT_and_Eigen_classes
 * \sa CVectorDynamic, CMatrixFixed, CVectorFixed
 * \ingroup mrpt_math_grp
 */
template <class T>
class CVectorDynamic
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
	mrpt::aligned_std_vector<T> m_data;

	/** Internal use only: It reallocs the memory for the 2D matrix, maintaining
	 * the previous contents if posible.
	 */
	void realloc(const size_t new_len, bool newElementsToZero = false)
	{
		const auto old_len = m_data.size();
		if (new_len == old_len) return;
		m_data.resize(new_len);
		if (newElementsToZero && new_len > old_len)
		{
			if constexpr (std::is_trivial_v<T>)
			    ::memset(&m_data[old_len], 0, sizeof(T) * (new_len - old_len));
			else
			    for (size_t k = old_len; k < new_len; k++) m_data[k] = T();
		}
	}

   public:
	void fill(const T& val) { std::fill(m_data.begin(), m_data.end(), val); }

	void swap(CVectorDynamic<T>& o) { std::swap(m_data); }

	CVectorDynamic() = default;

	/** Initializes to a vector of "N" zeros */
	CVectorDynamic(size_t N, bool initZero = true) { realloc(N, initZero); }

	/** Copy (casting from if needed) from another matrix  */
	template <typename U>
	explicit CVectorDynamic(const CVectorDynamic<U>& m)
	{
		(*this) = m;
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
	template <
	    typename ARRAY, typename = std::enable_if_t<std::is_array_v<ARRAY>>>
	CVectorDynamic(const ARRAY& data)
	{
		const auto N = std::size(data);
		static_assert(N != 0, "Empty array!");
		realloc(N);
		for (size_t i = 0; i < N; i++) m_data[i] = static_cast<T>(data[i]);
	}

	/** Number of rows in the vector */
	inline std::size_t rows() const { return m_data.size(); }

	/** Number of columns in the matrix (always 1) */
	inline std::size_t cols() const { return 1; }

	/** Get a 2-vector with [NROWS NCOLS] (as in MATLAB command size(x)) */
	inline std::size_t size() const { return m_data.size(); }

	/** Changes the size of matrix, maintaining the previous contents. */
	void setSize(size_t row, size_t col, bool zeroNewElements = false)
	{
		ASSERT_(col == 1);
		realloc(row, zeroNewElements);
	}
	inline void resize(std::size_t N, bool zeroNewElements = false)
	{
		setSize(N, zeroNewElements);
	}

	/** Subscript operator to get/set individual elements
	 */
	inline T& operator()(size_t row, size_t col)
	{
#if defined(_DEBUG) || (MRPT_ALWAYS_CHECKS_DEBUG_MATRICES)
		if (row >= m_data.size() || col > 0)
			THROW_EXCEPTION(format(
			    "Indexes (%lu,%lu) out of range. Vector is %lux%lu",
				static_cast<unsigned long>(row),
				static_cast<unsigned long>(col),
			    static_cast<unsigned long>(m_data.size()),
			    static_cast<unsigned long>(1)));
#endif
		return m_data[row];
	}

	/** Subscript operator to get individual elements
	 */
	inline const T& operator()(size_t row, size_t col) const
	{
#if defined(_DEBUG) || (MRPT_ALWAYS_CHECKS_DEBUG_MATRICES)
		if (row >= m_data.size() || col > 0)
			THROW_EXCEPTION(format(
			    "Indexes (%lu,%lu) out of range. Vector is %lux%lu",
			    static_cast<unsigned long>(row),
			    static_cast<unsigned long>(col),
			    static_cast<unsigned long>(m_data.size()),
			    static_cast<unsigned long>(1)));
#endif
		return m_data[row];
	}

	/** Subscript operator to get/set an individual element from a row or column
	 * matrix.
	 * \exception std::exception If the object is not a column or row matrix.
	 */
	inline T& operator[](size_t ith)
	{
#if defined(_DEBUG) || (MRPT_ALWAYS_CHECKS_DEBUG_MATRICES)
		if (ith >= m_data.size())
			THROW_EXCEPTION_FMT(
			    "Index %u out of range!", static_cast<unsigned>(ith));
#endif
		return m_data[ith];
	}

	/// \overload
	inline const T& operator[](size_t ith) const
	{
#if defined(_DEBUG) || (MRPT_ALWAYS_CHECKS_DEBUG_MATRICES)
		if (ith >= m_data.size())
			THROW_EXCEPTION_FMT(
			    "Index %u out of range!", static_cast<unsigned>(ith));
#endif
		return m_data[ith];
	}

	/** Get as an Eigen-compatible Eigen::Map object  */
	template <
	    typename EIGEN_VECTOR = Eigen::Matrix<T, -1, 1, 0, -1, 1>,
	    typename EIGEN_MAP = Eigen::Map<
	        EIGEN_VECTOR, MRPT_MAX_ALIGN_BYTES, Eigen::InnerStride<1>>>
	EIGEN_MAP asEigen()
	{
		return EIGEN_MAP(&m_data[0], m_data.size());
	}
	/** \overload (const version) */
	template <
	    typename EIGEN_VECTOR = Eigen::Matrix<T, -1, 1, 0, -1, 1>,
	    typename EIGEN_MAP = Eigen::Map<
	        const EIGEN_VECTOR, MRPT_MAX_ALIGN_BYTES, Eigen::InnerStride<1>>>
	EIGEN_MAP asEigen() const
	{
		return EIGEN_MAP(&m_data[0], m_data.size());
	}
};

using CVectorFloat = CVectorDynamic<float>;
using CVectorDouble = CVectorDynamic<double>;

mrpt::serialization::CArchive& operator<<(
    mrpt::serialization::CArchive& s, const CVectorFloat& a);
mrpt::serialization::CArchive& operator<<(
    mrpt::serialization::CArchive& s, const CVectorDouble& a);
mrpt::serialization::CArchive& operator>>(
    mrpt::serialization::CArchive& in, CVectorDouble& a);
mrpt::serialization::CArchive& operator>>(
    mrpt::serialization::CArchive& in, CVectorFloat& a);

}  // namespace mrpt::math

namespace mrpt::typemeta
{
// Extensions to mrpt::typemeta::TTypeName for matrices:
template <typename T>
struct TTypeName<mrpt::math::CVectorDynamic<T>>
{
	static auto get()
	{
		return literal("CVectorDynamic<") + TTypeName<T>::get() + literal(">");
	}
};
}  // namespace mrpt::typemeta