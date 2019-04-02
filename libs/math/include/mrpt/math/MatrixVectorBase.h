/* +------------------------------------------------------------------------+
   |                     Mobile Robot Programming Toolkit (MRPT)            |
   |                          https://www.mrpt.org/                         |
   |                                                                        |
   | Copyright (c) 2005-2019, Individual contributors, see AUTHORS file     |
   | See: https://www.mrpt.org/Authors - All rights reserved.               |
   | Released under BSD License. See: https://www.mrpt.org/License          |
   +------------------------------------------------------------------------+ */
#pragma once

#include <algorithm>  // fill()
#include <cstddef>  // size_t

namespace mrpt::math
{
/**  Base CRTP class for all MRPT vectors and matrices.
 *
 * \sa CMatrixTemplateNumeric
 * \ingroup mrpt_math_grp
 */
template <typename Scalar, class Derived>
class MatrixVectorBase
{
   public:
	Derived& mvbDerived() { return static_cast<Derived&>(*this); }
	const Derived& mvbDerived() const
	{
		return static_cast<const Derived&>(*this);
	}

	/*! Fill all the elements with a given value (Note: named "fillAll" since
	 * "fill" will be used by child classes) */
	void fill(const Scalar& val)
	{
		std::fill(mvbDerived().begin(), mvbDerived().end(), val);
	}

	static Derived Zero()
	{
		static_assert(
		    Derived::RowsAtCompileTime > 0 && Derived::ColsAtCompileTime > 0,
		    "Zero() without arguments can be used only for fixed-size "
		    "matrices/vectors");
		Derived m;
		m.setZero();
		return m;
	}
	static Derived Identity()
	{
		static_assert(
		    Derived::RowsAtCompileTime > 0 && Derived::ColsAtCompileTime > 0,
		    "Zero() without arguments can be used only for fixed-size "
		    "matrices/vectors");
		Derived m;
		m.setIdentity();
		return m;
	}

	inline void setZero() { fill(0); }
	inline void setZero(size_t nrows, size_t ncols)
	{
		mvbDerived().resize(nrows, ncols);
		fill(0);
	}
	inline void assign(const std::size_t N, const Scalar value)
	{
		mvbDerived().resize(N);
		fill(value);
	}
	void setDiagonal(const std::size_t N, const Scalar value)
	{
		mvbDerived().resize(N, N);
		for (std::size_t r = 0; r < mvbDerived().rows(); r++)
			for (std::size_t c = 0; c < mvbDerived().cols(); c++)
				mvbDerived()(r, c) = (r == c) ? value : 0;
	}
	void setDiagonal(const Scalar value)
	{
		ASSERT_EQUAL_(mvbDerived().cols(), mvbDerived().rows());
		setDiagonal(mvbDerived().cols(), value);
	}
	void setIdentity()
	{
		ASSERT_EQUAL_(mvbDerived().rows(), mvbDerived().cols());
		setDiagonal(mvbDerived().cols(), 1);
	}
	void setIdentity(const std::size_t N) { setDiagonal(N, 1); }

	bool isSquare() const { return mvbDerived().cols() == mvbDerived().rows; }

	Scalar det() const { return mvbDerived().det(); }

	Scalar minCoeff() const { return mvbDerived().minCoeff(); }
	Scalar maxCoeff() const { return mvbDerived().maxCoeff(); }

	/** return: H<sup>T</sup> * C * H */
	template <class MAT_C>
	Scalar multiply_HCHt_scalar(const MAT_C& C) const
	{
		return (mvbDerived().asEigen() * C * mvbDerived().asEigen().transpose())
		    .eval()(0, 0);
	}

	/** return:  H<sup>T</sup> * C * H */
	template <typename MAT_C>
	Scalar multiply_HtCH_scalar(const MAT_C& C) const
	{
		return (mvbDerived().asEigen().transpose() * C * mvbDerived().asEigen())
		    .eval()(0, 0);
	}

	template <int BLOCK_ROWS, int BLOCK_COLS>
	auto block(int start_row = 0, int start_col = 0)
	{
		return mvbDerived().asEigen().template block<BLOCK_ROWS, BLOCK_COLS>(
		    start_row, start_col);
	}
	template <int BLOCK_ROWS, int BLOCK_COLS>
	auto block(int start_row = 0, int start_col = 0) const
	{
		return mvbDerived().asEigen().template block<BLOCK_ROWS, BLOCK_COLS>(
		    start_row, start_col);
	}

	template <int NUM_ELEMENTS>
	auto tail()
	{
		return mvbDerived().asEigen().template tail<NUM_ELEMENTS>();
	}
	template <int NUM_ELEMENTS>
	auto tail() const
	{
		return mvbDerived().asEigen().template tail<NUM_ELEMENTS>();
	}
	template <int NUM_ELEMENTS>
	auto head()
	{
		return mvbDerived().asEigen().template head<NUM_ELEMENTS>();
	}
	template <int NUM_ELEMENTS>
	auto head() const
	{
		return mvbDerived().asEigen().template head<NUM_ELEMENTS>();
	}
	auto col(int colIdx) { return mvbDerived().asEigen().col(colIdx); }
	auto col(int colIdx) const { return mvbDerived().asEigen().col(colIdx); }

	auto row(int rowIdx) { return mvbDerived().asEigen().row(rowIdx); }
	auto row(int rowIdx) const { return mvbDerived().asEigen().row(rowIdx); }

	auto transpose() { return mvbDerived().asEigen().transpose(); }
	auto transpose() const { return mvbDerived().asEigen().transpose(); }

	auto operator-() const { return -mvbDerived().asEigen(); }

	void operator+=(Scalar s) { mvbDerived().asEigen() += s; }
	void operator-=(Scalar s) { mvbDerived().asEigen() -= s; }
	void operator*=(Scalar s) { mvbDerived().asEigen() *= s; }

	template <typename S2, class D2>
	auto operator+(const MatrixVectorBase<S2, D2>& m2) const
	{
		return mvbDerived().asEigen() + m2.mvbDerived().asEigen();
	}
	template <typename S2, class D2>
	void operator+=(const MatrixVectorBase<S2, D2>& m2) const
	{
		mvbDerived().asEigen() += m2.mvbDerived().asEigen();
	}

	template <typename S2, class D2>
	auto operator-(const MatrixVectorBase<S2, D2>& m2) const
	{
		return mvbDerived().asEigen() - m2.mvbDerived().asEigen();
	}
	template <typename S2, class D2>
	void operator-=(const MatrixVectorBase<S2, D2>& m2) const
	{
		mvbDerived().asEigen() -= m2.mvbDerived().asEigen();
	}

	template <typename S2, class D2>
	auto operator*(const MatrixVectorBase<S2, D2>& m2) const
	{
		return mvbDerived().asEigen() * m2.mvbDerived().asEigen();
	}
	auto operator*(const Scalar s) const { return mvbDerived().asEigen() * s; }
};

}  // namespace mrpt::math
