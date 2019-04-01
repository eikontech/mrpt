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
   protected:
	Derived& derived() { return static_cast<Derived&>(*this); }
	const Derived& derived() const
	{
		return static_cast<const Derived&>(*this);
	}

   public:
	/*! Fill all the elements with a given value (Note: named "fillAll" since
	 * "fill" will be used by child classes) */
	void fill(const Scalar& val)
	{
		std::fill(derived().begin(), derived().end(), val);
	}

	template <
		typename = std::enable_if_t<
			(Derived::RowsAtCompileTime > 0 && Derived::ColsAtCompileTime > 0)>>
	static Derived Zero()
	{
		Derived m;
		m.setZero();
		return m;
	}

	inline void setZero() { fill(0); }
	inline void setZero(size_t nrows, size_t ncols)
	{
		derived().resize(nrows, ncols);
		fill(0);
	}

	void setDiagonal(const std::size_t N, const Scalar value)
	{
		derived().resize(N, N);
		for (std::size_t r = 0; r < derived().rows(); r++)
			for (std::size_t c = 0; c < derived().cols(); c++)
				derived()(r, c) = (r == c) ? value : 0;
	}
	void setDiagonal(const Scalar value)
	{
		ASSERT_EQUAL_(derived().cols(), derived().rows());
		setDiagonal(derived().cols(), value);
	}
	void setIdentity()
	{
		ASSERT_EQUAL_(derived().rows(), derived().cols());
		setDiagonal(derived().cols(), 1);
	}
	void setIdentity(const std::size_t N) { setDiagonal(N, 1); }

	bool isSquare() const { return derived().cols() == derived().rows; }
};
}  // namespace mrpt::math
