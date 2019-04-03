/* +------------------------------------------------------------------------+
   |                     Mobile Robot Programming Toolkit (MRPT)            |
   |                          https://www.mrpt.org/                         |
   |                                                                        |
   | Copyright (c) 2005-2019, Individual contributors, see AUTHORS file     |
   | See: https://www.mrpt.org/Authors - All rights reserved.               |
   | Released under BSD License. See: https://www.mrpt.org/License          |
   +------------------------------------------------------------------------+ */
#pragma once

#include <mrpt/math/MatrixVectorBase.h>

namespace mrpt::math
{
/**  Base CRTP class for all MRPT matrices.
 *
 * See MatrixVectorBase
 *
 * \sa CMatrixTemplateNumeric
 * \ingroup mrpt_math_grp
 */
template <typename Scalar, class Derived>
class MatrixBase : public MatrixVectorBase<Scalar, Derived>
{
   public:
	Derived& mbDerived() { return static_cast<Derived&>(*this); }
	const Derived& mbDerived() const
	{
		return static_cast<const Derived&>(*this);
	}

	void setDiagonal(const std::size_t N, const Scalar value)
	{
		mbDerived().resize(N, N);
		for (typename Derived::Index r = 0; r < mbDerived().rows(); r++)
			for (typename Derived::Index c = 0; c < mbDerived().cols(); c++)
				mbDerived()(r, c) = (r == c) ? value : 0;
	}
	void setDiagonal(const Scalar value)
	{
		ASSERT_EQUAL_(mbDerived().cols(), mbDerived().rows());
		setDiagonal(mbDerived().cols(), value);
	}
	void setDiagonal(const std::vector<Scalar>& diags)
	{
		const std::size_t N = diags.size();
		mbDerived().setZero(N, N);
		for (std::size_t i = 0; i < N; i++) mbDerived()(i, i) = diags[i];
	}
	void setIdentity()
	{
		ASSERT_EQUAL_(mbDerived().rows(), mbDerived().cols());
		setDiagonal(mbDerived().cols(), 1);
	}
	void setIdentity(const std::size_t N) { setDiagonal(N, 1); }

	static Derived Identity()
	{
		ASSERTMSG_(
		    Derived::RowsAtCompileTime > 0 && Derived::ColsAtCompileTime > 0,
		    "Identity() without arguments can be used only for fixed-size "
		    "matrices/vectors");
		Derived m;
		m.setIdentity();
		return m;
	}

	/** this = A*B, with A & B of the same type of this.
	 * For products of different matrix types, use the regular * operator (which
	 * requires the `<Eigen/Dense>` header) */
	void matProductOf(const Derived& A, const Derived& B);

	/** @name Operations that DO require `#include <Eigen/Dense>` in user code
	 * @{ */
	/** return: H<sup>T</sup> * C * H */
	template <class MAT_C>
	Scalar multiply_HCHt_scalar(const MAT_C& C) const
	{
		internalAssertEigenDefined<Derived>();
		return (mbDerived().asEigen() * C.asEigen() *
		        mbDerived().asEigen().transpose())
			.eval()(0, 0);
	}

	/** return:  H<sup>T</sup> * C * H */
	template <typename MAT_C>
	Scalar multiply_HtCH_scalar(const MAT_C& C) const
	{
		internalAssertEigenDefined<Derived>();
		return (mbDerived().asEigen().transpose() * C.asEigen() *
		        mbDerived().asEigen())
			.eval()(0, 0);
	}

	/** Evaluates: R = this * C * this<sup>T</sup>, or the same
	 * with `R+=` if `accumInOutput=true`. */
	template <typename MAT_C, typename MAT_R>
	void multiply_HCHt(const MAT_C& C, MAT_R& R, bool accumInOutput = false)
	{
		internalAssertEigenDefined<Derived>();
		auto res =
		    (mbDerived().asEigen() * C.asEigen() *
		     mbDerived().asEigen().transpose());
		if (accumInOutput)
			R.asEigen() += res.eval();
		else
		{
			R.resize(C.rows(), C.cols());
			R.asEigen() = res.eval();
		}
	}

	/** this = A * A<sup>T</sup> */
	template <typename MAT_A>
	void multiply_AAt(const MAT_A& A)
	{
		internalAssertEigenDefined<Derived>();
		mbDerived().resize(A.rows(), A.rows());
		mbDerived().asEigen() = (A.asEigen() * A.asEigen().transpose()).eval();
	}

	auto col(int colIdx)
	{
		internalAssertEigenDefined<Derived>();
		return mbDerived().asEigen().col(colIdx);
	}
	auto col(int colIdx) const
	{
		internalAssertEigenDefined<Derived>();
		return mbDerived().asEigen().col(colIdx);
	}

	auto row(int rowIdx)
	{
		internalAssertEigenDefined<Derived>();
		return mbDerived().asEigen().row(rowIdx);
	}
	auto row(int rowIdx) const
	{
		internalAssertEigenDefined<Derived>();
		return mbDerived().asEigen().row(rowIdx);
	}
	/** @} */

	/** @name Standalone operations (do NOT require `#include <Eigen/Dense>`)
	 * @{ */
	/** Determinant of matrix. */
	Scalar det() const;

	/** Returns the inverse of a general matrix using LU */
	Derived inverse() const;

	/** Returns the inverse of a symmetric matrix using LLt */
	Derived inverse_LLt() const;

	/** Finds the rank of the matrix via LU decomposition.
	 * Uses Eigen's default threshold unless `threshold>0`. */
	int rank(Scalar threshold = 0) const;

	/** Cholesky M=U<sup>T</sup> * U decomposition for symmetric matrix
	 * (upper-half of the matrix is actually ignored.
	 * \return false if Cholesky fails
	 */
	bool chol(Derived& U) const;

	/** Computes the eigenvectors and eigenvalues for a square, general matrix.
	 * Use eig_symmetric() for symmetric matrices for better accuracy and
	 * performance.
	 * Eigenvectors are the columns of the returned matrix, and their order
	 * matches that of returned eigenvalues.
	 * \param[in] sorted If true, eigenvalues (and eigenvectors) will be sorted
	 * in ascending order.
	 * \param[out] eVecs The container where eigenvectors will be stored.
	 * \param[out] eVals The container where eigenvalues will be stored.
	 * \return false if eigenvalues could not be determined.
	 */
	bool eig(
		Derived& eVecs, std::vector<Scalar>& eVals, bool sorted = true) const;

	/** Read: eig() */
	bool eig_symmetric(
		Derived& eVecs, std::vector<Scalar>& eVals, bool sorted = true) const;
	/** Removes columns of the matrix.
	 * This "unsafe" version assumes indices sorted in ascending order. */
	void unsafeRemoveColumns(const std::vector<std::size_t>& idxs);

	/** Removes columns of the matrix. Indices may be unsorted and duplicated */
	void removeColumns(const std::vector<std::size_t>& idxsToRemove);

	/** Removes rows of the matrix.
	 * This "unsafe" version assumes indices sorted in ascending order. */
	void unsafeRemoveRows(const std::vector<std::size_t>& idxs);

	/** Removes rows of the matrix. Indices may be unsorted and duplicated */
	void removeRows(const std::vector<std::size_t>& idxsToRemove);

	/** Copies the given input submatrix/vector into this matrix/vector,
	 * starting at the given top-left coordinates. */
	template <typename OTHERMATVEC>
	void insertMatrix(
		const int row_start, const int col_start, const OTHERMATVEC& submat)
	{
		ASSERT_BELOW_(row_start + submat.rows(), mbDerived().rows());
		ASSERT_BELOW_(col_start + submat.cols(), mbDerived().cols());
		for (int r = 0; r < submat.rows(); r++)
			for (int c = 0; c < submat.cols(); c++)
				mbDerived()(r + row_start, c + col_start) = submat(r, c);
	}

	/** @} */
};

}  // namespace mrpt::math
