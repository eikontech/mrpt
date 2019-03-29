/* +------------------------------------------------------------------------+
   |                     Mobile Robot Programming Toolkit (MRPT)            |
   |                          https://www.mrpt.org/                         |
   |                                                                        |
   | Copyright (c) 2005-2019, Individual contributors, see AUTHORS file     |
   | See: https://www.mrpt.org/Authors - All rights reserved.               |
   | Released under BSD License. See: https://www.mrpt.org/License          |
   +------------------------------------------------------------------------+ */
#pragma once

#include <mrpt/math/eigen_frwds.h>  // forward declarations
#include <mrpt/math/math_frwds.h>  // forward declarations

#include <mrpt/math/CMatrixDynamic.h>
#include <mrpt/math/CMatrixFixed.h>
#include <mrpt/math/CVectorDynamic.h>

#include <mrpt/math/ops_containers.h>  // Many generic operations

/** \file ops_matrices.h
 * This file implements miscelaneous matrix and matrix/vector operations, and
 * internal functions in mrpt::math::detail
 */
namespace mrpt
{
namespace math
{
/** \addtogroup container_ops_grp
 *  @{ */

/** Transpose operator for matrices */
template <class Derived>
inline const typename Eigen::MatrixBase<Derived>::AdjointReturnType operator~(
	const Eigen::MatrixBase<Derived>& m)
{
	return m.adjoint();
}

/** Unary inversion operator. */
template <class Derived>
inline typename Eigen::MatrixBase<Derived>::PlainObject operator!(
	const Eigen::MatrixBase<Derived>& m)
{
	return m.inv();
}

/** @} */  // end MRPT matrices operators

/** R = H * C * H^t (with C symmetric) */
template <typename MAT_H, typename MAT_C, typename MAT_R>
inline void multiply_HCHt(
	const MAT_H& H, const MAT_C& C, MAT_R& R,
	bool accumResultInOutput)  //	bool allow_submatrix_mult)
{
	if (accumResultInOutput)
		R += ((H * C.template selfadjointView<Eigen::Lower>()).eval() *
			  H.adjoint())
				 .eval()
				 .template selfadjointView<Eigen::Lower>();
	else
		R = ((H * C.template selfadjointView<Eigen::Lower>()).eval() *
			 H.adjoint())
				.eval()
				.template selfadjointView<Eigen::Lower>();
}

/** r (a scalar) = H * C * H^t (with a vector H and a symmetric matrix C) */
template <typename VECTOR_H, typename MAT_C>
typename MAT_C::Scalar multiply_HCHt_scalar(const VECTOR_H& H, const MAT_C& C)
{
	return (H.matrix().adjoint() * C * H.matrix()).eval()(0, 0);
}

/** R = H^t * C * H  (with C symmetric) */
template <typename MAT_H, typename MAT_C, typename MAT_R>
void multiply_HtCH(
	const MAT_H& H, const MAT_C& C, MAT_R& R,
	bool accumResultInOutput)  // bool allow_submatrix_mult)
{
	if (accumResultInOutput)
		R +=
			((H.adjoint() * C.template selfadjointView<Eigen::Lower>()).eval() *
			 H)
				.eval()
				.template selfadjointView<Eigen::Lower>();
	else
		R = ((H.adjoint() * C.template selfadjointView<Eigen::Lower>()).eval() *
			 H)
				.eval()
				.template selfadjointView<Eigen::Lower>();
}

/** Computes the mean vector and covariance from a list of samples in an NxM
 * matrix, where each row is a sample, so the covariance is MxM.
 * \param v The set of data as a NxM matrix, of types: CMatrixDynamic
 * or CMatrixFixed
 * \param out_mean The output M-vector for the estimated mean.
 * \param out_cov The output MxM matrix for the estimated covariance matrix,
 * this can also be either a fixed-size of dynamic size matrix.
 * \sa mrpt::math::meanAndCovVec, math::mean,math::stddev, math::cov
 */
template <class MAT_IN, class VECTOR, class MAT_OUT>
void meanAndCovMat(const MAT_IN& v, VECTOR& out_mean, MAT_OUT& out_cov)
{
	const size_t N = v.rows();
	ASSERTMSG_(N > 0, "The input matrix contains no elements");
	const double N_inv = 1.0 / N;

	const size_t M = v.cols();
	ASSERTMSG_(M > 0, "The input matrix contains rows of length 0");

	// First: Compute the mean
	out_mean.assign(M, 0);
	for (size_t i = 0; i < N; i++)
		for (size_t j = 0; j < M; j++) out_mean[j] += v.coeff(i, j);
	out_mean *= N_inv;

	// Second: Compute the covariance
	//  Save only the above-diagonal part, then after averaging
	//  duplicate that part to the other half.
	out_cov.setZero(M, M);
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < M; j++)
			out_cov(j, j) += square(v(i, j) - out_mean[j]);

		for (size_t j = 0; j < M; j++)
			for (size_t k = j + 1; k < M; k++)
				out_cov(j, k) +=
					(v(i, j) - out_mean[j]) * (v(i, k) - out_mean[k]);
	}
	for (size_t j = 0; j < M; j++)
		for (size_t k = j + 1; k < M; k++) out_cov(k, j) = out_cov(j, k);
	out_cov *= N_inv;
}

/** Computes the covariance matrix from a list of samples in an NxM matrix,
 * where each row is a sample, so the covariance is MxM.
 * \param v The set of data, as a NxM matrix.
 * \param out_cov The output MxM matrix for the estimated covariance matrix.
 * \sa math::mean,math::stddev, math::cov
 */
template <class MATRIX>
inline Eigen::Matrix<
	typename MATRIX::Scalar, MATRIX::ColsAtCompileTime,
    MATRIX::ColsAtCompileTime, 0, MATRIX::ColsAtCompileTime, 1>
	cov(const MATRIX& v)
{
	CVectorDynamic<double> m;
	Eigen::Matrix<
		typename MATRIX::Scalar, MATRIX::ColsAtCompileTime,
	    MATRIX::ColsAtCompileTime, 1, MATRIX::ColsAtCompileTime,
	    MATRIX::ColsAtCompileTime>
		C;
	meanAndCovMat(v, m, C);
	return C;
}

/** A useful macro for saving matrixes to a file while debugging. */
#define SAVE_MATRIX(M) M.saveToTextFile(#M ".txt");

/** Only for vectors/arrays "v" of length3, compute out = A * Skew(v), where
 * Skew(v) is the skew symmetric matric generated from \a v (see
 * mrpt::math::skew_symmetric3)
 */
template <class MAT_A, class SKEW_3VECTOR, class MAT_OUT>
void multiply_A_skew3(const MAT_A& A, const SKEW_3VECTOR& v, MAT_OUT& out)
{
	MRPT_START
	ASSERT_EQUAL_(A.cols(), 3);
	ASSERT_EQUAL_(v.size(), 3);
	const size_t N = A.rows();
	out.setSize(N, 3);
	for (size_t i = 0; i < N; i++)
	{
		out.set_unsafe(i, 0, A(i, 1) * v[2] - A(i, 2) * v[1]);
		out.set_unsafe(i, 1, -A(i, 0) * v[2] + A(i, 2) * v[0]);
		out.set_unsafe(i, 2, A(i, 0) * v[1] - A(i, 1) * v[0]);
	}
	MRPT_END
}

/** Only for vectors/arrays "v" of length3, compute out = Skew(v) * A, where
 * Skew(v) is the skew symmetric matric generated from \a v (see
 * mrpt::math::skew_symmetric3)
 */
template <class SKEW_3VECTOR, class MAT_A, class MAT_OUT>
void multiply_skew3_A(const SKEW_3VECTOR& v, const MAT_A& A, MAT_OUT& out)
{
	MRPT_START
	ASSERT_EQUAL_(A.rows(), 3);
	ASSERT_EQUAL_(v.size(), 3);
	const size_t N = A.cols();
	out.setSize(3, N);
	for (size_t i = 0; i < N; i++)
	{
		out.set_unsafe(0, i, -A(1, i) * v[2] + A(2, i) * v[1]);
		out.set_unsafe(1, i, A(0, i) * v[2] - A(2, i) * v[0]);
		out.set_unsafe(2, i, -A(0, i) * v[1] + A(1, i) * v[0]);
	}
	MRPT_END
}

/**  @} */  // end of grouping

}  // namespace math
}  // namespace mrpt
