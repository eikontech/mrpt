/* +------------------------------------------------------------------------+
   |                     Mobile Robot Programming Toolkit (MRPT)            |
   |                          https://www.mrpt.org/                         |
   |                                                                        |
   | Copyright (c) 2005-2019, Individual contributors, see AUTHORS file     |
   | See: https://www.mrpt.org/Authors - All rights reserved.               |
   | Released under BSD License. See: https://www.mrpt.org/License          |
   +------------------------------------------------------------------------+ */
#pragma once

#include <mrpt/core/exceptions.h>
#include <mrpt/core/is_defined.h>
#include <mrpt/core/optional_ref.h>
#include <algorithm>  // fill()
#include <cstddef>  // size_t
#include <iosfwd>
#include <string>
#include <vector>

namespace mrpt::math
{
template <class T>
class CMatrixDynamic;
template <typename T, std::size_t ROWS, std::size_t COLS>
class CMatrixFixed;

#if defined _WIN32 || defined __CYGWIN__
#define DLL_PUBLIC  //__attribute__((dllexport))
#else
#define DLL_PUBLIC __attribute__((visibility("default")))
#endif

/*! Selection of the number format in MatrixVectorBase::saveToTextFile()
 * \ingroup mrpt_math_grp */
enum TMatrixTextFileFormat
{
	/** engineering format '%e' */
	MATRIX_FORMAT_ENG = 0,
	/** fixed floating point '%f' */
	MATRIX_FORMAT_FIXED = 1,
	/** intergers '%i' */
	MATRIX_FORMAT_INT = 2
};

/**  Base CRTP class for all MRPT vectors and matrices.
 *
 * Template methods whose implementation is not in this header file are
 * explicitly instantiated and exported for these derived types:
 * - CMatrixDynamic and CVectorDynamic, for `Scalar`: `float`, `double`
 * - CMatrixFixed and CVectorFixed, for `Scalar`: `float`, `double`, and sizes
 * from 2x2 to 6x6 and 2 to 6, respectively.
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

	/** @name Initialization methods
	 * @{ */

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
		for (typename Derived::Index r = 0; r < mvbDerived().rows(); r++)
			for (typename Derived::Index c = 0; c < mvbDerived().cols(); c++)
				mvbDerived()(r, c) = (r == c) ? value : 0;
	}
	void setDiagonal(const Scalar value)
	{
		ASSERT_EQUAL_(mvbDerived().cols(), mvbDerived().rows());
		setDiagonal(mvbDerived().cols(), value);
	}
	void setDiagonal(const std::vector<Scalar>& diags)
	{
		const std::size_t N = diags.size();
		mvbDerived().setZero(N, N);
		for (std::size_t i = 0; i < N; i++) mvbDerived()(i, i) = diags[i];
	}
	void setIdentity()
	{
		ASSERT_EQUAL_(mvbDerived().rows(), mvbDerived().cols());
		setDiagonal(mvbDerived().cols(), 1);
	}
	void setIdentity(const std::size_t N) { setDiagonal(N, 1); }

	/** this = A*B, with A & B of the same type of this.
	 * For products of different matrix types, use the regular * operator (which
	 * requires the `<Eigen/Dense>` header) */
	void matProductOf(const Derived& A, const Derived& B);
	/** @} */

	/** @name Operations that DO require `#include <Eigen/Dense>` in user code
	 * @{ */

	/** return: H<sup>T</sup> * C * H */
	template <class MAT_C>
	Scalar multiply_HCHt_scalar(const MAT_C& C) const
	{
		internalAssertEigenDefined<Derived>();
		return (mvbDerived().asEigen() * C.asEigen() *
				mvbDerived().asEigen().transpose())
			.eval()(0, 0);
	}

	/** return:  H<sup>T</sup> * C * H */
	template <typename MAT_C>
	Scalar multiply_HtCH_scalar(const MAT_C& C) const
	{
		internalAssertEigenDefined<Derived>();
		return (mvbDerived().asEigen().transpose() * C.asEigen() *
				mvbDerived().asEigen())
			.eval()(0, 0);
	}

	/** Evaluates: R = this * C * this<sup>T</sup>, or the same
	 * with `R+=` if `accumInOutput=true`. */
	template <typename MAT_C, typename MAT_R>
	void multiply_HCHt(const MAT_C& C, MAT_R& R, bool accumInOutput = false)
	{
		internalAssertEigenDefined<Derived>();
		auto res =
			(mvbDerived().asEigen() * C.asEigen() *
			 mvbDerived().asEigen().transpose());
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
		mvbDerived().resize(A.rows(), A.rows());
		mvbDerived().asEigen() = (A.asEigen() * A.asEigen().transpose()).eval();
	}

	template <int BLOCK_ROWS, int BLOCK_COLS>
	auto block(int start_row = 0, int start_col = 0)
	{
		internalAssertEigenDefined<Derived>();
		return mvbDerived().asEigen().template block<BLOCK_ROWS, BLOCK_COLS>(
			start_row, start_col);
	}
	template <int BLOCK_ROWS, int BLOCK_COLS>
	auto block(int start_row = 0, int start_col = 0) const
	{
		internalAssertEigenDefined<Derived>();
		return mvbDerived().asEigen().template block<BLOCK_ROWS, BLOCK_COLS>(
			start_row, start_col);
	}

	template <int NUM_ELEMENTS>
	auto tail()
	{
		internalAssertEigenDefined<Derived>();
		return mvbDerived().asEigen().template tail<NUM_ELEMENTS>();
	}
	template <int NUM_ELEMENTS>
	auto tail() const
	{
		internalAssertEigenDefined<Derived>();
		return mvbDerived().asEigen().template tail<NUM_ELEMENTS>();
	}
	template <int NUM_ELEMENTS>
	auto head()
	{
		internalAssertEigenDefined<Derived>();
		return mvbDerived().asEigen().template head<NUM_ELEMENTS>();
	}
	template <int NUM_ELEMENTS>
	auto head() const
	{
		internalAssertEigenDefined<Derived>();
		return mvbDerived().asEigen().template head<NUM_ELEMENTS>();
	}
	auto col(int colIdx)
	{
		internalAssertEigenDefined<Derived>();
		return mvbDerived().asEigen().col(colIdx);
	}
	auto col(int colIdx) const
	{
		internalAssertEigenDefined<Derived>();
		return mvbDerived().asEigen().col(colIdx);
	}

	auto row(int rowIdx)
	{
		internalAssertEigenDefined<Derived>();
		return mvbDerived().asEigen().row(rowIdx);
	}
	auto row(int rowIdx) const
	{
		internalAssertEigenDefined<Derived>();
		return mvbDerived().asEigen().row(rowIdx);
	}

	auto transpose()
	{
		internalAssertEigenDefined<Derived>();
		return mvbDerived().asEigen().transpose();
	}
	auto transpose() const
	{
		internalAssertEigenDefined<Derived>();
		return mvbDerived().asEigen().transpose();
	}

	auto array()
	{
		internalAssertEigenDefined<Derived>();
		return mvbDerived().asEigen().array();
	}
	auto array() const
	{
		internalAssertEigenDefined<Derived>();
		return mvbDerived().asEigen().array();
	}

	auto operator-() const
	{
		internalAssertEigenDefined<Derived>();
		return -mvbDerived().asEigen();
	}

	template <typename S2, class D2>
	auto operator+(const MatrixVectorBase<S2, D2>& m2) const
	{
		internalAssertEigenDefined<Derived>();
		return mvbDerived().asEigen() + m2.mvbDerived().asEigen();
	}
	template <typename S2, class D2>
	void operator+=(const MatrixVectorBase<S2, D2>& m2)
	{
		internalAssertEigenDefined<Derived>();
		mvbDerived().asEigen() += m2.mvbDerived().asEigen();
	}

	template <typename S2, class D2>
	auto operator-(const MatrixVectorBase<S2, D2>& m2) const
	{
		internalAssertEigenDefined<Derived>();
		return mvbDerived().asEigen() - m2.mvbDerived().asEigen();
	}
	template <typename S2, class D2>
	void operator-=(const MatrixVectorBase<S2, D2>& m2)
	{
		internalAssertEigenDefined<Derived>();
		mvbDerived().asEigen() -= m2.mvbDerived().asEigen();
	}

	template <typename S2, class D2>
	auto operator*(const MatrixVectorBase<S2, D2>& m2) const
	{
		internalAssertEigenDefined<Derived>();
		return mvbDerived().asEigen() * m2.mvbDerived().asEigen();
	}
	auto operator*(const Scalar s) const
	{
		internalAssertEigenDefined<Derived>();
		return mvbDerived().asEigen() * s;
	}

	/** @} */

	/** @name Standalone operations (do NOT require `#include <Eigen/Dense>`)
	 * @{ */

	Scalar& coeffRef(int r, int c) { return mvbDerived()(r, c); }
	const Scalar& coeff(int r, int c) const { return mvbDerived()(r, c); }

	/** Minimum value in the matrix/vector */
	DLL_PUBLIC Scalar minCoeff() const;

	/** Maximum value in the matrix/vector */
	DLL_PUBLIC Scalar maxCoeff() const;

	/** returns true if matrix is NxN */
	bool isSquare() const { return mvbDerived().cols() == mvbDerived().rows(); }

	/** returns true if matrix/vector has size=0 */
	bool empty() const
	{
		return mvbDerived().cols() == 0 && mvbDerived().rows == 0;
	}

	DLL_PUBLIC void operator+=(Scalar s);
	DLL_PUBLIC void operator-=(Scalar s);
	DLL_PUBLIC void operator*=(Scalar s);

	/** Determinant of matrix. */
	DLL_PUBLIC Scalar det() const;

	/** Sum of all elements in matrix/vector. */
	DLL_PUBLIC Scalar sum() const;

	/** Returns the inverse of a general matrix using LU */
	DLL_PUBLIC Derived inverse() const;

	/** Returns the inverse of a symmetric matrix using LLt */
	DLL_PUBLIC Derived inverse_LLt() const;

	/** Finds the rank of the matrix via LU decomposition.
	 * Uses Eigen's default threshold unless `threshold>0`. */
	DLL_PUBLIC int rank(Scalar threshold = 0) const;

	/** Cholesky M=U<sup>T</sup> * U decomposition for symmetric matrix
	 * (upper-half of the matrix is actually ignored.
	 * \return false if Cholesky fails
	 */
	DLL_PUBLIC bool chol(Derived& U) const;

	/** Returns a string representation of the vector/matrix, using Eigen's
	 * default settings. */
	DLL_PUBLIC std::string asStr() const;

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
	DLL_PUBLIC bool eig(
		Derived& eVecs, std::vector<Scalar>& eVals, bool sorted = true) const;

	/** Read: eig() */
	DLL_PUBLIC bool eig_symmetric(
		Derived& eVecs, std::vector<Scalar>& eVals, bool sorted = true) const;

	/** Reads a matrix from a string in Matlab-like format, for example:
	 *  "[1 0 2; 0 4 -1]"
	 * The string must start with '[' and end with ']'. Rows are separated by
	 * semicolons ';' and columns in each row by one or more whitespaces
	 * ' ', tabs '\t' or commas ','.
	 *
	 * This format is also used for CConfigFile::read_matrix.
	 *
	 * \return true on success. false if the string is malformed, and then the
	 * matrix will be resized to 0x0.
	 * \sa inMatlabFormat, CConfigFile::read_matrix
	 */
	DLL_PUBLIC bool fromMatlabStringFormat(
		const std::string& s,
		mrpt::optional_ref<std::ostream> dump_errors_here = std::nullopt);

	/** Exports the matrix as a string compatible with Matlab/Octave.
	 * \sa fromMatlabStringFormat()
	 */
	DLL_PUBLIC std::string inMatlabFormat(
		const std::size_t decimal_digits = 6) const;

	/** Saves the vector/matrix to a file compatible with MATLAB/Octave
	 * text format.
	 * \param file The target filename.
	 * \param fileFormat See TMatrixTextFileFormat. The format of the numbers in
	 * the text file.
	 * \param appendMRPTHeader Insert this header to the file "% File generated
	 * by MRPT. Load with MATLAB with: VAR=load(FILENAME);"
	 * \param userHeader Additional text to be written at the head of the file.
	 * Typically MALAB comments "% This file blah blah". Final end-of-line is
	 * not needed. \sa loadFromTextFile, CMatrixDynamic::inMatlabFormat,
	 * SAVE_MATRIX
	 */
	DLL_PUBLIC void saveToTextFile(
		const std::string& file,
		mrpt::math::TMatrixTextFileFormat fileFormat =
			mrpt::math::MATRIX_FORMAT_ENG,
		bool appendMRPTHeader = false,
		const std::string& userHeader = std::string()) const;

	/** Loads a vector/matrix from a text file, compatible with MATLAB text
	 * format. Lines starting with '%' or '#' are interpreted as comments and
	 * ignored.
	 * \exception std::runtime_error On format error.
	 * \sa saveToTextFile, fromMatlabStringFormat
	 */
	DLL_PUBLIC void loadFromTextFile(std::istream& f);

	/// \overload
	DLL_PUBLIC void loadFromTextFile(const std::string& file);

	/** Removes columns of the matrix.
	 * This "unsafe" version assumes indices sorted in ascending order. */
	DLL_PUBLIC void unsafeRemoveColumns(const std::vector<std::size_t>& idxs);

	/** Removes columns of the matrix. Indices may be unsorted and duplicated */
	DLL_PUBLIC void removeColumns(const std::vector<std::size_t>& idxsToRemove);

	/** Removes rows of the matrix.
	 * This "unsafe" version assumes indices sorted in ascending order. */
	DLL_PUBLIC void unsafeRemoveRows(const std::vector<std::size_t>& idxs);

	/** Removes rows of the matrix. Indices may be unsorted and duplicated */
	DLL_PUBLIC void removeRows(const std::vector<std::size_t>& idxsToRemove);

	/** Copies the given input submatrix/vector into this matrix/vector,
	 * starting at the given top-left coordinates. */
	template <typename OTHERMATVEC>
	void insertMatrix(
		const int row_start, const int col_start, const OTHERMATVEC& submat)
	{
		ASSERT_BELOW_(row_start + submat.rows(), mvbDerived().rows());
		ASSERT_BELOW_(col_start + submat.cols(), mvbDerived().cols());
		for (int r = 0; r < submat.rows(); r++)
			for (int c = 0; c < submat.cols(); c++)
				mvbDerived()(r + row_start, c + col_start) = submat(r, c);
	}

	template <typename OTHERMATVEC>
	bool operator==(const OTHERMATVEC& o) const
	{
		const auto& d = mvbDerived();
		if (d.cols() != o.cols() || d.rows() != o.rows()) return false;
		for (typename OTHERMATVEC::Index r = 0; r < d.rows(); r++)
			for (typename OTHERMATVEC::Index c = 0; c < d.cols(); c++)
				if (d(r, c) != o(r, c)) return false;
		return true;
	}
	template <typename OTHERMATVEC>
	bool operator!=(const OTHERMATVEC& o) const
	{
		return !(*this == o);
	}

	/** @} */

   protected:
	/** Issues a static_assert() error if trying to compile a method that
	 * requires Eigen headers, without including them. */
	template <typename DER>
	void internalAssertEigenDefined() const
	{
		if constexpr (!mrpt::is_defined_v<typename DER::eigen_t>)
		{
			static_assert(
				mrpt::is_defined_v<typename DER::eigen_t>,
				"Using this method requires including `<Eigen/Dense>` in the "
				"calling C++ file");
		}
	}
};

template <typename Scalar, class Derived>
std::ostream& operator<<(
	std::ostream& o, const MatrixVectorBase<Scalar, Derived>& m)
{
	return o << m.asStr();
}

}  // namespace mrpt::math
