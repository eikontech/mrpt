/* +------------------------------------------------------------------------+
   |                     Mobile Robot Programming Toolkit (MRPT)            |
   |                          https://www.mrpt.org/                         |
   |                                                                        |
   | Copyright (c) 2005-2019, Individual contributors, see AUTHORS file     |
   | See: https://www.mrpt.org/Authors - All rights reserved.               |
   | Released under BSD License. See: https://www.mrpt.org/License          |
   +------------------------------------------------------------------------+ */

#pragma once

#include <mrpt/core/optional_ref.h>
#include <mrpt/math/math_frwds.h>
#include <algorithm>
#include <cstdio>  // fopen(),...
#include <sstream>  // stringstream
#include <stdexcept>
#include <vector>

/* Utility functions to ease the use of Eigen matrices */

namespace mrpt::math
{
/*! Selection of the number format in mrpt::math::saveToTextFile() */
enum TMatrixTextFileFormat
{
	/** engineering format '%e' */
	MATRIX_FORMAT_ENG = 0,
	/** fixed floating point '%f' */
	MATRIX_FORMAT_FIXED = 1,
	/** intergers '%i' */
	MATRIX_FORMAT_INT = 2
};

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

// Generic version for all kind of matrices:
template <int R, int C>
struct MatOrVecResizer
{
	template <typename S, int Opt, int MaxR, int MaxC>
	static inline void doit(
		Eigen::Matrix<S, R, C, Opt, MaxR, MaxC>& mat, size_t new_rows,
		size_t new_cols)
	{
		::mrpt::math::detail::TAuxResizer<
			Eigen::Matrix<S, R, C, Opt, MaxR, MaxC>,
			Eigen::Matrix<S, R, C, Opt, MaxR, MaxC>::SizeAtCompileTime>::
			internal_resize(mat, new_rows, new_cols);
	}
};
// Specialization for column matrices:
template <int R>
struct MatOrVecResizer<R, 1>
{
	template <typename S, int Opt, int MaxR, int MaxC>
	static inline void doit(
		Eigen::Matrix<S, R, 1, Opt, MaxR, MaxC>& mat, size_t new_rows, size_t)
	{
		::mrpt::math::detail::TAuxResizer<
			Eigen::Matrix<S, R, 1, Opt, MaxR, MaxC>,
			Eigen::Matrix<S, R, 1, Opt, MaxR, MaxC>::SizeAtCompileTime>::
			internal_resize(mat, new_rows);
	}
};
// Specialization for row matrices:
template <int C>
struct MatOrVecResizer<1, C>
{
	template <typename S, int Opt, int MaxR, int MaxC>
	static inline void doit(
		Eigen::Matrix<S, 1, C, Opt, MaxR, MaxC>& mat, size_t, size_t new_cols)
	{
		::mrpt::math::detail::TAuxResizer<
			Eigen::Matrix<S, 1, C, Opt, MaxR, MaxC>,
			Eigen::Matrix<S, 1, C, Opt, MaxR, MaxC>::SizeAtCompileTime>::
			internal_resize(mat, new_cols);
	}
};
template <>
struct MatOrVecResizer<1, 1>
{
	template <typename S, int Opt, int MaxR, int MaxC>
	static inline void doit(
		Eigen::Matrix<S, 1, 1, Opt, MaxR, MaxC>& mat, size_t, size_t new_cols)
	{
		::mrpt::math::detail::TAuxResizer<
			Eigen::Matrix<S, 1, 1, Opt, MaxR, MaxC>,
			Eigen::Matrix<S, 1, 1, Opt, MaxR, MaxC>::SizeAtCompileTime>::
			internal_resize(mat, new_cols);
	}
};

/** Compute the eigenvectors and eigenvalues, both returned as matrices:
 * eigenvectors are the columns, and eigenvalues
 */
template <class Derived, class MATRIX1, class VECTOR1>
bool eigenVectorsVec(
	const Eigen::MatrixBase<Derived>& m, MATRIX1& eVecs, VECTOR1& eVals)
{
	Eigen::EigenSolver<Derived> es(m, true);
	// Keep only the real part of complex matrix
	if (es.info() != Eigen::Success) return false;
	eVecs = es.eigenvectors().real();
	eVals = es.eigenvalues().real();

	// Sort by ascending eigenvalues:
	std::vector<std::pair<Scalar, Index>> D;
	D.reserve(eVals.size());
	for (Index i = 0; i < eVals.size(); i++)
		D.push_back(std::pair<Scalar, Index>(eVals.coeff(i, 0), i));
	std::sort(D.begin(), D.end());
	MATRIX1 sortedEigs;
	sortedEigs.resizeLike(eVecs);
	for (std::size_t i = 0; i < eVals.size(); i++)
	{
		eVals.coeffRef(i, 0) = D[i].first;
		sortedEigs.col(i) = eVecs.col(D[i].second);
	}
	eVecs = std::move(sortedEigs);
	return true;
}

/** Compute the eigenvectors and eigenvalues, both returned as matrices:
 * eigenvectors are the columns, and eigenvalues. \return false on error.
 */
template <class Derived, class MATRIX1, class MATRIX2>
bool eigenVectors(
	const Eigen::MatrixBase<Derived>& m, MATRIX1& eVecs, MATRIX2& eVals)
{
	Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> evals;
	if (!eigenVectorsVec(m, eVecs, evals)) return false;
	eVals.resize(evals.size(), evals.size());
	eVals.setZero();
	eVals.diagonal() = evals;
	return true;
}

/** Compute the eigenvectors and eigenvalues, both returned as matrices:
 * eigenvectors are the columns, and eigenvalues
 */
template <class Derived, class MATRIX1, class VECTOR1>
void eigenVectorsSymmetricVec(
	const Eigen::MatrixBase<Derived>& m, MATRIX1& eVecs, VECTOR1& eVals) const
{
	// This solver returns the eigenvectors already sorted.
	Eigen::SelfAdjointEigenSolver<Derived> eigensolver(m);
	eVecs = eigensolver.eigenvectors();
	eVals = eigensolver.eigenvalues();
}
/** Compute the eigenvectors and eigenvalues, both returned as matrices:
 * eigenvectors are the columns, and eigenvalues
 */
template <class Derived, class MATRIX1, class MATRIX2>
void eigenVectorsSymmetric(
	const Eigen::MatrixBase<Derived>& m, MATRIX1& eVecs, MATRIX2& eVals)
{
	Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> evals;
	eigenVectorsSymmetricVec(m, eVecs, evals);
	eVals.resize(evals.size(), evals.size());
	eVals.setZero();
	eVals.diagonal() = evals;
}

template <class Derived>
bool fromMatlabStringFormat(
	Eigen::MatrixBase<Derived>& m, const std::string& s,
	mrpt::optional_ref<std::ostream> dump_errors_here = std::nullopt)
{
	using Scalar = typename Derived::Scalar;
	// Start with a (0,0) matrix:
	if (Derived::RowsAtCompileTime == Eigen::Dynamic) m = Derived();

	// Look for starting "[".
	size_t ini = s.find_first_not_of(" \t\r\n");
	if (ini == std::string::npos || s[ini] != '[')
	{
		return false;
	}

	size_t end = s.find_last_not_of(" \t\r\n");
	if (end == std::string::npos || s[end] != ']') return false;

	if (ini > end) return false;

	std::vector<Scalar> lstElements;

	size_t i = ini + 1;
	size_t nRow = 0;

	while (i < end)
	{
		// Extract one row:
		size_t end_row = s.find_first_of(";]", i);
		if (end_row == std::string::npos)
		{
			return false;
		}

		// We have one row in s[ i : (end_row-1) ]
		std::stringstream ss(s.substr(i, end_row - i));
		lstElements.clear();
		try
		{
			while (!ss.eof())
			{
				Scalar val;
				ss >> val;
				if (ss.bad() || ss.fail()) break;
				lstElements.push_back(val);
			}
		}
		catch (...)
		{
		}  // end of line

		// Empty row? Only for the first row, then this is an empty matrix:
		if (lstElements.empty())
		{
			if (nRow > 0)
				return false;
			else
			{
				// Else, this may be an empty matrix... if there is no next row,
				// we'll return with a (0,0) matrix
				if (Derived::RowsAtCompileTime == Eigen::Dynamic) m = Derived();
			}
		}
		else
		{
			const size_t N = lstElements.size();

			// Check valid width: All rows must have the same width
			if ((nRow > 0 && size_t(cols()) != N) ||
				(nRow == 0 && Derived::ColsAtCompileTime != Eigen::Dynamic &&
				 Derived::ColsAtCompileTime != int(N)))
			{
				if (dump_errors_here)
					dump_errors_here->get()
						<< "[fromMatlabStringFormat] Row " << nRow + 1
						<< " has invalid number of columns.\n";
				return false;
			}

			// Append to the matrix:
			if (Derived::RowsAtCompileTime == Eigen::Dynamic ||
				Derived::ColsAtCompileTime == Eigen::Dynamic)
				internal::MatOrVecResizer<
					Derived::RowsAtCompileTime,
					Derived::ColsAtCompileTime>::doit(m, nRow + 1, N);
			else if (
				Derived::RowsAtCompileTime != Eigen::Dynamic &&
				int(nRow) >= Derived::RowsAtCompileTime)
			{
				if (dump_errors_here)
					dump_errors_here->get()
						<< "[fromMatlabStringFormat] Read more "
						   "rows than the capacity of the "
						   "fixed sized matrix.\n";
				return false;
			}
			for (size_t q = 0; q < N; q++) m.coeffRef(nRow, q) = lstElements[q];
			// Go for the next row:
			nRow++;
		}
		i = end_row + 1;
	}
	// For fixed sized matrices, check size:
	if (Derived::RowsAtCompileTime != Eigen::Dynamic &&
		int(nRow) != Derived::RowsAtCompileTime)
	{
		if (dump_errors_here)
			dump_errors_here->get()
				<< "[fromMatlabStringFormat] Read less rows "
				   "than the capacity of the fixed sized "
				   "matrix.\n";
		return false;
	}
	return true;  // Ok
}

template <class Derived>
std::string inMatlabFormat(
	Eigen::MatrixBase<Derived>& m, const size_t decimal_digits = 6)
{
	using Index = typename Derived::Index;
	std::stringstream s;
	s << "[" << std::scientific;
	s.precision(decimal_digits);
	for (Index i = 0; i < m.rows(); i++)
	{
		for (Index j = 0; j < m.cols(); j++) s << m.coeff(i, j) << " ";
		if (i < n.rows() - 1) s << ";";
	}
	s << "]";
	return s.str();
}

/** Save matrix to a text file, compatible with MATLAB text format (see also the
 * methods of matrix classes themselves).
 * \param theMatrix It can be a CMatrixTemplate or a CMatrixFixedNumeric.
 * \param file The target filename.
 * \param fileFormat See TMatrixTextFileFormat. The format of the numbers in
 * the text file.
 * \param appendMRPTHeader Insert this header to the file "% File generated
 * by MRPT. Load with MATLAB with: VAR=load(FILENAME);"
 * \param userHeader Additional text to be written at the head of the file.
 * Typically MALAB comments "% This file blah blah". Final end-of-line is not
 * needed.
 * \sa loadFromTextFile, CMatrixTemplate::inMatlabFormat, SAVE_MATRIX
 */
template <class Derived>
void saveToTextFile(
	const Eigen::MatrixBase<Derived>& m, const std::string& file,
	mrpt::math::TMatrixTextFileFormat fileFormat =
		mrpt::math::MATRIX_FORMAT_ENG,
	bool appendMRPTHeader = false,
	const std::string& userHeader = std::string())
{
	using Index = typename Derived::Index;
	// Use a secure version in Visual Studio 2005+
#if defined(_MSC_VER) && (_MSC_VER >= 1400)
	FILE* f;
	if (0 != ::fopen_s(&f, file.c_str(), "wt")) f = nullptr;
#else
	FILE* f = ::fopen(file.c_str(), "wt");
#endif
	if (!f)
		throw std::runtime_error(
			std::string("saveToTextFile: Error opening file ") + file +
			std::string("' for writing a matrix as text."));

	if (!userHeader.empty()) fprintf(f, "%s", userHeader.c_str());

	if (appendMRPTHeader)
	{
		time_t rawtime;
		::time(&rawtime);
		// Use a secure version in Visual Studio 2005+
#if defined(_MSC_VER) && (_MSC_VER >= 1400)
		struct tm timeinfo_data;
		struct tm* timeinfo;
		if (0 != ::localtime_s(&timeinfo_data, &rawtime))
			timeinfo = nullptr;
		else
			timeinfo = &timeinfo_data;
#else
		struct tm* timeinfo = ::localtime(&rawtime);
#endif

#if defined(_MSC_VER) && (_MSC_VER >= 1400)
		char strTimeBuf[100];
		if (0 != asctime_s(strTimeBuf, sizeof(strTimeBuf), timeinfo))
			strTimeBuf[0] = '\0';
		char* strTime = &strTimeBuf[0];
#else
		char* strTime = asctime(timeinfo);
#endif
		fprintf(
			f,
			"%% File generated with mrpt-math at %s\n"
			"%%------------------------------------\n",
			strTime);
	}

	for (Index i = 0; i < m.rows(); i++)
	{
		for (Index j = 0; j < m.cols(); j++)
		{
			switch (fileFormat)
			{
				case mrpt::math::MATRIX_FORMAT_ENG:
					::fprintf(f, "%.16e", static_cast<double>(m(i, j)));
					break;
				case mrpt::math::MATRIX_FORMAT_FIXED:
					::fprintf(f, "%.16f", static_cast<double>(m(i, j)));
					break;
				case mrpt::math::MATRIX_FORMAT_INT:
					::fprintf(f, "%i", static_cast<int>(m(i, j)));
					break;
				default:
					throw std::runtime_error(
						"Unsupported value for the parameter 'fileFormat'!");
			};
			// Separating blank space
			if (j < (m.cols() - 1)) ::fprintf(f, " ");
		}
		::fprintf(f, "\n");
	}
	::fclose(f);
}

/** Load matrix from a text file, compatible with MATLAB text format.
 *  Lines starting with '%' or '#' are interpreted as comments and ignored.
 * \sa saveToTextFile, fromMatlabStringFormat
 */
template <class Derived>
void loadFromTextFile(Eigen::MatrixBase<Derived>& m, std::istream& f)
{
	using Index = typename Derived::Index;
	std::string str;
	std::vector<double> fil(512);
	size_t nRows = 0;
	while (!f.eof() && !f.fail())
	{
		std::getline(f, str);
		if (str.size() && str[0] != '#' && str[0] != '%')
		{
			// Parse row to floats:
			const char* ptr = str.c_str();
			char* ptrEnd = nullptr;
			size_t i = 0;
			// Process each number in this row:
			while (ptr[0] && ptr != ptrEnd)
			{
				// Find next number: (non white-space character):
				while (ptr[0] &&
					   (ptr[0] == ' ' || ptr[0] == ',' || ptr[0] == '\t' ||
						ptr[0] == '\r' || ptr[0] == '\n'))
					ptr++;
				if (fil.size() <= i) fil.resize(fil.size() + (fil.size() >> 1));
				// Convert to "double":
				fil[i] = strtod(ptr, &ptrEnd);
				// A valid conversion has been done?
				if (ptr != ptrEnd)
				{
					i++;  // Yes
					ptr = ptrEnd;
					ptrEnd = nullptr;
				}
			};  // end while procesing this row

			// "i": # of columns:
			if ((Derived::ColsAtCompileTime != Eigen::Dynamic &&
				 Index(i) != Derived::ColsAtCompileTime))
				throw std::runtime_error(
					"loadFromTextFile: The matrix in the text file does not "
					"match fixed matrix size");
			if (Derived::ColsAtCompileTime == Eigen::Dynamic && nRows > 0 &&
				Index(i) != cols())
				throw std::runtime_error(
					"loadFromTextFile: The matrix in the text file does not "
					"have the same number of columns in all rows");

			// Append to the matrix:
			if (Derived::RowsAtCompileTime == Eigen::Dynamic ||
				Derived::ColsAtCompileTime == Eigen::Dynamic)
			{
				if (m.rows() < static_cast<int>(nRows + 1) ||
					m.cols() < static_cast<int>(i))
				{
					const size_t extra_rows =
						std::max(static_cast<size_t>(1), nRows >> 1);
					internal::MatOrVecResizer<
						Derived::RowsAtCompileTime,
						Derived::ColsAtCompileTime>::
						doit(m, nRows + extra_rows, i);
				}
			}
			else if (
				Derived::RowsAtCompileTime != Eigen::Dynamic &&
				int(nRows) >= Derived::RowsAtCompileTime)
				throw std::runtime_error(
					"loadFromTextFile: Read more rows than the capacity of the "
					"fixed sized matrix.");

			for (size_t q = 0; q < i; q++) coeffRef(nRows, q) = Scalar(fil[q]);

			nRows++;
		}  // end if fgets
	}  // end while not feof

	// Final resize to the real size (in case we allocated space in advance):
	if (Derived::RowsAtCompileTime == Eigen::Dynamic ||
		Derived::ColsAtCompileTime == Eigen::Dynamic)
		internal::MatOrVecResizer<
			Derived::RowsAtCompileTime,
			Derived::ColsAtCompileTime>::doit(m, nRows, cols());

	// Report error as exception
	if (!nRows)
		throw std::runtime_error(
			"loadFromTextFile: Error loading from text file");
}

/// \overload
template <class Derived>
void loadFromTextFile(Eigen::MatrixBase<Derived>& m, const std::string& file)
{
	std::ifstream f(file.c_str());
	if (f.fail())
		throw std::runtime_error(
			std::string("loadFromTextFile: can't open file:") + file);
	loadFromTextFile(m, f);
}

}  // namespace mrpt::math
