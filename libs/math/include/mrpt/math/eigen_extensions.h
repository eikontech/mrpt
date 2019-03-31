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
#include <Eigen/Eigenvalues>  // EigenSolver
#include <algorithm>
#include <cstdio>  // fopen(),...
#include <ctime>  // time(),...
#include <fstream>  // ifstream
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

// Generic version for all kind of matrices:
template <int R, int C>
struct MatOrVecResizer
{
	template <typename MAT>
	static inline void doit(MAT& mat, size_t new_rows, size_t new_cols)
	{
		::mrpt::math::detail::TAuxResizer<MAT, MAT::SizeAtCompileTime>::
			internal_resize(mat, new_rows, new_cols);
	}
};
// Specialization for column matrices:
template <int R>
struct MatOrVecResizer<R, 1>
{
	template <typename MAT>
	static inline void doit(MAT& mat, size_t new_rows, size_t)
	{
		::mrpt::math::detail::TAuxResizer<
			MAT, MAT::SizeAtCompileTime>::internal_resize(mat, new_rows);
	}
};
// Specialization for row matrices:
template <int C>
struct MatOrVecResizer<1, C>
{
	template <typename MAT>
	static inline void doit(MAT& mat, size_t, size_t new_cols)
	{
		::mrpt::math::detail::TAuxResizer<
			MAT, MAT::SizeAtCompileTime>::internal_resize(mat, new_cols);
	}
};
template <>
struct MatOrVecResizer<1, 1>
{
	template <typename MAT>
	static inline void doit(MAT& mat, size_t, size_t new_cols)
	{
		::mrpt::math::detail::TAuxResizer<
			MAT, MAT::SizeAtCompileTime>::internal_resize(mat, new_cols);
	}
};
}  // namespace detail

template <class MATRIX>
bool fromMatlabStringFormat(
	MATRIX& m, const std::string& s,
	mrpt::optional_ref<std::ostream> dump_errors_here = std::nullopt)
{
	using Scalar = typename MATRIX::Scalar;
	// Start with a (0,0) matrix:
	if (MATRIX::RowsAtCompileTime == -1 /*Eigen::Dynamic*/) m = MATRIX();

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
				if (MATRIX::RowsAtCompileTime == -1 /*Eigen::Dynamic*/)
					m = MATRIX();
			}
		}
		else
		{
			const size_t N = lstElements.size();

			// Check valid width: All rows must have the same width
			if ((nRow > 0 && size_t(m.cols()) != N) ||
				(nRow == 0 &&
				 MATRIX::ColsAtCompileTime != -1 /*Eigen::Dynamic*/ &&
				 MATRIX::ColsAtCompileTime != int(N)))
			{
				if (dump_errors_here)
					dump_errors_here->get()
						<< "[fromMatlabStringFormat] Row " << nRow + 1
						<< " has invalid number of columns.\n";
				return false;
			}

			// Append to the matrix:
			if (MATRIX::RowsAtCompileTime == -1 /*Eigen::Dynamic*/ ||
				MATRIX::ColsAtCompileTime == -1 /*Eigen::Dynamic*/)
				detail::MatOrVecResizer<
					MATRIX::RowsAtCompileTime,
					MATRIX::ColsAtCompileTime>::doit(m, nRow + 1, N);
			else if (
				MATRIX::RowsAtCompileTime != -1 /*Eigen::Dynamic*/ &&
				int(nRow) >= MATRIX::RowsAtCompileTime)
			{
				if (dump_errors_here)
					dump_errors_here->get()
						<< "[fromMatlabStringFormat] Read more "
						   "rows than the capacity of the "
						   "fixed sized matrix.\n";
				return false;
			}
			for (size_t q = 0; q < N; q++) m(nRow, q) = lstElements[q];
			// Go for the next row:
			nRow++;
		}
		i = end_row + 1;
	}
	// For fixed sized matrices, check size:
	if (MATRIX::RowsAtCompileTime != -1 /*Eigen::Dynamic*/ &&
		int(nRow) != MATRIX::RowsAtCompileTime)
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
	const Eigen::MatrixBase<Derived>& m, const size_t decimal_digits = 6)
{
	using Index = typename Derived::Index;
	std::stringstream s;
	s << "[" << std::scientific;
	s.precision(decimal_digits);
	for (Index i = 0; i < m.rows(); i++)
	{
		for (Index j = 0; j < m.cols(); j++) s << m.coeff(i, j) << " ";
		if (i < m.rows() - 1) s << ";";
	}
	s << "]";
	return s.str();
}

/** Save matrix to a text file, compatible with MATLAB text format (see also the
 * methods of matrix classes themselves).
 * \param theMatrix It can be a CMatrixDynamic or a CMatrixFixed.
 * \param file The target filename.
 * \param fileFormat See TMatrixTextFileFormat. The format of the numbers in
 * the text file.
 * \param appendMRPTHeader Insert this header to the file "% File generated
 * by MRPT. Load with MATLAB with: VAR=load(FILENAME);"
 * \param userHeader Additional text to be written at the head of the file.
 * Typically MALAB comments "% This file blah blah". Final end-of-line is not
 * needed.
 * \sa loadFromTextFile, CMatrixDynamic::inMatlabFormat, SAVE_MATRIX
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
	using Scalar = typename Derived::Scalar;
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
				Index(i) != m.cols())
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
					detail::MatOrVecResizer<
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

			for (size_t q = 0; q < i; q++) m(nRows, q) = Scalar(fil[q]);

			nRows++;
		}  // end if fgets
	}  // end while not feof

	// Final resize to the real size (in case we allocated space in advance):
	if (Derived::RowsAtCompileTime == Eigen::Dynamic ||
		Derived::ColsAtCompileTime == Eigen::Dynamic)
		detail::MatOrVecResizer<
			Derived::RowsAtCompileTime,
			Derived::ColsAtCompileTime>::doit(m, nRows, m.cols());

	// Report error as exception
	if (!nRows)
		throw std::runtime_error(
			"loadFromTextFile: Error loading from text file");
}

/// \overload
template <class MATRIX>
void loadFromTextFile(MATRIX& m, const std::string& file)
{
	std::ifstream f(file.c_str());
	if (f.fail())
		throw std::runtime_error(
			std::string("loadFromTextFile: can't open file:") + file);
	loadFromTextFile(m, f);
}

/** Remove columns of the matrix. The unsafe version assumes that, the indices
 * are sorted in ascending order. */
template <class MATRIX>
void unsafeRemoveColumns(MATRIX& m, const std::vector<std::size_t>& idxs)
{
	std::size_t k = 1;
	const auto nR = m.rows();
	for (auto it = idxs.rbegin(); it != idxs.rend(); ++it, ++k)
	{
		const auto nC = m.cols() - *it - k;
		if (nC > 0)
			m.asEigen().block(0, *it, nR, nC) =
				m.asEigen().block(0, *it + 1, nR, nC).eval();
	}
	m.setSize(nR, m.cols() - idxs.size());
}

/** Remove columns of the matrix. indices may be unsorted and duplicated. */
template <class MATRIX>
void removeColumns(MATRIX& m, const std::vector<std::size_t>& idxsToRemove)
{
	std::vector<std::size_t> idxs = idxsToRemove;
	std::sort(idxs.begin(), idxs.end());
	auto itEnd = std::unique(idxs.begin(), idxs.end());
	idxs.resize(itEnd - idxs.begin());
	unsafeRemoveColumns(m, idxs);
}

/** Remove rows of the matrix. The unsafe version assumes that, the indices are
 * sorted in ascending order. */
template <class MATRIX>
void unsafeRemoveRows(MATRIX& m, const std::vector<size_t>& idxs)
{
	std::size_t k = 1;
	const auto nC = m.cols();
	for (auto it = idxs.rbegin(); it != idxs.rend(); ++it, ++k)
	{
		const auto nR = m.rows() - *it - k;
		if (nR > 0)
			m.asEigen().block(*it, 0, nR, nC) =
				m.asEigen().block(*it + 1, 0, nR, nC).eval();
	}
	m.setSize(m.rows() - idxs.size(), nC);
}

/** Remove rows of the matrix. */
template <class MATRIX>
void removeRows(MATRIX& m, const std::vector<size_t>& idxsToRemove)
{
	std::vector<std::size_t> idxs = idxsToRemove;
	std::sort(idxs.begin(), idxs.end());
	auto itEnd = std::unique(idxs.begin(), idxs.end());
	idxs.resize(itEnd - idxs.begin());
	unsafeRemoveRows(m, idxs);
}

}  // namespace mrpt::math
