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
#include <mrpt/math/MatrixVectorBase.h>
#include <Eigen/Eigenvalues>  // EigenSolver
#include <cstdint>
#include <cstdio>  // fopen(),...
#include <ctime>  // time(),...
#include <fstream>  // ifstream
#include <sstream>  // stringstream
#include <stdexcept>
#include <vector>

namespace mrpt::math
{
template <typename Scalar, class Derived>
bool MatrixVectorBase<Scalar, Derived>::fromMatlabStringFormat(
	const std::string& s, mrpt::optional_ref<std::ostream> dump_errors_here)
{
	// Start with a (0,0) matrix:
	if (Derived::RowsAtCompileTime == Eigen::Dynamic) mvbDerived() = Derived();

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
				if (Derived::RowsAtCompileTime == Eigen::Dynamic)
					mvbDerived() = Derived();
			}
		}
		else
		{
			const size_t N = lstElements.size();

			// Check valid width: All rows must have the same width
			if ((nRow > 0 && size_t(mvbDerived().cols()) != N) ||
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
				mvbDerived().resize(nRow + 1, N);
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
			for (size_t q = 0; q < N; q++)
				mvbDerived()(nRow, q) = lstElements[q];
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

template <typename Scalar, class Derived>
std::string MatrixVectorBase<Scalar, Derived>::inMatlabFormat(
	const std::size_t decimal_digits) const
{
	using Index = typename Derived::Index;
	std::stringstream s;
	s << "[" << std::scientific;
	s.precision(decimal_digits);
	for (Index i = 0; i < mvbDerived().rows(); i++)
	{
		for (Index j = 0; j < mvbDerived().cols(); j++)
			s << mvbDerived().coeff(i, j) << " ";
		if (i < mvbDerived().rows() - 1) s << ";";
	}
	s << "]";
	return s.str();
}

template <typename Scalar, class Derived>
void MatrixVectorBase<Scalar, Derived>::saveToTextFile(
	const std::string& file, mrpt::math::TMatrixTextFileFormat fileFormat,
	bool appendMRPTHeader, const std::string& userHeader) const
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

	for (Index i = 0; i < mvbDerived().rows(); i++)
	{
		for (Index j = 0; j < mvbDerived().cols(); j++)
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
			if (j < (mvbDerived().cols() - 1)) ::fprintf(f, " ");
		}
		::fprintf(f, "\n");
	}
	::fclose(f);
}

template <typename Scalar, class Derived>
void MatrixVectorBase<Scalar, Derived>::loadFromTextFile(std::istream& f)
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
				Index(i) != mvbDerived().cols())
				throw std::runtime_error(
					"loadFromTextFile: The matrix in the text file does not "
					"have the same number of columns in all rows");

			// Append to the matrix:
			if (Derived::RowsAtCompileTime == Eigen::Dynamic ||
				Derived::ColsAtCompileTime == Eigen::Dynamic)
			{
				if (mvbDerived().rows() < static_cast<int>(nRows + 1) ||
					mvbDerived().cols() < static_cast<int>(i))
				{
					const size_t extra_rows =
						std::max(static_cast<size_t>(1), nRows >> 1);
					mvbDerived().resize(nRows + extra_rows, i);
				}
			}
			else if (
				Derived::RowsAtCompileTime != Eigen::Dynamic &&
				int(nRows) >= Derived::RowsAtCompileTime)
				throw std::runtime_error(
					"loadFromTextFile: Read more rows than the capacity of the "
					"fixed sized matrix.");

			for (size_t q = 0; q < i; q++)
				mvbDerived()(nRows, q) = Scalar(fil[q]);

			nRows++;
		}  // end if fgets
	}  // end while not feof

	// Final resize to the real size (in case we allocated space in advance):
	if (Derived::RowsAtCompileTime == Eigen::Dynamic ||
		Derived::ColsAtCompileTime == Eigen::Dynamic)
		mvbDerived().resize(nRows, mvbDerived().cols());

	// Report error as exception
	if (!nRows)
		throw std::runtime_error(
			"loadFromTextFile: Error loading from text file");
}

template <typename Scalar, class Derived>
void MatrixVectorBase<Scalar, Derived>::loadFromTextFile(
	const std::string& file)
{
	std::ifstream f(file.c_str());
	if (f.fail())
		throw std::runtime_error(
			std::string("loadFromTextFile: can't open file:") + file);
	loadFromTextFile(f);
}

template <typename Scalar, class Derived>
void MatrixVectorBase<Scalar, Derived>::unsafeRemoveColumns(
	const std::vector<std::size_t>& idxs)
{
	std::size_t k = 1;
	const auto nR = mvbDerived().rows();
	for (auto it = idxs.rbegin(); it != idxs.rend(); ++it, ++k)
	{
		const auto nC = mvbDerived().cols() - *it - k;
		if (nC > 0)
			mvbDerived().asEigen().block(0, *it, nR, nC) =
				mvbDerived().asEigen().block(0, *it + 1, nR, nC).eval();
	}
	mvbDerived().setSize(nR, mvbDerived().cols() - idxs.size());
}

template <typename Scalar, class Derived>
void MatrixVectorBase<Scalar, Derived>::removeColumns(
	const std::vector<std::size_t>& idxsToRemove)
{
	std::vector<std::size_t> idxs = idxsToRemove;
	std::sort(idxs.begin(), idxs.end());
	auto itEnd = std::unique(idxs.begin(), idxs.end());
	idxs.resize(itEnd - idxs.begin());
	unsafeRemoveColumns(idxs);
}

template <typename Scalar, class Derived>
void MatrixVectorBase<Scalar, Derived>::unsafeRemoveRows(
	const std::vector<size_t>& idxs)
{
	std::size_t k = 1;
	const auto nC = mvbDerived().cols();
	for (auto it = idxs.rbegin(); it != idxs.rend(); ++it, ++k)
	{
		const auto nR = mvbDerived().rows() - *it - k;
		if (nR > 0)
			mvbDerived().asEigen().block(*it, 0, nR, nC) =
				mvbDerived().asEigen().block(*it + 1, 0, nR, nC).eval();
	}
	mvbDerived().setSize(mvbDerived().rows() - idxs.size(), nC);
}

template <typename Scalar, class Derived>
void MatrixVectorBase<Scalar, Derived>::removeRows(
	const std::vector<size_t>& idxsToRemove)
{
	std::vector<std::size_t> idxs = idxsToRemove;
	std::sort(idxs.begin(), idxs.end());
	auto itEnd = std::unique(idxs.begin(), idxs.end());
	idxs.resize(itEnd - idxs.begin());
	unsafeRemoveRows(idxs);
}

template <typename Scalar, class Derived>
std::string MatrixVectorBase<Scalar, Derived>::asStr() const
{
	std::stringstream ss;
	ss << mvbDerived();
	return ss.str();
}

template <typename Scalar, class Derived>
Scalar MatrixVectorBase<Scalar, Derived>::det() const
{
	return mvbDerived().asEigen().det();
}

template <typename Scalar, class Derived>
Scalar MatrixVectorBase<Scalar, Derived>::sum() const
{
	return mvbDerived().asEigen().array().det();
}

template <typename Scalar, class Derived>
Scalar MatrixVectorBase<Scalar, Derived>::minCoeff() const
{
	return mvbDerived().asEigen().minCoeff();
}

template <typename Scalar, class Derived>
Scalar MatrixVectorBase<Scalar, Derived>::maxCoeff() const
{
	return mvbDerived().asEigen().maxCoeff();
}

template <typename Scalar, class Derived>
void MatrixVectorBase<Scalar, Derived>::operator+=(Scalar s)
{
	mvbDerived().asEigen().array() += s;
}

template <typename Scalar, class Derived>
void MatrixVectorBase<Scalar, Derived>::operator-=(Scalar s)
{
	mvbDerived().asEigen().array() -= s;
}

template <typename Scalar, class Derived>
void MatrixVectorBase<Scalar, Derived>::operator*=(Scalar s)
{
	mvbDerived().asEigen().array() *= s;
}

namespace detail
{
// Aux func to sort by ascending eigenvalues:
template <typename Scalar, typename VEC1, typename MATRIX1, typename MATRIX2>
void sortEigResults(
	const VEC1& eVals, const MATRIX1& eVecs, std::vector<Scalar>& sorted_eVals,
	MATRIX2& sorted_eVecs)
{
	const int64_t N = static_cast<int64_t>(eVals.size());
	std::vector<std::pair<Scalar, int64_t>> D;
	D.reserve(N);
	for (int64_t i = 0; i < N; i++) D.emplace_back(eVals[i], i);
	std::sort(D.begin(), D.end());

	// store:
	sorted_eVecs.resize(eVecs.rows(), eVecs.cols());
	sorted_eVals.resize(N);
	for (int64_t i = 0; i < N; i++)
	{
		sorted_eVals[i] = D[i].first;
		sorted_eVecs.col(i) = eVecs.col(D[i].second);
	}
}
}  // namespace detail

template <typename Scalar, class Derived>
bool MatrixVectorBase<Scalar, Derived>::eig(
	Derived& eVecs, std::vector<Scalar>& eVals, bool sorted) const
{
	Eigen::EigenSolver<typename Derived::eigen_t> es(mvbDerived());
	if (es.info() != Eigen::Success) return false;
	const auto eigenVal = es.eigenvalues();
	ASSERT_EQUAL_(eigenVal.rows(), mvbDerived().rows());
	const auto N = eigenVal.rows();

	if (sorted)
	{
		detail::sortEigResults(eigenVal, es.eigenvectors(), eVals, eVecs);
	}
	else
	{
		eVals.resize(N);
		eVecs = es.eigenvectors();
		for (int i = 0; i < N; i++) eVals[i] = eigenVal[i];
	}
	return true;
}

template <typename Scalar, class Derived>
bool MatrixVectorBase<Scalar, Derived>::eig_symmetric(
	Derived& eVecs, std::vector<Scalar>& eVals, bool sorted) const
{
	Eigen::SelfAdjointEigenSolver<typename Derived::eigen_t> es(mvbDerived());
	if (es.info() != Eigen::Success) return false;
	const auto eigenVal = es.eigenvalues();
	ASSERT_EQUAL_(eigenVal.rows(), mvbDerived().rows());
	const auto N = eigenVal.rows();

	if (sorted)
	{
		detail::sortEigResults(eigenVal, es.eigenvectors(), eVals, eVecs);
	}
	else
	{
		eVals.resize(N);
		eVecs = es.eigenvectors();
		for (int i = 0; i < N; i++) eVals[i] = eigenVal[i];
	}
	return true;
}

template <typename Scalar, class Derived>
int MatrixVectorBase<Scalar, Derived>::rank(Scalar threshold) const
{
	Eigen::FullPivLU<typename Derived::eigen_t> lu(mvbDerived());
	if (threshold > 0) lu.setThreshold(threshold);
	return lu.rank();
}

template <typename Scalar, class Derived>
bool MatrixVectorBase<Scalar, Derived>::chol(Derived& U) const
{
	Eigen::LLT<typename Derived::PlainObject> Chol =
		mvbDerived().template selfadjointView<Eigen::Lower>().llt();
	if (Chol.info() == Eigen::NoConvergence) return false;
	U = Derived(Chol.matrixU());
	return true;
}

template <typename Scalar, class Derived>
void MatrixVectorBase<Scalar, Derived>::matProductOf(
	const Derived& A, const Derived& B)
{
	*this = (A.asEigen() * B.asEigen()).eval();
}

template <typename Scalar, class Derived>
Derived MatrixVectorBase<Scalar, Derived>::inverse() const
{
	ASSERT_EQUAL_(mvbDerived().cols(), mvbDerived().rows());
	const auto N = mvbDerived().cols();
	const auto I = Derived::eigen_t::Identity(N, N);
	Derived inv(mrpt::math::UNINITIALIZED_MATRIX);
	inv.asEigen() = mvbDerived().asEigen().llt().solve(I).eval();
	return inv;
}

template <typename Scalar, class Derived>
Derived MatrixVectorBase<Scalar, Derived>::inverse_LLt() const
{
	ASSERT_EQUAL_(mvbDerived().cols(), mvbDerived().rows());
	const auto N = mvbDerived().cols();
	const auto I = Derived::eigen_t::Identity(N, N);
	Derived inv(mrpt::math::UNINITIALIZED_MATRIX);
	inv.asEigen() = mvbDerived().asEigen().llt().solve(I).eval();
	return inv;
}

}  // namespace mrpt::math
