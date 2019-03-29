/* +------------------------------------------------------------------------+
   |                     Mobile Robot Programming Toolkit (MRPT)            |
   |                          https://www.mrpt.org/                         |
   |                                                                        |
   | Copyright (c) 2005-2019, Individual contributors, see AUTHORS file     |
   | See: https://www.mrpt.org/Authors - All rights reserved.               |
   | Released under BSD License. See: https://www.mrpt.org/License          |
   +------------------------------------------------------------------------+ */
#pragma once

// Eigen forward declarations:
#include <mrpt/config.h>
#include <cstddef>  // size_t

// Minimum Eigen fwrd-decls:
namespace Eigen
{
template <typename Derived>
struct EigenBase;
template <typename Derived>
class MatrixBase;
}  // namespace Eigen

namespace mrpt
{
namespace math
{
/** ContainerType<T>::element_t exposes the value of any STL or Eigen container
 */
template <typename CONTAINER>
struct ContainerType;
/** Specialization for Eigen containers */
template <typename Derived>
struct ContainerType<Eigen::EigenBase<Derived>>
{
	using element_t = typename Derived::Scalar;
};

// Dynamic size:
// template <class T>
// class CMatrixDynamic;
// using CMatrixFloat = CMatrixDynamic<float>;
// using CMatrixDouble = CMatrixDynamic<double>;
// template <typename T>
// class dynamic_vector;
// using CVectorFloat = dynamic_vector<float>;
// using CVectorDouble = dynamic_vector<double>;

// Fixed size:
template <typename T, size_t NROWS, size_t NCOLS>
class CMatrixFixed;

/** @name Typedefs for common sizes
	@{ */
using CMatrixDouble22 = CMatrixFixed<double, 2, 2>;
using CMatrixDouble23 = CMatrixFixed<double, 2, 3>;
using CMatrixDouble32 = CMatrixFixed<double, 3, 2>;
using CMatrixDouble33 = CMatrixFixed<double, 3, 3>;
using CMatrixDouble44 = CMatrixFixed<double, 4, 4>;
using CMatrixDouble66 = CMatrixFixed<double, 6, 6>;
using CMatrixDouble77 = CMatrixFixed<double, 7, 7>;
using CMatrixDouble13 = CMatrixFixed<double, 1, 3>;
using CMatrixDouble31 = CMatrixFixed<double, 3, 1>;
using CMatrixDouble12 = CMatrixFixed<double, 1, 2>;
using CMatrixDouble21 = CMatrixFixed<double, 2, 1>;
using CMatrixDouble61 = CMatrixFixed<double, 6, 1>;
using CMatrixDouble16 = CMatrixFixed<double, 1, 6>;
using CMatrixDouble71 = CMatrixFixed<double, 7, 1>;
using CMatrixDouble17 = CMatrixFixed<double, 1, 7>;
using CMatrixDouble51 = CMatrixFixed<double, 5, 1>;
using CMatrixDouble15 = CMatrixFixed<double, 1, 5>;
using CMatrixDouble41 = CMatrixFixed<double, 4, 1>;
using CMatrixDouble6_12 = CMatrixFixed<double, 6, 12>;
using CMatrixDouble12_6 = CMatrixFixed<double, 12, 6>;
using CMatrixDouble39 = CMatrixFixed<double, 3, 9>;
using CMatrixDouble93 = CMatrixFixed<double, 9, 3>;

using CMatrixFloat22 = CMatrixFixed<float, 2, 2>;
using CMatrixFloat23 = CMatrixFixed<float, 2, 3>;
using CMatrixFloat32 = CMatrixFixed<float, 3, 2>;
using CMatrixFloat33 = CMatrixFixed<float, 3, 3>;
using CMatrixFloat44 = CMatrixFixed<float, 4, 4>;
using CMatrixFloat66 = CMatrixFixed<float, 6, 6>;
using CMatrixFloat77 = CMatrixFixed<float, 7, 7>;
using CMatrixFloat13 = CMatrixFixed<float, 1, 3>;
using CMatrixFloat31 = CMatrixFixed<float, 3, 1>;
using CMatrixFloat12 = CMatrixFixed<float, 1, 2>;
using CMatrixFloat21 = CMatrixFixed<float, 2, 1>;
using CMatrixFloat61 = CMatrixFixed<float, 6, 1>;
using CMatrixFloat16 = CMatrixFixed<float, 1, 6>;
using CMatrixFloat71 = CMatrixFixed<float, 7, 1>;
using CMatrixFloat17 = CMatrixFixed<float, 1, 7>;
using CMatrixFloat51 = CMatrixFixed<float, 5, 1>;
using CMatrixFloat15 = CMatrixFixed<float, 1, 5>;
/**  @} */
}  // namespace math
}  // namespace mrpt
