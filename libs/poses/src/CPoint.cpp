/* +------------------------------------------------------------------------+
   |                     Mobile Robot Programming Toolkit (MRPT)            |
   |                          https://www.mrpt.org/                         |
   |                                                                        |
   | Copyright (c) 2005-2019, Individual contributors, see AUTHORS file     |
   | See: https://www.mrpt.org/Authors - All rights reserved.               |
   | Released under BSD License. See: https://www.mrpt.org/License          |
   +------------------------------------------------------------------------+ */

#include "poses-precomp.h"  // Precompiled headers

#include <mrpt/poses/CPoint.h>
#include <mrpt/poses/CPoint2D.h>
#include <mrpt/poses/CPoint3D.h>
#include <iostream>

using namespace mrpt::poses;

template <class DERIVEDCLASS>
void CPoint<DERIVEDCLASS>::fromString(const std::string& s)
{
	mrpt::math::CMatrixDouble m;
	if (!m.fromMatlabStringFormat(s))
		THROW_EXCEPTION("Malformed expression in ::fromString");
	ASSERT_EQUAL_(m.rows(), 1);
	ASSERT_EQUAL_(m.cols(), DERIVEDCLASS::static_size);
	for (int i = 0; i < DERIVEDCLASS::static_size; i++)
		derived().m_coords[i] = m(0, i);
}

template <class DERIVEDCLASS>
void CPoint<DERIVEDCLASS>::asString(std::string& s) const
{
	s = (!DERIVEDCLASS::is3DPoseOrPoint())
			? mrpt::format(
				  "[%f %f]", static_cast<const DERIVEDCLASS*>(this)->x(),
				  static_cast<const DERIVEDCLASS*>(this)->y())
			: mrpt::format(
				  "[%f %f %f]", static_cast<const DERIVEDCLASS*>(this)->x(),
				  static_cast<const DERIVEDCLASS*>(this)->y(),
				  static_cast<const DERIVEDCLASS*>(this)->m_coords[2]);
}

// Explicit instantiations:
template class CPoint<CPoint2D>;
template class CPoint<CPoint3D>;

template <class DERIVEDCLASS>
std::ostream& operator<<(std::ostream& o, const CPoint<DERIVEDCLASS>& p)
{
	o << "(" << p[0] << "," << p[1];
	if (p.is3DPoseOrPoint()) o << "," << p[2];
	o << ")";
	return o;
}
