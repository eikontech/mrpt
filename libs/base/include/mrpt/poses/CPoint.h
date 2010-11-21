/* +---------------------------------------------------------------------------+
   |          The Mobile Robot Programming Toolkit (MRPT) C++ library          |
   |                                                                           |
   |                   http://mrpt.sourceforge.net/                            |
   |                                                                           |
   |   Copyright (C) 2005-2010  University of Malaga                           |
   |                                                                           |
   |    This software was written by the Machine Perception and Intelligent    |
   |      Robotics Lab, University of Malaga (Spain).                          |
   |    Contact: Jose-Luis Blanco  <jlblanco@ctima.uma.es>                     |
   |                                                                           |
   |  This file is part of the MRPT project.                                   |
   |                                                                           |
   |     MRPT is free software: you can redistribute it and/or modify          |
   |     it under the terms of the GNU General Public License as published by  |
   |     the Free Software Foundation, either version 3 of the License, or     |
   |     (at your option) any later version.                                   |
   |                                                                           |
   |   MRPT is distributed in the hope that it will be useful,                 |
   |     but WITHOUT ANY WARRANTY; without even the implied warranty of        |
   |     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
   |     GNU General Public License for more details.                          |
   |                                                                           |
   |     You should have received a copy of the GNU General Public License     |
   |     along with MRPT.  If not, see <http://www.gnu.org/licenses/>.         |
   |                                                                           |
   +---------------------------------------------------------------------------+ */
#ifndef CPOINT_H
#define CPOINT_H

#include <mrpt/poses/CPoseOrPoint.h>

namespace mrpt
{
	namespace poses
	{
		/** A base class for representing a point in 2D or 3D.
		 *   For more information refer to the <a href="http://www.mrpt.org/2D_3D_Geometry">2D/3D Geometry tutorial</a> online.
		 * \note This class is based on the CRTP design pattern
		 * \sa CPoseOrPoint, CPose
		 */
		template <class DERIVEDCLASS>
		class CPoint : public CPoseOrPoint<DERIVEDCLASS>
		{
		public:
			/** @name Methods common to all 2D or 3D points
			    @{ */

			/** Scalar addition of all coordinates.
			  * This is diferent from poses/point composition, which is implemented as "+" operators in classes derived from "CPose"
			  */
			template <class OTHERCLASS>
			inline void AddComponents(const OTHERCLASS &b)
			{
				const int dims = std::min( size_t(DERIVEDCLASS::static_size), size_t(OTHERCLASS::is3DPoseOrPoint() ? 3:2));
				for (int i=0;i<dims;i++)
					static_cast<DERIVEDCLASS*>(this)->m_coords[i]+= static_cast<const OTHERCLASS*>(&b)->m_coords[i];
			}

			/** Scalar multiplication. */
			inline void operator *=(const double s)
			{
				for (int i=0;i<DERIVEDCLASS::static_size;i++)
					static_cast<DERIVEDCLASS*>(this)->m_coords[i] *= s;
			}

			/** Return the pose or point as a 1x2 or 1x3 vector [x y] or [x y z] */
			inline void getAsVector(vector_double &v) const
			{
				v.resize(DERIVEDCLASS::static_size);
				for (int i=0;i<DERIVEDCLASS::static_size;i++)
					v[i] = static_cast<const DERIVEDCLASS*>(this)->m_coords[i];
			}
			//! \overload
			inline vector_double getAsVector() const { vector_double v; getAsVector(v); return v; }

			/** Returns the corresponding 4x4 homogeneous transformation matrix for the point(translation) or pose (translation+orientation).
			  * \sa getInverseHomogeneousMatrix
			  */
			void getHomogeneousMatrix(CMatrixDouble44 & out_HM ) const
			{
				out_HM.unit(4,1.0);
				out_HM.get_unsafe(0,3)= static_cast<const DERIVEDCLASS*>(this)->x();
				out_HM.get_unsafe(1,3)= static_cast<const DERIVEDCLASS*>(this)->y();
				if (DERIVEDCLASS::is3DPoseOrPoint())
					out_HM.get_unsafe(2,3)= static_cast<const DERIVEDCLASS*>(this)->m_coords[2];
			}

			/** Returns a human-readable textual representation of the object (eg: "[0.02 1.04]" )
			* \sa fromString
			*/
			void asString(std::string &s) const
			{
				s = (DERIVEDCLASS::is3DPoseOrPoint()) ?
					mrpt::format("[%f %f]", static_cast<const DERIVEDCLASS*>(this)->x(), static_cast<const DERIVEDCLASS*>(this)->y()) :
					mrpt::format("[%f %f %f]",static_cast<const DERIVEDCLASS*>(this)->x(), static_cast<const DERIVEDCLASS*>(this)->y(), static_cast<const DERIVEDCLASS*>(this)->m_coords[2]);
			}
			inline std::string asString() const { std::string s; asString(s); return s; }

			/** Set the current object value from a string generated by 'asString' (eg: "[0.02 1.04]" )
			* \sa asString
			* \exception std::exception On invalid format
			*/
			void fromString(const std::string &s)
			{
				CMatrixDouble  m;
				if (!m.fromMatlabStringFormat(s)) THROW_EXCEPTION("Malformed expression in ::fromString");
				ASSERT_EQUAL_(mrpt::math::size(m,1),1)
				ASSERT_EQUAL_(mrpt::math::size(m,2),DERIVEDCLASS::static_size)
				for (int i=0;i<DERIVEDCLASS::static_size;i++)
					static_cast<DERIVEDCLASS*>(this)->m_coords[i] = m.get_unsafe(0,i);
			}

			inline const double &operator[](unsigned int i) const { return static_cast<const DERIVEDCLASS*>(this)->m_coords[i]; }
			inline       double &operator[](unsigned int i)       { return static_cast<DERIVEDCLASS*>(this)->m_coords[i]; }

			/** @} */

		}; // End of class def.

	/** Dumps a point as a string [x,y] or [x,y,z]  */
	template <class DERIVEDCLASS>
	std::ostream &operator << (std::ostream& o, const CPoint<DERIVEDCLASS>& p)
	{
		o << "(" << p[0] << "," << p[1];
		if (p.is3DPoseOrPoint())	o << "," << p[2];
		o <<")";
		return o;
	}

	/** Used by STL algorithms */
	template <class DERIVEDCLASS>
	bool operator < (const CPoint<DERIVEDCLASS> &a, const CPoint<DERIVEDCLASS> &b)
	{
		if (a.x()<b.x()) return true;
		else
		{
			if (!a.is3DPoseOrPoint())
				return a.y()<b.y();
			else  if (a.y()<b.y())
				return true;
			else return a[2]<b[2];
		}
	}

	template <class DERIVEDCLASS>
	bool operator==(const CPoint<DERIVEDCLASS> &p1,const CPoint<DERIVEDCLASS> &p2)
	{
		for (int i=0;i<DERIVEDCLASS::static_size;i++)
			if (p1[i]!=p2[i])	return false;
		return true;
	}

	template <class DERIVEDCLASS>
	bool operator!=(const CPoint<DERIVEDCLASS> &p1,const CPoint<DERIVEDCLASS> &p2)
	{
		for (int i=0;i<DERIVEDCLASS::static_size;i++)
			if (p1[i]!=p2[i])	return true;
		return false;
	}



	} // End of namespace
} // End of namespace

#endif
