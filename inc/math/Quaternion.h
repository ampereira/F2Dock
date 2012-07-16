/*
  Copyright 2011 The University of Texas at Austin

        Authors: Rezaul Alam Chowdhury <shaikat@cs.utexas.edu>, Vinay Siddavanahalli <skvinay@cs.utexas.edu>
        Advisor: Chandrajit Bajaj <bajaj@cs.utexas.edu>

  This file is part of F2Dock.

  F2Dock is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License version 2.1 as published by the Free Software Foundation.

  F2Dock is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
*/

// Quaternion.h: interface for the Quaternion class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_QUATERNION_H__4A5485F3_5ADE_437D_A2C9_6D864A63C23B__INCLUDED_)
#define AFX_QUATERNION_H__4A5485F3_5ADE_437D_A2C9_6D864A63C23B__INCLUDED_

#include "Tuple.h"

namespace CCVOpenGLMath {

class Vector;
class Matrix;
class Ray;

class Quaternion : public Tuple  
{
public:
	Quaternion();
	virtual ~Quaternion();
	Quaternion(const Quaternion& copy);
	Quaternion& operator=(const Quaternion& copy);
	Quaternion(float w, float x, float y, float z);

	Quaternion& set(float w, float x, float y, float z);
	Quaternion& set(float* array);
	Quaternion& set(const Quaternion& copy);

	Quaternion operator*(const Quaternion& quat) const;
	Quaternion operator*(float scalar) const;
	Quaternion& operator*=(float scalar);
	Quaternion operator/(float scalar) const;
	Quaternion& operator/=(float scalar);

	Quaternion& preMultiply(const Quaternion& quat);
	Quaternion& postMultiply(const Quaternion& quat);
	Quaternion& rotate(float angle, float x, float y, float z);
	Quaternion& normalize();

	Quaternion conjugate() const;
	Quaternion inverse() const;
	float norm() const;

	Vector applyRotation(const Vector& vec) const;
	Ray applyRotation(const Ray& ray) const;
	Matrix buildMatrix() const;
	Quaternion power(double scalar);

        Quaternion slerp(const Quaternion& rhs, const float t) {
                Quaternion lhs = *this;
                Quaternion rv = lhs * ((lhs.inverse() * rhs).power(t));
                return rv;
        }
        // Quaternion operator+(const Quaternion& rhs);
        // Quaternion operator-(const Quaternion& rhs);
        // Quaternion slerp(const Quaternion& rhs, const float t);

        static Quaternion rotation(float angle, float x, float y, float z);
	static Quaternion rotation(float angle, const Vector& axis);
protected:
	explicit Quaternion(const Vector& vec);

};

};

#endif // !defined(AFX_QUATERNION_H__4A5485F3_5ADE_437D_A2C9_6D864A63C23B__INCLUDED_)
