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

// Vector.h: interface for the Vector class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_VECTOR_H__C73C6BDB_2D75_4770_B86A_E3329C07A817__INCLUDED_)
#define AFX_VECTOR_H__C73C6BDB_2D75_4770_B86A_E3329C07A817__INCLUDED_

#include "Tuple.h"

namespace CCVOpenGLMath {

class Vector : public Tuple  
{
public:
	Vector(float x, float y, float z, float w);
	Vector(float* array);
	Vector();
	virtual ~Vector();
	Vector(const Vector& copy);
	Vector& operator=(const Vector& copy);

	Vector& set(float x, float y, float z, float w);
	Vector& set(float* array);
	Vector& set(const Vector& copy);

	Vector cross(const Vector& vec) const;
	Vector& crossEquals(const Vector& vec);
	float dot(const Vector& vec) const;

	Vector operator+(const Vector vec) const;
	Vector& operator+=(const Vector vec);
	Vector operator-(const Vector vec) const;
	Vector& operator-=(const Vector vec);

	Vector operator*(float scalar) const;
	Vector& operator*=(float scalar);

	Vector operator-() const;

	Vector& normalize();
	float norm() const;

	bool isBad();

	static Vector badVector();
	static bool getCorners(double* min, double* max, CCVOpenGLMath::Vector* vCorner);
	virtual Vector* clone() const;
};

};

#endif // !defined(AFX_VECTOR_H__C73C6BDB_2D75_4770_B86A_E3329C07A817__INCLUDED_)
