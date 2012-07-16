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

// Vector.cpp: implementation of the Vector class.
//
//////////////////////////////////////////////////////////////////////

#include "Vector.h"
#include <math.h>

using CCVOpenGLMath::Tuple;
using CCVOpenGLMath::Vector;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

const float EPS = 0.00001f;

Vector::Vector(float x, float y, float z, float w) : Tuple(x,y,z,w)
{
}

Vector::Vector() : Tuple()
{
}

Vector::Vector(float* array)
{
	set(array);
}

Vector::~Vector()
{

}

Vector::Vector(const Vector& copy): Tuple(copy)
{
}

Vector& Vector::operator=(const Vector& copy)
{
	if (this!=&copy) {
		set(copy);
	}
	return *this;
}


Vector& Vector::set(float x, float y, float z, float w)
{
	Tuple::set(x,y,z,w);
	return *this;
}

Vector& Vector::set(float* array)
{
	Tuple::set(array);
	return *this;
}

Vector& Vector::set(const Vector& copy)
{
	Tuple::set(copy);
	return *this;
}

Vector Vector::cross(const Vector& vec) const
{
	return Vector(
		p[1]*vec[2] - p[2]*vec[1],
		p[2]*vec[0] - p[0]*vec[2],
		p[0]*vec[1] - p[1]*vec[0],		
		0.0f		
		);
}

Vector& Vector::crossEquals(const Vector& vec)
{
	return set(
		p[1]*vec[2] - p[2]*vec[1],
		p[2]*vec[0] - p[0]*vec[2],
		p[0]*vec[1] - p[1]*vec[0],		
		0.0f		
		);
}

float Vector::dot(const Vector& vec) const
{
	return p[0]*vec[0] + p[1]*vec[1] + p[2]*vec[2] + p[3]*vec[3]; 
}


Vector Vector::operator+(const Vector vec) const
{
	return Vector(
		p[0]+vec[0],
		p[1]+vec[1],
		p[2]+vec[2],
		p[3]+vec[3]);
}

Vector& Vector::operator+=(const Vector vec)
{
	return set(
		p[0]+vec[0],
		p[1]+vec[1],
		p[2]+vec[2],
		p[3]+vec[3]);
}

Vector Vector::operator-(const Vector vec) const
{
	return Vector(
		p[0]-vec[0],
		p[1]-vec[1],
		p[2]-vec[2],
		p[3]-vec[3]);
}

Vector& Vector::operator-=(const Vector vec)
{
	return set(
		p[0]-vec[0],
		p[1]-vec[1],
		p[2]-vec[2],
		p[3]-vec[3]);
}

Vector Vector::operator*(float scalar) const
{
	return Vector(p[0]*scalar, p[1]*scalar, p[2]*scalar, p[3]);
}

Vector& Vector::operator*=(float scalar)
{
	return set(p[0]*scalar, p[1]*scalar, p[2]*scalar, p[3]);
}

Vector Vector::operator-() const
{
	return Vector(-p[0], -p[1], -p[2], p[3]);
}

Vector& Vector::normalize()
{
	if ((float)fabs(p[3])<=EPS) {
		float length = (float)sqrt((double)(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]));
		return set(p[0]/length,p[1]/length,p[2]/length,0.0f);
	}
	else {
		return set(p[0]/p[3], p[1]/p[3], p[2]/p[3], 1.0f);
	}
}

float Vector::norm() const
{
	return (float)sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
}

bool Vector::isBad()
{
	return (p[0]==0.0f && p[1]==0.0f && p[2]==0.0f && p[3]==0.0f);
}

Vector Vector::badVector()
{
	return Vector(0.0f, 0.0f, 0.0f, 0.0f);
}

Vector* Vector::clone() const
{
	return new Vector(*this);
}

bool Vector::getCorners(double* min, double* max, CCVOpenGLMath::Vector* vCorner)
{
	if( !min || !max || !vCorner ) return false;

	vCorner[0].set((float)min[0], (float)min[1], (float)min[2], 1.f);
	vCorner[1].set((float)max[0], (float)min[1], (float)min[2], 1.f);
	vCorner[2].set((float)min[0], (float)max[1], (float)min[2], 1.f);
	vCorner[3].set((float)max[0], (float)max[1], (float)min[2], 1.f);
	vCorner[4].set((float)min[0], (float)min[1], (float)max[2], 1.f);
	vCorner[5].set((float)max[0], (float)min[1], (float)max[2], 1.f);
	vCorner[6].set((float)min[0], (float)max[1], (float)max[2], 1.f);
	vCorner[7].set((float)max[0], (float)max[1], (float)max[2], 1.f);

	return true;
}
