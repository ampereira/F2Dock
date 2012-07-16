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
// Ray.cpp: implementation of the Ray class.
//
//////////////////////////////////////////////////////////////////////

#include "Ray.h"

#include <math.h>
#include <stdio.h>

using CCVOpenGLMath::Vector; 
using CCVOpenGLMath::Ray; 

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Ray::Ray() : m_Origin(0.0f, 0.0f, 0.0f, 1.0f), m_Dir(0.0f, 0.0f, 1.0f, 0.0f)
{
}

Ray::Ray(const Vector& origin, const Vector& dir) : m_Origin(origin), m_Dir(dir)
{
	
}

Ray::~Ray()
{

}

Vector Ray::getPointOnRay(float t) const
{
	return m_Origin+m_Dir*t;
}

float Ray::nearestTOnXAxis(Vector Origin) const
{
	Origin[3] = 0;
	Ray ray(m_Origin-Origin, m_Dir);
	float distance = ray.distanceToXAxis(Origin);
	float t = -(ray.m_Origin[1]*ray.m_Dir[1] + ray.m_Origin[2]*ray.m_Dir[2])/
		((ray.m_Dir[1]*ray.m_Dir[1]+ray.m_Dir[2]*ray.m_Dir[2]) );
	return t;
}

float Ray::nearestTOnYAxis(Vector Origin) const
{
	Origin[3] = 0;
	Ray ray(m_Origin-Origin, m_Dir);
	float distance = ray.distanceToYAxis(Origin);
	float t = -(ray.m_Origin[0]*ray.m_Dir[0] + ray.m_Origin[2]*ray.m_Dir[2])/
		((ray.m_Dir[0]*ray.m_Dir[0]+ray.m_Dir[2]*ray.m_Dir[2]));
	return t;
}

float Ray::nearestTOnZAxis(Vector Origin) const
{
	Origin[3] = 0;
	Ray ray(m_Origin-Origin, m_Dir);
	float distance = ray.distanceToZAxis(Origin);
	float t = -(ray.m_Origin[0]*ray.m_Dir[0] + ray.m_Origin[1]*ray.m_Dir[1])/
		((ray.m_Dir[1]*ray.m_Dir[1]+ray.m_Dir[0]*ray.m_Dir[0]));
	return t;
}

Vector Ray::nearestPointOnXAxis(Vector Origin) const
{
	Origin[3] = 0;
	float t = nearestTOnXAxis(Origin);
	Vector result = getPointOnRay(t);
	// project to axis
	result[1] = Origin[1];
	result[2] = Origin[2];
	//result+=Origin;
	return result;
}

Vector Ray::nearestPointOnYAxis(Vector Origin) const
{
	Origin[3] = 0;
	float t = nearestTOnYAxis(Origin);
	Vector result = getPointOnRay(t);
	// project to axis
	result[0] = Origin[0];
	result[2] = Origin[2];
	//result+=Origin;
	return result;
}

Vector Ray::nearestPointOnZAxis(Vector Origin) const
{
	Origin[3] = 0;
	float t = nearestTOnZAxis(Origin);
	Vector result = getPointOnRay(t);
	// project to axis
	result[0] = Origin[0];
	result[1] = Origin[1];
	return result;
}

float Ray::distanceToXAxis(Vector Origin) const
{
	Origin[3] = 0;
	Ray ray(m_Origin-Origin, m_Dir);
	return (float)fabs( 
		( ray.m_Origin[2]*ray.m_Dir[1]-ray.m_Origin[1]*m_Dir[2] ) /
		(float)sqrt( ray.m_Dir[2]*ray.m_Dir[2] + ray.m_Dir[1]*ray.m_Dir[1] )
		);
}

float Ray::distanceToYAxis(Vector Origin) const
{
	Origin[3] = 0;
	Ray ray(m_Origin-Origin, m_Dir);
	return (float)fabs( 
		( ray.m_Origin[2]*ray.m_Dir[0]-ray.m_Origin[0]*m_Dir[2] ) /
		(float)sqrt( ray.m_Dir[2]*ray.m_Dir[2] + ray.m_Dir[0]*ray.m_Dir[0] )
		);
}

float Ray::distanceToZAxis(Vector Origin) const
{
	Origin[3] = 0;
	Ray ray(m_Origin-Origin, m_Dir);
	return (float)fabs( 
		( ray.m_Origin[0]*ray.m_Dir[1]-ray.m_Origin[1]*m_Dir[0] ) /
		(float)sqrt( ray.m_Dir[0]*ray.m_Dir[0] + ray.m_Dir[1]*ray.m_Dir[1] )
		);
}

/*********************************************************/
/*                                                       */
/*  Returns false if there is no intersection.           */
/*  Else, it returns both the points and the values of   */
/*  the parameter where the intersections took place.    */
/*                                                       */
/*********************************************************/
bool Ray::intersectSphere( Vector center, float radius, Vector *point1, Vector* point2, float *distance1, float* distance2 )
{
	if( !point1 || !point2 ) return false;
	if( radius <= 0 ) return false;


	/// solve quadratic equation /////
	////	A = Xd^2 + Yd^2 + Zd^2
	////	B = 2 * (Xd * (X0 - Xc) + Yd * (Y0 - Yc) + Zd * (Z0 - Zc))
	////	C = (X0 - Xc)^2 + (Y0 - Yc)^2 + (Z0 - Zc)^2 - Sr^2
	///////////////////////////////////


	float A =	m_Dir[0]*m_Dir[0] + 
				m_Dir[1]*m_Dir[1] + 
				m_Dir[2]*m_Dir[2];
	  
	float B = 2* (	m_Dir[0] * (m_Origin[0] - center[0]) +
					m_Dir[1] * (m_Origin[1] - center[1]) +
					m_Dir[2] * (m_Origin[2] - center[2]) );

	float C = (m_Origin[0] - center[0])*(m_Origin[0] - center[0]) +
			  (m_Origin[1] - center[1])*(m_Origin[1] - center[1]) +
			  (m_Origin[2] - center[2])*(m_Origin[2] - center[2]) -
			  radius*radius;

	float discriminant = B*B - 4*A*C;
	if( discriminant < 0 ) return false;

	*distance1 = (float)(( -B - sqrt( discriminant ) ) / ( 4.0 * A * C ));
	*distance2 = (float)(( -B + sqrt( discriminant ) ) / ( 4.0 * A * C ));

	*point1 = m_Origin + m_Dir * (*distance1);
	*point2 = m_Origin + m_Dir * (*distance2);

	return true;
}

void Ray::print()
{
	printf("[%5.3lf %5.3lf %5.3lf %5.3lf ] [%5.3lf %5.3lf %5.3lf %5.3lf ]\n", 
		m_Origin[0], m_Origin[1], m_Origin[2], m_Origin[3], 
		m_Dir[0], m_Dir[1], m_Dir[2], m_Dir[3] );
}
