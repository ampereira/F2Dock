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

// Tuple.cpp: implementation of the Tuple class.
//
//////////////////////////////////////////////////////////////////////

#include "Tuple.h"
#include <stdio.h>

using CCVOpenGLMath::Tuple;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Tuple::Tuple(float x, float y, float z, float w)
{
	set(x,y,z,w);
}

Tuple::Tuple()
{
	set(0.0, 0.0, 0.0, 0.0);
}

Tuple::~Tuple()
{

}

Tuple::Tuple(const Tuple& copy)
{
	set(copy);
}

Tuple& Tuple::operator=(const Tuple& copy)
{
	return set(copy);
}

Tuple& Tuple::set(float x, float y, float z, float w)
{
	p[0] = x;
	p[1] = y;
	p[2] = z;
	p[3] = w;
	return *this;
}

Tuple& Tuple::set(float* array)
{
	p[0] = array[0];
	p[1] = array[1];
	p[2] = array[2];
	p[3] = array[3];
	return *this;
}


Tuple& Tuple::set(const Tuple& copy)
{
	if (this!=&copy) {
		p[0] = copy.p[0];
		p[1] = copy.p[1];
		p[2] = copy.p[2];
		p[3] = copy.p[3];
	}
	return *this;
}



float& Tuple::operator[](unsigned int i)
{
	return p[i];
}


const float& Tuple::operator[](unsigned int i) const
{
	return p[i];	
}

void Tuple::print() const
{
	printf("[%5.3lf %5.3lf %5.3lf %5.3lf ]\n", p[0], p[1], p[2], p[3] );
}

