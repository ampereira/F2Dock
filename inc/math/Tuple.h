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
// Tuple.h: interface for the Tuple class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_TUPLE_H__5AD4C604_B71A_4924_941A_15A0955C4E4E__INCLUDED_)
#define AFX_TUPLE_H__5AD4C604_B71A_4924_941A_15A0955C4E4E__INCLUDED_

namespace CCVOpenGLMath {

class Tuple  
{
public:
	Tuple(float x, float y, float z, float w);
	Tuple();
	virtual ~Tuple();
	Tuple(const Tuple& copy);
	Tuple& operator=(const Tuple& copy);

	void print() const;

	Tuple& set(float x, float y, float z, float w);
	Tuple& set(float* array);
	Tuple& set(const Tuple& copy);
	
	float& operator[](unsigned int i);
	const float& operator[](unsigned int i) const;

protected:
	float p[4];

};

};

#endif // !defined(AFX_TUPLE_H__5AD4C604_B71A_4924_941A_15A0955C4E4E__INCLUDED_)
