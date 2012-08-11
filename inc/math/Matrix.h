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

// Matrix.h: interface for the Matrix class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MATRIX_H__AB7171AD_869D_4D29_A455_5919551F0B67__INCLUDED_)
#define AFX_MATRIX_H__AB7171AD_869D_4D29_A455_5919551F0B67__INCLUDED_

namespace CCVOpenGLMath {

class Quaternion;
class Vector;
class Ray;

class Matrix  
{
public:
	Matrix();
	Matrix(
		float m00, float m01, float m02, float m03,
		float m10, float m11, float m12, float m13,
		float m20, float m21, float m22, float m23,
		float m30, float m31, float m32, float m33
		);
	Matrix(const Quaternion& quat);
	virtual ~Matrix();
	Matrix(const Matrix& copy);
	Matrix& operator=(const Matrix& copy);

	void print() const;

	Matrix& set (
		float m00, float m01, float m02, float m03,
		float m10, float m11, float m12, float m13,
		float m20, float m21, float m22, float m23,
		float m30, float m31, float m32, float m33
		);
	Matrix& set(const Matrix& copy);
	Matrix& reset();
	inline float get(int row, int column) const;
	inline void set(int row, int column, float value);	
	const float* getMatrix() const;

	Vector operator*(const Vector& vec) const;
	Ray operator*(const Ray& ray) const;
	Matrix operator*(const Matrix& mat) const;
	Matrix& preMultiplication(const Matrix& mat);
	Matrix& postMultiplication(const Matrix& mat);

	Matrix inverse() const;
	Matrix inverseTranspose() const;
	Matrix transpose() const;
	
        int isAlmostEqual(const Matrix& mat);
        Matrix sqrt( );

	float determinant() const;

	static Matrix rotationX(float angle);
	static Matrix rotationY(float angle);
	static Matrix rotationZ(float angle);
	static Matrix translation(float x, float y, float z);
	static Matrix translation(const Vector& vec);
	static Matrix scale(float x, float y, float z);
        
protected:
	float m[16];
        
};
        
        inline float Matrix::get(int row, int column) const
        {
                return m[row + column*4];
        }
        
        inline void Matrix::set(int row, int column, float value)
        {
                m[ row + column * 4 ] = value;
        }
                
};

#endif // !defined(AFX_MATRIX_H__AB7171AD_869D_4D29_A455_5919551F0B67__INCLUDED_)
