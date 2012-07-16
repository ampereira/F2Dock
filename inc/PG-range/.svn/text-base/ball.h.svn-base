/*
  Copyright 2011 The University of Texas at Austin

        Authors: Muhibur Rasheed <muhibur@cs.utexas.edu>
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


#ifndef INC_BALL_H
#define INC_BALL_H
#include "Point.h"
//#include "curve_and_face.h"
#include <vector>
#include <slist.h>
#include <hash_map.h>
#include <string>
using namespace std;

struct eqcurve
{
  bool operator()(const char* c1, const char* c2) const
  {
    return strcmp(c1,c2)==0;
  }
};
/*
namespace __gnu_cxx  
{                                                                                             
  template<> struct hash<const boundary_curve*>                                                       
    {                                                                                           
      size_t operator()(const boundary_curve* x) const                                           
      {                                                                                         
	return hash<const boundary_curve*>()(x);                                              
      }                                                                                         
    };                                                                                          
}   
*/
struct ball {

  //Arrangement_2 arr_up, arr_down, arr_front, arr_rear, arr_right, arr_left;
  //hash_map<const char*, bool, hash<const char*>, eqcurve> lil_circle_list;
  float semi_side;
  Point center;
  float radius;
/*  bool exposed[6];
  slist<exposed_face*> top;
  slist<exposed_face*> bottom;
  slist<exposed_face*> front;
  slist<exposed_face*> rear;
  slist<exposed_face*> right;
  slist<exposed_face*> left;*/

  ball(Point p, float r) {
    center = p;
    radius = r;
    semi_side = r/1.73205;
 //   exposed[0] = exposed[1] = exposed[2] = exposed[3] = exposed[4] = exposed[5] = true;

  }
  bool equals(ball* b) {
    return (radius == b->radius && center == b->center);
  }
/*  bool isBuried(void) {
    return !(exposed[0] || exposed[1] || exposed[2] || exposed[3] || exposed[4] || exposed[5]);
    
  }*/
};
#endif
