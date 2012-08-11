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


#ifndef INC_CELL_H
#define INC_CELL_H
#include "ball.h"
#include <vector>

#include <cmath>

using namespace std;

struct Point {
  float x;
  float y;
  float z;
  Point() {}
  Point(float a, float b, float c) {
    x = a;
    y = b;
    z = c;
  }
  bool operator==(const Point& p) {
    return (x==p.x && y==p.y && z==p.z);
  }
 
  float distsq(Point a, Point b) {
  float res;
  double dx, dy, dz;
  dx = a.x - b.x;
  dy = a.y - b.y;
  dz = a.z - b.z;  
  res = dx*dx + dy*dy + dz*dz;
  return res;
}

};

struct cellID {
  int x;
  int y;
  int z;
};

struct lineID {
  int y;
  int z;
};

struct planeID {
  int z;
};

class atom;

struct gridcell {
  cellID ID;
  vector<atom*> balls;
  gridcell(int a, int b, int c) {
    ID.x = a;
    ID.y = b;
    ID.z = c;
  }
};
#endif
