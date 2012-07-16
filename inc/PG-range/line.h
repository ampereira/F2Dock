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


#ifndef INC_LINE_H
#define INC_LINE_H
#include "cell.h"
#include "1D-IntegerRange.h"
//#include "ddriver.h"
using namespace std;

struct gridcell;

struct line {
  lineID ID;
  IntegerRange<gridcell*> RR;
  line(int b, int c) {
    ID.y = b;
    ID.z = c;
  }
};
#endif
