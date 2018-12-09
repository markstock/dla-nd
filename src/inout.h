/*************************************************************
 *
 *  inout.h - input/output subroutines for dla-nd
 *
 *  Copyright (C) 2000-14  Mark J. Stock, mstock@umich.edu
 *
 *  This file is part of dla-nd.
 *
 *  dla-nd is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  dla-nd is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with dla-nd; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 * 
 *********************************************************** */

#pragma once

#include "structs.h"

extern int write_output(sim_ptr,cell_ptr);
extern int read_png (char*, int*, int*, int, FLOAT***, FLOAT, FLOAT, FLOAT***, FLOAT, FLOAT, FLOAT***, FLOAT, FLOAT);
extern png_byte** allocate_2d_array_pb (int,int,int);
