/*************************************************************
 *
 *  density.h - subroutines for the density reconstruction
 *
 *  Copyright (C) 2001,2003-4  Mark J. Stock, mstock@umich.edu
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

#if (DIM > 2)
int create_density_field_3d(cell_ptr,cell_ptr,cell_ptr,field3_ptr);
#endif
int create_density_field_2d(cell_ptr,cell_ptr,cell_ptr,field2_ptr);

