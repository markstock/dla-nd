/*************************************************************
 *
 *  density.c - subroutines for the density reconstruction
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


#include "structs.h"

#if DIM>2
int add_particle_to_field_3d(cell_ptr,particle_ptr,field3_ptr);
#endif
int add_particle_to_field_2d(cell_ptr,particle_ptr,cell_ptr,field2_ptr);


#if DIM>2
/*
 *  Fill the cell_count array in sim with an accurate count of the cells
 */
int create_density_field_3d(cell_ptr top,cell_ptr curr_cell,field3_ptr ff) {

   int i,j,k;
   particle_ptr curr;

   // if this level is top, zero the field values and begin
   if (curr_cell == top) {
      for (i=0;i<ff->n[0];i++)
         for (j=0;j<ff->n[1];j++)
            for (k=0;k<ff->n[2];k++)
               ff->rho[i][j][k] = 0.0;
   }

   // add all particles in this cell
   curr = curr_cell->first;
   while (curr) {
      add_particle_to_field_3d(top,curr,ff);
      curr = curr->next;
   }

   // or add all particles from sublevels
   if (curr_cell->has_subcells)
      for (i=0;i<NCHILD;i++)
         create_density_field_3d(top,curr_cell->s[i],ff);

   return(0);
}


/*
 *  Add the effect of the particle to the density field
 */
int add_particle_to_field_3d(cell_ptr top,particle_ptr curr,field3_ptr ff) {

   int i,j,k;
   int start[3],end[3];
   FLOAT rad = 2.0*curr->rad;
   FLOAT radsq = rad*rad;
   FLOAT factor = 4.0 / (1000000.0*pow(curr->rad,3)*M_PI);
   FLOAT twopi = M_PI*2.0;
   FLOAT distsqx,distsqy,distsq,dist;
   FLOAT dummy = 0.0;

   // fprintf(stdout,"rad is %g\n",rad);

   // make them darker
   factor *= pow(top->max[0]-top->min[0],2);
   // rad *= 2.;
   // radsq = rad*rad;
   // factor = 4.0 / (1000000.0*pow(curr->rad,3)*M_PI);

   // find min and max cells affected
   for (i=0; i<3; i++) {
      start[i] = 0.5+(curr->x[i]-rad-top->min[i])/ff->d[i];
      end[i] = (curr->x[i]+rad-top->min[i])/ff->d[i] - 0.5;
      if (start[i] < 0) start[i] = 0;
      if (end[i] >= ff->n[i]) end[i] = ff->n[i]-1;
   }
   // fprintf(stdout,"Adding value in the ranges %d:%d %d:%d %d:%d\n",start[0],end[0],start[1],end[1],start[2],end[2]);

   // apply spherical kernel across that range
   for (i=start[0]; i<=end[0]; i++) {
      distsqx = pow( curr->x[0] - (top->min[0] + ff->d[0]*(0.5+(FLOAT)i)) , 2);
      // fprintf(stdout,"  i %d, dx %g\n",i,sqrt(distsqx));
      for (j=start[1]; j<=end[1]; j++) {
         distsqy = distsqx + pow( curr->x[1] - (top->min[1] + ff->d[1]*(0.5+(FLOAT)j)) , 2);
         for (k=start[2]; k<=end[2]; k++) {
            distsq = distsqy + pow( curr->x[2] - (top->min[2] + ff->d[2]*(0.5+(FLOAT)k)) , 2);;
            if (distsq < radsq) {
               dist = sqrt(distsq);
               // fprintf(stdout,"%d %d %d is at rad %g, within %g\n",i,j,k,dist,rad);
               // fprintf(stdout,"      adding %g\n",factor*(1.0+cos(M_PI*dist/rad)));
               ff->rho[i][j][k] += factor * (1.0+cos(M_PI*dist/rad));
               // dummy += factor * (1.0+cos(M_PI*dist/rad));
            }
         }
      }
   }
   // fprintf(stdout,"Added %g to the field\n",dummy);

   return(0);
}
#endif


/*
 *  Fill the cell_count array in sim with an accurate count of the cells
 */
int create_density_field_2d(cell_ptr top, cell_ptr curr_cell,
                            cell_ptr plotzone, field2_ptr ff) {

   int i,j,k;
   particle_ptr curr;

   // if this level is top, zero the field values and begin
   if (curr_cell == top) {
      for (i=0;i<ff->n[0];i++)
         for (j=0;j<ff->n[1];j++)
            ff->rho[i][j] = 0.0;
   }

   // add all particles in this cell
   curr = curr_cell->first;
   while (curr) {
      add_particle_to_field_2d(top,curr,plotzone,ff);
      curr = curr->next;
   }

   // or add all particles from sublevels
   if (curr_cell->has_subcells)
      for (i=0;i<NCHILD;i++)
         create_density_field_2d(top,curr_cell->s[i],plotzone,ff);

   return(0);
}


/*
 *  Add the effect of the particle to the density field
 */
int add_particle_to_field_2d(cell_ptr top, particle_ptr curr,
                             cell_ptr plotzone, field2_ptr ff) {

   int i,j,k;
   int lx = 0;		// sets the local (planar) axes to be 0..2, or x..z
   int ly = 1;
   int ld = 2;		// sets local depth direction --- not used anymore
   int start[DIM],end[DIM];
   FLOAT rad,radsq,factor;
   // FLOAT rad = 2.0*curr->rad;
   // FLOAT radsq = rad*rad;
   // FLOAT factor = 4.0 / (1000000.0*pow(curr->rad,3)*M_PI);
   // FLOAT factor = 4.0 * M_PI / pow(curr->rad,2);
   FLOAT twopi = M_PI*2.0;
   FLOAT distsqx,distsqy,distsq,dist;

   // scale factor due to ???
   // factor *= pow(top->max[0]-top->min[0],2);

   // make sure radius is large enough!
   if (curr->rad < ff->d[0]) rad = 2.0*ff->d[0];
   else rad = 2.0*curr->rad;
   // fprintf(stdout,"rad is %g\n",rad);
   radsq = rad*rad;

   // then properly set the factor
   // factor = 4.0 * M_PI / pow(rad,2);
   // factor = 0.5;
   factor = 1.0/rad;

   // find min and max cells affected

   // define the start and end bounds for each dimension
   start[lx] = 0.5+(curr->x[lx] - rad - plotzone->min[lx])/ff->d[0];
   end[lx] = (curr->x[lx] + rad - plotzone->min[lx])/ff->d[0] - 0.5;

   // make sure that these are in bounds!
   if (start[lx] < 0) start[lx] = 0;
   if (end[lx] >= ff->n[0]) end[lx] = ff->n[0];

   // do same for other axis
   start[ly] = 0.5+(curr->x[ly] - rad - plotzone->min[ly])/ff->d[1];
   end[ly] = (curr->x[ly] + rad - plotzone->min[ly])/ff->d[1] - 0.5;
   if (start[ly] < 0) start[ly] = 0;
   if (end[ly] >= ff->n[1]) end[ly] = ff->n[1];

   // fprintf(stdout,"Adding value in the ranges %d:%d %d:%d\n",start[lx],end[lx],start[ly],end[ly]);
   // fprintf(stdout,"  ff->d is %d:%d\n",ff->d[lx],ff->d[ly]);

   // apply spherical kernel across that range
   for (i=start[lx]; i<=end[lx]; i++) {
      // distsqx = pow( curr->x[lx] - (top->min[lx] + ff->d[lx]*(0.5+(FLOAT)i)) , 2);
      distsqx = pow( curr->x[lx] - (plotzone->min[lx] + ff->d[lx]*(0.5+(FLOAT)i)) , 2);
      // fprintf(stdout,"  i %d, dx %g\n",i,sqrt(distsqx));
      for (j=start[ly]; j<=end[ly]; j++) {
         distsq = distsqx + pow( curr->x[ly] -
                   (plotzone->min[ly]+ff->d[ly]*(0.5+(FLOAT)j)) , 2);
                   // (top->min[ly]+ff->d[ly]*(0.5+(FLOAT)j)) , 2);
//       distsqy = distsqx + pow( curr->x[ly] - (top->min[ly]+ff->d[ly]*(0.5+(FLOAT)j)) , 2);
//       for (k=start[ld]; k<=end[ld]; k++) {
//          distsq = distsqy + pow( curr->x[ld] - (top->min[ld]+ff->d[ld]*(0.5+(FLOAT)k)) , 2);
            if (distsq < radsq) {
               dist = sqrt(distsq);
               // fprintf(stdout,"%d %d %d is at rad %g, within %g\n",i,j,k,dist,rad);
               // fprintf(stdout,"      adding %g\n",factor*(1.0+cos(M_PI*dist/rad)));
               ff->rho[i][j] += factor * (1.0+cos(M_PI*dist/rad));
            }
//       }
      }
   }
   // fprintf(stdout,"Added %g to the field\n",dummy);

   return(0);
}

