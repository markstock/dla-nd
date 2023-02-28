/*************************************************************
 *
 *  ndtree.c - subroutines controlling the ndtree structure
 *
 *  Copyright (C) 2000-2005  Mark J. Stock, mstock@umich.edu
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

cell_ptr build_new_tree(sim_ptr,cell_ptr);
int remake_box(sim_ptr,cell_ptr,cell_ptr);
int realloc_box(sim_ptr,cell_ptr);
int find_particle_bounds(cell_ptr,FLOAT*,FLOAT*);
int add_particle_to_cell(sim_ptr,particle_ptr,cell_ptr);
int split_cell(sim_ptr,cell_ptr);
int count_cells(sim_ptr,cell_ptr);
int replace_all_particles(sim_ptr,cell_ptr,cell_ptr);
int clean_up_all_cells(sim_ptr,cell_ptr);
int clean_up_cells(sim_ptr,cell_ptr,int);
int merge_cell(cell_ptr);
int find_all_cells_cm(sim_ptr,cell_ptr);
int find_cell_cm(cell_ptr,int);


/*
 *  Build a completely a new tree from an old tree
 */
cell_ptr build_new_tree(sim_ptr sim,cell_ptr oldtop) {

   int i;
   FLOAT min[DIM+1],max[DIM+1],maxsize,size;
   //particle_ptr curr,last;
   cell_ptr newtop = (CELL*)malloc(sizeof(CELL));

   // First, find the min and max bounds of the entire tree
   for (i=0;i<DIM+1;i++) min[i] = 9.9e+99;
   for (i=0;i<DIM+1;i++) max[i] = -9.9e+99;
   (void) find_particle_bounds(oldtop,min,max);

   // set some parameters for the new top cell
   newtop->level = 0;
   newtop->has_subcells = FALSE;
   newtop->first = NULL;
   for (i=0;i<NCHILD;i++)
      newtop->s[i] = NULL;
   newtop->num = 0;
   newtop->mass = 0.0;

   // then, set bounds for the new tree
   maxsize = 0.0;
   for (i=0;i<DIM;i++) {
      size = 0.66*(max[i]-min[i])+2.0*sim->new_part_rad;
      if (size > maxsize) maxsize = size;
   }
   for (i=0;i<DIM;i++) {
      newtop->mid[i] = 0.5*(min[i]+max[i]);
      newtop->min[i] = newtop->mid[i] - maxsize;
      newtop->max[i] = newtop->mid[i] + maxsize;
      //fprintf(stdout,"    %g %g\n",newtop->min[i],newtop->max[i]);
   }
   for (i=0; i<2; i++) sim->ff2->d[i] = (newtop->max[i]-newtop->min[i])/sim->ff2->n[i];
   for (i=0; i<2; i++) sim->out->d[i] = (newtop->max[i]-newtop->min[i])/sim->out->n[i];
   fprintf(stdout,"  making new tree with x bounds %g %g\n",newtop->min[0],newtop->max[0]);

   // lastly, march through all particles in the old tree and place
   // them in the new tree
   remake_box(sim,oldtop,newtop);

   // we hope this worked
   return(newtop);
}


/*
 * Determine the bounds of the entire particle cloud
 */
int remake_box(sim_ptr sim,cell_ptr curr_cell,cell_ptr newtop) {

   int i;
   particle_ptr curr,next;

   if (curr_cell->has_subcells) {

      for (i=0;i<NCHILD;i++) {
         (void) remake_box(sim,curr_cell->s[i],newtop);
         free(curr_cell->s[i]);
      }

   } else {

      curr = curr_cell->first;
      while (curr) {

         // set the pointer to the next particle in the current cell
         next = curr->next;

         // add this particle to the new top level box
         //fprintf(stdout,"in remake_box\n");
         add_particle_to_cell(sim,curr,newtop);

         // and reset the next pointer
         curr = next;
      }
   }

   return(0);
}


/*
 * Reallocate all particles into contiguous memory
 *
 * Using malloc like this sucks. Just make arrays for each cell.
 * Duh.
 */
int realloc_box (sim_ptr sim,cell_ptr curr_cell) {

   int i;
   particle_ptr curr,next,startp;

   if (curr_cell->has_subcells) {

      for (i=0;i<NCHILD;i++) {
         (void) realloc_box(sim,curr_cell->s[i]);
      }

   } else {

      // allocate particles in this cell to a contiguous spot in memory

      // first, count particles and malloc the space
      curr = curr_cell->first;
      curr_cell->num = 0;
      while (curr) {
         next = curr->next;
         curr_cell->num++;
         curr = next;
      }
      //fprintf(stderr,"%d particles in this box\n",curr_cell->num);

      // malloc the space
      startp = (PARTICLE*)malloc(curr_cell->num*sizeof(PARTICLE));

      // then, copy the particles
      curr = curr_cell->first;
      curr_cell->first = startp;

      while (curr) {

         // set the pointer to the next particle in the current cell
         next = curr->next;

         // copy this particle's attributes to the new location
         // still need to do this
         startp->index = curr->index;
         startp->rad = curr->rad;
         startp->mass = curr->mass;
         startp->bitfield = curr->bitfield;
         startp->counter = curr->counter;
         for (i=0;i<DIM;i++)
            startp->x[i] = curr->x[i];
         // the "next" pointer: somewhat complicated
         if (next==NULL) startp->next = NULL;
         else startp->next = startp+1;
         // the rootward particle pointer: very complicated
         startp->root = curr->root;

#if DIM==2
         fprintf(stderr,"oldp at %g %g\n",curr->x[0],curr->x[1]);
         fprintf(stderr,"newp at %g %g\n",startp->x[0],startp->x[1]);
#else
         fprintf(stderr,"oldp at %g %g %g\n",curr->x[0],curr->x[1],curr->x[2]);
         fprintf(stderr,"newp at %g %g %g\n",startp->x[0],startp->x[1],startp->x[2]);
#endif
         // segfaults here
         fprintf(stderr,"oldproot is index %d\n",curr->root->index);
         fprintf(stderr,"newproot is index %d\n",startp->root->index);

         // free the memory at the old location
         free(curr);

         // and reset the next pointer
         curr = next;
         startp++;
      }


      curr = curr_cell->first;
      while (curr) {
#if DIM==2
         fprintf(stderr,"newp at %g %g\n",curr->x[0],curr->x[1]);
#else
         fprintf(stderr,"newp at %g %g %g\n",curr->x[0],curr->x[1],curr->x[2]);
#endif
         curr = curr->next;
      }

      // reset last particle pointer to point to NULL
      //startp--;
      //startp->next = NULL;
   }

   return(0);
}


/*
 * Determine the bounds of the entire particle cloud
 */
int find_particle_bounds(cell_ptr curr_cell,FLOAT* min,FLOAT* max) {

   int i;
   particle_ptr curr;

   if (curr_cell->has_subcells) {

      for (i=0;i<NCHILD;i++)
         (void) find_particle_bounds(curr_cell->s[i],min,max);

   } else {

      curr = curr_cell->first;
      while (curr) {
         for (i=0;i<DIM;i++) {
            if (curr->x[i] < min[i]) min[i] = curr->x[i];
            if (curr->x[i] > max[i]) max[i] = curr->x[i];
         }
         curr = curr->next;
      }
   }

   return(0);
}


/*
 *  Add the particle to the given cell
 *
 *  returns -1 if particle is outside of top's bounds
 */
int add_particle_to_cell(sim_ptr sim,particle_ptr curr,cell_ptr cell){

   int i,d;
   int retval = 0;
   particle_ptr oldfirst;

   if (cell->has_subcells) {
      // if the cell has subcells, then the particle must go in one of those
      i=0;
      for (d=0;d<DIM;d++) {
         if (curr->x[d] > cell->mid[d]) i+=pow(2,d);
      }
      // fprintf(stdout,"  putting particle %d in cell %d\n",curr->index,i);
      retval = add_particle_to_cell(sim,curr,cell->s[i]);

   } else {
      // otherwise the particle must go in this cell
      oldfirst = cell->first;
      curr->next = oldfirst;
      cell->first = curr;
      cell->num++;
      //fprintf(stdout,"added particle to cell at level %d, now has %d particles\n",cell->level,cell->num);

      // is the particle actually inside of the bounds of this cell?
      for (d=0;d<DIM;d++) {
         if (curr->x[d] > cell->max[d]) return(-1);
         if (curr->x[d] < cell->min[d]) return(-1);
      }

      if (cell->num > sim->max_parts_in_cell && cell->level+1 < sim->max_levels) {
         // now there are too many particles in this cell, split the cell into subcells
         // fprintf(stdout,"  cell at level %d has %d particles, splitting\n",cell->level,cell->num);
         split_cell(sim,cell);
      }
   }

   return(retval);
}


/*
 *  Split this cell into NCHILD subcells, and distribute all particles
 */
int split_cell(sim_ptr sim,cell_ptr parent){

   int i,j,d;
   FLOAT start[DIM],end[DIM];
   particle_ptr curr,thenext;

   // create the NCHILD new cells
   for (i=0;i<NCHILD;i++) {
      for (d=0;d<DIM;d++) {
         // if (mod(i,pow(2,d+1)) == 1) {
         if (mod(i,pow(2,d+1)) < pow(2,d)) {
            start[d] = parent->min[d];
            end[d] = parent->min[d] + (parent->max[d]-parent->min[d])/2.0;
         } else {
            start[d] = parent->min[d] + (parent->max[d]-parent->min[d])/2.0;
            end[d] = parent->max[d];
         }
      }
      parent->s[i] = (CELL*)malloc(sizeof(CELL));
      parent->s[i]->level = parent->level+1;
      parent->s[i]->has_subcells = FALSE;
      parent->s[i]->first = NULL;
      parent->s[i]->num = 0;
      for (j=0;j<DIM;j++) {
         parent->s[i]->min[j] = start[j];
         parent->s[i]->mid[j] = (start[j]+end[j])/2.0;
         parent->s[i]->max[j] = end[j];
      }
      /* fprintf(stdout,"    new cell %d %d %d has coords %g %g  %g %g  %g %g\n",i,j,k,xstart,xend,ystart,yend,zstart,zend); */
   }

   /* then, iterate thru all particles and distribute them */
   /* fprintf(stdout,"    parent cell had %d particles\n",parent->num); */
   curr = parent->first;
   while (curr) {
      thenext = curr->next;
      i = 0;
      for (d=0;d<DIM;d++) {
         if (curr->x[d] > parent->mid[d]) i+=pow(2,d);
      }
      //fprintf(stdout,"in split_cell\n");
      add_particle_to_cell(sim,curr,parent->s[i]);
      curr = thenext;
   }
   parent->first = NULL;
   parent->num = 0;

   /*
   for (i=0;i<NCHILD;i++)
      fprintf(stdout,"    child cell %d now has %d particles\n",i,parent->s[i]->num);
      if (parent->s[i]->first)
         fprintf(stdout,"      one is at %g %g %g\n",parent->s[i]->first->x[0],
                                                     parent->s[i]->first->x[1],
                                                     parent->s[i]->first->x[2]);
      }
   fprintf(stdout,"    parent cell now has %d particles\n",parent->num);
   */

   /* set the flag and we're ready to go */
   parent->has_subcells = TRUE;

   return(0);
}


/*
 *  Fill the cell_count array in sim with an accurate count of the cells
 */
int count_cells(sim_ptr sim,cell_ptr curr_cell) {

   int i;

   /* count this level */
   sim->cell_count[(int)curr_cell->level]++;

   /* and count all sublevels */
   if (curr_cell->has_subcells)
      for (i=0;i<NCHILD;i++)
         count_cells(sim,curr_cell->s[i]);

   return(0);
}


/*
 *  Replace the particles by completely re-allocating them
 */
int replace_all_particles(sim_ptr sim,cell_ptr top,cell_ptr curr_cell){

   int i,retval,out_of_bounds;
   particle_ptr curr,last,moving_particle;

   /* for each particle, see if it needs to be moved, if so,
    * move it back up to the top cell where it will reallocate
    * downward to an appropriate cell */

   if (curr_cell->has_subcells) {

      for (i=0;i<NCHILD;i++) {
         retval = replace_all_particles(sim,top,curr_cell->s[i]);
         // if a particle is completely OOB, remake the entire tree
         if (retval < 0) return(retval);
      }

   } else {

      curr = curr_cell->first;
      last = curr;
      while (curr) {
         out_of_bounds = FALSE;

         /* check to see if the particle is no longer in the cell's bounds */
         /* if the acc memory is freed, one of these x[i]'s can't be seen */
         for (i=0;i<DIM;i++)
            if (curr->x[i] < curr_cell->min[i] || curr->x[i] > curr_cell->max[i]) {
               out_of_bounds = TRUE;
               break;
            }

         if (out_of_bounds) {

            moving_particle = curr;

            // see if it's out-of-bounds of the top-level cell first
            for (i=0;i<DIM;i++)
               if (moving_particle->x[i] < top->min[i] || moving_particle->x[i] > top->max[i]) {
                  return(-1);
               }

            // if not, then:
            /* remove this particle from the current cell */
            /* if (curr->next) fprintf(stdout,"    %d is out of bounds %lf %lf %lf last is %d, next is %d\n",curr->index,curr_cell->min[i],curr->x[i],curr_cell->max[i],last->index,curr->next->index); */
            /* else fprintf(stdout,"    %d is out of bounds %lf %lf %lf last is %d, no next\n",curr->index,curr_cell->min[i],curr->x[i],curr_cell->max[i],last->index); */
            
            if (curr == curr_cell->first) {
               curr_cell->first = curr->next;
               last = curr;
            } else {
               /* fprintf(stdout,"    last is %d, last->next is %d\n",last->index,last->next->index); */
               /* fprintf(stdout,"    curr is %d, curr->next is %d\n",curr->index,curr->next->index); */
               last->next = curr->next;
            }
            curr = last->next;
            curr_cell->num--;

            /* and add it to the top-level cell, but only if it not completely out-of-bounds */
            // out_of_bounds = FALSE;
            // for (i=0;i<DIM;i++)
               // if (moving_particle->x[i] < top->min[i] || moving_particle->x[i] > top->max[i]) {
                  // out_of_bounds = TRUE;
                  // i+=DIM;
               // }
            // if (out_of_bounds) {
               // fprintf(stdout,"  particle removed because loc is %g %g %g\n",moving_particle->x[0],moving_particle->x[1],moving_particle->x[2]);
               // free(moving_particle);
            // } else {
               /* it's OK to add it to the top-level cell */
               //fprintf(stdout,"in replace_all_particles\n");
               add_particle_to_cell(sim,moving_particle,top);
            // }


         /* if it's not out of bounds */
         } else {

            /* if (curr->next) fprintf(stdout,"    %d is OK, last is %d, next is %d\n",curr->index,last->index,curr->next->index); */
            /* else fprintf(stdout,"    %d is OK, last is %d, no next\n",curr->index,last->index); */

            /* jump to the next one */
            last = curr;
            curr = curr->next;
         }

      }

   }

   return(0);
}


/*
 *  Then, march backwards from the deepest level and remove
 *  cells with too few particles
 */
int clean_up_all_cells(sim_ptr sim,cell_ptr top){

   int i;

   /* fprintf(stdout,"\n  doing clean_up_all_cells\n"); */

   /* if max_levels is 10, iterate from 9 to 0 */
   for (i=sim->max_levels-1; i>=0; i--) {
      /* fprintf(stdout,"  level %d\n",i); */
      clean_up_cells(sim,top,i);
   }

   /* return 0 if all went well */
   return(0);
}


/*
 *  For this level, either count particles or check to see if the
 *  NCHILD child subcells could be merged together
 */
int clean_up_cells(sim_ptr sim,cell_ptr curr_cell,int do_level){

   int i,cnt;
   int any_sub_subcells;
   particle_ptr curr;

   /* if there are subcells, then there are no particles attached to this level */
   if (curr_cell->has_subcells) {

      /* if we're to compute this level, then all levels beneath are done */
      if (curr_cell->level == do_level) {

         /* count the particles in the NCHILD subcells, AND confirm that all
          * subcells have no subcells of their own */
         cnt = 0;
         any_sub_subcells = FALSE;
         for (i=0;i<NCHILD;i++) {
            cnt += curr_cell->s[i]->num;
            if (curr_cell->s[i]->has_subcells) any_sub_subcells = TRUE;
         }

         /* make sure particle count was consistent */
         /* if (cnt != curr_cell->num) {
            fprintf(stderr,"WARNING (clean_up_cells): particle count inconsistent\n");
            fprintf(stderr,"    Count is %d, previous saved count was %d\n",cnt,curr_cell->num);
         } */
         curr_cell->num = cnt;

         /* if there are no sub-subcells AND the particle count is too low, merge */
         if (!any_sub_subcells && cnt<sim->max_parts_in_cell) {
            /* fprintf(stdout,"    merging subcells at level %d, only %d particles\n",curr_cell->level+1,cnt); */
            merge_cell(curr_cell);
         }

      /* otherwise, run this subroutine on all of the subcells */
      } else {
         for (i=0;i<NCHILD;i++)
            clean_up_cells(sim,curr_cell->s[i],do_level);
      }

   /* and if there are no subcells, doublecheck the particle count */
   } else {

      /* cycle thru all particles, make sure the particle count is accurate */
      curr = curr_cell->first;
      cnt = 0;
      while (curr) {
         cnt++;
         curr = curr->next;
      }

      /* compare this count to the cell's count */
      /* if (cnt != curr_cell->num) {
         fprintf(stderr,"WARNING (clean_up_cells): particle count inconsistent\n");
         fprintf(stderr,"    Count is %d, previous saved count was %d\n",cnt,curr_cell->num);
      } */
      curr_cell->num = cnt;
   }

   return(0);
}


/*
 *  Merge this cell's NCHILD children back into the cell
 */
int merge_cell(cell_ptr curr_cell){

   int i;
   particle_ptr scurr,snext,oldfirst;

   /* continuity check on parent cell */
   if (curr_cell->first) {
      fprintf(stderr,"WARNING (merge_cell): parent cell had particles attached to it!\n");
   }

   /* reset this to zero, and add the particles from the subcells */
   curr_cell->num = 0;

   /* then, move all of the NCHILD subcells' particles to the parent's list */
   for (i=0;i<NCHILD;i++) {

      /* loop through the Subcell's CURRent particle */
      scurr = curr_cell->s[i]->first;
      while (scurr) {

         /* set the Subcell's NEXT particle now */
         snext = scurr->next;

         /* remove the particle from the subcell */
         curr_cell->s[i]->first = snext;

         /* add the particle to the parent cell */
         oldfirst = curr_cell->first;
         scurr->next = oldfirst;
         curr_cell->first = scurr;
         curr_cell->num++;

         scurr = snext;
      }
   }

   /* set the appropriate pointers */
   curr_cell->has_subcells = FALSE;

   /* free the memory from the NCHILD subcells */
   for (i=0;i<NCHILD;i++) {
      free(curr_cell->s[i]);
      curr_cell->s[i] = NULL;
   }

   return(0);
}


/*
 * find the centers of mass of all cells
 */
int find_all_cells_cm(sim_ptr sim,cell_ptr top){

   int i;

   /* if max_levels is 10, iterate from 9 to 0 */
   for (i=sim->max_levels-1; i>=0; i--) {
      /* fprintf(stdout,"\n  doing level %d\n",i); */
      find_cell_cm(top,i);
   }

   /* return 0 if all went well */
   return(0);
}


/*
 * find the centers of mass of all cells
 */
int find_cell_cm(cell_ptr curr_cell,int do_level){

   int do_print = FALSE;
   int i,l;
   FLOAT com[DIM];
   FLOAT mtot;
   particle_ptr curr;

   for (i=0;i<DIM;i++) com[i] = 0.0;
   mtot = 0.0;

   /* if there are subcells, then there are no particles attached to this level */
   if (curr_cell->has_subcells) {

      /* if we're to compute this level, then all levels beneath are done */
      if (curr_cell->level == do_level) {
         curr_cell->num = 0;
         for (i=0;i<NCHILD;i++) {
            mtot += curr_cell->s[i]->mass;
            for (l=0;l<DIM;l++) com[l] += curr_cell->s[i]->mass * curr_cell->s[i]->cm[l];
            curr_cell->num += curr_cell->s[i]->num;
         }

         /* if, for some reason, total mass is zero, set com to origin */
         if (fabs(mtot) > 1e-6) {
            for (i=0;i<DIM;i++) com[i] = com[i]/mtot;
            if (do_print) fprintf(stdout,"  level %d cell has mass %g, com %g %g %g\n",curr_cell->level,mtot,com[0],com[1],com[2]);
         } else {
            for (i=0;i<DIM;i++) com[i] = 0.0;
            if (do_print) fprintf(stdout,"  level %d cell has mass %g\n",curr_cell->level,mtot);
         }

         /* set sums in cell structure */
         curr_cell->mass = mtot;
         for (i=0;i<DIM;i++) curr_cell->cm[i] = com[i];

      /* otherwise, run this subroutine on all of the subcells */
      } else {
         for (i=0;i<NCHILD;i++)
            find_cell_cm(curr_cell->s[i],do_level);
      }
   } else {

      /* cycle thru all particles, sum mass and center of mass */
      curr = curr_cell->first;
      while (curr) {
         mtot += curr->mass;
         for (i=0;i<DIM;i++) com[i] += curr->x[i] * curr->mass;
         curr = curr->next;
      }

      /* if, for some reason, the mass is too low, set com to origin */
      if (fabs(mtot) > 1e-6) {
         for (i=0;i<DIM;i++) com[i] = com[i]/mtot;
         if (do_print) fprintf(stdout,"  level %d particles have mass %g, com %g %g %g\n",curr_cell->level,mtot,com[0],com[1],com[2]);
      } else {
         for (i=0;i<DIM;i++) com[i] = 0.0;
         if (do_print) fprintf(stdout,"  level %d particles have mass %g\n",curr_cell->level,mtot);
      }

      /* set the sums in the cell structure */
      curr_cell->mass = mtot;
      for (i=0;i<DIM;i++) curr_cell->cm[i] = com[i];

   }

   return(0);
}

