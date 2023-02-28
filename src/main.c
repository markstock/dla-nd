/*************************************************************
 *
 *  main.c - Arbitrary-dimensional diffusion-limited aggregation
 *
 *  Copyright (C) 2000-18  Mark J. Stock, mstock@umich.edu
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
#include "inout.h"

int run_sim (sim_ptr,cell_ptr);
int diffuse_new_particle (sim_ptr,cell_ptr,FLOAT*);
void pick_point_on_sphere (FLOAT*,FLOAT*,FLOAT,BULKVEL,FLOAT*);
FLOAT find_dist_to_closest (cell_ptr,FLOAT*,FLOAT,particle_ptr,particle_ptr*);
FLOAT find_dist_to_closest2 (cell_ptr,FLOAT*,FLOAT,particle_ptr,particle_ptr*);
FLOAT find_dist_to_farthest (cell_ptr,FLOAT*,FLOAT);
void sum_tipward_mass (cell_ptr);
void add_to_rootward_masses (particle_ptr);
void zero_all_masses (cell_ptr);
void increment_counters (cell_ptr);
void zero_counters (cell_ptr);
FLOAT vec_dot (FLOAT*,FLOAT*);
FLOAT vec_length (FLOAT*);
FLOAT vec_length_sq (FLOAT*);
void gaussian_rand (FLOAT*,FLOAT*);
void gaussian_rand2 (FLOAT*,FLOAT*);
void gaussian_rand3 (FLOAT*,FLOAT*);
FLOAT solve_quartic (FLOAT*);
int correct_for_chiral (FLOAT*, particle_ptr, FLOAT, FLOAT);
int find_start_point (sim_ptr,cell_ptr,FLOAT*,FLOAT*,FLOAT*);


int main(int argc,char **argv) {

   SIMP sim_props;		// initialize struct for simulation props
   sim_ptr sim = &sim_props;
   CELL top_cell;		// initialize struct for ndtree
   cell_ptr top = &top_cell;

   // get things ready
   set_defaults (sim,top);

   // Parse command-line args, includes reading input files
   parse_args (argc,argv,sim,top);

   // Initialize the top-level cell and participating particles
   initialize_system (sim,top);

   // Do the time integration
   run_sim (sim,top);

   fprintf(stderr,"\nDone.\n");
   exit(0);
}



/*
 *  Run the time integration of the flowfield simulation
 */
int run_sim(sim_ptr sim,cell_ptr top) {

   int d,retval,next_update;
   FLOAT this_step_time,total_time,vel[DIM];//,start[DIM];
   static int wrote_stats_header_line = FALSE;

   next_update = 100;

   sim->next_output_step = sim->output_step;

   // open the stats file if necessary
   if (sim->write_stats) {
      sim->statsfp = fopen(sim->stats_file,"w");
      if (sim->statsfp==NULL) {
         fprintf(stderr,"Could not open output file %s\n",sim->stats_file);
         fflush(stderr);
         exit(0);
      }
   }

   // make a good tree
   top = build_new_tree(sim,top);

   // start the timer
   sim->last_tics = clock();
   total_time = 0.0;

   // for Nittmann and Stanley's fake surface tension (noise reduced
   //   DLA), zero counters once, here
   if (sim->use_stubborn) zero_counters(top);

   // set all masses (if we loaded files beforehand)
   sum_tipward_mass(top);

   // loop over max number of particles
   for (sim->step=1;sim->step<sim->particle_cnt+1;sim->step++) {

      // Write some stuff
      //if (sim->step%1000 == 0) {
      if (sim->step == next_update) {
         fprintf(stdout,"Step %d, next print %g, radius %g\n",sim->step,sim->next_output_step,sim->overall_rad);
         //zero_all_masses(top);
         //sum_tipward_mass(top);
         //fprintf(stdout,"Step %d, next print %g, %d particles within %g radius\n",sim->step,sim->next_output_step,(int)top->mass,sim->overall_rad);
         fflush(stdout);
         if (sim->step < 2000) {
            next_update += 100;
         } else {
            next_update = (int)(1.05*(float)next_update);
         }
      }
/*
      // write particle count and cell utilization
      if (sim->step%10000 == 0) {
         fprintf(stdout,"  %d particles\n",top->num);
         for (i=0; i<sim->max_levels; i++) sim->cell_count[i] = 0;
         count_cells(sim,top);
         fprintf(stdout,"  ");
         for (i=0; i<sim->max_levels; i++) fprintf(stdout,"%d ",sim->cell_count[i]);
         fprintf(stdout,"cells\n  utilization ");
         for (i=0; i<sim->max_levels; i++) fprintf(stdout,"%g ",sim->cell_count[i]/pow(8,i));
         fprintf(stdout,"\n");
      }
      if (sim->step%100 == 0) {
         fprintf(stdout,"."); fflush(stdout);
      }
*/

      // check current step vs. next_output_step, write output if necessary
      if (sim->step == sim->next_output_step) {
         // stop the timer, save the time
         sim->tics = clock();
         this_step_time = ((FLOAT)sim->tics-(FLOAT)sim->last_tics)/CLOCKS_PER_SEC;
         if (this_step_time < 0.0) this_step_time += 4.2949673e+09/(FLOAT)CLOCKS_PER_SEC;
         total_time += this_step_time;

         // write the output
         fprintf(stdout,"\n  writing output step %d\n",sim->next_output_index);
         write_output(sim,top);

         // write a line to the stats file
         if (sim->write_stats) {
            if (!wrote_stats_header_line) {
               wrote_stats_header_line = TRUE;
               fprintf(sim->statsfp,"# %s %s %s %s %s %s %s\n",
                    "Step",
                    "NumNew",
                    "Radius",
                    "InverseSuccessRatio",
                    "JumpsPerParticle",
                    "SecPerParticle",
                    "TotalRunTime");
            }
            fprintf(sim->statsfp,"%d %d %g %g %g %g %g\n",
                    sim->step,
                    sim->output_step,
                    sim->overall_rad,
                    (FLOAT)sim->total_particles_tried/(FLOAT)sim->output_step,
                    (FLOAT)sim->total_steps_taken/(FLOAT)sim->output_step,
                    (FLOAT)this_step_time/(FLOAT)sim->output_step,
                    total_time);
            fflush(sim->statsfp);
         } else {
            fprintf(stdout,"step time %g, total time %g\n",this_step_time,total_time);
         }

         // reset the counters
         sim->total_particles_tried = 0;
         sim->total_steps_taken = 0;
         sim->next_output_step += sim->output_step;
         sim->next_output_index++;

         // start a new timer
         sim->last_tics = clock();
      }

      // place a new particle
      //for (d=0;d<DIM;d++) start[d]=0.0;

      if (sim->use_bulk_vel) {
         // fprintf(stdout,"Using bulk vel, stick=%g, pen=%g\n",sim->stickiness,sim->penetration);
         retval = diffuse_new_particle(sim,top,sim->bulk_vel);
      } else {
         //fprintf(stdout,"Step %d, no bulk vel, stick=%g, pen=%g\n",sim->step,sim->stickiness,sim->penetration);
         for (d=0;d<DIM;d++) vel[d]=0.0;
         retval = diffuse_new_particle(sim,top,vel);
      }

      /* replace the particles that are not in the correct cell */
      /* it doesn't recalculate the proper number of particles per cell */
      if (retval < 0) {
         top = build_new_tree(sim,top);
         //realloc_box(sim,top);
         //fprintf(stderr,"got here\n");
      }

      /* remove any cells with no particles, and combine cells with
       * too few particles, it also recalculates the proper number of
       * particles per cell */
      // In DLA, we don't need to do this!
      // clean_up_all_cells(sim,top);
   }

   // close the stats file
   if (sim->write_stats) {
      fclose(sim->statsfp);
   }

   // return code 1 - ran out of particles
   fprintf(stdout,"Simulation completed.\n");
   return(1);
}


/*
 *  Run the random-walk method until a hit occurs
 */
int diffuse_new_particle(sim_ptr sim,cell_ptr top,FLOAT* vel) {

   int d,trying_particles,still_walking,is_outside_volume,ix,iy;
   int retval = 0;
   int num_particles_tried = 0;	// how many particles did we try?
   int num_steps_taken = 0;	// how many steps did it take to stick?
   FLOAT loc[DIM],center[DIM],rad,initvel[DIM];
   FLOAT overlap,moveaway;
   particle_ptr newpart;	// pointer to the new particle
   particle_ptr closepart;	// pointer to particle that newpart stuck to

   // overlap is non-dimensional
   overlap = 2.0*(1.0-sim->penetration);

   // this is dimensional
   moveaway = (2.0+sim->penetration)*sim->new_part_rad;

   // create the new part right here
   for (d=0;d<DIM;d++) center[d]=0.0;
   newpart = new_stationary_particle (sim->step,sim->new_part_rad,center);

   trying_particles = TRUE;
   while (trying_particles) {

      num_particles_tried++;

      // find particle start point
      if (sim->particle_source == inf) {
         for (d=0;d<DIM;d++) center[d]=0.0;
         (void) find_start_point (sim,top,center,loc,initvel);
      } else {
         for (d=0;d<DIM;d++) loc[d] = sim->source_vec[d];
         for (d=0;d<DIM;d++) initvel[d] = vel[d];
      }
      //fprintf(stderr,"%d, start at %g %g %g\n",sim->step,loc[0],loc[1],loc[2]); fflush(stdout);

      // loop until particle hits a static particle
      still_walking = TRUE;
      while (still_walking) {

         // check for contact between this part and a static one
         closepart = NULL;
         //rad = find_dist_to_closest2 (top,loc,LARGE,newpart,&closepart);
         rad = find_dist_to_closest (top,loc,LARGE,newpart,&closepart);

         num_steps_taken++;
         //fprintf(stderr,"  nearest is %g\n",rad);
         //if (rad > 9.e+99) exit(1);

         // if it's close enough to stick...
         if (rad < 2.0*sim->new_part_rad) {

            // but if the rootward particle is stubborn
            if (sim->use_stubborn && closepart->counter < sim->stubborn) {
               // then don't stick, but make rootward particle less stubborn
               closepart->counter++;
               // it didn't stick, push it away slightly
               //for (d=0;d<DIM;d++) loc[d] = closepart->x[d] + moveaway*(loc[d]-closepart->x[d])/rad;
               // no, kill the particle and start a new one
               //   say noise-reduced model of Tang 85, Kertesz & Vicsek 85
               still_walking = FALSE;
               rad = 2.0*sim->new_part_rad + 1.;
               //break;

            // or if the rootward particle isn't stubborn anymore...
            // see if stickiness allows particle to stick...
            } else if (sim->use_stickiness) {
               if ((double)(rand())/(double)(RAND_MAX) < sim->stickiness) {
                  // it sticks!
                  still_walking = FALSE;
                  // fprintf(stderr,"point %g %g stuck to particle at %g %g\n",
                  //   loc[0],loc[1],closepart->x[0],closepart->x[1]);

               // not this time, pushback and keep walking
               } else {
                  // fprintf(stderr,"point %g %g didn't stick to particle at %g %g, ",
                  //   loc[0],loc[1],closepart->x[0],closepart->x[1]);
                  // for (d=0;d<DIM;d++) loc[d] += moveaway*(loc[d]-closepart->x[d])/rad;
                  // it didn't stick, push it away slightly
                  for (d=0;d<DIM;d++)
                     loc[d] = closepart->x[d] + moveaway*(loc[d]-closepart->x[d])/rad;
                  // fprintf(stderr,"pushed back to %g %g\n",loc[0],loc[1]);
               }
            // we always stick
            } else {
               still_walking = FALSE;
            }
         }

         // if we strayed too far away, quit anyways
         // this started at 50, then it was 10...it keeps getting lower
         if (rad > 3.*sim->overall_rad && rad > 10.) still_walking = FALSE;

         // if we were stubborn, but successful, reset the closepart's counter
         if (sim->use_stubborn && !still_walking) {
            // we need to do this to kill the particle and start a new one
         }

         // break out if we should stop now
         if (!still_walking) break;

         // move the particle to a sphere of that radius
         for (d=0;d<DIM;d++) center[d]=loc[d];
         // fprintf(stderr,"  center at %g %g %g\n",center[0],center[1],center[2]);
         pick_point_on_sphere(loc,center,rad-overlap*sim->new_part_rad,sim->use_bulk_vel,vel);
         //fprintf(stderr,"  point at %g %g\n",loc[0],loc[1]);

      }

      // if we quit because we contacted another particle...
      if (rad < 2.0*sim->new_part_rad) {

         // apply the final touches

         // if grip is on, slide it closer
         if (sim->use_grip) {
            // particle "grips" and slides toward the contacted particle
            for (d=0;d<DIM;d++) loc[d] -= sim->grip*(loc[d]-closepart->x[d]);
         }

         // if chiral is on, rotate it around the parent particle
         if (sim->use_chiral) {
            // particle rotates around the contacted particle in x-y plane only
            // first, find the vector from the closepart to the particle
            //fprintf(stderr,"loc starts %g %g\n",loc[0],loc[1]);
            correct_for_chiral(loc,closepart,sim->chiral_angle,sim->chiral_power);
            //fprintf(stderr,"  and ends %g %g\n",loc[0],loc[1]);
         }

         // test against a fixed volume
         is_outside_volume = FALSE;
         if (sim->trim_method == cube) {
            // keep inside a cube
            for (d=0;d<DIM;d++) {
               if (fabs(loc[d]) > sim->trim_rad) {
                  is_outside_volume = TRUE;
               }
            }

         } else if (sim->trim_method == plane) {
            // make a flat plane
            if (fabs(loc[DIM-1]) > sim->trim_rad) {
               is_outside_volume = TRUE;
            }

         } else if (sim->trim_method == sphere) {
            // keep inside a sphere
            if (vec_length_sq(loc) > pow(sim->trim_rad,2)) {
               is_outside_volume = TRUE;
            }

         } else if (sim->trim_method == png) {
            // keep inside a mask (i.e not outside the volume)
            is_outside_volume = TRUE;
            ix = (int)(0.5*(loc[0]+1.0)*sim->trim_nx);
            iy = (int)(0.5*(loc[1]+1.0)*sim->trim_nx);
            if (ix>-1 && iy>-1 && ix<sim->trim_nx && iy<sim->trim_ny) {
               //if (sim->trim_mask[ix][iy] > 0.5 && fabs(loc[2]) < 0.07) {
#if DIM>2
               if (sim->trim_mask[ix][iy] > 0.5 && loc[2] > 0.0) {
#else
               if (sim->trim_mask[ix][iy] > 0.5) {
#endif
                  // if outside mask and above zero, keep it
                  is_outside_volume = FALSE;
               }
#if DIM>2
               if (loc[2] > 0.29) {
                  // if over the height of the bldngs, keep it
                  is_outside_volume = FALSE;
               }
            } else {
               // if outside of the mask area entirely, but above ground, keep it
               if (loc[2] > 0.0) {
                  is_outside_volume = FALSE;
               }
#endif
            }
         }

         // sometimes allow the segment even if we shouldn't
         //if (is_outside_volume && rand() < RAND_MAX/40) {
         //   is_outside_volume = FALSE;
         //}

         // finally, check the above criteria and quit before making the particle
         if (is_outside_volume) {
            trying_particles = FALSE;
            break;
         }

         // only now locate the particle
         for (d=0;d<DIM;d++) newpart->x[d] = loc[d];

         // put the particle in the tree, and return a negative number
         //   if the new part is out of bounds
         //fprintf(stdout,"in diffuse_new_particle\n");
         retval = add_particle_to_cell(sim,newpart,top);

         // define the rootward particle
         newpart->root = closepart;

         // add to list of root's tipward particles
         newpart->next_tip = newpart->root->tip_head;
         newpart->root->tip_head = newpart;

         // set this particles mass, and add 1.0 to all rootward masses
         add_to_rootward_masses (newpart);

         // if junction flow is on, adjust the position of the contacted particle!
         if (sim->use_junction_flow) {
            static int dbg_jctn = FALSE;
            if (dbg_jctn) fprintf(stderr,"\nmarching toward root\n");
            // march down the branches to the root, adjusting every node as you go
            particle_ptr nextroot = newpart->root;
            while (nextroot) {
               // the rootward node of this one must exist, also!
               if (nextroot->root) {

                  // find outflow and sum of all inflows
                  FLOAT in[DIM];
                  for (d=0; d<DIM; d++) in[d] = 0.0;
                  FLOAT inmass = 0.0;
                  particle_ptr tip = nextroot->tip_head;
                  while (tip) {
                     for (d=0; d<DIM; d++) in[d] += tip->mass * (tip->x[d] - nextroot->x[d]);
                     inmass += tip->mass;
                     tip = tip->next_tip;
                  }
                  // should we normalize the vector, or just divide by weight?
                  if (TRUE) {
                     inmass = 0.0;
                     for (d=0; d<DIM; d++) inmass += in[d]*in[d];
                     inmass = sqrt(inmass) / sim->new_part_rad;
                  }
                  for (d=0; d<DIM; d++) in[d] /= inmass;
                  if (dbg_jctn) fprintf(stderr,"in mass %g and vector %g %g %g\n",inmass,in[0],in[1],in[2]);

                  // subtract one from rootward mass so that masses balance on each side
                  FLOAT out[DIM];
                  FLOAT outmass = nextroot->mass - 1.0;
                  for (d=0; d<DIM; d++) out[d] = outmass * (nextroot->root->x[d] - nextroot->x[d]);
                  if (TRUE) {
                     outmass = 0.0;
                     for (d=0; d<DIM; d++) outmass += out[d]*out[d];
                     outmass = sqrt(outmass) / sim->new_part_rad;
                  }
                  for (d=0; d<DIM; d++) out[d] /= outmass;
                  if (dbg_jctn) fprintf(stderr,"out mass %g and vector %g %g %g\n",outmass,out[0],out[1],out[2]);

                  // finally, shift the point a little
                  for (d=0; d<DIM; d++) nextroot->x[d] += sim->junction_coeff*(in[d]+out[d])
                                                         / (1.0 + sim->junction_coeff);;
               }
               nextroot = nextroot->root;
            }
         }

         //fprintf(stderr,"  planting particle at %g %g, rad is %g, npr is %g\n\n",loc[0],loc[1],rad,sim->new_part_rad);
         trying_particles = FALSE;

         // and if we were successful, even with stubbornness, reset counter
         if (sim->use_stubborn) {
            closepart->counter = 0;
         }

      } else {
         // or, because we were too far away
         // fprintf(stderr,"  particle at %g %g %g, rad is %g, is too far\n",loc[0],loc[1],loc[2],rad);
      }
   }

   // set statistics
   sim->total_particles_tried += num_particles_tried;
   sim->total_steps_taken += num_steps_taken;

   return(retval);
}


/*
 * Determine an appropriate particle start point (for now, on a sphere
 * the size of the whole assembly)
 *
 * Inputs are sim, top, center
 * Outputs are startloc, startvel
 */
int find_start_point (sim_ptr sim, cell_ptr top, FLOAT *center,
                      FLOAT *startloc, FLOAT *startvel) {

   int d;
   FLOAT rad;

   // take one step to start
   rad = find_dist_to_farthest(top,center,0.0) + sim->new_part_rad;
   // rad = find_dist_to_closest(top,startloc,LARGE,newpart,&closepart) 
   //       -2.0*sim->new_part_rad;

   // set overall size, in case we need it
   sim->overall_rad = rad;

   // fprintf(stderr,"  maxrad is %g\n",rad);
   // fprintf(stderr,"starting at %g %g\n",startloc[0],startloc[1]);
   //for (d=0;d<DIM;d++) center[d]=startloc[d];

   // use negative velocity for initial particle placement?
   // for (j=0;j<DIM;j++) startvel[j]=-1.0*vel[j];
   for (d=0;d<DIM;d++) startvel[d]=0.0;

   (void) pick_point_on_sphere(startloc, center, rad+sim->new_part_rad,
                               nobulk, startvel);
   // fprintf(stderr,"  first move to %g %g\n",startloc[0],startloc[1]);

   // testing Lichtenberg figures
   //   choose a random point inside of a fixed volume!
   //for (d=0;d<DIM;d++) startloc[d] = 
   //startloc[0] = -0.5 + rand()/(RAND_MAX+1.0);
   //startloc[1] = rand()/(RAND_MAX+1.0);
   //startloc[2] = -0.05 + 0.1*rand()/(RAND_MAX+1.0);
   //fprintf(stderr,"starting at %g %g %g\n",startloc[0],startloc[1],startloc[2]);

   return(0);
}


/*
 * find a random point on a given sphere shell
 *
 * inputs
 *   center[DIM]  starting location
 *
 * outputs
 *   loc[DIM]  new point
 */
void pick_point_on_sphere (FLOAT* loc, FLOAT* center, FLOAT rad,
                           BULKVEL use_vel, FLOAT* vel) {

   // int rotex = TRUE;
   int d,td;
   FLOAT len = 2.;
   FLOAT dt = rad*rad;
   FLOAT a[5],sdev[DIM+1],r1,r2;
   FLOAT thisvel[DIM];
   for (d=0;d<DIM;d++) thisvel[d] = 0.0;

   // first, account for bulk velocity
   if (use_vel == rotex) {
      // This was the old rotex growth thing, but all it does is influence
      //   the bulk velocity, this should be incorported differently, maybe
      //   as a list of potential flow primitives (freestream, vortex, source)
      // point velocity around counterclockwise
      // dt = vec_length_sq(center)*vec_length(center);
      float dist_from_origin = vec_length(center);
      #if DIM==2
         // 2d flow, use vel[0] as vortex strength
         thisvel[0] = 10. * (vel[0] * -center[1] - center[0]);// / dist_from_origin;
         thisvel[1] = 10. * (vel[0] *  center[0] - center[1]);// / dist_from_origin;
      #endif
      #if DIM==3
         // 3d flow, use vel[0:2] as vectorial vortex strength
         thisvel[0] = (vel[1]*center[2] - vel[2]*center[1]) / dist_from_origin;
         thisvel[1] = (vel[2]*center[0] - vel[0]*center[2]) / dist_from_origin;
         thisvel[2] = (vel[0]*center[1] - vel[1]*center[0]) / dist_from_origin;
      #endif
      // don't take steps so long that particles don't curve
      dt = 0.25*dist_from_origin;///vel[0];
      if (rad > dt) rad = dt;

   } else if (use_vel == straight) {
      // fprintf(stderr,"bulk vel is %g %g\n",vel[0],vel[1]);
      for (d=0;d<DIM;d++) thisvel[d] = vel[d];
   }

   if (use_vel != nobulk) {
      // first, find the number of standard deviations in each direction
      for (td=0;td<(DIM+1)/2;td++) {
         gaussian_rand(&r1,&r2);
         sdev[2*td] = r1;
         sdev[2*td+1] = r2;
      }
      // fprintf(stderr,"std devs are %g %g\n",sdev[0],sdev[1]);
      // then, set up the quartic equation
      a[4] = vec_length_sq(thisvel);
      a[3] = 2.0*vec_dot(vel,sdev);
      a[2] = vec_length_sq(sdev);
      a[1] = 0.0;
      a[0] = -1.0*rad*rad;
      // solve it (for sqrt(dt))
      dt = solve_quartic(a);
      // and account for the random motion
      for (d=0;d<DIM;d++) loc[d] = center[d]+dt*sdev[d];
      // determine the real dt
      dt *= dt;
      // move the particle due to convection
      for (d=0;d<DIM;d++) loc[d] += dt*thisvel[d];
      // now we're done!
      return;
#ifdef OLD
      // determine speed of convection
      len = vec_length(vel);
      // solve for dt such that the particle never leaves a sphere of radius rad
      // WARNING: this is still not quite accurate; how would one determine how
      // long a *specific* particle took to diffuse a specific distance?
      rad = (sqrt(1.0+4.0*rad*len)-1.0)/(2.0*len);
      dt = pow(rad,2);
      // move the particle due to convection
      for (d=0;d<DIM;d++) center[d] += dt*vel[d];
#endif
   }

   // then move it due to random motion
#if DIM==2
#ifdef USE_SINE_SPHERE
   // second method uses no iteration at all, because sin/cos is accelerated
   len = 6.2831853071795864 * rand()/(RAND_MAX+1.0);
   loc[0] = cos(len);
   loc[1] = sin(len);
   for (d=0;d<DIM;d++) loc[d] = center[d]+rad*loc[d];
#else
   // new method needs no sqrt
   FLOAT temp[2];
   len = 2.;
   while (len > 1.) {
      for (d=0;d<2;d++) temp[d] = 2.*rand()/(RAND_MAX+1.0) - 1.;
      // len = vec_length(loc);
      len = temp[0]*temp[0]+temp[1]*temp[1];
   }
   // for (d=0;d<DIM;d++) loc[d] = center[d]+rad*loc[d]/len;
   loc[0] = (temp[0]*temp[0]-temp[1]*temp[1])/len;
   loc[1] = 2.*temp[0]*temp[1]/len;
   // printf("pt on sphere %g %g\n",loc[0],loc[1]);
   for (d=0;d<DIM;d++) loc[d] = center[d]+rad*loc[d];
   // printf("scaled pt %g %g\n",loc[0],loc[1]);
#endif
#elif DIM==3
   // new method needs no sqrt
   FLOAT temp[2];
   len = 2.;
   while (len > 1.) {
      // for (d=0;d<DIM;d++) loc[d] = 2.*rand()/(RAND_MAX+1.0) - 1.;
      for (d=0;d<2;d++) temp[d] = 2.*rand()/(RAND_MAX+1.0) - 1.;
      // len = vec_length(loc);
      len = temp[0]*temp[0]+temp[1]*temp[1];
   }
   // printf("pt on circle %g %g\n",temp[0],temp[1]);
   // for (d=0;d<DIM;d++) loc[d] = center[d]+rad*loc[d]/len;
   // loc[0] = (temp[0]*temp[0]-temp[1]*temp[1])/len;
   // loc[1] = 2.*temp[0]*temp[1]/len;
   loc[2] = 2.*len -1.;
   len = 2.*sqrt(1.-len);
   loc[0] = temp[0]*len;
   loc[1] = temp[1]*len;
   // printf("pt on sphere %g %g %g\n",loc[0],loc[1],loc[2]);
   for (d=0;d<DIM;d++) loc[d] = center[d]+rad*loc[d];
   // printf("new pt %g %g %g\n",loc[0],loc[1],loc[2]);
#else
   // else use the standard method
   len = 2.;
   while (len > 1.) {
      for (d=0;d<DIM;d++) loc[d] = 2.*rand()/(RAND_MAX+1.0) - 1.;
      len = vec_length_sq(loc);
   }
   len = sqrt(len);
   for (d=0;d<DIM;d++) loc[d] = center[d]+rad*loc[d]/len;
#endif

   return;
}


/*
 *  Find the distance to the closest particle
 *  based on is_this_space_open2
 *
 *  cell is cell to search within
 *  loc is point to search against
 *  nearest is current value for nearest
 */
FLOAT find_dist_to_closest(cell_ptr cell,FLOAT* loc,FLOAT neardist,
                           particle_ptr newpart,particle_ptr *closest){

   int i,d,index;
   //int s[DIM];
   // int need_to_do_sqrt = FALSE;
   FLOAT temp[DIM],dist,nsq;
   particle_ptr curr;

   nsq = neardist*neardist;

   // are all corners of this cell 
   for (d=0;d<DIM;d++) {
      if (loc[d] < cell->min[d]-neardist) return(neardist);
      if (loc[d] > cell->max[d]+neardist) return(neardist);
   }

   if (cell->has_subcells) {
      // fprintf(stdout,"    cell has %d subcells\n",NCHILD);
      // must check all of the subcells, but check closest first
      index = 0;
      for (d=0;d<DIM;d++)
        if (loc[d] > cell->mid[d])
          index += (1<<d);
      for (i=index;i<index+NCHILD;i++) {
      // the slow way does in 27 sec what the fast way does in 19
      // for (i=0;i<NCHILD;i++) {
         dist = find_dist_to_closest(cell->s[i%NCHILD],
                   loc,neardist,newpart,closest);
         if (dist < neardist) neardist = dist;
         // fprintf(stderr,"    closest is %g\n",neardist);
      }
   } else {
      // fprintf(stdout,"    cell has %d particles\n",cell->num);
      // check vs. all of the particles in this cell
      // fprintf(stderr,"    checking vs %d particles in box at level %d\n",cell->num,cell->level);
      curr = cell->first;
      while (curr) {

         // it DOES NOT PAY to test this sum at every addition!
         for (d=0;d<DIM;d++) temp[d] = loc[d]-curr->x[d];
         dist = vec_length_sq(temp);

         // fprintf(stderr,"    dist to this particle is %g\n",dist);
         if (dist < nsq) {
            nsq = dist;
            // it DOES NOT PAY to postpone the sqrt until later!
            neardist = sqrt(dist);
            // need_to_do_sqrt = TRUE;
            (*closest) = curr;
            newpart->root = curr;
         }
         curr = curr->next;
      }
   }

   // if (need_to_do_sqrt) neardist = sqrt(dist);

   return(neardist);
}


/*
 *  Find the distance to the closest particle - take 2
 *
 *  Using speedup ideas from Andy Lomas, whereby only the closest
 *  is traversed first. On the return steps, any other cells can
 *  be checked
 *
 *  cell is cell to search within
 *  loc is point to search against
 *  nearest is current value for nearest
 */
FLOAT find_dist_to_closest2 (cell_ptr cell,FLOAT* loc,FLOAT neardist,
                             particle_ptr newpart,particle_ptr *closest){

   int i,d;
   //int s[DIM];
   FLOAT temp[DIM],dist,nsq;
   particle_ptr curr;
   unsigned char index,adj_index;//,mask;
   //mask = 0;

   // are all corners of this cell 
   for (d=0;d<DIM;d++) {
      if (loc[d] < cell->min[d]-neardist) return(neardist);
      if (loc[d] > cell->max[d]+neardist) return(neardist);
   }

   // first, find the closest subcell and check it first
   if (cell->has_subcells) {
      // fprintf(stdout,"    cell has %d subcells\n",NCHILD);
      // must check all of the subcells, but check closest first
      // compute cell index using bitwise AND
      index = 0;
      for (d=0;d<DIM;d++)
        if (loc[d] > cell->mid[d])
          index += (1<<d);

      // first, check the closest subcell
      dist = find_dist_to_closest2 (cell->s[index],
                loc,neardist,newpart,closest);
      if (dist < neardist) neardist = dist;

/*
      // now, make sure no other subcells could possibly hold a closer one
      //for (d=DIM-1;d>-1;d--) {
      for (d=0;d<DIM;d++) {
        if (abs(loc[d]-cell->mid[d]) < neardist) {
          // HOW THE HELL DO I DO THIS? THE SPEEDUP IS PROMISING!
          mask = 1<<d;
          adj_index = index ^ mask;
          dist = find_dist_to_closest2 (cell->s[adj_index],
                    loc,neardist,newpart,closest);
          if (dist < neardist) neardist = dist;
        }
      }
*/

      // dumb way: just check the remaining 2^DIM - 1 cells
      for (adj_index=0; adj_index<NCHILD; adj_index++) {
        if (adj_index!=index) {
          dist = find_dist_to_closest2 (cell->s[adj_index],
                    loc,neardist,newpart,closest);
          if (dist < neardist) neardist = dist;
        }
      }

      // really dumb way: don't bother to check the other cells,
      //   just return

   // otherwise, check vs. all of the particles in this cell
   } else {
      // fprintf(stdout,"    cell has %d particles\n",cell->num);
      // fprintf(stderr,"    checking vs %d particles in box at level %d\n",cell->num,cell->level);
      nsq = neardist*neardist;
      curr = cell->first;
      while (curr) {

         // it DOES NOT PAY to test this sum at every addition!
         for (d=0;d<DIM;d++) temp[d] = loc[d]-curr->x[d];
         //dist = vec_length_sq(temp);
         dist = 0.0;
         for (i=0;i<DIM;i++) dist += temp[i]*temp[i];

         // fprintf(stderr,"    dist to this particle is %g\n",dist);
         if (dist < nsq) {
            nsq = dist;
            // it DOES NOT PAY to postpone the sqrt until later!
            neardist = sqrt(dist);
            (*closest) = curr;
            newpart->root = curr;
         }
         curr = curr->next;
      }
   }

   return(neardist);
}


/*
 *  Find the distance to the farthest particle
 *  based on is_this_space_open2
 */
// FLOAT find_dist_to_farthest(cell_ptr cell,FLOAT loc[DIM],FLOAT farthest){
FLOAT find_dist_to_farthest(cell_ptr cell,FLOAT* loc,FLOAT farthest){

   int i,d;
   FLOAT temp[DIM],dsq,fsq,dmin,dmax,dist;
   particle_ptr curr;

   // if (cell->level == 0) farthest = 0.;
   // farthest = 0.;
   fsq = farthest*farthest;

   // must see if whole cell is inside of sphere radius "farthest"
   // so, find the distance to the farthest corner on the cube
   dsq = 0.;
   for (d=0;d<DIM;d++) {
      dmin = fabs(loc[d]-cell->min[d]);
      dmax = fabs(loc[d]-cell->max[d]);
      if (dmin > dmax) dsq += dmin*dmin;
      else dsq += dmax*dmax;
   }
   if (dsq < fsq) return(farthest);

   if (cell->has_subcells) {
      // must check all of the subcells
      for (i=0;i<NCHILD;i++) {
         dist = find_dist_to_farthest(cell->s[i],loc,farthest);
         if (dist > farthest) farthest = dist;
         // fprintf(stderr,"    farthest is %g\n",farthest);
      }
   } else {
      // check vs. all of the particles in this cell
      // fprintf(stderr,"    checking vs %d particles in box at level %d\n",cell->num,cell->level);
      curr = cell->first;
      while (curr) {
         for (d=0;d<DIM;d++) temp[d] = loc[d]-curr->x[d];
         //dist = vec_length_sq(temp);
         dist = 0.0;
         for (i=0;i<DIM;i++) dist += temp[i]*temp[i];
         // fprintf(stderr,"    dist to this particle is %g\n",dist);
         if (dist > fsq) {
            fsq = dist;
            farthest = sqrt(dist);
         }
         curr = curr->next;
      }
   }

   return(farthest);
}


/*
 * Account for chiral growth based on weight of rootward particle
 *
 * overwrites location with the new location
 */
int correct_for_chiral (FLOAT* loc, particle_ptr closest, FLOAT angle, FLOAT power) {

   int i;
   particle_ptr oldest;
   FLOAT rvec[2];
   FLOAT weight, omw;

/*
   // take 1: pure chiral, only chiral if tip growth
   if (closest->mass < 1.1) {
      // move new particle location to follow arc exactly

      // find the list of particles, oldest -> closest -> this (loc)
      if (closest->root) {
         oldest = closest->root;

         //fprintf(stderr,"old %g %g\n",

         // find the arc from oldest to closest
         for (i=0; i<2; i++) rvec[i] = closest->x[i] - oldest->x[i];

         // rotate that vector and add it to closest to get loc
         loc[0] = closest->x[0] + rvec[0]*cos(angle) - rvec[1]*sin(angle);
         loc[1] = closest->x[1] + rvec[0]*sin(angle) + rvec[1]*cos(angle);

      } else {
         // if there's no root, do nothing
      }

   } else {
      // do nothing - keep impact point as is
   }
*/

   // take 2, weighted chiral
   if (closest->root) {

      oldest = closest->root;

      // find the arc from oldest to closest
      for (i=0; i<2; i++) rvec[i] = closest->x[i] - oldest->x[i];

      // now, find the weight of the true arc
      weight = pow(1./closest->mass, power);
      omw = 1. - weight;

      // rotate that vector and add it to closest to get loc
      loc[0] = omw*loc[0] + weight*(closest->x[0] + rvec[0]*cos(angle) - rvec[1]*sin(angle));
      loc[1] = omw*loc[1] + weight*(closest->x[1] + rvec[0]*sin(angle) + rvec[1]*cos(angle));

   }

   return(0);
}


/*
 * Starting with each particle, march toward the center
 * and add 1.0=mass to every node along the way
 */
void sum_tipward_mass(cell_ptr cell) {

   int i;
   particle_ptr curr,march;

   if (cell->has_subcells) {
      // must check all of the subcells
      for (i=0;i<NCHILD;i++)
         sum_tipward_mass(cell->s[i]);
   } else {
      // check vs. all of the particles in this cell
      curr = cell->first;
      while (curr) {
         curr->mass += 1.0;
         march = curr->root;
         while (march) {
            march->mass += 1.0;
            march = march->root;
         }
         curr = curr->next;
      }
   }

   return;
}


/*
 * Add this particle's mass to all rootward particles
 */
void add_to_rootward_masses (particle_ptr this) {

   particle_ptr march;

   this->mass += 1.0;
   march = this->root;
   while (march) {
      march->mass += 1.0;
      // double-check that this is the real root and that the sub is working
      //fprintf(stderr,"Root weights %g now\n",march->mass);
      march = march->root;
   }


   return;
}


/*
 * zero all masses
 */
void zero_all_masses(cell_ptr cell) {

   int i;
   particle_ptr curr;

   if (cell->has_subcells) {
      // must check all of the subcells
      for (i=0;i<NCHILD;i++)
         zero_all_masses(cell->s[i]);
   } else {
      // check vs. all of the particles in this cell
      curr = cell->first;
      while (curr) {
         curr->mass = 0.0;
         curr = curr->next;
      }
   }

   return;
}


/*
 * zero all counters - this is not normally needed
 */
void zero_counters(cell_ptr cell) {

   int i;
   particle_ptr curr;

   if (cell->has_subcells) {
      // must check all of the subcells
      for (i=0;i<NCHILD;i++)
         zero_counters(cell->s[i]);
   } else {
      // check vs. all of the particles in this cell
      curr = cell->first;
      while (curr) {
         curr->counter = 0;
         curr = curr->next;
      }
   }

   return;
}


/*
 * increment all counters - this is not normally needed
 */
void increment_counters(cell_ptr cell) {

   int i;
   particle_ptr curr;

   if (cell->has_subcells) {
      // must check all of the subcells
      for (i=0;i<NCHILD;i++)
         increment_counters(cell->s[i]);
   } else {
      // check vs. all of the particles in this cell
      curr = cell->first;
      while (curr) {
         curr->counter++;
         curr = curr->next;
      }
   }

   return;
}


/*
 * find the 2-norm of a d-length vector
 */
FLOAT vec_dot(FLOAT* x,FLOAT* y) {

   int i;
   FLOAT sum = 0.0;

   for (i=0;i<DIM;i++) sum += x[i]*y[i];

   return(sum);
}


/*
 * find the 2-norm of a d-length vector
 */
FLOAT vec_length(FLOAT* x) {

   int i;
   FLOAT sum = 0.0;

   for (i=0;i<DIM;i++) sum += x[i]*x[i];
   sum = sqrt(sum);

   return(sum);
}


/*
 * find the 2-norm squared of a d-length vector
 */
FLOAT vec_length_sq(FLOAT* x) {

   int i;
   FLOAT sum = 0.0;

   for (i=0;i<DIM;i++) sum += x[i]*x[i];

   return(sum);
}


/*
 * return two random numbers from the normal distribution, mean=0, sigma=1
 *
 * uses the Box-Mueller transformation
 */
void gaussian_rand(FLOAT* r1,FLOAT* r2) {

   double t1,t2;

   t1 = (double)(rand())/(double)(RAND_MAX);
   t2 = (double)(rand())/(double)(RAND_MAX);

   *r1 = (FLOAT)(sqrt(-2.*log(t1))*cos(2.*M_PI*t2));
   *r2 = (FLOAT)(sqrt(-2.*log(t1))*sin(2.*M_PI*t2));

   return;
}

inline void gaussian_rand2(FLOAT* r1,FLOAT* r2) {

   const double t1 = (double)(rand())/(double)(RAND_MAX);
   const double t2 = 2.0*M_PI*(double)(rand())/(double)(RAND_MAX);

   const double sntlt = sqrt(-2.0*log(t1));
   double sinval, cosval;
   (void)sincos(t2, &sinval, &cosval);

   *r1 = (FLOAT)(sntlt*cosval);
   *r2 = (FLOAT)(sntlt*sinval);

   return;
}

inline void gaussian_rand3(FLOAT* r1,FLOAT* r2) {

   const float t1 = (float)(rand())/(float)(RAND_MAX);
   const float t2 = 2.f*M_PI*(float)(rand())/(float)(RAND_MAX);

   const float sntlt = sqrtf(-2.f*logf(t1));
   float sinval, cosval;
   (void)sincosf(t2, &sinval, &cosval);

   *r1 = (FLOAT)(sntlt*cosval);
   *r2 = (FLOAT)(sntlt*sinval);

   return;
}


/*
 * solve the quartic equation numerically for the real positive root;
 * equation is a[4] x^4 + a[3] x^3 + a[2] x^2 + a[1] x + a[0] = 0
 */
FLOAT solve_quartic(FLOAT* a) {

   int i,it,maxit = 100;
   FLOAT xmult;
   FLOAT bottomx;//,bottomval;
   FLOAT topx,topval;
   FLOAT newx,newval;

   // evaluate the function at the bottom bound (always negative in our case)
   bottomx = 0.0;
   //bottomval = a[0];
   // fprintf(stderr,"initial bottom bound f(%g)=%g\n",bottomx,bottomval);

   // evaluate the function at the top bound
   topx = fabs((a[3]+a[0])/2.0*a[4]);
   topval = -1.0;
   while (topval < 0.0) {
      topx *= 2.0;
      topval=a[0];
      xmult=1.;
      for (i=1;i<5;i++) {
         xmult *= topx;
         topval += a[i]*xmult;
      }
      // fprintf(stderr,"initial top bound f(%g)=%g\n",topx,topval);
      // we hope the top bound is positive
   }

   newval = 1.0;
   it = 0;
   while (fabs(newval) > 1.0e-6 && it < maxit) {
      // take a new guess for x
      newx = (topx+bottomx)/2.0;
      // evaluate the expression
      newval=a[0];
      xmult=1.;
      for (i=1;i<5;i++) {
         xmult *= newx;
         newval += a[i]*xmult;
      }
      // fprintf(stderr,"   interim value f(%g)=%g\n",newx,newval);
      if (newval < 0.0) {
         // replace the bottom value if negative
         bottomx = newx;
      } else {
         // replace the top value if positive
         topx = newx;
      }
      it++;
   }

   return newx;
}

