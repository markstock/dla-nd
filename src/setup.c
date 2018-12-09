/*************************************************************
 *
 *  setup.c - subroutines for setting up simulations
 *
 *  Copyright (C) 2000-8,14  Mark J. Stock, mstock@umich.edu
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

#include <ctype.h>


/*
 *  Initialize the system - make static and free element structures, etc.
 */
int set_defaults (sim_ptr sim, cell_ptr top) {

   int i,d;

   // -------------------------------------------------------------------------

   /* set default values for variables */
   // strcpy(sim->out_fn_root,"output_");
   strcpy(sim->out_fn_root,"");
   sim->write_pgm = FALSE;
   sim->write_png = TRUE;
   sim->write_dot = FALSE;
   sim->write_rad = FALSE;
   sim->write_obj = FALSE;
   sim->write_seg = FALSE;
   sim->write_part = FALSE;
   sim->image_depth = 8;
   sim->out_img_size = 512;

   sim->write_stats = TRUE;
   strcpy(sim->stats_file,"stats.dat");

   // particle masks
   sim->trim_method = notrim;
   sim->trim_rad = 0.5;
   strcpy(sim->trim_file,"trim.png");

   // particle sources
   sim->particle_source = inf;
   for (d=0;d<DIM;d++) sim->source_vec[d] = 0.0;
   sim->use_source_mask = FALSE;
   strcpy(sim->source_file,"source.png");

   sim->seed = 1;
   sim->step = 0;
   sim->output_step = 100;
   sim->max_levels = 50;
   sim->max_parts_in_cell = 40;

   sim->next_output_step = 0;
   sim->next_output_index = 1;

   sim->penetration = 0.1;

   sim->use_stickiness = FALSE;
   sim->stickiness = 1.0;

   sim->use_grip = FALSE;
   sim->grip = 0.0;

   sim->use_stubborn = FALSE;
   sim->stubborn = 0;

   sim->use_bulk_vel = nobulk;
   for (d=0;d<DIM;d++) {
      sim->bulk_vel[d] = 0.0;
   }

   sim->use_chiral = FALSE;
   sim->chiral_diameter = 500.;
   sim->chiral_angle = 2.*M_PI/sim->chiral_diameter;
   sim->chiral_power = 2.;

   sim->use_junction_flow = FALSE;
   sim->junction_coeff = 0.0;

   for (d=0;d<DIM;d++) {
      sim->bdry[d][0] = WALL;
      sim->bdry[d][1] = WALL;
   }

   sim->particle_cnt = 1000;
   sim->new_part_rad = 0.005;
   sim->new_part_mass = 1.0;

   sim->num_blocks = 0;
   sim->num_spheres = 0;
   sim->num_indiv_parts = 0;
   sim->num_read_part = 0;
   sim->num_read_stat = 0;

   sim->write_2d_dens = FALSE;
   sim->write_3d_dens = FALSE;
   sim->use_dens_zone = FALSE;

   sim->out = &sim->out_field;
   sim->ff2 = &sim->density_field2;
   sim->ff3 = &sim->density_field3;

   /* assign default values for the top-level cell */
   top->level = 0;
   top->has_subcells = FALSE;
   top->first = NULL;
   for (i=0;i<NCHILD;i++)
      top->s[i] = NULL;
   for (d=0;d<DIM;d++) {
      top->min[d] = 0.0;
      top->max[d] = 1.0;
   }
   top->num = 0;
   top->mass = 0.0;

   /* return 0, everything OK */
   return(0);
}


/*
 * Parse the command-line arguments
 */
int parse_args (int argc,char **argv,sim_ptr sim,cell_ptr top) {

   int i,j,d,ret_val;

   // copy executable name
   (void) strcpy(sim->exectuable_fn,argv[0]);

   if (argc < 2) {
      // if no command-line arg, then input must be stdin
      sim->use_input_file = FALSE;
   }

   for (i=1; i<argc; i++) {

      // if the next arg begins with a dash, it is not a filename
      //if (strncmp(argv[i], "-",1) == 0) {
      if (argv[i][0] == '-') {

         // if command-line arg is "help", write usage
         if (strncmp(argv[i], "--help", 3) == 0) {
            (void) Usage(sim->exectuable_fn,0);
         } else if (strncmp(argv[i], "-end", 4) == 0) {
            sim->particle_cnt = atoi(argv[++i]);
         } else if (strncmp(argv[i], "-every", 6) == 0) {
            sim->output_step = atoi(argv[++i]);
         } else if (strncmp(argv[i], "-tm", 3) == 0) {
            sim->trim_method = png;
            if (i == argc-1) Usage(sim->exectuable_fn,1);
            strcpy(sim->trim_file,argv[++i]);
         } else if (strncmp(argv[i], "-cube", 3) == 0) {
            sim->trim_method = cube;
            if (i == argc-1) Usage(sim->exectuable_fn,1);
            sim->trim_rad = atof(argv[++i]);
         } else if (strncmp(argv[i], "-sphere", 4) == 0) {
            sim->trim_method = sphere;
            if (i == argc-1) Usage(sim->exectuable_fn,1);
            sim->trim_rad = atof(argv[++i]);
         } else if (strncmp(argv[i], "-sp", 3) == 0) {
            sim->particle_source = point;
            if (i == argc-DIM) Usage(sim->exectuable_fn,1);
            for (d=0;d<DIM;d++) sim->source_vec[d] = atof(argv[++i]);
         } else if (strncmp(argv[i], "-pr", 3) == 0) {
            sim->new_part_rad = atof(argv[++i]);
         } else if (strncmp(argv[i], "-pen", 4) == 0) {
            sim->penetration = atof(argv[++i]);
         } else if (strncmp(argv[i], "-vel", 4) == 0) {
            sim->use_bulk_vel = straight;
            if (i == argc-DIM) Usage(sim->exectuable_fn,1);
            for (d=0;d<DIM;d++) sim->bulk_vel[d] = atof(argv[++i]);
         } else if (strncmp(argv[i], "-stick", 6) == 0) {
            sim->use_stickiness = TRUE;
            if (i == argc-1) Usage(sim->exectuable_fn,1);
            sim->stickiness = atof(argv[++i]);
         } else if (strncmp(argv[i], "-grip", 5) == 0) {
            sim->use_grip = TRUE;
            if (i == argc-1) Usage(sim->exectuable_fn,1);
            sim->grip = atof(argv[++i]);
         } else if (strncmp(argv[i], "-stub", 5) == 0) {
            sim->use_stubborn = TRUE;
            if (i == argc-1) Usage(sim->exectuable_fn,1);
            sim->stubborn = atoi(argv[++i]);
         } else if (strncmp(argv[i], "-rotex", 4) == 0) {
            sim->use_bulk_vel = rotex;
            if (i == argc-DIM) Usage(sim->exectuable_fn,1);
            for (d=0;d<DIM;d++) sim->bulk_vel[d] = atof(argv[++i]);
         } else if (strncmp(argv[i], "-junction", 5) == 0) {
            sim->use_junction_flow = TRUE;
            if (i == argc-1) Usage(sim->exectuable_fn,1);
            sim->junction_coeff = atof(argv[++i]);
         } else if (strncmp(argv[i], "-rad", 4) == 0) {
            sim->write_rad = TRUE;
         } else if (strncmp(argv[i], "-obj", 4) == 0) {
            sim->write_obj = TRUE;
         } else if (strncmp(argv[i], "-seg", 4) == 0) {
            sim->write_seg = TRUE;
         } else if (strncmp(argv[i], "-res", 4) == 0) {
            if (i == argc-1) Usage(sim->exectuable_fn,1);
            sim->out_img_size = atoi(argv[++i]);
         } else if (strncmp(argv[i], "-dot", 4) == 0) {
            sim->write_dot = TRUE;
         } else if (strncmp(argv[i], "-dens", 4) == 0) {
            sim->write_2d_dens = TRUE;
         } else if (strncmp(argv[i], "-slices", 4) == 0) {
            sim->write_3d_dens = TRUE;
            sim->out_img_size = 64;
         } else if (strncmp(argv[i], "-zone", 4) == 0) {
            sim->use_dens_zone = TRUE;
            // do we have enough args?
            if (i > argc-2*DIM-1) Usage(sim->exectuable_fn,1);
            // load them in
            for (d=0;d<DIM;d++) {
              sim->plotzone_dens.min[d] = atof(argv[++i]);
              sim->plotzone_dens.max[d] = atof(argv[++i]);
            }
         } else if (strncmp(argv[i], "-p", 2) == 0) {
            sim->indiv_part[sim->num_indiv_parts][0] = 1;
            j = sim->num_indiv_parts;
            if (i == argc-DIM) Usage(sim->exectuable_fn,1);
            for (d=0;d<DIM;d++) sim->indiv_part[j][d+1] = atof(argv[++i]);
            for (d=0;d<DIM;d++) sim->indiv_part[j][d+1+DIM] = 0.0;
            sim->num_indiv_parts++;
         } else if (strncmp(argv[i], "-help", 2) == 0) {
            (void) Usage(sim->exectuable_fn,0);

         } else {
            fprintf(stderr,"Unrecognized argument (%s)\n",argv[i]);
            (void) Usage(sim->exectuable_fn,0);
         }

      } else {

         // if not, though, save it as the input file name
         (void) strcpy(sim->input_fn,argv[i]);
         sim->use_input_file = TRUE;

         // now, if the file listed on the command line exists, read it
         ret_val = read_input_file(sim,top);

         if (ret_val == 1) {
            fprintf(stdout,"ERROR: cannot read input file, please check syntax\n");
            exit(1);
         }
      }
   }

   return(0);
}


/*
 *  Initialize the system - make static and free element structures, etc.
 */
int initialize_system (sim_ptr sim,cell_ptr top) {

   int i,j,d;
   FLOAT loc[DIM],vel[DIM];
   particle_ptr newpart;

   // -------------------------------------------------------------------------

   // "seed" the random number generator
   srand(sim->seed);

   // cell centers are for tree-traversing
   for (d=0;d<DIM;d++) top->mid[d] = (top->min[d]+top->max[d])/2.0;

   // I think we only need these for the self-gravitation thing
   for (d=0;d<DIM;d++) top->cm[d] = (top->min[d]+top->max[d])/2.0;

   // -------------------------------------------------------------------------
   // now, place the seed particle(s)

   for (i=0; i<sim->num_read_stat; i++)
      add_particles_from_file(sim,top,i,TRUE);

   // then, place specific particles
   for (i=0; i<sim->num_indiv_parts; i++) {
      for (d=0;d<DIM;d++) loc[d] = sim->indiv_part[i][d+1];
      for (d=0;d<DIM;d++) vel[d] = sim->indiv_part[i][d+DIM+1];
      newpart = new_particle(i,sim->indiv_part[i][0],sim->new_part_rad,loc,vel);
      add_particle_to_cell(sim,newpart,top);
   }

   // -------------------------------------------------------------------------
   // load the masks

   if (sim->trim_method == png) {
     read_png(sim->trim_file,&sim->trim_nx,&sim->trim_ny,FALSE,
              &sim->trim_mask,0.,1.,NULL,0.,1.,NULL,0.,1.);
     fprintf(stderr,"# Mask file (%s) is %d x %d\n",sim->trim_file,sim->trim_nx,sim->trim_ny);
     //fprintf(stderr,"# outside %g\n",mask[0][0]);
     //fprintf(stderr,"# %g %g %g %g\n",mask[0][0],mask[50][50],mask[100][100],mask[150][150]);
   }

   // -------------------------------------------------------------------------

   if (sim->write_2d_dens) {
      // set the cell sizes
      for (i=0; i<2; i++) sim->ff2->n[i] = sim->out_img_size;
      for (i=0; i<2; i++) sim->ff2->d[i] = (top->max[i]-top->min[i])/sim->ff2->n[i];
      // fprintf(stdout,"ff2 has %d %d cells, each %g by %g\n",sim->ff2->n[0],sim->ff2->n[1],sim->ff2->d[0],sim->ff2->d[1]);

      // allocate memory for the array
      // ff->rho = allocate_3d_array_F(ff->n[0],ff->n[1],ff->n[2]);
      // sim->ff2->rho = allocate_2d_array_F(sim->ff2->n[0],sim->ff2->n[2]);
      sim->ff2->rho = allocate_2d_array_F(sim->ff2->n[0],sim->ff2->n[1]);

      // zero the array
      for (i=0;i<sim->ff2->n[0];i++)
         for (j=0;j<sim->ff2->n[1];j++)
            sim->ff2->rho[i][j] = 0.0;
            // for (k=0;k<sim->ff2->n[2];k++)
               // ff->rho[i][j][k] = 0.0;
   }

   if (sim->write_3d_dens) {
#if DIM>2
      // set the cell sizes
      for (i=0; i<3; i++) sim->ff3->n[i] = sim->out_img_size;
      for (i=0; i<3; i++) sim->ff3->d[i] = (top->max[i]-top->min[i])/sim->ff3->n[i];

      // allocate memory for the array
      sim->ff3->rho = allocate_3d_array_F(sim->ff3->n[0],sim->ff3->n[1],sim->ff3->n[2]);

      // zero the array
      for (i=0;i<sim->ff3->n[0];i++)
         for (j=0;j<sim->ff3->n[1];j++)
            for (int k=0;k<sim->ff3->n[2];k++)
               sim->ff3->rho[i][j][k] = 0.0;
#endif
   }

   // fprintf(stdout,"\nStarting with %d particles",top->num);

   // and the output file (just dots)
   for (i=0; i<2; i++) sim->out->d[i] = (top->max[i]-top->min[i])/sim->out->n[i];
   if (sim->write_dot || sim->write_png)
      sim->out->rho = allocate_2d_array_F(sim->out_img_size,sim->out_img_size);

   // malloc the space for png writing, if necessary
   if (sim->write_png)
      sim->image = allocate_2d_array_pb(sim->out_img_size,sim->out_img_size,sim->image_depth);

   /* return 0, everything OK */
   return(0);
}


/*
 *  Read particle locations from a file, include them in the simulation
 *  Do this for both static and non-static particles
 */
int add_particles_from_file(sim_ptr sim,cell_ptr top,int i,char is_stat) {

   int d;
   int cnt = 0;
   FLOAT loc[DIM],rad,vel[DIM];
   char token[7][16];
   char twochar[2];
   char sbuf[512];
   particle_ptr newpart;
   FILE *infile;

   /* open file for reading */
   if (is_stat) {
      infile = fopen(sim->read_stat[i],"r");
      if (infile==NULL) {
         fprintf(stderr,"Could not open particle file %s\n",sim->read_stat[i]);
         fflush(stderr);
         return(1);
      }
      fprintf(stdout,"Opening file %s\n",sim->read_stat[i]);
      fflush(stdout);
   } else {
      infile = fopen(sim->read_part[i],"r");
      if (infile==NULL) {
         fprintf(stderr,"Could not open particle file %s\n",sim->read_part[i]);
         fflush(stderr);
         return(1);
      }
      fprintf(stdout,"Opening file %s\n",sim->read_part[i]);
      fflush(stdout);
   }

   /* read a line from the input file */
   while (fscanf(infile,"%[^\n]",sbuf) != EOF) {

      /* fprintf(stdout,"%s\n",sbuf); */

      /* grab the line */
      sscanf(sbuf,"%s %s %s %s %s %s %s",token[0],token[1],token[2],token[3],token[4],token[5],token[6]);
      /* fprintf(stdout,"   first token is %s\n",token[0]); */

      if (token[0][0] == '#') {
         // read a comment line, or some other descriptor
         fscanf(infile,"%[\n]",twochar);    // read up to newline
         // fprintf(stdout,"%s\n",sbuf);        // write comment

      } else {
         /* read the particle parameters */
         /* fprintf(stdout,"   second token is %s\n",token[1]); */

         for (d=0;d<DIM;d++) loc[d] = atof(token[d]);
         rad = atof(token[DIM]);

         for (d=0;d<DIM;d++) {
            loc[d] += sim->stat_xform[i][d];
            vel[d] = 0.0;
         }
         newpart = new_particle(cnt,sim->new_part_mass,rad,loc,vel);

         add_particle_to_cell(sim,newpart,top);

         cnt++;
         if (cnt/1000 != (cnt+1)/1000) {
            fprintf(stdout,".");
            fflush(stdout);
         }

         fscanf(infile,"%[\n]",twochar);    /* read newline */
      }
   }

   fprintf(stdout,"\n");
   fprintf(stdout,"Placed %d particles\n",cnt);
   fflush(stdout);

   fclose(infile);

   /* ret_val is 0 is file is OK and read, 1 if not */
   return(0);
}


/*
 *  Create a new particle
 */
particle_ptr new_particle(int i,FLOAT m,FLOAT r,FLOAT* x,FLOAT* u){

   particle_ptr newpart = new_stationary_particle(i,r,x);
   newpart->bitfield = 0;	// 0 means not stationary
   // fprintf(stderr,"Added particle %d with mass %g and radius %g\n",newpart->index,newpart->mass,newpart->rad);
   return(newpart);
}


/*
 *  Create a new stationary particle
 */
particle_ptr new_stationary_particle(int i,FLOAT r,FLOAT* x){

   int d;
   particle_ptr newpart;

   newpart = (PARTICLE*)malloc(sizeof(PARTICLE));
   newpart->index = i;
   newpart->mass = 0.0;
   newpart->rad = r;
   for (d=0;d<DIM;d++) newpart->x[d] = x[d];
   newpart->counter = 0;
   newpart->bitfield = 1;	// 1 means stationary
   newpart->next = NULL;
   newpart->root = NULL;
   newpart->tip_head = NULL;
   newpart->next_tip = NULL;

   return(newpart);
}


/*
 * Read the input file and set the listed parameters
 */
int read_input_file(sim_ptr sim,cell_ptr top) {

   int d;
   char token[13][32];
   char twochar[2];
   char sbuf[512];
   FILE *infile;

   if (sim->use_input_file) {

      // open file for reading
      infile = fopen(sim->input_fn,"r");
      if (infile==NULL) {
         fprintf(stderr,"Could not open input file %s\n",sim->input_fn);
         fflush(stderr);
         return(1);
      }
      fprintf(stdout,"Opening file %s\n",sim->input_fn);
      fflush(stdout);

   } else {

      // open stdin for reading
      infile = stdin;
      fprintf(stdout,"Opening stdin for reading...\n");
      fflush(stdout);
   }


   /* read a line from the input file */
   while (fscanf(infile,"%[^\n]",sbuf) != EOF) {

      // fprintf(stdout,"%s\n",sbuf);

      /* grab the first word */
      sscanf(sbuf,"%s %s %s %s %s %s %s %s %s %s %s %s %s",token[0],token[1],token[2],token[3],token[4],token[5],token[6],token[7],token[8],token[9],token[10],token[11],token[12]);
      /* fprintf(stdout,"   first token is %s\n",token[0]); */

      if (token[0][0] == '#') {
         // read a comment line, or some other descriptor
         fscanf(infile,"%[\n]",twochar);    // read up to newline
         // fprintf(stdout,"%s\n",sbuf);        // write comment

      } else {
         /* first word is a parameter name */
         /* read the parameter value */
         /* fprintf(stdout,"   second token is %s\n",token[1]); */

         if (strncmp(token[0],"outfile_root",12) == 0) {
            strcpy(sim->out_fn_root,token[1]);
         } else if (strncmp(token[0],"particles",9) == 0) {
            sim->particle_cnt = atoi(token[1]);
         } else if (strncmp(token[0],"write_pgm",9) == 0) {
            if (strncmp(token[1],"yes",1) == 0)
               sim->write_pgm = TRUE;
            else
               sim->write_pgm = FALSE;
         } else if (strncmp(token[0],"write_png",9) == 0) {
            if (strncmp(token[1],"yes",1) == 0)
               sim->write_png = TRUE;
            else
               sim->write_png = FALSE;
         } else if (strncmp(token[0],"write_dot",9) == 0) {
            if (strncmp(token[1],"yes",1) == 0)
               sim->write_dot = TRUE;
            else
               sim->write_dot = FALSE;
         } else if (strncmp(token[0],"write_obj",9) == 0) {
            if (strncmp(token[1],"yes",1) == 0)
               sim->write_obj = TRUE;
            else
               sim->write_obj = FALSE;
         } else if (strncmp(token[0],"write_seg",9) == 0) {
            if (strncmp(token[1],"yes",1) == 0)
               sim->write_seg = TRUE;
            else
               sim->write_seg = FALSE;
         } else if (strncmp(token[0],"write_rad",9) == 0) {
            if (strncmp(token[1],"yes",1) == 0)
               sim->write_rad = TRUE;
            else
               sim->write_rad = FALSE;
         } else if (strncmp(token[0],"write_part",10) == 0) {
            if (strncmp(token[1],"yes",1) == 0)
               sim->write_part = TRUE;
            else
               sim->write_part = FALSE;
         } else if (strncmp(token[0],"write_stats",11) == 0) {
            sim->write_stats = TRUE;
            if (isalnum(token[1][0])) {
               strcpy(sim->stats_file,token[1]);
            }
         } else if (strncmp(token[0],"output_step",9) == 0) {
            sim->output_step = atoi(token[1]);
         } else if (strncmp(token[0],"max_levels",10) == 0) {
            sim->max_levels = atoi(token[1]);
         } else if (strncmp(token[0],"max_ppc",7) == 0) {
            sim->max_parts_in_cell = atoi(token[1]);
         } else if (strncmp(token[0],"seed",4) == 0) {
            sim->seed = atoi(token[1]);
         } else if (strncmp(token[0],"penetration",5) == 0) {
            sim->penetration = atof(token[1]);
         } else if (strncmp(token[0],"stickiness",5) == 0) {
            sim->use_stickiness = TRUE;
            sim->stickiness = atof(token[1]);
         } else if (strncmp(token[0],"comp_domain",11) == 0) {
            /* set initial size of first cell */
            for (d=0;d<DIM;d++) {
               top->min[d] = atof(token[2*d+1]);
               top->max[d] = atof(token[2*d+2]);
            }
         } else if (strncmp(token[0],"boundary",8) == 0) {
            if (strncmp(token[1],"wall",4) == 0) {
               for (d=0;d<DIM;d++) {
                  sim->bdry[d][0] = WALL;
                  sim->bdry[d][1] = WALL;
               }
            } else if (strncmp(token[1],"open",4) == 0) {
               for (d=0;d<DIM;d++) {
                  sim->bdry[d][0] = OPEN;
                  sim->bdry[d][1] = OPEN;
               }
            }
         } else if (strncmp(token[0],"image_size",10) == 0) {
            sim->out_img_size = atoi(token[1]);
            for (d=0;d<2;d++) sim->out->n[d] = atoi(token[1]);
            // fprintf(stdout,"img size is %d\n",sim->out_img_size);
         } else if (strncmp(token[0],"image_depth",11) == 0) {
            sim->image_depth = atoi(token[1]);
         } else if (strncmp(token[0],"rad",3) == 0) {
            sim->new_part_rad = atof(token[1]);
         //} else if (strncmp(token[0],"mass",4) == 0) {
         //   sim->new_part_mass = atof(token[1]);
         } else if (strncmp(token[0],"add_part",8) == 0) {
            sim->indiv_part[sim->num_indiv_parts][0] = atof(token[1]);
            for (d=1;d<2*DIM+1;d++) {
               sim->indiv_part[sim->num_indiv_parts][d] = atof(token[d+1]);
            }
            sim->num_indiv_parts++;
         //} else if (strncmp(token[0],"read_part",9) == 0) {
         //   strcpy(sim->read_part[sim->num_read_part],token[1]);
         //   for (d=0;d<2*DIM;d++)
         //      sim->part_xform[sim->num_read_part][d] = atof(token[d+2]);
         //   sim->num_read_part++;
         } else if (strncmp(token[0],"read_stat",9) == 0) {
            strcpy(sim->read_stat[sim->num_read_stat],token[1]);
            for (d=0;d<DIM;d++)
               sim->stat_xform[sim->num_read_stat][d] = atof(token[d+2]);
            sim->num_read_stat++;
         } else if (strncmp(token[0],"write_dens",10) == 0) {
            // sim->ff2->n[0] = atoi(token[1]);
            // sim->ff2->n[1] = atoi(token[1]);
            // sim->ff2->n[2] = atoi(token[1]);
            // sim->ff3->n[0] = atoi(token[1]);
            // sim->ff3->n[1] = atoi(token[1]);
            // sim->ff3->n[2] = atoi(token[1]);
            // sim->write_2d_dens = TRUE;
            if (strncmp(token[1],"yes",1) == 0)
               sim->write_2d_dens = TRUE;
            else
               sim->write_2d_dens = FALSE;
         } else if (strncmp(token[0],"grip",4) == 0) {
            sim->use_grip = TRUE;
            sim->grip = atof(token[1]);
         } else if (strncmp(token[0],"stub",4) == 0) {
            sim->use_stubborn = TRUE;
            sim->stubborn = atoi(token[1]);
         } else if (strncmp(token[0],"bulk_vel",8) == 0) {
            sim->use_bulk_vel = straight;
            for (d=0;d<DIM;d++) sim->bulk_vel[d] = atof(token[d+1]);
         } else if (strncmp(token[0],"chiral",6) == 0) {
            sim->use_chiral = TRUE;
            if (isalnum(token[1][0])) {
               sim->chiral_diameter = atof(token[1]);
               if (abs(sim->chiral_diameter) < 1.) sim->chiral_diameter = 100.;
               sim->chiral_angle = 2.*M_PI/sim->chiral_diameter;
               if (isalnum(token[2][0])) {
                  sim->chiral_power = atof(token[2]);
               }
            }
         }

         fscanf(infile,"%[\n]",twochar);    /* read newline */
      }
   }

   fclose(infile);

   /* ret_val is 0 is file is OK and read, 1 if not */
   return(0);
}


/*
 * This function writes basic usage information to stderr,
 * and then quits. Too bad.
 */
int Usage(char progname[80],int status) {

   /* Usage for grav3d */
   static char **cpp, *help_message[] =
   {
       "  The input file should be a raw ASCII file containing 2 columns, the first",
       "  with a parameter name (dt,particles,G,write_rad) and the second with either",
       "  an integer or floating-point value or with a \'yes\' or \'no\'.",
       " ",
       "  There is running output to stdout, but PGM and Radiance data are written",
       "  to disk.",
       " ",
       NULL
   };

   fprintf(stderr,"\n  Usage:  %s infile\n\n", progname);
   for (cpp = help_message; *cpp; cpp++) {
      fprintf(stderr, "%s\n", *cpp);
      fflush(stderr);
   }
   exit(status);
   return(0);
}


/*
 * allocate memory for a three-dimensional array of FLOATs
 */
FLOAT*** allocate_3d_array_F(int nx, int ny, int nz) {

   int i,j;
   FLOAT ***array = (FLOAT ***)malloc(nx * sizeof(FLOAT **));

   array[0] = (FLOAT **)malloc(nx * ny * sizeof(FLOAT *));
   array[0][0] = (FLOAT *)malloc(nx * ny * nz * sizeof(FLOAT));

   for (i=1; i<nx; i++)
      array[i] = array[0] + i * ny;

   for (i=0; i<nx; i++) {
      if (i!=0)
         array[i][0] = array[0][0] + i * ny * nz;
      for (j=1; j<ny; j++)
         array[i][j] = array[i][0] + j * nz;
   }

   return(array);
}


/*
 * allocate memory for a two-dimensional array of FLOATs
 */
FLOAT** allocate_2d_array_F(int nx, int ny) {

   int i;
   FLOAT **array = (FLOAT **)malloc(nx * sizeof(FLOAT *));

   array[0] = (FLOAT *)malloc(nx * ny * sizeof(FLOAT));

   for (i=1; i<nx; i++)
      array[i] = array[0] + i * ny;

   return(array);
}

