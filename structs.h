/***********************************************************
 *
 *  structs.h - data structures supporting dla-nd
 *
 *  Copyright (C) 2000-8,14-15  Mark J. Stock
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
 ***********************************************************/


/* 
 * Here are the includes
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "png.h"


/*
 * Here are some commonly-used defines:
 */

// use a (potentially) faster spherical random number generator (2D, 3D only)
#define USE_SINE_SPHERE

// all floating-point numbers are double
#define FLOAT double

// number of spatial dimensions --- must set NCHILD manually (sorry!)
//#define NCHILD (pow(2,DIM))	// doesn't work

// you won't need to change any other defines
#define TRUE 1
#define FALSE 0
#define PI 3.1415926535
#define TWOPI 6.28318531
#define PIOT 1.57079633
#define OPEN 0
#define WALL 1
#define PERIODIC 2
#define SQRT2 1.41421356
#define LARGE 9.9e+99
#define mod(a,b) ((a)-(b)*floor((a)/(b)))
#define FNLEN 255


/*
 * Pointers to the two types of data structures below
 */
typedef struct particle_record *particle_ptr;
typedef struct cell_record *cell_ptr;


/*
 * A data structure definition of a participating particle.
 * The 1st dimension in x and u can be extended to 2 if a 
 * predictor-corrector method is to be used in the advection routine.
 */
typedef struct particle_record {
   FLOAT x[DIM];		/* node's location in space */
   FLOAT rad;			/* particle's radius */
   FLOAT mass;			/* particle's mass */
   unsigned int index;		/* index of the node */
   unsigned short int bitfield;	/* bitfield, uses masks */
   unsigned short int counter;	/* counter for stubbornness */
   //char stationary;		/* is the particle static (immobile) */
   //int flag;			/* generic flag */
   // pointers for tree
   particle_ptr next;		/* pointer to the next node in the list */

   // pointer to rootward
   particle_ptr root;		/* pointer to the rootward node */

   // pointers to tipward
   particle_ptr tip_head;	/* pointer to the head of the list of tipward nodes */
   particle_ptr next_tip;	/* pointer to the next node in the list of tipward nodes */
} PARTICLE;

//unsigned short int stationary_mask = 1;
//unsigned short int flag_mask = 2;


/*
 * A data structure definition of a cell in the octree
 */
typedef struct cell_record {
   char level;			/* recursion level */
   char has_subcells;		/* generic flag */

   /* for now, include both sets of pointers */
   particle_ptr first;		/* pointer to the first element, CCW */
   /* or */
   cell_ptr s[NCHILD];		/* pointers to the 2^DIM subcells */

   FLOAT min[DIM];		/* minimum coordinate of cube */
   FLOAT mid[DIM];		/* midpoint of cube */
   FLOAT max[DIM];		/* minimum coordinate of cube */

   unsigned int num;		/* number of particles in this cell */
   FLOAT mass;			/* total mass of particles in this cell */
   FLOAT cm[DIM];		/* physical center of mass */
} CELL;


/*
 * Field structure and pointer
 */
typedef struct field3_record *field3_ptr;
typedef struct field3_record {
   int n[3];                    /* integer dimensions of the field values */
   FLOAT d[3];                  /* size of each cell */
   FLOAT ***rho;                /* density array */
} FIELD3;

typedef struct field2_record *field2_ptr;
typedef struct field2_record {
   int n[2];                    /* integer dimensions of the field values */
   FLOAT d[2];                  /* size of each cell */
   FLOAT **rho;                 /* density array */
} FIELD2;

typedef enum diffuse_source {
   inf,		// diffuse from infinity
   point,	// always from a single point
   mask		// from a point selected from a mask image
} SOURCE;

typedef enum geometry_trimmer {
   none,	// do not trim to geometry
   cube,	// trim to origin-centered cube
   plane,	// trim to origin-centered thick planar band (in z)
   sphere,	// trim to origin-centered sphere
   png		// trim to input image mask
} TRIMMER;

typedef struct simulation_properties *sim_ptr;
typedef struct simulation_properties {

   int use_input_file;		// TRUE if input file is given on command line
   char input_fn[FNLEN];	/* name of the input file, if not stdin */
   char exectuable_fn[FNLEN];	/* name of the executable */
   char out_fn_root[FNLEN];	/* root name of the output files */
   int write_gif;		/* write a 2D GIF file of the scene */
   int write_pgm;		/* write a 2D ASCII PGM file of the scene */
   int write_png;		/* write PNG files */
   int write_dot;		/* write a simple 1-pixel-per-particle image */
   int write_rad;		/* write a 3D Radiance description of the scene */
   int write_obj;		/* write a scene description in obj format */
   int write_seg;		/* write a scene description in seg format */
   int write_part;		/* write a part3d-native description of the scene */
   int out_img_size;		/* the height and width of the raster output */
   FIELD2 out_field;		/* initialize struct for output array, re-usable, now */
   field2_ptr out;
   int image_depth;		/* bits per pixel in the output image, 8 or 16 */
   png_byte **image;		// memory for the png image
   int write_stats;		// write a statistics datafile?
   char stats_file[FNLEN];	// name of the stats datafile
   FILE *statsfp;		// the file pointer
   int write_temp;		// write a statistics datafile?
   char temp_file[FNLEN];	// name of the stats datafile

   TRIMMER trim_method;		// which, if any, method to trim new segments
   FLOAT trim_rad;		// if sphere or cube, distance to trim
   char trim_file[FNLEN];	// name of the png mask image
   int trim_nx,trim_ny;
   FLOAT **trim_mask;

   SOURCE particle_source;	// diffusing particle start location
   FLOAT source_vec[DIM];
   int use_source_mask;
   char source_file[FNLEN];	// name of the stats datafile

   int seed;			/* random seed */
   int step;			/* index of current particle being placed */
   int output_step;		/* num particles between output */
   double next_output_step;	/* step at which the next set of output file will be written */
   int next_output_index;	/* integer index for next set of output files */
   int max_levels;		/* maximum number of levels in octree */
   int max_parts_in_cell;	/* maximum number of particles in a cell */
   int cell_count[100];		/* count of number of cells at various levels */
   int particle_cnt;		/* initial particle count */
   FLOAT overall_rad;		/* radius of entire structure */
   FLOAT new_part_rad;		/* default particle radius */
   FLOAT new_part_mass;		/* default particle mass */

   FLOAT penetration;		// fraction of radius req'd to contact

   int use_stickiness;
   FLOAT stickiness;		// probability of sticking per contact

   int use_bulk_vel;
   FLOAT bulk_vel[DIM];		// add bulk velocity to random walk

   int use_grip;
   FLOAT grip;			// slide particle this fraction toward root

   int use_stubborn;
   int stubborn;		// how many hits until a particle will allow
 				//   new particles to stick to it

   int use_chiral;
   FLOAT chiral_diameter;	// how many particles fit across a perfect circle
   FLOAT chiral_angle;		// the corresponding angle per particle
   FLOAT chiral_power;		// exponent on arc-following weight

   int use_junction_flow;	// adjust contacted node toward mean flow vector
   FLOAT junction_coeff;	// amount of adjustment (1.0 = average old and new positions)

   char bdry[DIM][2];		/* boundary types, 0=OPEN, 1=WALL, 2=PERIODIC */

   char use_density_field;      /* flags use of density field calculations */
   char construct_dens_surface; /* flags construction of an isosurface of density */
   FIELD2 density_field2;	/* initialize struct for density field */
   field2_ptr ff2;
   FIELD3 density_field3;	/* initialize struct for density field */
   field3_ptr ff3;

   /* construction variables */
   FLOAT block[100][2*DIM+1];	/* store data for placing blocks of particles */
   int num_blocks;
   FLOAT niblock[100][2*DIM+1];	/* same, but for non-intersecting blocks */
   int num_niblocks;

   FLOAT sphere[100][DIM+2];	/* store data for placing blocks of particles */
   int num_spheres;
   FLOAT nisphere[100][DIM+2];	/* same, but for non-intersecting spheres */
   int num_nispheres;

   FLOAT indiv_part[100][2*DIM+1];/* store data for one individual particle */
   int num_indiv_parts;

   char read_part[10][80];	/* part files to be read in and included */
   FLOAT part_xform[10][2*DIM];	/* translation and velocity of particles to be
 				 * read in from said file */
   int num_read_part;

   char read_stat[10][80];	/* part files to be read in and included as stationary particles */
   FLOAT stat_xform[10][DIM];	/* translation of particles to be read in from said file */
   int num_read_stat;

   // stats variables
   int total_particles_tried;
   int total_steps_taken;

   // timing variables
   unsigned long int tics;      /* number of tics (CPU time) used by program so far */
   unsigned long int last_tics; /* number of tics at the end of the last time step */
} SIMP;

/*
 * A timer, flags when stuff should be done
 */
typedef struct event_timer {
   double start_time;		/* time at the start of the simulation */
   double dt;			/* time step for what the process is timing */
   double next_time;		/* next flagged time */
   int next_index;		/* integer index for next flagged time */
} TIMER;


// from ndtree.c
extern cell_ptr build_new_tree(sim_ptr,cell_ptr);
extern int remake_box(sim_ptr,cell_ptr,cell_ptr);
extern int find_particle_bounds(cell_ptr,FLOAT*,FLOAT*);
extern int add_particle_to_cell(sim_ptr,particle_ptr,cell_ptr);
extern int split_cell(sim_ptr,cell_ptr);
extern int count_cells(sim_ptr,cell_ptr);
extern int replace_all_particles(sim_ptr,cell_ptr,cell_ptr);
extern int clean_up_all_cells(sim_ptr,cell_ptr);
extern int clean_up_cells(sim_ptr,cell_ptr,int);
extern int merge_cell(cell_ptr);
extern int find_all_cells_cm(sim_ptr,cell_ptr);
extern int find_cell_cm(cell_ptr,int);
extern int realloc_box(sim_ptr,cell_ptr);

// from setup.c
extern int set_defaults (sim_ptr,cell_ptr);
extern int initialize_system (sim_ptr,cell_ptr);
extern int add_particles_from_file (sim_ptr,cell_ptr,int,char);
extern particle_ptr new_particle (int,FLOAT,FLOAT,FLOAT*,FLOAT*);
extern particle_ptr new_stationary_particle (int,FLOAT,FLOAT*);
extern int read_input_file (sim_ptr,cell_ptr);
extern int parse_args (int,char**,sim_ptr,cell_ptr);
extern int Usage (char*,int);
extern FLOAT*** allocate_3d_array_F (int,int,int);
extern FLOAT** allocate_2d_array_F (int,int);

// from writeout.c
extern int write_output(sim_ptr,cell_ptr);
extern int write_particle_count(cell_ptr);
extern int write_2d_dots(sim_ptr,cell_ptr,int);
//extern int write_pgm_density_3d(cell_ptr,field3_ptr,char*);
//extern int write_2d_density(sim_ptr,cell_ptr,field2_ptr,int);
extern int put_part_into_array(cell_ptr,cell_ptr,FLOAT**,int,int);
extern int write_rad(sim_ptr,cell_ptr,char*);
extern int write_rad_cell(cell_ptr,FILE*,int);
extern int write_part(cell_ptr,char*);
extern int write_part_cell(cell_ptr,FILE*);

// from density.c
//extern int create_density_field_3d(cell_ptr,cell_ptr,field3_ptr);
//extern int add_particle_to_field_3d(cell_ptr,particle_ptr,field3_ptr);
//extern int create_density_field_2d(cell_ptr,cell_ptr,cell_ptr,field2_ptr);
//extern int add_particle_to_field_2d(cell_ptr,particle_ptr,cell_ptr,field2_ptr);

// from readin.c
extern int read_png (char*, int*, int*, int, FLOAT***, FLOAT, FLOAT, FLOAT***, FLOAT, FLOAT, FLOAT***, FLOAT, FLOAT);
extern png_byte** allocate_2d_array_pb (int,int,int);
extern png_byte** allocate_2d_rgb_array_pb (int,int,int);
extern int free_2d_array_pb (png_byte**);

/* end */

