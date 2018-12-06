/*************************************************************
 *
 *  inout.c - input/output subroutines for dla-nd
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

#include "structs.h"

#define png_infopp_NULL (png_infopp)NULL
#define int_p_NULL (int*)NULL

int push_bounds_out(cell_ptr, cell_ptr);
int write_obj(sim_ptr,cell_ptr,char*);
int write_obj_nodes(cell_ptr,FILE*,int);
int write_obj_segs(cell_ptr,FILE*,int);
int write_seg(sim_ptr,cell_ptr,char*);
int write_seg_nodes(cell_ptr,FILE*,int);
int write_seg_segs(cell_ptr,FILE*,int);
int write_png(char*,png_byte**,int,int,int);


/*
 * read a PNG, write it to 1 or 3 channels
 */
int read_png (char *infile, int *nx, int *ny,
   int expect_three_channel,
   FLOAT ***red, FLOAT redmin, FLOAT redrange,
   FLOAT ***grn, FLOAT grnmin, FLOAT grnrange,
   FLOAT ***blu, FLOAT blumin, FLOAT blurange) {

   //int autorange = FALSE;
   int high_depth;
   int three_channel;
   int i,j; //,printval,bit_depth,color_type,interlace_type;
   //FLOAT newminrange,newmaxrange;
   FILE *fp;
   unsigned char header[8];
   png_uint_32 height,width;
   int bit_depth,color_type,interlace_type;
   png_structp png_ptr;
   png_infop info_ptr;
   png_byte **img;

   // check the file
   fp = fopen(infile,"rb");
   if (fp==NULL) {
      fprintf(stderr,"Could not open input file %s\n",infile);
      fflush(stderr);
      exit(0);
   }

   // check to see that it's a PNG
   fread (&header, 1, 8, fp);
   if (png_sig_cmp(header, 0, 8)) {
      fprintf(stderr,"File %s is not a PNG\n",infile);
      fflush(stderr);
      exit(0);
   }

   /* Create and initialize the png_struct with the desired error handler
    * functions.  If you want to use the default stderr and longjump method,
    * you can supply NULL for the last three parameters.  We also supply the
    * the compiler header file version, so that we know if the application
    * was compiled with a compatible version of the library.  REQUIRED
    */
   png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING,
      NULL, NULL, NULL);

   /* Allocate/initialize the memory for image information.  REQUIRED. */
   info_ptr = png_create_info_struct(png_ptr);
   if (info_ptr == NULL) {
      fclose(fp);
      png_destroy_read_struct(&png_ptr, png_infopp_NULL, png_infopp_NULL);
      exit(0);
   }

   /* Set error handling if you are using the setjmp/longjmp method (this is
    * the normal method of doing things with libpng).  REQUIRED unless you
    * set up your own error handlers in the png_create_read_struct() earlier.  */
   if (setjmp(png_jmpbuf(png_ptr))) {
      /* Free all of the memory associated with the png_ptr and info_ptr */
      png_destroy_read_struct(&png_ptr, &info_ptr, png_infopp_NULL);
      fclose(fp);
      /* If we get here, we had a problem reading the file */
      exit(0);
   }

   /* One of the following I/O initialization methods is REQUIRED */
   /* Set up the input control if you are using standard C streams */
   png_init_io(png_ptr, fp);

   /* If we have already read some of the signature */
   png_set_sig_bytes(png_ptr, 8);

   /* The call to png_read_info() gives us all of the information from the
    * PNG file before the first IDAT (image data chunk).  REQUIRED */
   png_read_info(png_ptr, info_ptr);

   png_get_IHDR(png_ptr, info_ptr, &width, &height, &bit_depth, &color_type,
       &interlace_type, int_p_NULL, int_p_NULL);

   /* Set up the data transformations you want.  Note that these are all
    * optional.  Only call them if you want/need them.  Many of the
    * transformations only work on specific types of images, and many
    * are mutually exclusive.  */

   /* tell libpng to strip 16 bit/color files down to 8 bits/color */
   //png_set_strip_16(png_ptr);

   /* Extract multiple pixels with bit depths of 1, 2, and 4 from a single
    * byte into separate bytes (useful for paletted and grayscale images).  */
   png_set_packing(png_ptr);

   /* Expand paletted colors into true RGB triplets */
   //if (color_type == PNG_COLOR_TYPE_PALETTE)
   //   png_set_palette_rgb(png_ptr);

   /* Expand grayscale images to the full 8 bits from 1, 2, or 4 bits/pixel */
   if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
      png_set_expand_gray_1_2_4_to_8(png_ptr);

   /* Optional call to gamma correct and add the background to the palette
    * and update info structure.  REQUIRED if you are expecting libpng to
    * update the palette for you (ie you selected such a transform above).
    */
   //png_read_update_info(png_ptr, info_ptr);

   // check image type for applicability
   if (bit_depth != 8 && bit_depth != 16) {
     fprintf(stderr,"INCOMPLETE: read_png expect 8-bit or 16-bit images\n");
     exit(0);
   }
   if (color_type != PNG_COLOR_TYPE_GRAY && color_type != PNG_COLOR_TYPE_RGB) {
     fprintf(stderr,"INCOMPLETE: read_png expect grayscale or RGB images\n");
     exit(0);
   }

   // set channels
   if (color_type == PNG_COLOR_TYPE_GRAY) three_channel = FALSE;
   else three_channel = TRUE;

   if (expect_three_channel && !three_channel) {
     fprintf(stderr,"ERROR: expecting 3-channel PNG, but input is 1-channel\n");
     fprintf(stderr,"  file (%s)",infile);
     fprintf(stderr,"  Convert file to color and try again.\n");
     exit(0);
   }

   if (!expect_three_channel && three_channel) {
     fprintf(stderr,"ERROR: not expecting 3-channel PNG, but input is 3-channel\n");
     fprintf(stderr,"  file (%s)",infile);
     fprintf(stderr,"  Convert file to grayscale and try again.\n");
     exit(0);
   }

   // set specific bit depth
   if (bit_depth == 16) high_depth = TRUE;
   else high_depth = FALSE;

   // make sure pixel sizes match exactly!
   //if (ny != height || nx != width) {
   //  fprintf(stderr,"INCOMPLETE: read_png expects image resolution to match\n");
   //  fprintf(stderr,"  the simulation resolution.");
   //  fprintf(stderr,"  simulation %d x %d",nx,ny);
   //  fprintf(stderr,"  image %d x %d",width,height);
   //  fprintf(stderr,"  file (%s)",infile);
   //  exit(0);
   //}

   // set the sizes so that we can understand them
   (*ny) = height;
   (*nx) = width;

   // allocate the space for the image array
   if (three_channel) {
      img = allocate_2d_rgb_array_pb((*nx),(*ny),bit_depth);
   } else {
      img = allocate_2d_array_pb((*nx),(*ny),bit_depth);
   }

   // allocate space for the mask channel (red)
   (*red) = allocate_2d_array_F(width,height);

   /* Now it's time to read the image.  One of these methods is REQUIRED */
   png_read_image(png_ptr, img);

   /* read rest of file, and get additional chunks in info_ptr - REQUIRED */
   png_read_end(png_ptr, info_ptr);

   /* At this point you have read the entire image */

   /* clean up after the read, and free any memory allocated - REQUIRED */
   png_destroy_read_struct(&png_ptr, &info_ptr, png_infopp_NULL);

   /* close the file */
   fclose(fp);


   // now convert the data to stuff we can use
   if (three_channel) {

     // no scaling, 16-bit per channel, RGB
     if (high_depth) {
         for (j=(*ny)-1; j>=0; j--) {
           for (i=0; i<(*nx); i++) {
             (*red)[i][j] = redmin+redrange*(img[(*ny)-1-j][6*i]*256+img[(*ny)-1-j][6*i+1])/65535.;
             (*grn)[i][j] = grnmin+grnrange*(img[(*ny)-1-j][6*i+2]*256+img[(*ny)-1-j][6*i+3])/65535.;
             (*blu)[i][j] = blumin+blurange*(img[(*ny)-1-j][6*i+4]*256+img[(*ny)-1-j][6*i+5])/65535.;
           }
         }

     // no scaling, 8-bit per channel, RGB
     } else {
         for (j=(*ny)-1; j>=0; j--) {
           for (i=0; i<(*nx); i++) {
             (*red)[i][j] = redmin+redrange*img[(*ny)-1-j][3*i]/255.;
             (*grn)[i][j] = grnmin+grnrange*img[(*ny)-1-j][3*i+1]/255.;
             (*blu)[i][j] = blumin+blurange*img[(*ny)-1-j][3*i+2]/255.;
           }
         }
     }

   // monochrome image, read data from red array
   } else {

     // no scaling, 16-bit per channel
     if (high_depth) {
         for (j=(*ny)-1; j>=0; j--) {
           for (i=0; i<(*nx); i++) {
             (*red)[i][j] = redmin+redrange*(img[(*ny)-1-j][2*i]*256+img[(*ny)-1-j][2*i+1])/65535.;
           }
         }

     // no scaling, 8-bit per channel
     } else {
         for (j=(*ny)-1; j>=0; j--) {
           for (i=0; i<(*nx); i++) {
             (*red)[i][j] = redmin+redrange*img[(*ny)-1-j][i]/255.;
             //fprintf(stderr,"%d %d %g\n",j,i,(*red)[i][j]);
           }
         }
     }

   }

   // free the data array
   free_2d_array_pb(img);

   return(0);
}


/*
 * allocate memory for a two-dimensional array of png_byte
 */
/*
png_byte** allocate_2d_array_pb(int nx, int ny, int depth) {

   int i,bytesperpixel;
   png_byte **array;

   if (depth <= 8) bytesperpixel = 1;
   else bytesperpixel = 2;
   array = (png_byte **)malloc(ny * sizeof(png_byte *));
   array[0] = (png_byte *)malloc(bytesperpixel * nx * ny * sizeof(png_byte));

   for (i=1; i<ny; i++)
      array[i] = array[0] + i * bytesperpixel * nx;

   return(array);
}
*/

png_byte** allocate_2d_rgb_array_pb(int nx, int ny, int depth) {

   int i,bytesperpixel;
   png_byte **array;

   if (depth <= 8) bytesperpixel = 3;
   else bytesperpixel = 6;
   array = (png_byte **)malloc(ny * sizeof(png_byte *));
   array[0] = (png_byte *)malloc(bytesperpixel * nx * ny * sizeof(png_byte));

   for (i=1; i<ny; i++)
      array[i] = array[0] + i * bytesperpixel * nx;

   return(array);
}

int free_2d_array_pb(png_byte** array){
   free(array[0]);
   free(array);
   return(0);
}


/*
 * Controls writing of output, once decision to write has been made
 */
int write_output(sim_ptr sim,cell_ptr top) {

   char outfile[511];
   //cell_ptr plotzone = (CELL*)malloc(sizeof(CELL));

   /* actually write the output here */
   if (sim->write_dot) {
      write_2d_dots(sim,top,sim->next_output_index);
   }
   if (sim->write_rad) {
      //zero_all_masses(top);
      //sum_tipward_mass(top);
      sprintf(outfile,"%spart%04d.rad",sim->out_fn_root,sim->next_output_index);
      write_rad(sim,top,outfile);
   }
   if (sim->write_part) {
      sprintf(outfile,"%spart%04d.part",sim->out_fn_root,sim->next_output_index);
      write_part(top,outfile);
   }
   if (sim->write_obj) {
      //zero_all_masses(top);
      //sum_tipward_mass(top);
      sprintf(outfile,"%spart%04d.obj",sim->out_fn_root,sim->next_output_index);
      write_obj(sim,top,outfile);
   }
   if (sim->write_seg) {
      sprintf(outfile,"%spart%04d.seg",sim->out_fn_root,sim->next_output_index);
      write_seg(sim,top,outfile);
   }

   // if we are to write a temporary, replaceable file; but only if it
   //   wouldn't be written normally
   if (sim->write_temp && !sim->write_part) {
      write_part(top,sim->temp_file);
   }

   // create the density field, hopefully more later
/*
   if (sim->use_density_field) {
      // fprintf(stderr,"to write density field\n");
      // first, determine the plot zone (sufficient to hold all particles)
      define_plotzone(top,plotzone);
      // then, set the size of each cell in the output image
      sim->ff2->d[0] = (plotzone->max[0]-plotzone->min[0])/sim->ff2->n[0];
      sim->ff2->d[1] = (plotzone->max[1]-plotzone->min[1])/sim->ff2->n[1];
      fprintf(stdout,"plotzone is [%g %g] [%g %g]\n",plotzone->min[0],
              plotzone->max[0],plotzone->min[1],plotzone->max[1]);
      fprintf(stdout,"ff2 has n [%d %d] and d [%g %g]\n",sim->ff2->n[0],
              sim->ff2->n[1],sim->ff2->d[0],sim->ff2->d[1]);
      // then, fill the ff2 field
      create_density_field_2d(top,top,plotzone,sim->ff2);
      // finally, scale and write the image
      write_2d_density(top,sim->ff2,sim->next_output_index);
   }
*/

   return(0);
}


/*
 * write a heirarchical description of the particles in cells
 */
int write_particle_count(cell_ptr cell){

   int i;

   if (cell->level == 0) fprintf(stdout,"\n");

   if (cell->num > 0) {
      for (i=-1;i<cell->level;i++) fprintf(stdout,"  ");
      fprintf(stdout,"cell at level %d has %d particles\n",cell->level,cell->num);
   }

   if (cell->has_subcells)
      for (i=0;i<NCHILD;i++)
         write_particle_count(cell->s[i]);

   return(0);
}


/*
 * Determine the overall bounds of the simulation for purposes of 
 * plotting
 */
int define_plotzone(cell_ptr top, cell_ptr plotzone) {

   int i;
   int always_center_on_origin = FALSE;
   FLOAT halfwidth,maxhalfwidth;

   // first, set the plot zone boundaries
   for (i=0;i<DIM;i++) {
      plotzone->min[i] = 9.9e+9;
      plotzone->max[i] = -9.9e+9;
   }

   // then, recurse through the entire tree, pushing the boundaries
   //    out to their maximum
   push_bounds_out(top,plotzone);

   if (always_center_on_origin) {
      // Center at origin always, scale such that all particle are inside
      maxhalfwidth = 0.;
      for (i=0;i<DIM;i++) {
         // plotzone->mid[i] = 0.5 * (plotzone->min[i]+plotzone->max[i]);
         plotzone->mid[i] = 0.;
         halfwidth = plotzone->max[i]-plotzone->mid[i];
         if (halfwidth > maxhalfwidth) maxhalfwidth = halfwidth;
         halfwidth = plotzone->mid[i]-plotzone->min[i];
         if (halfwidth > maxhalfwidth) maxhalfwidth = halfwidth;
      }
      maxhalfwidth *= 1.1;
      for (i=0;i<DIM;i++) {
         plotzone->min[i] = plotzone->mid[i] - maxhalfwidth;
         plotzone->max[i] = plotzone->mid[i] + maxhalfwidth;
      }

   } else {
      // Center the structure in the middle of the image
      maxhalfwidth = 0.;
      for (i=0;i<DIM;i++) {
         plotzone->mid[i] = 0.5 * (plotzone->min[i]+plotzone->max[i]);
         halfwidth = plotzone->max[i]-plotzone->mid[i];
         if (halfwidth > maxhalfwidth) maxhalfwidth = halfwidth;
         halfwidth = plotzone->mid[i]-plotzone->min[i];
         if (halfwidth > maxhalfwidth) maxhalfwidth = halfwidth;
      }
      maxhalfwidth *= 1.1;
      for (i=0;i<DIM;i++) {
         plotzone->min[i] = plotzone->mid[i] - maxhalfwidth;
         plotzone->max[i] = plotzone->mid[i] + maxhalfwidth;
      }

   }

   return(0);
}


/*
 * Determine the overall bounds of the particles
 */
int push_bounds_out(cell_ptr cell, cell_ptr bounds) {

   int i;
   particle_ptr curr;

   if (cell->has_subcells)
      for (i=0;i<NCHILD;i++)
         push_bounds_out(cell->s[i],bounds);

   curr = cell->first;
   while (curr) {
      for (i=0;i<DIM;i++) {
         if (curr->x[i] < bounds->min[i]) bounds->min[i] = curr->x[i];
         if (curr->x[i] > bounds->max[i]) bounds->max[i] = curr->x[i];
      }
      curr = curr->next;
   }

   return(0);
}


/*
 * Write a PGM image of the points projected onto the xy-plane
 */
int write_2d_dots (sim_ptr sim,cell_ptr cell,int index){

   int do_fireworks = FALSE;
   int i,j;
   int nx = sim->out_img_size;
   int ny = sim->out_img_size;
   char filename[511];
   FLOAT scale = 1.0;
   //FIELD2 my_array;		// initialize struct for density field
   //field2_ptr array = &my_array;
   int printval;
   FILE *outfile;

   // fprintf(stderr,"I'm here %d\n",sim->out_img_size); fflush(stderr);

   // make the full filename, default to png
   if (sim->write_pgm) {
      sprintf(filename,"%sdots_%04d.pgm",sim->out_fn_root,index);
   } else {
      sprintf(filename,"%sdots_%04d.png",sim->out_fn_root,index);
   }

   /* zero the array */
   /* Not zeroing this array must use the pervious numbers, and this
      thing just adds the new locations to the old, looks like the
      burst of a firework! */
   if (do_fireworks) {
      /* do not clear memory */
   } else {
      for (j=ny-1; j>=0; j--) {
         for (i=0; i<nx; i++) {
            sim->out->rho[i][j] = 0.0;
         }
      }
   }

   /* fill in the array */
   put_part_into_array(cell,cell,sim->out->rho,nx,ny);

   /* if it's a PGM, write it here, else, write a PNG */
   if (sim->write_pgm) {

      /* open file for writing */
      outfile = fopen(filename,"w");
      if (outfile==NULL) {
         fprintf(stderr,"Could not open output file %s\n",filename);
         fflush(stderr);
         exit(0);
      }

      /* plot a z-plane */
      if (sim->image_depth == 8) {
         fprintf(outfile,"P2\n%d %d\n%d\n",nx,ny,255);
         for (j=ny-1; j>=0; j--) {
            for (i=0; i<nx; i++) {
               printval = (int)(256*sim->out->rho[i][j]*scale);
               if (printval<0) printval = 0;
               if (printval>255) printval = 255;
               fprintf(outfile,"%d\n",printval);
            }
         }
      } else {
         fprintf(outfile,"P2\n%d %d\n%d\n",nx,ny,65535);
         for (j=ny-1; j>=0; j--) {
            for (i=0; i<nx; i++) {
               printval = (int)(65536*sim->out->rho[i][j]*scale);
               if (printval<0) printval = 0;
               if (printval>65535) printval = 65535;
               fprintf(outfile,"%d\n",printval);
            }
         }
      }

      /* close file */
      fclose(outfile);

   } else {

      // convert the floating pt "rho" into "png_byte"
      if (sim->image_depth == 8) {
         for (j=ny-1; j>=0; j--) {
            for (i=0; i<nx; i++) {
               printval = (int)(256*sim->out->rho[i][j]*scale);
               if (printval<0) printval = 0;
               if (printval>255) printval = 255;
               sim->image[ny-1-j][i] = (png_byte)printval;
            }
         }
         write_png(filename,sim->image,nx,ny,8);
      } else {
         for (j=ny-1; j>=0; j--) {
            for (i=0; i<nx; i++) {
               printval = (int)(65536*sim->out->rho[i][j]*scale);
               if (printval<0) printval = 0;
               if (printval>65535) printval = 65535;
               sim->image[ny-1-j][2*i] = (png_byte)(printval/256);
               sim->image[ny-1-j][2*i+1] = (png_byte)(printval%256);
            }
         }
         write_png(filename,sim->image,nx,ny,16);
      }
   }

   /* if all went well, return zero */
   return(0);
}


/*
 * Write a PGM image of the density field
 */
int write_pgm_density_3d(cell_ptr cell,field3_ptr ff,char *filename){

   int do_middle_slice = FALSE;
   int i,j,k;
   int kindex = ff->n[1]/2;
   int nx = ff->n[0];
   int ny = ff->n[2];
   FLOAT scale = 2.0;
   FLOAT **array = (FLOAT **)malloc(nx * sizeof(FLOAT *));
   int printval;
   FILE *outfile;

   /* allocate space for 2D array */
   array[0] = (FLOAT *)malloc(nx * ny * sizeof(FLOAT));
   for (i=1; i<nx; i++)
      array[i] = array[0] + i * ny;

   /* zero the array */
   for (j=ny-1; j>=0; j--)
      for (i=0; i<nx; i++)
         array[i][j] = 0.0;

   /* fill the array */
   for (i=0; i<nx; i++)
      for (j=0; j<ny; j++) {
         if (do_middle_slice) array[i][j] = ff->rho[i][kindex][j];
         else {
            for (k=0; k<ff->n[1]; k++)
               array[i][j] += ff->rho[i][k][j];
            array[i][j] *= 0.0025;
         }
      }

   // scale the array?


   /* open file for writing */
   outfile = fopen(filename,"w");
   if (outfile==NULL) {
      fprintf(stderr,"Could not open output file %s\n",filename);
      fflush(stderr);
      exit(0);
   }

   /* plot a y-plane */
   fprintf(outfile,"P2\n%d %d\n%d\n",nx,ny,255);
   for (j=ny-1; j>=0; j--) {
      for (i=0; i<nx; i++) {
         printval = (int)(256.0*array[i][j]*scale);
         if (printval<0) printval = 0;
         if (printval>255) printval = 255;
         fprintf(outfile,"%d\n",printval);
      }
   }

   /* close file */
   fclose(outfile);

   /* free memory from 2D array! */
   free(array[0]);
   free(array);

   /* if all went well, return zero */
   return(0);
}


/*
 * Write a PGM image of the density field
 */
int write_2d_density (sim_ptr sim,cell_ptr cell,field2_ptr ff,int index){

   int contrast_enhance = TRUE;
   int i,j;
   int nx = ff->n[0];
   int ny = ff->n[1];
   static FLOAT maxval = -1.0;
   //FLOAT scale = 1.0;
   int printval;
   char filename[511];
   FILE *outfile;

   // make the full filename, default to png
   if (sim->write_pgm) {
      sprintf(filename,"%sdens_%04d.pgm",sim->out_fn_root,index);
   } else {
      sprintf(filename,"%sdens_%04d.png",sim->out_fn_root,index);
   }

   // scale ff->rho to unit-maximum
   if (contrast_enhance || maxval < 0.0) {
      for (j=ny-1; j>=0; j--) {
         for (i=0; i<nx; i++) {
             if (ff->rho[i][j] > maxval) maxval = ff->rho[i][j];
         }
      }
      // allow larger areas to be white
      if (DIM==2) maxval *= 0.8;	// this one works well
      else if (DIM==3) maxval *= 0.5;	// this one...not so well.
      else maxval *= 0.2;
   }
   // fprintf(stdout,"maxval is %g\n",maxval);


   /* if it's a PGM, write it here, else, write a PNG */
   if (sim->write_pgm) {

      /* open file for writing */
      outfile = fopen(filename,"w");
      if (outfile==NULL) {
         fprintf(stderr,"Could not open output file %s\n",filename);
         fflush(stderr);
         exit(0);
      }

      /* plot a y-plane */
      if (sim->image_depth == 8) {
         fprintf(outfile,"P2\n%d %d\n%d\n",nx,ny,255);
         for (j=ny-1; j>=0; j--) {
            for (i=0; i<nx; i++) {
               printval = (int)(256.0*ff->rho[i][j]/maxval);
               if (printval<0) printval = 0;
               if (printval>255) printval = 255;
               fprintf(outfile,"%d\n",printval);
            }
         }
      } else {
         fprintf(outfile,"P2\n%d %d\n%d\n",nx,ny,65535);
         for (j=ny-1; j>=0; j--) {
            for (i=0; i<nx; i++) {
               printval = (int)(65536.0*ff->rho[i][j]/maxval);
               if (printval<0) printval = 0;
               if (printval>65535) printval = 65535;
               fprintf(outfile,"%d\n",printval);
            }
         }
      }

      /* close file */
      fclose(outfile);

   } else {

      // convert the floating pt "rho" into "png_byte"
      if (sim->image_depth == 8) {
         for (j=ny-1; j>=0; j--) {
            for (i=0; i<nx; i++) {
               printval = (int)(256*ff->rho[i][j]/maxval);
               if (printval<0) printval = 0;
               if (printval>255) printval = 255;
               sim->image[ny-1-j][i] = (png_byte)printval;
            }
         }
         write_png(filename,sim->image,nx,ny,8);
      } else {
         for (j=ny-1; j>=0; j--) {
            for (i=0; i<nx; i++) {
               printval = (int)(65536*ff->rho[i][j]/maxval);
               if (printval<0) printval = 0;
               if (printval>65535) printval = 65535;
               sim->image[ny-1-j][2*i] = (png_byte)(printval/256);
               sim->image[ny-1-j][2*i+1] = (png_byte)(printval%256);
            }
         }
         write_png(filename,sim->image,nx,ny,16);
      }
   }

   /* if all went well, return zero */
   return(0);
}


/*
 * Put all of the particles in this cell onto the array
 */
int put_part_into_array(cell_ptr top,cell_ptr cell,FLOAT** array,int nx,int ny){

   int i,j;
   int xdim,ydim;
   FLOAT xstart,xsize,mult;//ystart,ysize;
   particle_ptr curr;

   xdim = 0;	/* use x-values for x-dim of image */
   // ydim = 2;	/* use z-values for y-dim of image */
   ydim = 1;	/* use z-values for y-dim of image */

   if (cell->has_subcells) {
      for (i=0;i<NCHILD;i++)
         /* put_part_into_array(top,cell,&array[0][0],nx,ny); */
         put_part_into_array(top,cell->s[i],array,nx,ny);
   } else {
      /* run through the list of elems and place each on the grid */
      // mult = nx*ny*0.05;
      mult = 255.0*sqrt(nx*ny);
      xstart = top->min[xdim];
      xsize = top->max[xdim]-xstart;
      //ystart = top->min[ydim];
      //ysize = top->max[ydim]-ystart;
      curr = cell->first;
      while (curr) {
         i = (int)(nx*(curr->x[xdim]-xstart)/xsize);
         j = (int)(ny*(curr->x[ydim]-xstart)/xsize);
         if (i>-1 && j>-1 && i<nx && j<ny)
            array[i][j] += mult * curr->rad;
            // array[i][j] += mult * curr->mass;
         curr = curr->next;
      }
   }
   return(0);
}


/*
 * Write a description of the nodes and segments in modified .obj file
 */
int write_obj(sim_ptr sim,cell_ptr top,char *filename){

   int numnodes,numsegs;
   FILE *outfile;

   /* open file for writing */
   outfile = fopen(filename,"w");
   // outfile = fopen(filename,"a");
   if (outfile==NULL) {
      fprintf(stderr,"Could not open output file %s\n",filename);
      fflush(stderr);
      exit(0);
   }

   // header
   fprintf(outfile,"# scene description of dla-nd run\n");
   // each node will have 4 floating pt values, 4th is radius
   fprintf(outfile,"d %d\n",DIM+1);

   // march through nodes, writing, and setting indicies
   numnodes = write_obj_nodes(top,outfile,0);

   // march through segments and write them
   numsegs = write_obj_segs(top,outfile,0);

   /* close file */
   fclose(outfile);

   fprintf(stdout,"Wrote %d nodes and %d segments.\n",numnodes,numsegs);

   /* if all went well, return zero */
   return(0);
}


/*
 * Write the nodes' locations, and set their indexes, starting at 1
 */
int write_obj_nodes(cell_ptr c,FILE* outfile,int index){

   int i;//,cnt;
   FLOAT rad;
   particle_ptr curr;

   // fprintf(outfile,"\n# cell at level %d\n",c->level);

   if (c->has_subcells) {
      /* recurse for the subcells */
      for (i=0;i<NCHILD;i++)
         index = write_obj_nodes(c->s[i],outfile,index);
   } else {
      /* write out the nodes */
      curr = c->first;
      //nt = 0;
      while (curr) {
         // Does Murray's Law work in 2D, 3D or what?
         // rad = 0.05*curr->rad*pow(curr->mass,1./3.);
         rad = 0.05*curr->rad*pow(curr->mass,1./DIM);
         fprintf(outfile,"v ");
         for (i=0;i<DIM;i++) fprintf(outfile,"%g ",curr->x[i]);
         fprintf(outfile,"%g\n",rad);
         // fprintf(outfile,"v %g %g %g %g\n",curr->x[0],curr->x[1],curr->x[2],rad);
         curr->index = index++;
         // index++;
         curr = curr->next;
      }
   }

   // return total count
   return(index);
}


/*
 * Write the segments' connectivity
 */
int write_obj_segs(cell_ptr c,FILE* outfile,int index){

   int i;//,cnt;
   particle_ptr curr;

   // fprintf(outfile,"\n# cell at level %d\n",c->level);

   if (c->has_subcells) {
      /* recurse for the subcells */
      for (i=0;i<NCHILD;i++)
         index = write_obj_segs(c->s[i],outfile,index);
   } else {
      /* write out the segments */
      curr = c->first;
      //cnt = 0;
      while (curr) {
         if (curr->root) {
            fprintf(outfile,"s %d %d\n",curr->index,curr->root->index);
            index++;
         }
         curr = curr->next;
      }
   }

   // return total count
   return(index);
}


/*
 * Write a description of the nodes and segments in modified .obj file
 */
int write_seg (sim_ptr sim,cell_ptr top,char *filename){

   int numnodes,numsegs;
   FILE *outfile;

   /* open file for writing */
   outfile = fopen(filename,"w");
   if (outfile==NULL) {
      fprintf(stderr,"Could not open output file %s\n",filename);
      fflush(stderr);
      exit(0);
   }

   // header
   fprintf(outfile,"# scene description of dla-nd run\n");
   fprintf(outfile,"d %d\n",DIM);
   fprintf(outfile,"gr %d\n",-1);

   // march through nodes, writing, and setting indicies
   numnodes = write_seg_nodes(top,outfile,0);

   // march through segments and write them
   numsegs = write_seg_segs(top,outfile,0);

   /* close file */
   fclose(outfile);

   fprintf(stdout,"Wrote %d nodes and %d segments.\n",numnodes,numsegs);

   /* if all went well, return zero */
   return(0);
}


/*
 * Write the nodes' locations, and set their indexes, starting at 1
 */
int write_seg_nodes(cell_ptr c,FILE* outfile,int index){

   //int i;//,cnt;
   //FLOAT rad;
   //particle_ptr curr;

   // fprintf(outfile,"\n# cell at level %d\n",c->level);

   if (c->has_subcells) {
      /* recurse for the subcells */
      for (int i=0;i<NCHILD;i++)
         index = write_seg_nodes(c->s[i],outfile,index);
   } else {
      /* write out the nodes */
      particle_ptr curr = c->first;
      //cnt = 0;
      while (curr) {
         // New: do not write radii here!
         // Does Murray's Law work in 2D, 3D or what?
         //rad = 0.05*curr->rad*pow(curr->mass,1./DIM);
         fprintf(outfile,"v");
         for (int i=0;i<DIM;i++) fprintf(outfile," %g",curr->x[i]);
         fprintf(outfile,"\n");
         //fprintf(outfile,"\nvr %g\n",rad);
         // fprintf(outfile,"v %g %g %g %g\n",curr->x[0],curr->x[1],curr->x[2],rad);
         curr->index = ++index;
         // index++;
         curr = curr->next;
      }
   }

   // return total count
   return(index);
}


/*
 * Write the segments' connectivity
 */
int write_seg_segs (cell_ptr c, FILE* outfile, int index){

   // fprintf(outfile,"\n# cell at level %d\n",c->level);

   if (c->has_subcells) {
      /* recurse for the subcells */
      for (int i=0;i<NCHILD;i++)
         index = write_seg_segs(c->s[i],outfile,index);
   } else {
      /* write out the segments */
      particle_ptr curr = c->first;
      while (curr) {
         if (curr->root) {
            // write two unique radii first
            FLOAT rad1 = 0.05*curr->rad*pow(curr->mass,    1./DIM);
            FLOAT rad2 = 0.05*curr->rad*pow(curr->mass+1.0,1./DIM);
            fprintf(outfile,"vr %g\n",rad1);
            fprintf(outfile,"vr %g\n",rad2);
            // then write the segment
            fprintf(outfile,"s %d/%d %d/%d\n", curr->index, 2*index+1,
                    curr->root->index, 2*index+2);
            index++;
         }
         curr = curr->next;
      }
   }

   // return total count
   return(index);
}


/*
 * Write a Radiance description of the particles and cell edges
 */
int write_rad(sim_ptr sim,cell_ptr top,char *filename){

   int type = 3;
   FILE *outfile;

   /* open file for writing */
   outfile = fopen(filename,"w");
   // outfile = fopen(filename,"a");
   if (outfile==NULL) {
      fprintf(stderr,"Could not open output file %s\n",filename);
      fflush(stderr);
      exit(0);
   }

   fprintf(outfile,"# Radiance scene description of dla-nd run\n");
   if (type == 0) {
      /* write all cell edges, and particles as lights */
      fprintf(outfile,"void plastic edgec 0 0 5 0.2 0.2 0.2 0.0 0.0\n");
      fprintf(outfile,"void light partc 0 0 3 1000 1000 1000\n");
   } else if (type == 1) {
      /* write stars as plastic objects, and a light source overhead */
      /* fprintf(outfile,"void plastic wallc 0 0 5  0.2 0.2 0.2 0.0 0.0\n"); */
      /* fprintf(outfile,"wallc polygon floor 0 0 12  0 0 0  1 0 0  1 1 0  0 1 0\n"); */
      fprintf(outfile,"void light lightc 0 0 3  10 10 10\n");
      fprintf(outfile,"lightc polygon l1 0 0 12  0 0 2  0 1 2  1 1 2  1 0 2\n");
      fprintf(outfile,"void plastic partc 0 0 5  0.5 0.5 0.5 0.0 0.0\n");
      if (sim->bdry[2][0] == WALL) fprintf(outfile,"partc polygon floor 0 0 12 %g %g %g %g %g %g %g %g %g %g %g %g\n", top->min[0],top->min[1],top->min[2],top->max[0],top->min[1],top->min[2], top->max[0],top->max[1],top->min[2],top->min[0],top->max[1],top->min[2]);
   } else if (type == 2) {
      /* write stars as glow objects */
      fprintf(outfile,"void glow partc 0 0 4  1.0 1.0 1.0 0.0\n");
   }

   write_rad_cell(top,outfile,type);

   /* close file */
   fclose(outfile);

   /* if all went well, return zero */
   return(0);
}


/*
 * Write a Radiance description of just this cell
 */
int write_rad_cell(cell_ptr c,FILE* outfile,int type){

   int i,cnt;
   FLOAT zp1,zp2,rad;
   particle_ptr curr;
   zp2 = 0.0;

   // all branchy-like:
   //FLOAT scale = 0.05;
   //FLOAT exponent = 1./3.;
   // straight-up uniform
   FLOAT scale = 0.4;
   FLOAT exponent = 0.;

   // fprintf(outfile,"\n# cell at level %d\n",c->level);

   if (c->has_subcells) {
      /* recurse for the subcells */
      for (i=0;i<NCHILD;i++)
         write_rad_cell(c->s[i],outfile,type);
   } else {
      if (type == 0) {
         /* write out the edges */
         rad = 0.001;
#if DIM==2
         fprintf(outfile,"edgec cylinder e1 0 0 7 %g %g %g %g %g %g %g\n",c->min[0],c->min[1],0.0,c->max[0],c->min[1],0.0,rad);
         fprintf(outfile,"edgec cylinder e2 0 0 7 %g %g %g %g %g %g %g\n",c->min[0],c->min[1],0.0,c->max[0],c->min[1],0.0,rad);
         fprintf(outfile,"edgec cylinder e3 0 0 7 %g %g %g %g %g %g %g\n",c->min[0],c->max[1],0.0,c->max[0],c->max[1],0.0,rad);
         fprintf(outfile,"edgec cylinder e4 0 0 7 %g %g %g %g %g %g %g\n",c->min[0],c->max[1],0.0,c->max[0],c->max[1],0.0,rad);
         fprintf(outfile,"edgec cylinder e5 0 0 7 %g %g %g %g %g %g %g\n",c->min[0],c->min[1],0.0,c->min[0],c->max[1],0.0,rad);
         fprintf(outfile,"edgec cylinder e6 0 0 7 %g %g %g %g %g %g %g\n",c->min[0],c->min[1],0.0,c->min[0],c->max[1],0.0,rad);
         fprintf(outfile,"edgec cylinder e7 0 0 7 %g %g %g %g %g %g %g\n",c->max[0],c->min[1],0.0,c->max[0],c->max[1],0.0,rad);
         fprintf(outfile,"edgec cylinder e8 0 0 7 %g %g %g %g %g %g %g\n",c->max[0],c->min[1],0.0,c->max[0],c->max[1],0.0,rad);
         fprintf(outfile,"edgec cylinder e9 0 0 7 %g %g %g %g %g %g %g\n",c->min[0],c->min[1],0.0,c->min[0],c->min[1],0.0,rad);
         fprintf(outfile,"edgec cylinder e10 0 0 7 %g %g %g %g %g %g %g\n",c->min[0],c->max[1],0.0,c->min[0],c->max[1],0.0,rad);
         fprintf(outfile,"edgec cylinder e11 0 0 7 %g %g %g %g %g %g %g\n",c->max[0],c->min[1],0.0,c->max[0],c->min[1],0.0,rad);
         fprintf(outfile,"edgec cylinder e12 0 0 7 %g %g %g %g %g %g %g\n",c->max[0],c->max[1],0.0,c->max[0],c->max[1],0.0,rad);
#else
         fprintf(outfile,"edgec cylinder e1 0 0 7 %g %g %g %g %g %g %g\n",c->min[0],c->min[1],c->min[2],c->max[0],c->min[1],c->min[2],rad);
         fprintf(outfile,"edgec cylinder e2 0 0 7 %g %g %g %g %g %g %g\n",c->min[0],c->min[1],c->max[2],c->max[0],c->min[1],c->max[2],rad);
         fprintf(outfile,"edgec cylinder e3 0 0 7 %g %g %g %g %g %g %g\n",c->min[0],c->max[1],c->min[2],c->max[0],c->max[1],c->min[2],rad);
         fprintf(outfile,"edgec cylinder e4 0 0 7 %g %g %g %g %g %g %g\n",c->min[0],c->max[1],c->max[2],c->max[0],c->max[1],c->max[2],rad);
         fprintf(outfile,"edgec cylinder e5 0 0 7 %g %g %g %g %g %g %g\n",c->min[0],c->min[1],c->min[2],c->min[0],c->max[1],c->min[2],rad);
         fprintf(outfile,"edgec cylinder e6 0 0 7 %g %g %g %g %g %g %g\n",c->min[0],c->min[1],c->max[2],c->min[0],c->max[1],c->max[2],rad);
         fprintf(outfile,"edgec cylinder e7 0 0 7 %g %g %g %g %g %g %g\n",c->max[0],c->min[1],c->min[2],c->max[0],c->max[1],c->min[2],rad);
         fprintf(outfile,"edgec cylinder e8 0 0 7 %g %g %g %g %g %g %g\n",c->max[0],c->min[1],c->max[2],c->max[0],c->max[1],c->max[2],rad);
         fprintf(outfile,"edgec cylinder e9 0 0 7 %g %g %g %g %g %g %g\n",c->min[0],c->min[1],c->min[2],c->min[0],c->min[1],c->max[2],rad);
         fprintf(outfile,"edgec cylinder e10 0 0 7 %g %g %g %g %g %g %g\n",c->min[0],c->max[1],c->min[2],c->min[0],c->max[1],c->max[2],rad);
         fprintf(outfile,"edgec cylinder e11 0 0 7 %g %g %g %g %g %g %g\n",c->max[0],c->min[1],c->min[2],c->max[0],c->min[1],c->max[2],rad);
         fprintf(outfile,"edgec cylinder e12 0 0 7 %g %g %g %g %g %g %g\n",c->max[0],c->max[1],c->min[2],c->max[0],c->max[1],c->max[2],rad);
#endif
         fflush(outfile);
      } else if (type >= 1) {
         /* don't write out the edges */
      }

      /* write out the particles */
      curr = c->first;
      cnt = 0;
      while (curr) {
         /* if (type == 0) rad = 0.001*sqrt(curr->mass); */
         /* else if (type == 1) rad = 0.5*sqrt(curr->mass); */
         /* else if (type >= 1) rad = pow(curr->mass,1./3.)/50.; */
         /* rad = pow(curr->mass,1./3.)/50.; */
         // fprintf(outfile,"partc sphere p%d 0 0 4 %g %g %g %g\n",cnt,curr->x[0],curr->x[1],curr->x[2],curr->rad);

         // allow printing radiance files for 2-dimensional objects
         if (DIM < 3) {
            zp1 = 0.0;
            zp2 = 0.0;
         } else {
            zp1 = curr->x[2];
            if (curr->root) zp2 = curr->root->x[2];
         }

         // if we're using cylinders instead of spheres:
         // for 50,000 particles+, use 0.05; for less, try 0.1
         if (type==3) {
            rad = scale*curr->rad*pow(curr->mass,exponent);
            // fprintf(outfile,"partc sphere p%d 0 0 4 %g %g %g %g\n",cnt,curr->x[0],curr->x[1],curr->x[2],rad);
            fprintf(outfile,"partc sphere p%d 0 0 4 %g %g %g %g\n",cnt,curr->x[0],curr->x[1],zp1,rad);
            if (curr->root) {
               // fprintf(outfile,"partc cylinder p%d 0 0 7 %g %g %g %g %g %g %g\n",cnt,curr->x[0],curr->x[1],curr->x[2],curr->root->x[0],curr->root->x[1],curr->root->x[2],rad);
               fprintf(outfile,"partc cylinder p%d 0 0 7 %g %g %g %g %g %g %g\n",cnt,curr->x[0],curr->x[1],zp1,curr->root->x[0],curr->root->x[1],zp2,rad);
            }
         } else {
            // just write the spheres
            // fprintf(outfile,"partc sphere p%d 0 0 4 %g %g %g %g\n",cnt,curr->x[0],curr->x[1],curr->x[2],curr->rad);
            fprintf(outfile,"partc sphere p%d 0 0 4 %g %g %g %g\n",cnt,curr->x[0],curr->x[1],zp1,curr->rad);
         }

         // fflush(outfile);
         cnt++;
         // moved this outside of the loop bc of changes made 2002-09-08
         curr = curr->next;
      }
   }

   /* if all went well, return zero */
   return(0);
}


/*
 * Write a re-readable native particle3d description of the particles
 */
int write_part(cell_ptr top,char *filename){

   FILE *outfile;

   /* open file for writing */
   outfile = fopen(filename,"w");
   if (outfile==NULL) {
      fprintf(stderr,"Could not open output file %s\n",filename);
      fflush(stderr);
      exit(0);
   }

   fprintf(outfile,"# scene description of dla-nd run\n");

   write_part_cell(top,outfile);

   /* close file */
   fclose(outfile);

   /* if all went well, return zero */
   return(0);
}


/*
 * Write a native particle3d description of just this cell
 */
int write_part_cell(cell_ptr c,FILE* outfile){

   int i,cnt;
   particle_ptr curr;

   /* fprintf(outfile,"\n# cell at level %d\n",c->level); */

   if (c->has_subcells) {
      /* recurse for the subcells */
      for (i=0;i<NCHILD;i++)
         write_part_cell(c->s[i],outfile);
   } else {
      /* write out the particles */
      curr = c->first;
      cnt = 0;
      while (curr) {
         // fprintf(outfile,"%g %g %g %g %g %g %g\n",curr->x[0],curr->x[1],curr->x[2],curr->rad,curr->u[0],curr->u[1],curr->u[2]);
         fprintf(outfile,"%g %g %g %g\n",curr->x[0],curr->x[1],curr->x[2],curr->rad);
         /* fflush(outfile); */
         cnt++;
         curr = curr->next;
      }
   }

   /* if all went well, return zero */
   return(0);
}


/*
 * write a png file
 */
int write_png(char *file_name,png_byte** image,int xres,int yres,int depth) {

   png_uint_32 height,width;
   // png_byte image[height][width*bytes_per_pixel];
   // png_bytep row_pointers[yres];
   //int bytes_per_pixel=1;
   FILE *fp;
   png_structp png_ptr;
   png_infop info_ptr;
   //png_colorp palette;
   //png_voidp user_error_ptr;
   // char *file_name = "out.png";

   height=yres;
   width=xres;

   /* open the file */
   fp = fopen(file_name, "wb");
   // fp = stdout;
   if (fp == NULL)
      return (-1);

   /* Create and initialize the png_struct with the desired error handler
    * functions.  If you want to use the default stderr and longjump method,
    * you can supply NULL for the last three parameters.  We also check that
    * the library version is compatible with the one used at compile time,
    * in case we are using dynamically linked libraries.  REQUIRED.
    */
   png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
      NULL, NULL, NULL);
      // user_error_ptr, user_error_fn, user_warning_fn);

   if (png_ptr == NULL)
   {
      fclose(fp);
      return (-1);
   }

   /* Allocate/initialize the image information data.  REQUIRED */
   info_ptr = png_create_info_struct(png_ptr);
   if (info_ptr == NULL)
   {
      fclose(fp);
      png_destroy_write_struct(&png_ptr,(png_infopp)NULL);
      return (-1);
   }

   /* Set error handling.  REQUIRED if you aren't supplying your own
    * error handling functions in the png_create_write_struct() call.
    */
   if (setjmp(png_jmpbuf(png_ptr)))
   {
      /* If we get here, we had a problem reading the file */
      fclose(fp);
      png_destroy_write_struct(&png_ptr, &info_ptr);
      return (-1);
   }

   /* set up the output control if you are using standard C streams */
   png_init_io(png_ptr, fp);

   /* Set the image information here.  Width and height are up to 2^31,
    * bit_depth is one of 1, 2, 4, 8, or 16, but valid values also depend on
    * the color_type selected. color_type is one of PNG_COLOR_TYPE_GRAY,
    * PNG_COLOR_TYPE_GRAY_ALPHA, PNG_COLOR_TYPE_PALETTE, PNG_COLOR_TYPE_RGB,
    * or PNG_COLOR_TYPE_RGB_ALPHA.  interlace is either PNG_INTERLACE_NONE or
    * PNG_INTERLACE_ADAM7, and the compression_type and filter_type MUST
    * currently be PNG_COMPRESSION_TYPE_BASE and PNG_FILTER_TYPE_BASE. REQUIRED
    */
   png_set_IHDR(png_ptr, info_ptr, height, width, depth, PNG_COLOR_TYPE_GRAY,
      PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

   /* Optional gamma chunk is strongly suggested if you have any guess
    * as to the correct gamma of the image.
    */
   png_set_gAMA(png_ptr, info_ptr, 2.2);

   /* Write the file header information.  REQUIRED */
   png_write_info(png_ptr, info_ptr);

   /* The easiest way to write the image (you may have a different memory
    * layout, however, so choose what fits your needs best).  You need to
    * use the first method if you aren't handling interlacing yourself.
    */
   // for (k = 0; k < height; k++)
     // row_pointers[k] = image + k*width*bytes_per_pixel;

   /* One of the following output methods is REQUIRED */
   // png_write_image(png_ptr, row_pointers);
   png_write_image(png_ptr, image);

   /* It is REQUIRED to call this to finish writing the rest of the file */
   png_write_end(png_ptr, info_ptr);

   /* clean up after the write, and free any memory allocated */
   png_destroy_write_struct(&png_ptr, &info_ptr);

   /* close the file */
   fclose(fp);

   /* that's it */
   return (0);
}

