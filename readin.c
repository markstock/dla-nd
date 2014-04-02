/*
 * readin.c - routines for png input (from vic2d)
 *
 * Copyright 2004-8 Mark J. Stock mstock@umich.edu
 */

#include "structs.h"


/*
 * read a PNG, write it to 1 or 3 channels
 */
int read_png (char *infile, int *nx, int *ny,
   int expect_three_channel,
   FLOAT ***red, FLOAT redmin, FLOAT redrange,
   FLOAT ***grn, FLOAT grnmin, FLOAT grnrange,
   FLOAT ***blu, FLOAT blumin, FLOAT blurange) {

   int autorange = FALSE;
   int high_depth;
   int three_channel;
   int i,j,printval; //,bit_depth,color_type,interlace_type;
   FLOAT newminrange,newmaxrange;
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
      png_set_gray_1_2_4_to_8(png_ptr);

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
