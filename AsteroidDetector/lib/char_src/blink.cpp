#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "watdefs.h"
#include "findstar.h"
#include "charon.h"
#include "miscell.h"

#ifdef DEBUG_MEM
#include "checkmem.h"
#endif

IMAGE *get_blinkable_images( const IMAGE *ref_image, int *n_found,
           const double xfactor, const double yfactor);     /* blink.h */
int dark_subtract_and_flat_field( IMAGE *image, const IMAGE *dark,
                           const IMAGE *flat);           /* blink.cpp */
int add_wcs_header_data( const IMAGE *img);              /* blink.cpp */

#define DPT struct dpt

DPT
   {
   double x, y;
   };

PIXEL *get_blink_image( const char *iline, const IMAGE *ref_image,
           int *err_code, const int check_only,
           const double xfactor, const double yfactor)
{
   int i, loc = 0, nbytes, rval = 0;
   double corners[8], xmax, xmin, ymax, ymin;
   unsigned low_end, high_end;
   int format, xsize, ysize;
   long data_offset;
   const char *filename;
   IMAGE temp_image;
   extern int print_to_screen;

   for( i = 0; i < 4 && !rval; i++)
      {
      double ra, dec, xy[2];

      sscanf( iline + loc, "%lf %lf%n", &ra, &dec, &nbytes);
      loc += nbytes;
      ra_dec_to_pixel( xy, ra, dec, ref_image->xform);
      xy[0] /= xfactor;
      xy[1] /= yfactor;
                       /* Ensure truncation works correctly by adding one: */
//    xy[0]--;
//    xy[1]--;
      corners[  i + i  ] = xy[0];
      corners[i + i + 1] = xy[1];
      if( !i || xmax < xy[0])
         xmax = xy[0];
      if( !i || xmin > xy[0])
         xmin = xy[0];
      if( !i || ymax < xy[1])
         ymax = xy[1];
      if( !i || ymin > xy[1])
         ymin = xy[1];
      }

                  /* Check against rectangle for a quick reject: */
   if( !rval)
      {
      if( ymin > (double)ref_image->ysize || ymax < 0.)
         rval = -2;
      if( xmin > (double)ref_image->xsize || xmax < 0.)
         rval = -2;
      }

   sscanf( iline + loc, "%u %u %d %d %ld %d%n", &low_end,
               &high_end, &xsize, &ysize, &data_offset,
               &format, &nbytes);
   if( ref_image->xsize != xsize || ref_image->ysize != ysize)
      rval = -6;
   loc += nbytes;
   while( iline[loc] == ' ' && iline[loc])
      loc++;
   if( iline[loc] == '!')         /* from CD drive */
      rval = -4;

   filename = iline + loc;
   if( !strcmp( filename, ref_image->filename))
      rval = -9;        /* avoid "own goals" */

   if( !rval && check_only)
      {
      FILE *ifile = fopen( filename, "rb");

      if( !ifile)
         rval = -5;
      else
         fclose( ifile);
      }

   if( rval || check_only)   /* image doesn't work,  or don't wanna load it */
      {
      *err_code = rval;
      return( NULL);
      }

   memcpy( &temp_image, ref_image, sizeof( IMAGE));
   print_to_screen = 0;
   if( !rval && load_image( filename, &temp_image))
      {
      *err_code = -7;
      return( NULL);
      }
   print_to_screen = 1;

   if( temp_image.is_inverted)
      flip_image( &temp_image);

   *err_code = rval;
   free( temp_image.filename);
   return( temp_image.img);
}

IMAGE *get_blinkable_images( const IMAGE *ref_image, int *n_found,
           const double xfactor, const double yfactor)
{
   FILE *ifile = fopen( "images.dat", "rb");
   IMAGE *rval = NULL;
   char ibuff[200];
   int curr_line_no = 0, last_linear_line = -999;
   double xform[N_FIND_PARAMS];

   *n_found = 0;
   if( !ifile)
      return( NULL);
   while( fgets( ibuff, 200, ifile))
      {
      curr_line_no++;
      if( *ibuff == ' ' && curr_line_no == last_linear_line + 1)
         {
         int i, err_code;
         PIXEL *ret_img;

         for( i = 0; ibuff[i] != 10 && ibuff[i] != 13; i++)
            ;
         ibuff[i] = 0;
         ret_img = get_blink_image( ibuff, ref_image, &err_code, 0,
                                       xfactor, yfactor);
         if( ret_img && !err_code)
            {
            IMAGE *rval2 = (IMAGE *)calloc( *n_found + 1, sizeof( IMAGE));

            if( rval)
               {
               memcpy( rval2, rval, *n_found * sizeof( IMAGE));
               free( rval);
               }
            rval = rval2;
                     /* Start out with images the same: */
            memcpy( rval + *n_found, ref_image, sizeof( IMAGE));

            rval[*n_found].filename = (char *)malloc( strlen( ibuff + 114) + 1);
            strcpy( rval[*n_found].filename, ibuff + 114);
            rval[*n_found].img = ret_img;
            memcpy( rval[*n_found].xform, xform, N_FIND_PARAMS * sizeof( double));
            (*n_found)++;
            }
         }
      if( *ibuff == 'L')
         {
         int i, loc = 1, bytes_scanned;
         const int linear_terms[9] = { 3, 4, 5, 6, 7, 5 + N_FIND_COEFFS,
                         6 + N_FIND_COEFFS, 7 + N_FIND_COEFFS, 1 };

         last_linear_line = curr_line_no;
         for( i = 0; i < N_FIND_PARAMS; i++)
            xform[i] = 0.;
         for( i = 0; i < 9; i++, loc += bytes_scanned)
            {
            sscanf( ibuff + loc, "%lf%n", xform + linear_terms[i],
                                                 &bytes_scanned);
            if( i >= 2)
               xform[linear_terms[i]] *= 1.e-6;
            }
         xform[0] = 1.;
         }

      if( *ibuff == 'Q')
         {
         int i, j, loc = 1, bytes_scanned;

         last_linear_line = curr_line_no;
         for( j = 0; j <= N_FIND_COEFFS; j += N_FIND_COEFFS)
            for( i = 8; i < 11; i++, loc += bytes_scanned)
               {
               sscanf( ibuff + loc, "%lf%n", xform + i + j, &bytes_scanned);
               xform[i + j] *= 1.e-9;
               }
         xform[0] = 2.;
         }

      if( *ibuff == 'C')
         {
         int i, j, loc = 1, bytes_scanned;

         last_linear_line = curr_line_no;
         for( j = 0; j <= N_FIND_COEFFS; j += N_FIND_COEFFS)
            for( i = 11; i < 15; i++, loc += bytes_scanned)
               {
               sscanf( ibuff + loc, "%lf%n", xform + i + j, &bytes_scanned);
               xform[i + j] *= 1.e-12;
               }
         xform[0] = 3.;
         }
      }
   fclose( ifile);
   return( rval);
}

static void max_min_diff( const IMAGE *i1, const IMAGE *i2,
            long *max_diff, long *min_diff)
{
   long i, n_pixels = i1->xsize * i1->ysize, delta;

   for( i = 0; i < n_pixels; i++)
      {
      if( i2)
         delta = (long)i1->img[i] - (long)i2->img[i];
      else
         delta = (long)i1->img[i];
      if( max_diff)
         if( !i || *max_diff < delta)
            *max_diff = delta;
      if( min_diff)
         if( !i || *min_diff > delta)
            *min_diff = delta;
      }
}

int dark_subtract_and_flat_field( IMAGE *image, const IMAGE *dark,
                           const IMAGE *flat)
{
   long max_image_minus_dark, min_image_minus_dark;
   long max_image;
   long i, n_pixels = image->xsize * image->ysize;

   max_min_diff( image, dark, &max_image_minus_dark, &min_image_minus_dark);
   printf( "image-dark runs from %ld to %ld\n", min_image_minus_dark,
                                                max_image_minus_dark);
// getch( );
          /* The following may not be true right now,  but it will be after */
          /* the dark is subtracted:                                        */
   max_image = max_image_minus_dark - min_image_minus_dark;
          /* There is some risk,  though not much,  that the subtraction    */
          /* will require losing the least significant bit.                 */
   if( max_image < 0xffff)
      {
      for( i = 0; i < n_pixels; i++)
         image->img[i] -= dark->img[i] + min_image_minus_dark;
      }
   else           /* guess we gotta lose that bit... */
      {
      max_image >>= 1;
      for( i = 0; i < n_pixels; i++)
         {
         unsigned long tval = (unsigned long)
            ((long)image->img[i] - (long)dark->img[i] - min_image_minus_dark);

         image->img[i] = (PIXEL)( tval >> 1);
         }
      }

   if( flat)
      {
      long max_flat_minus_dark;

      max_min_diff( flat, dark, &max_flat_minus_dark, NULL);
      for( i = 0; i < n_pixels; i++)
         {
         long flat_minus_dark = (long)flat->img[i] - (long)dark->img[i];

                  /* The 'flat' pixel is supposed to be brighter than the */
                  /* 'dark' one.  If it isn't... run in circles...        */
         if( flat_minus_dark < 1)
            image->img[i] = 0;
         else
            {
            unsigned long tval = (unsigned long)image->img[i] *
                                 (unsigned long)max_flat_minus_dark /
                                 (unsigned long)flat_minus_dark;

            image->img[i] = (PIXEL)tval;
            }
         }
      }
   return( 0);
}

static PIXEL convolve_at_edge( const IMAGE *image, const int x, const int y,
             const int kernel_size, const double *kernel, const double normal)
{
   int i, j, khalf = kernel_size / 2;
   PIXEL rval, *tptr = image->img + (x - khalf)
                                  + (y - khalf) * image->xsize;
   double normal2 = 0., sum = 0.;

   for( j = y - khalf; j <= y + khalf; j++, tptr += image->xsize - kernel_size)
      for( i = x - khalf; i <= x + khalf; i++, kernel++, tptr++)
         if( j >= 0 && j < image->ysize && i >= 0 && i < image->xsize)
            {
            sum += *kernel * (double)*tptr;
            normal2 += *kernel;
            }
   sum *= normal / normal2;
   if( sum < 0.)
      rval = 0;
   else if( sum >= (double)MAX_PIXEL)
      rval = MAX_PIXEL;
   else
      rval = (PIXEL)sum;
   return( rval);
}

int apply_convolution( IMAGE *image, const int kernel_size, const double *kernel)
{
   int k2 = kernel_size * kernel_size, khalf = kernel_size / 2;
   int lines_buffered = (kernel_size + 1) / 2;
   PIXEL *tbuff = (PIXEL *)calloc( lines_buffered * image->xsize, sizeof( PIXEL));
   double normal = 0.;
   int x, y, i, j;

   for( i = 0; i < k2; i++)
      normal += kernel[i];

   for( y = 0; y < image->ysize; y++)
      {
      PIXEL *tptr = tbuff + (y % lines_buffered) * image->xsize;

      if( y >= lines_buffered)         /* copy results back into image */
         memcpy( image->img + (y - lines_buffered) * image->xsize, tptr,
                                          image->xsize * sizeof( PIXEL));


                  /* The first few lines and the last few will be partly */
                  /* "off edge",  and need to use convolve_at_edge:      */
      if( y < khalf || y >= image->ysize - khalf)
         for( x = 0; x < image->xsize; x++)
            tptr[x] = convolve_at_edge( image, x, y, kernel_size, kernel, normal);
      else           /* "normal",  internal line: */
         {
         PIXEL *image_ptr0 = image->img + (y - khalf) * image->xsize - khalf;

         for( x = 0; x < khalf; x++)
            tptr[x] = convolve_at_edge( image, x, y, kernel_size, kernel, normal);
         for( x = image->xsize - khalf; x < image->xsize; x++)
            tptr[x] = convolve_at_edge( image, x, y, kernel_size, kernel, normal);
         for( x = khalf; x < image->xsize - khalf; x++)
            {
            double sum = 0.;
            PIXEL *image_ptr = image_ptr0 + x;
            const double *kernel_ptr = kernel;

            for( j = kernel_size; j; j--, image_ptr += image->xsize - kernel_size)
               for( i = kernel_size; i; i--, kernel_ptr++, image_ptr++)
//                if( *kernel_ptr)
                     sum += *kernel_ptr * (double)*image_ptr;
            if( sum < 0)
               tptr[x] = 0;
            else if( sum >= MAX_PIXEL)
               tptr[x] = MAX_PIXEL;
            else
               tptr[x] = (PIXEL)sum;
            }
         }
      }

                  /* copy final results: */
   for( y = image->ysize; y < image->ysize + lines_buffered; y++)
      {
      PIXEL *tptr = tbuff + (y % lines_buffered) * image->xsize;

      memcpy( image->img + (y - lines_buffered), tptr,
                                             image->xsize * sizeof( PIXEL));
      }

   free( tbuff);
   return( 0);
}

int load_kernel_file( const char *filename, double *kernel)
{
   FILE *ifile = fopen( filename, "rb");
   int rval, i;
   char buff[80];
   double scale;

   if( !ifile)
      return( -1);
   fgets( buff, sizeof( buff), ifile);
   if( memcmp( buff, "CONV NORM", 9))
      {
      fclose( ifile);
      return( -2);
      }
   scale = atof( buff + 9);
   if( !scale)
      scale = 1.;
   fgets( buff, sizeof( buff), ifile);
   rval = atoi( buff + 1);
   if( kernel)
      for( i = 0; i < rval * rval; i++)
         {
         fscanf( ifile, "%lf", kernel + i);
         kernel[i] /= scale;
         }
   fclose( ifile);
   return( rval);
}

static void fill_image_polygon( IMAGE *image, const int npts,
                const int *xy, const PIXEL color)
{
   int y, i;

   for( y = 0; y < image->ysize; y++)
      {
      int x1, y1, x2 = xy[npts * 2 - 2], y2 = xy[npts * 2 - 1];
      int intercepts[20], n_intercepts = 0;

      for( i = 0; i < npts; i++)
         {
         x1 = xy[i + i];
         y1 = xy[i + i + 1];
         if( y1 <= y && y2 > y || y2 <= y && y1 > y)
            {
            int temp = x1 + (int)
                  ( (long)(x2 - x1) * (long)(y - y1) / (long)(y2 - y1));

            if( temp < 0)
               temp = 0;
            if( temp > image->xsize)
               temp = image->xsize;
            intercepts[n_intercepts++] = temp;
            }
         x2 = x1;
         y2 = y1;
         }
      for( i = 0; i < n_intercepts - 1; i++)
         if( intercepts[i] > intercepts[i + 1])
            {
            int temp = intercepts[i + 1];

            intercepts[i + 1] = intercepts[i];
            intercepts[i] = temp;
            if( i)
               i -= 2;
            }
      for( i = 0; i < n_intercepts; i += 2)
         {
         int j;
         PIXEL *tptr = image->img + y * image->xsize + intercepts[i];

         for( j = intercepts[i + 1] - intercepts[i]; j; j--)
            *tptr++ = color;
         }
      }
}

int apply_polygon_mask( IMAGE *image, const char *mask_filename)
{
   FILE *ifile = fopen( mask_filename, "rb");
   char buff[90];
   int rval = 0;

   if( !ifile)
      rval = -1;
   while( !rval && fgets( buff, sizeof( buff), ifile))
      switch( *buff)
         {
         case 'r':         /* fill a rectangle */
         case 'e':         /* fill an ellipse */
            {
            int xmin, xmax, ymin, ymax, x1, x2, y;

            sscanf( buff + 1, "%d %d %d %d", &xmin, &xmax, &ymin, &ymax);
            for( y = ymin; y < ymax; y++)
               {
               if( *buff == 'r')
                  {
                  x1 = xmin;
                  x2 = xmax;
                  }
               else if( *buff == 'e')
                  {
                  int dy1 = ymax - ymin, dy2 = ymax + ymin - 2 * y;
                  double fraction = (double)dy2 / (double)dy1;
                  double width = (double)( xmax - xmin)
                             * sqrt( 1. - fraction * fraction);

                  x1 = (xmax + xmin) / 2 - (int)( width * .5);
                  x2 = x1 + (int)width;
                  }
               if( x1 < 0)
                  x1 = 0;
               if( x2 > image->xsize)
                  x2 = image->xsize;
               memset( image->img + y * image->xsize + x1, 0,
                                 (x2 - x1) * sizeof( PIXEL));
               }
            }
            break;
         case 'p':
            {
            int i, npts = atoi( buff + 1), xy[50];

            for( i = 0; i < npts; i++)
               {
               fgets( buff, sizeof( buff), ifile);
               sscanf( buff, "%d %d", xy + i + i, xy + i + i + 1);
               }
            fill_image_polygon( image, npts, xy, 0);
            }
            break;
         case '#':         /* comment line */
            break;
         default:
            rval = -2;
            break;
         }

   if( ifile)
      fclose( ifile);
   return( rval);
}

double get_pixel_dimension( const double *params, const int true_if_height);

#define PI 3.1415926535897932384626433

extern FILE *log_file;

int add_wcs_header_data( const IMAGE *img)
{
   int rval = -1;

   if( img->image_type == 'f')
      {
      FILE *ifile = fopen( img->filename, "rb");
      int i, j, n_lines = img->data_offset / 80, orig_lines, lines_out;
      char *header = (char *)calloc( 80, n_lines + 72), *tptr;
      double crota = img->tilt_angle * 180. / PI, crpix[2];
      double pixel_xsize_in_degrees =
                        img->pixel_xsize / (img->focal_len * 1000.);
      double pixel_ysize_in_degrees =
                        img->pixel_ysize / (img->focal_len * 1000.);
#if 0
      double pixel_xsize_in_degrees =
                      get_pixel_dimension( img->xform, 0) * 180. / PI;
      double pixel_ysize_in_degrees =
                      get_pixel_dimension( img->xform, 1) * 180. / PI;
#endif
      double xfactor = img->pixel_xsize / (img->focal_len * 1000.);
      double yfactor = img->pixel_ysize / (img->focal_len * 1000.);

      if( !img->is_inverted)
         pixel_ysize_in_degrees *= -1.;
      pixel_xsize_in_degrees /= (3600. * get_pixel_dimension( img->xform, 0));
      pixel_ysize_in_degrees /= (3600. * get_pixel_dimension( img->xform, 0));
      log_printf( log_file, "%d lines\n", n_lines);
      if( !ifile)
         return( -2);
      if( crota < 0.)
         crota += 360.;
      fread( header, 80, n_lines, ifile);
      orig_lines = n_lines;
      for( i = 0; i < n_lines && memcmp( header + i * 80, "END ", 4); i++)
         ;
      n_lines = i;
      log_printf( log_file, "Header read; %d real lines\n", n_lines);
      for( i = 0; i < n_lines; i++)
         {
         int remove_line = 0;
         const char *reserved_strings[7] = {
                  "CDELT", "CROTA", "CRVAL", "CRPIX", "CTYPE",
                  "CUNIT", NULL };

         tptr = header + i * 80;
         for( j = 0; reserved_strings[j]; j++)
            if( !memcmp( tptr, reserved_strings[j],
                                        strlen( reserved_strings[j])))
               remove_line = 1;
         if( remove_line)
            {
            n_lines--;
            memmove( tptr, tptr + 80, (n_lines - i) * 80);
            i--;
            }
         }
      log_printf( log_file, "%d lines left\n", n_lines);
      tptr = header + n_lines * 80;
      memset( tptr, ' ', 36 * 80);
                     /* x/yfactors are in radians/pix; cvt to arcsec/pix */
      xfactor *= (180. / PI) * 3600.;
      yfactor *= (180. / PI) * 3600.;
      ra_dec_to_pixel( crpix, img->xform[3], img->xform[4], img->xform);
                  /* crpix is in arcseconds;  convert to actual pixels: */
      crpix[0] /= xfactor;
      crpix[1] /= yfactor;
      for( i = 1; i <= 2; i++)
         {
         sprintf( tptr, "CTYPE%d  = 'DEC--TAN'           / Coordinate type", i);
         log_printf( log_file, "%s\n", tptr);
         tptr += 80;
         sprintf( tptr, "CRPIX%d  = %20lf / Pixel coordinate of %c-reference",
                   i, crpix[i - 1],
                   (i == 2) ? 'y' : 'x');
         log_printf( log_file, "%s\n", tptr);
         tptr += 80;
         sprintf( tptr, "CRVAL%d  = %20lf / Coordinate at %c-reference",
                  i, img->xform[i + 2],
                   (i == 2) ? 'y' : 'x');
         log_printf( log_file, "%s\n", tptr);
         tptr += 80;
         sprintf( tptr, "CDELT%d  = %20.10lf / Image scale on %c-axis, deg per pixel",
               i, (i == 2 ? pixel_ysize_in_degrees : pixel_xsize_in_degrees),
               (i == 2) ? 'y' : 'x');
         log_printf( log_file, "%s\n", tptr);
         tptr += 80;
         sprintf( tptr, "CROTA%d  = %20lf / Rotation of coordinate",
                  i, 360. - crota,
                  (i == 2) ? 'y' : 'x');
         log_printf( log_file, "%s\n", tptr);
         tptr += 80;
         sprintf( tptr, "CUNIT%d  = 'deg     '           / Unit of coordinate",
                  i);
         log_printf( log_file, "%s\n", tptr);
         tptr += 80;
         n_lines += 6;
         }
      strcpy( tptr, "END");
      n_lines++;
      for( i = 0; i < n_lines * 80; i++)
         if( !header[i])
            header[i] = ' ';
      lines_out = ((n_lines + 35) / 36) * 36;
      log_printf( log_file, "%d lines out\n", lines_out);
      memset( header + n_lines * 80, ' ', (lines_out - n_lines) * 80);
      log_printf( log_file, "Memory set\n");
      if( (n_lines + 35) / 36 == (orig_lines + 35) / 36)
         {
         FILE *ofile;

         fclose( ifile);
         ofile = fopen( img->filename, "r+b");
         if( !ofile)
            return( -3);
         fwrite( header, lines_out, 80, ofile);
         fclose( ofile);
         }
      else
         {
         FILE *ofile = fopen( "tempfit.ick", "wb");
         int bytes_read;

         if( !ofile)
            return( -4);
         fwrite( header, lines_out, 80, ofile);
         fseek( ifile, img->data_offset, SEEK_SET);
         while( bytes_read = fread( header, 1, img->data_offset, ifile))
            fwrite( header, 1, bytes_read, ofile);
         fclose( ofile);
         fclose( ifile);
         }
      free( header);
      rval = 0;
      }
   return( rval);
}
