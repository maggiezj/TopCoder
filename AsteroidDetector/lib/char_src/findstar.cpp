#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "findstar.h"
#include "lsquare.h"

#ifdef DEBUG_MEM
#include "checkmem.h"
#endif

#define MIN_PIXELS_TO_FIND 6
#define N_ITER 10
#define MAX_CELL 15
#define SHELLSORT_GAP0 29524

int log_printf( FILE *log_file, const char *format, ...);   /* charon.cpp */

extern FILE *log_file;
extern int debug_level;
extern int saturation_point, max_background;
static int width, height;
static PIXEL *ibuff;

int psf_fit_star( FOUND_STAR *found, int cell_size);       /* findstar.c */

#define SWAP( A, B, TEMP)        { TEMP = A;    A = B;      B = TEMP;  }

static int fix_section( PIXEL *p, int n_vals)
{
   int low = 1, high = n_vals - 1;
   PIXEL pivot, tval;

   if( n_vals > 6)         /* may be worth adjusting pivot in this case... */
      {
      int mid = n_vals / 2, idx = 0;

      if( p[mid] > *p)
         {
         if( p[high] > p[mid])
            idx = mid;
         else if( p[high] > *p)
            idx = high;
         }
      else
         {
         if( p[high] < p[mid])
            idx = mid;
         else if( p[high] < *p)
            idx = high;
         }
      if( idx)
         SWAP( *p, p[idx], tval);
      }

   pivot = *p;
   while( low <= high)
      {
      while( p[low] <= pivot && low <= high)
         low++;
      while( p[high] > pivot && low <= high)
         high--;
      if( high > low)
         {
         SWAP( p[low], p[high], tval);
         low++;
         high--;
         }
      }
   SWAP( p[high], *p, tval);
   return( high);
}

static PIXEL find_nth_value( PIXEL *p, int n_vals, int desired_val)
{
   int i;

   while( n_vals > 1)
      {
      i = fix_section( p, n_vals);
      if( i == desired_val)
         return( p[desired_val]);
      else if( i > desired_val)
         n_vals = i;
      else
         {
         i++;
         p += i;
         n_vals -= i;
         desired_val -= i;
         }
      }
   return( p[desired_val]);
}

static int grab_star( int x0, int y0, FOUND_STAR *found)
{
   const int median_size = 4;
   PIXEL *tptr = ibuff + x0 + y0 * width;
   int i, j, dx, dy, total_above_median = 0;
   long xsum = 0L, ysum = 0L, bright = 0L;
   PIXEL backgnd[8 * 4], *median_ptr, median_val;

// if( tptr[1] >= *tptr || tptr[width] >= *tptr)
   if( tptr[1] >= *tptr || tptr[width] >= *tptr || *tptr >= saturation_point)
      return( -1);

   median_ptr = backgnd;
   for( i = -median_size; i < median_size; i++)
      {
      dx = median_size;
      dy = i;
      for( j = 4; j; j--)
         {
         int temp, x = x0 + dx, y = y0 + dy;

         if( x > 0 && x < width && y > 0 && y < height)
            {
            PIXEL zval = ibuff[x + y * width];

            if( zval > *tptr || !zval)
               return( -1);
            *median_ptr++ = zval;
            }
         temp = dx;
         dx = dy;
         dy = -temp;
         }
      }

   i = median_ptr - backgnd;
#ifdef OBSOLETE_CODE
   if( i < median_size * 2 + 5 || i > 8 * median_size || x0 + y0 == 4)
      printf( "x = %d   y = %d  i = %d\n", x0, y0, i);
#endif
   median_val = find_nth_value( backgnd, i, i / 3);

   for( dx = -2; dx < 3; dx++)
      for( dy = -2; dy < 3; dy++)
         {
         PIXEL zval = tptr[dx + dy * width];

         if( zval > *tptr || zval == (PIXEL)-1)
            return( -1);                  /* can't get a star here */
         if( zval == *tptr && (dx || dy))
            return( -1);                  /* can't get a star here */
         if( zval > median_val)
            {
            total_above_median++;
            zval -= median_val;
            xsum += (long)dx * (long)zval;
            ysum += (long)dy * (long)zval;
            bright += (long)zval;
            }
         }
   if( total_above_median < MIN_PIXELS_TO_FIND)
      return( -1);

   found->x = (double)x0 + (double)xsum / (double)bright;
   found->y = (double)y0 + (double)ysum / (double)bright;
   found->bright = bright;
   return( 0);
}

int psf_fit_star( FOUND_STAR *found, int cell_size)
{
   int iter, half_size = cell_size / 2;
   double i0, r0 = 2., a0;

   found->bright = -1L;
   for( iter = N_ITER; iter; iter--)
      {
      int x0 = (int)found->x + 1, y0 = (int)found->y + 1;
      int i, j, err, n_good_pixels;
      double diffs[5], dx[MAX_CELL], dy[MAX_CELL];
      double fx[MAX_CELL], fy[MAX_CELL], r0_squared;
      void *lsquare = lsquare_init( 5);
      PIXEL *tptr;

      if( x0 < half_size) x0 = half_size;
      if( y0 < half_size) y0 = half_size;
      if( x0 >= width - half_size) x0 = width - half_size - 1;
      if( y0 >= height - half_size) y0 = height - half_size - 1;
      tptr = ibuff + (x0 - half_size) + (y0 - half_size) * width;
      if( iter == N_ITER)           /* on first pass through,  provide... */
         {                          /* ...'decent' starting values */
         PIXEL zval = 0;

         i0 = 0.;          /* assume no backgnd intensity at first... */
         for( i = 0; i < cell_size; i++)
            for( j = 0; j < cell_size; j++)
               if( tptr[i + j * width] > zval)
                  zval = tptr[i + j * width];
         a0 = (double)zval - i0;
         }
      r0_squared = r0 * r0;
      for( i = 0; i < cell_size; i++)
         {
         dx[i] = (double)( x0 + i - half_size) - found->x;
         fx[i] = exp( -dx[i] * dx[i] / r0_squared);
         dy[i] = (double)( y0 + i - half_size) - found->y;
         fy[i] = exp( -dy[i] * dy[i] / r0_squared);
         }
      n_good_pixels = 0;
      for( i = 0; i < cell_size; i++)
         for( j = 0; j < cell_size; j++)
            {
            PIXEL zval = tptr[i + j * width];

            if( zval < saturation_point)
               {
               double exp_val = fx[i] * fy[j], slopes[5], residual;

               residual = i0 + a0 * exp_val - (double)zval;
               slopes[0] = 1.;                         /* slope in i0 */
               slopes[1] = exp_val;                    /* slope in a0 */
               slopes[2] = 2. * a0 * dx[i] * exp_val / r0_squared;
               slopes[3] = 2. * a0 * dy[j] * exp_val / r0_squared;
               slopes[4] = 2. * a0 * (dx[i] * dx[i] + dy[j] * dy[j]) * exp_val /
                                    (r0 * r0_squared);
               lsquare_add_observation( lsquare, -residual, 1., slopes);
               n_good_pixels++;
               }
            }
      err = lsquare_solve( lsquare, diffs);
      lsquare_free( lsquare);
                /* we need at least 6 unsaturated pixels to solve 5 params */
      if( n_good_pixels == cell_size * cell_size && !err)
         {
         i0 += diffs[0];
         a0 += diffs[1];
         for( i = 2; i < 4; i++)
            {
#if 0
            if( diffs[i] > 1.) diffs[i] = 1.;
            if( diffs[i] < -1.) diffs[i] = -1.;
#endif
            if( fabs( diffs[i]) > 1.)
               return( -1);
            }
         found->x += diffs[2];
         found->y += diffs[3];
         if( diffs[4] > r0 / 10.)
            diffs[4] = r0 / 10.;
         if( diffs[4] < -r0 / 10.)
            diffs[4] = -r0 / 10.;
         r0 += diffs[4];
         }
      else
         return( -1);
      }
   if( max_background && i0 > (double)max_background)
      return( -1);
   found->bright = (long)( a0 * 100.);
   return( 0);
}

int find_stars( PIXEL *image, int xsize, int ysize,
            PIXEL limit, FOUND_STAR *found, int n_max, int cell_size)
{
   int rval = 0, i, j, k, discard_it, gap, counter = 0;
   FOUND_STAR tstar;
   extern int psf_fitting_on;
   extern int ignored_north, ignored_south, ignored_east, ignored_west;
   int x1 = ignored_west + cell_size / 2;
   int x2 = xsize - (ignored_west + cell_size / 2);
   int y1 = ignored_north + cell_size / 2;
   int y2 = ysize - (ignored_south + cell_size / 2);

   ibuff = image;
   width = xsize;
   height = ysize;
   if( debug_level)
      log_printf( log_file, "find_stars: ");
                      /* We go through the entire image, looking for stars */
                      /* using the grab_star( ) function.  To evade checking */
                      /* the same star twice,  we require that the center   */
                      /* pixel passed to grab_star( ) be the brightest in   */
                      /* the cell.  By comparing the center pixel to four of */
                      /* its neighbors,  and by insisting that it be above  */
                      /* "limit",  we can immediately reject a lot of cases */
                      /* without even having to do a function call.  Since  */
                      /* grab_star( ) is slow,  this is a Good Thing.       */
   for( j = y1; j < y2; j++)
      {
      PIXEL *tptr = image + x1 + xsize * j;

      for( i = x1; i < x2; i++, tptr++)
         if( *tptr > limit)
            if( tptr[-xsize] < *tptr && tptr[-1] < *tptr &&
                tptr[ xsize] < *tptr && tptr[ 1] < *tptr)
               if( !grab_star( i, j, &tstar))
                  {
                  extern int image_pixel_separation;
                  double sep = (double)image_pixel_separation;

                  discard_it = 0;
                  for( k = 0; k < rval; k++)
                     if( fabs( tstar.x - found[k].x) < sep &&
                         fabs( tstar.y - found[k].y) < sep)
                        discard_it = 1;
                  if( !discard_it)
                     {
                     for( k = 0; k < rval && found[k].bright > tstar.bright; k++)
                        ;
                     memmove( found + k + 1, found + k,
                               (rval - k) * sizeof( FOUND_STAR));
                     if( k < n_max)
                        found[k] = tstar;
                     if( rval < n_max)
                        rval++;
                     }
                  }
      if( debug_level > 1)
         {
         counter += 40;
         while( counter >= y2 - y1)
            {
            log_printf( log_file, ".");
            counter -= y2 - y1;
            }
         }
      }

   if( debug_level)
      log_printf( log_file, "found, ");
                  /* If we're not just finding stars by centroid,  but */
                  /* are also PSF-fing them,  then we go through one by */
                  /* one and apply psf_fit_star() to them.  In some cases, */
                  /* this tells us the star isn't PSF-able,  in which case */
                  /* it's discarded from the list.                         */
   if( psf_fitting_on)
      for( i = 0; i < rval; i++)
         if( psf_fit_star( found + i, cell_size))
            {
            rval--;
            memmove( found + i, found + i + 1,
                          (rval - i) * sizeof( FOUND_STAR));
            }

                        /* Now ShellSort the stars found in order of mag: */
   for( gap = SHELLSORT_GAP0; gap; gap /= 3)
      for( i = 0; i < gap; i++)
         for( j = i; j + gap < rval; j += gap)
            if( found[j].bright < found[j + gap].bright)
               {
               FOUND_STAR temp = found[j];

               found[j] = found[j + gap];
               found[j + gap] = temp;
               if( j >= gap)
                  j -= (gap << 1);
               }
   if( debug_level)
      log_printf( log_file, "sorted, ");

   for( i = 0; i < rval; i++)          /* move stars to _center_ of pixel */
      {
      found[i].x += .5;
      found[i].y += .5;
      }
   if( debug_level)
      log_printf( log_file, "done\n");
   return( rval);
}
