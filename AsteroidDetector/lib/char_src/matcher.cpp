#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "findstar.h"
#include "lsquare.h"
#include "charon.h"
#include "miscell.h"

#ifdef __WATCOMC__
#include <conio.h>
#endif

#ifdef DEBUG_MEM
#include "checkmem.h"
#endif

#define THREEBYTE_CONVERSION ( 16777216. / 360.)
#define N_PASSES 6
#define PI 3.14159265358979323
#define STAR_PAIR struct star_pair

STAR_PAIR
   {
   int start, end;
   double dist;
   };

extern char *strings[];
extern int use_input_mags;
#ifdef INCLUDE_DEBUG_STATEMENTS
extern int debug_level;
static long debug_times[10];
#endif

static void deproject_star( REF_STAR *star, const double x0, const double y0);
double get_pixel_dimension( const double *params, const int true_if_height);

static int ref_star_compare( const void *entry1, const void *entry2)
{
   double diff = ((REF_STAR *)entry2)->ra - ((REF_STAR *)entry1)->ra;

   if( diff)
      return( (diff > 0.) ? -1 : 1);
   else
      return( 0);
}

static int pair_compare( const void *entry1, const void *entry2)
{
   double diff = ((STAR_PAIR *)entry2)->dist - ((STAR_PAIR *)entry1)->dist;

   if( diff)
      return( (diff > 0.) ? -1 : 1);
   else
      return( 0);
}

static void stuff_ang_into_buff( char *buff, const double val)
{
   long z;
   double dval;

   dval = fmod( (double)val + 180., 360.);
   z = (long)( dval * THREEBYTE_CONVERSION);
   memcpy( buff, &z, 3);
}

static void compute_powers( double *ret_arr, const double x, const double y,
                                             const int order)
{
   int i, j, k;

   for( i = 0; i <= order; i++)
      for( j = 0; j <= i; j++)
         {
         double coeff = 1.;

         for( k = 0; k < j; k++)
            coeff *= y;
         for( k = j; k < i; k++)
            coeff *= x;
         *ret_arr++ = coeff;
         }
               /* the following was great back when only linear,  quad,
                  and cubic fits were supported.   */
#ifdef OLDER_SIMPLER_CODE
   *ret_arr++ = 1.;
   *ret_arr++ = x;
   *ret_arr++ = y;
   if( order > 1)
      {
      *ret_arr++ = x * x;
      *ret_arr++ = x * y;
      *ret_arr++ = y * y;
      }
   if( order > 2)
      {
      *ret_arr++ = x * x * x;
      *ret_arr++ = x * x * y;
      *ret_arr++ = x * y * y;
      *ret_arr++ = y * y * y;
      }
#endif
}


         /* 66 coefficients are enough for a 10th-order fit: */
#define XPARAMS( n) params[(n)+5]
#define YPARAMS( n) params[(n)+5+N_FIND_COEFFS]

static void calc_xform( double *params,
           const FOUND_STAR *p1, const FOUND_STAR *p2,
           const REF_STAR *q1, const REF_STAR *q2)
{
   double dx = p2->x - p1->x, dy = p2->y - p1->y;
   double du = q2->ra - q1->ra, dv = q2->dec - q1->dec;
   double dist;

   dist = dx * dx + dy * dy;
   YPARAMS(2) = (dy * dv - dx * du) / dist;
   YPARAMS(1) = (dy * du + dx * dv) / dist;
   XPARAMS(1) = -YPARAMS(2);
   XPARAMS(2) =  YPARAMS(1);
   XPARAMS(0) = q1->ra - (p1->x * XPARAMS(1) + p1->y * XPARAMS(2));
   YPARAMS(0) = q1->dec - (p1->x * YPARAMS(1) + p1->y * YPARAMS(2));
   params[1] = 0.;         /* zero out magnitude factor */
}

static void calc_high_order_offset( const double *params, const FOUND_STAR *p,
                                     double *loc)
{
   double powers[N_FIND_COEFFS];
   int i, order = (int)( *params + .5);
   int n_coeffs = (order + 1) * (order + 2) / 2;

   compute_powers( powers, p->x, p->y, order);
   for( i = 3; i < n_coeffs; i++)         /* skip first three,  low-order */
      {
      loc[0] += XPARAMS(i) * powers[i];
      loc[1] += YPARAMS(i) * powers[i];
      }
}

static void xform_point( const double *params, const FOUND_STAR *p, REF_STAR *q)
{
   double loc[2];

   loc[0] = XPARAMS(0) + XPARAMS(1) * p->x + XPARAMS(2) * p->y;
   loc[1] = YPARAMS(0) + YPARAMS(1) * p->x + YPARAMS(2) * p->y;
   if( *params > 1.5)            /* higher-order terms */
      calc_high_order_offset( params, p, loc);
   q->ra = loc[0];
   q->dec = loc[1];
}

double get_pixel_dimension( const double *params, const int true_if_height)
{
   params += true_if_height;
   return( sqrt( XPARAMS(1)*XPARAMS(1)+YPARAMS(1)*YPARAMS(1)) * (PI / 180.));
}

static void de_xform_point( const double *params, const REF_STAR *q,
                                             FOUND_STAR *p)
{
   double x = q->ra - XPARAMS(0);
   double y = q->dec - YPARAMS(0);
   double determ = XPARAMS(1) * YPARAMS(2) - XPARAMS(2) * YPARAMS(1);
   int iteration;

   p->x = (x * YPARAMS(2) - y * XPARAMS(2)) / determ;
   p->y = (y * XPARAMS(1) - x * YPARAMS(1)) / determ;
   if( *params > 1.5)            /* higher-order terms */
      for( iteration = 0; iteration < 15; iteration++)
         {
         REF_STAR q1;
         double dx, dy;

         xform_point( params, p, &q1);
         dx = q1.ra - q->ra;
         dy = q1.dec - q->dec;
         p->x -= (dx * YPARAMS(2) - dy * XPARAMS(2)) / determ;
         p->y -= (dy * XPARAMS(1) - dx * YPARAMS(1)) / determ;
         }
}

static double find_closest( const int n_stars, const REF_STAR *stars,
                         const REF_STAR *search,
                         int *ret_idx, const double search_dist)
{
   int i, step, loc = 0, loc1;
   double dx, dy, best_dist = search_dist;
   static int step0;

   if( !step0)
      for( step0 = 1; (step0 << 1) < n_stars; step0 <<= 1)
         ;
   for( step = step0; step; step >>= 1)
      if( (loc1 = (loc | step)) < n_stars && stars[loc1].ra < search->ra)
         loc = loc1;

   for( i = loc; i >= 0 && (dx = search->ra - stars[i].ra) < best_dist; i--)
      {
      dy = search->dec - stars[i].dec;

      if( dy < best_dist && -dy < best_dist)
         {
         double dist;

         dist = sqrt( (double)dx * (double)dx + (double)dy * (double)dy);
         if( dist < best_dist)
            {
            best_dist = dist;
            *ret_idx = i;
            }
         }
      }

   for( i = loc+1; i < n_stars && (dx = stars[i].ra-search->ra) <best_dist; i++)
      {
      dy = search->dec - stars[i].dec;

      if( dy < best_dist && -dy < best_dist)
         {
         double dist;

         dist = sqrt( (double)dx * (double)dx + (double)dy * (double)dy);
         if( dist < best_dist)
            {
            best_dist = dist;
            *ret_idx = i;
            }
         }
      }
   return( best_dist);
}

static void project_star( REF_STAR *star, const double x0, const double y0)
{
   double dtheta = (star->ra - x0) * (PI / 180.);
   double phi = star->dec * (PI / 180.), phi0;
   double x, y, z, cos_phi0, sin_phi0, cos_phi, sin_phi;
   extern int is_drift_scan;

   phi0 = y0 * (PI / 180.);           /* now we're all in radians */
   cos_phi0 = cos( phi0);
   if( is_drift_scan)
      {
      star->ra = (star->ra - x0) * cos_phi0;
      star->dec -= y0;
      }
   else
      {
      sin_phi0 = sin( phi0);
      cos_phi = cos( phi);
      sin_phi = sin( phi);
      x = sin( dtheta) * cos_phi;
      y = cos_phi0 * sin_phi - cos_phi *  cos( dtheta) * sin_phi0;
      z = sqrt( 1. - x * x - y * y);
#ifdef STEREOGRAPHIC_INSTEAD_OF_TANGENT_PLANE
      z = 2. / (1. + z);
#else
      z = 1. / z;
#endif
      star->ra = x * z * (180. / PI);
      star->dec = y * z * (180. / PI);
      }
}

int ra_dec_to_pixel( double *pixel, const double ra, const double dec,
                                       const double *params)
{
   REF_STAR star;
   FOUND_STAR pixel_loc;

   star.ra = ra;
   star.dec = dec;
   project_star( &star, params[3], params[4]);
#if 1
   star.dec = -star.dec;
#endif
   de_xform_point( params, &star, &pixel_loc);
   *pixel++ = pixel_loc.x;
   *pixel++ = pixel_loc.y;
   return( 0);
}

int pixel_to_ra_dec( double *ra_dec, const double *pixel, const double *params)
{
   REF_STAR star;
   FOUND_STAR pixel_loc;

   pixel_loc.x = *pixel++;
   pixel_loc.y = *pixel++;
   xform_point( params, &pixel_loc, &star);
#if 1
   star.dec = -star.dec;
#endif
   deproject_star( &star, params[3], params[4]);
   *ra_dec++ = star.ra;
   *ra_dec++ = star.dec;
   return( 0);
}

static void deproject_star( REF_STAR *star, const double x0, const double y0)
{
   double x, y, z, y1, z1, d1;
   double cos_phi0, sin_phi0, phi0;
   extern int is_drift_scan;

   phi0 = y0 * (PI / 180.);           /* now we're all in radians */
   cos_phi0 = cos( phi0);
   x = (double)star->ra * (PI / 180.);
   y = (double)star->dec * (PI / 180.);
   if( is_drift_scan)
      {
      star->dec = y0 + star->dec;
      star->ra = x0 + star->ra / cos_phi0;
      }
   else
      {
      sin_phi0 = sin( phi0);
      d1 = x * x + y * y;
#ifdef STEREOGRAPHIC_INSTEAD_OF_TANGENT_PLANE
      d1 = 4. / (4. + d1);
#else
      d1 = 1. / sqrt( 1. + d1);
#endif
      x *= d1;
      y *= d1;
      z = sqrt( 1. - x * x - y * y);
      y1 = y * sin_phi0 - z * cos_phi0;
      z1 = y * cos_phi0 + z * sin_phi0;
      star->dec = asin( z1) * (180. / PI);
      star->ra = x0 + atan2( x, -y1) * (180. / PI);
      }
}

/*
   The following get_initial_match( ) function is the "core" part of the
algorithm that takes a list of stars found in an image,  a list of catalog
stars,  and does the pattern matching between them.

   First,  a copy of the list of catalog stars is made and converted from
(RA, dec) to (xi, eta) values.    The result is that we can now think in
terms of "ordinary" (non-spherical) geometry.  Next,  the list is sorted
in x;  we'll have to do a lot of searching through this list, and that
will be much faster with the sorted list.

   We also make a list of "pairs" among the image stars.  Our strategy
will be to make the initial match by comparing pairs of image stars with
pairs of catalog stars.  Saying that "the pair of stars A and B in the
image match stars X and Y in the catalog" defines exactly an affine
transformation between image and catalog space,  of the sort that can
convert an image pixel (x, y) value to a plate location (xi, eta):

xi  = xi0  + Ax + By
eta = eta0 + Bx - Ay

   This _could_ be done with the following,  "brute force" algorithm:

for( i = 0; i < n_gsc_stars; i++)
   for( j = i + 1; j < n_gsc_stars; j++)
      for( k = 0; k < n0; k++)
         for( l = 0; l < n0; l++)
            if( k != l)
               {
               consider the possibility that stars i and j in the
               catalog match k and l in the image;  compute the affine
               transformation that would do this;  if the image scale is
               within tolerances,  see how many of the image stars have
               catalog star counterparts suitably nearby
               }
Take the "possible match" that results in the maximum number of image
stars being within the tolerance distance of catalog stars.

   Essentially,  each quadruplet (i, j, k, l) defines a possible match
between the image and the catalog stars,  and we look for the quadruplet
that causes the best match between the two lists of stars.

   That was exactly what I did when I first tried this.  The main difficulty
is that as shown above,  we have to consider order (n_gsc_stars * n0)^2
quadruplets,  and that's extremely slow.  Most of the weird code in
get_initial_match( ) results from efforts to drastically reduce the number
of "possible matches" we need to consider.

   But once this is done,  we have a good affine match,  and can use the
method of least squares to refine plate constants.
*/

int get_initial_match( const FOUND_STAR *found, const int n_stars,
           double *params, const double *search_loc, const int n0,
           const REF_STAR *gsc_stars, const int n_gsc_stars,
           const int photo_band)
{
   int n_source = 0, i, j, order;
   int best_n_within_tolerance = 0;
   int n_pairs = n0 * (n0 - 1), step0;
   long t2;
   double search_dist;
   double max_dist = 0.;
   double best_rms_sum = 0., best_xform[N_FIND_PARAMS];
   double max_allowed_residual = search_loc[4] / 3600.;
   STAR_PAIR *pairs, *pair_ptr;
   REF_STAR *istars, *end_star_ptr;
   extern FILE *log_file;
   extern FILE *overlay_file;

                  /* We want to have a copy of the GSC (or other catalog)  */
                  /* stars,  that can be sorted in x.  This speeds up the  */
                  /* process of finding out "what star is closest to (x,y)"*/
                  /* immensely (see the function find_closest( ) above to  */
                  /* see how this helps.)                                  */
   n_source = n_gsc_stars;
   istars = (REF_STAR *)calloc( n_source, sizeof( REF_STAR));
   if( !istars)
      return( -1);
                        /* Copy in the catalog stars... */
   memcpy( istars, gsc_stars, n_source * sizeof( REF_STAR));
   end_star_ptr = istars + n_source;

                       /* Convert 'em from RA/dec to "plate" coordinates */
                       /* (basically,  eta/xi)                           */
   for( i = 0; i < n_source; i++)
      {
      project_star( istars + i, search_loc[0], search_loc[1]);
      istars[i].dec = -istars[i].dec;   /* get around an axis flipping problem */
      }
                       /* Sort the stars in x: */
   qsort( istars, n_source, sizeof( REF_STAR), ref_star_compare);

                       /* Now make a list of all pairs of image stars that */
                       /* we're going to look at.  For each pair,  we want */
                       /* to know the distance between them.               */
                       /* We also keep track of the maximum distance found */
                       /* in this manner.  The way this speeds matters up  */
                       /* will be seen later.                              */

   pairs = pair_ptr = (STAR_PAIR *)calloc( n_pairs, sizeof( STAR_PAIR));
   for( i = 0; i < n0; i++)
      for( j = 0; j < n0; j++)
         if( i != j)
            {
            double dx = found[i].x - found[j].x;
            double dy = found[i].y - found[j].y;

            pair_ptr->start = i;
            pair_ptr->end = j;
            pair_ptr->dist = dx * dx + dy * dy;
            if( max_dist < pair_ptr->dist)
               max_dist = pair_ptr->dist;
            pair_ptr++;
            }

   max_dist = sqrt( max_dist);

               /* Now sort the list of pairs,  in order of distance: */
   qsort( pairs, n_pairs, sizeof( STAR_PAIR), pair_compare);

   t2 = time( NULL);
   search_dist = max_dist / 3600.;
               /* We'll be doing some binary searching in the list of */
               /* pairs,  so it helps to compute step0.  Step0 is the */
               /* largest power of 2 that is less than n_pairs.       */
   for( step0 = 1; step0 * 2 < n_pairs; step0 <<= 1)
      ;
               /* Now we look at every catalog star,  and try to match */
               /* it to an image star.                                 */
   for( i = 0; i < n_source; i++)
      {
      double ymin, ymax, xmax;
      REF_STAR *tstar = istars + i + 1;

      xmax  = istars[i].ra + search_dist;
      ymax  = istars[i].dec + search_dist;
      ymin  = istars[i].dec - search_dist;
      if( t2 != time( NULL))
         {
         printf( "%4.1lf%% complete\r",
                               (double)i * 100. / (double)n_source);
         t2 = time( NULL);
#ifdef __WATCOMC__
         if( kbhit( ))
            if( getch( ) == 27)
               exit( 0);
#endif
         }

            /* Now we look from the 'current' catalog star to the end    */
            /* of the list of catalog stars.  There are some quick tests */
            /* we can use to reduce the number of possible matches that  */
            /* we have to check,  using a bounding box.  Also,  we can   */
            /* make excellent use of the fact that the catalog stars are */
            /* sorted in order of x;  once we reach tstar->x < xmax,     */
            /* it's not necessary to go any farther.                     */

      for( ; tstar < end_star_ptr && tstar->ra < xmax; tstar++)
         if( tstar->dec > ymin && tstar->dec < ymax)
            {
            double dist, minimum, maximum, xform[N_FIND_PARAMS];
            int k = 0, m, step;
            double dx = tstar->ra  - istars[i].ra;
            double dy = tstar->dec - istars[i].dec;

            *xform = 1.;                  /* use a linear xformation only */
            dist = dx * dx + dy * dy;
            dist *= 3600. * 3600.;
                          /* Search_loc[3] contains the 'scale tolerance'  */
                          /* parameter.  The two catalog stars in question */
                          /* are separated by 'dist';  therefore,  we have */
                          /* to look hard at any pair of image stars that  */
                          /* are separated by 'minimum' to 'maximum' units. */
            minimum = dist / search_loc[3];
            maximum = dist * search_loc[3];
                          /* We can immediately skip over any image star     */
                          /* pairs whose separation is less than 'minimum'.. */
            for( step = step0; step; step >>= 1)
               if( (m = (k | step)) < n_pairs && minimum > pairs[m].dist)
                  k = m;
                          /* ...and,  we can stop considering image star   */
                          /* pairs when their separation is > 'maximum'.   */
                          /* This is why we went through all the trouble   */
                          /* to make a sorted list of image star pairs.    */
            for( k++; k < n_pairs && maximum > pairs[k].dist; k++)
               {
               int n_within_tolerance = 2;
               int check_it = 1;
               double rms_sum = 0.;

                            /* We construct a linear,  orthogonal transform */
                           /* between image space and catalog stars,  based */
                          /* on the assumption that the two image stars    */
                         /* happen to be the same as the two catalog stars. */
               calc_xform( xform, found + pairs[k].start,
                                  found + pairs[k].end,
                                  istars + i, tstar);
                      /* 'n_within_tolerance = 2' reflects the fact that     */
                      /* we know full well that at least those two will      */
                      /* get matched.  But how about the other n0-2          */
                      /* image stars?  Do they have catalog counterparts?    */
                      /* (Check first to be sure that we have a valid        */
                      /* transformation,  not tilted at too great an angle.  */
                    /* The tilt tolerance is in search_loc[5].  The expected */
                    /* tilt is in search_loc[6];  if it's zero,  we can skip */
                    /* an atan2( ) call,  a pleasant speedup.)               */
               if( search_loc[5])
                  if( !search_loc[6])        /* "standard" N/S tilt expected */
                     {
                     double axis_tilt = fabs( xform[7] / xform[6]);

                     if( axis_tilt > fabs( search_loc[5]))
                        check_it = 0;     /* too much tilt to the N/S axis */
//  Why did I        if( search_loc[5] > 0.)   /* check north_up only */
//  have these          if( xform[6] > 0.)
//  three lines??          check_it = 0;
                     }
                  else
                     {
                     double tilt = atan2( xform[7], xform[6]) - search_loc[6];

                     while( tilt > PI)
                        tilt -= PI + PI;
                     while( tilt < -PI)
                        tilt += PI + PI;
                     if( tilt > search_loc[5] || tilt < -search_loc[5])
                        check_it = 0;
                     }
               if( check_it)
                  for( m = 0; m < n0; m++)
                     if( m != pairs[k].start && m != pairs[k].end)
                        {
                        REF_STAR loc;
                        double d;
                        int closest;
                        double x = found[m].x, y = found[m].y;
                        const double *yform = xform + N_FIND_COEFFS;

                             /* Transform the image star to catalog space. */
                        loc.ra = xform[5] + xform[6]*x + xform[7] * y;
                        loc.dec = yform[5] + yform[6]*x + yform[7] * y;
             /* Look for the result in the list of catalog stars. */
                        d = find_closest( n_source, istars, &loc, &closest,
                                       max_allowed_residual * 5.);
             /* If the result happens to be within tolerances,  add it in. */
                        if( d < max_allowed_residual)
                           {
                           n_within_tolerance++;
                           rms_sum += d * d;
                           }
                        }
                           /* Next,  ask:  did this 'possible match' get */
                           /* more stars within tolerance than any   */
                           /* other we've tried,  or at least reduce our */
                           /* RMS residuals?  i.e.,  is it 'the best yet'? */
               if( n_within_tolerance > best_n_within_tolerance ||
                            (n_within_tolerance == best_n_within_tolerance &&
                                      rms_sum < best_rms_sum))
                  {
                           /* If we've got a winner,  save the xformation */
                           /* and similar data. */
                  best_n_within_tolerance = n_within_tolerance;
                  best_rms_sum = rms_sum;
                  memcpy( best_xform, xform, N_FIND_PARAMS * sizeof( double));
                  }
               }
            }
      }


               /* OK,  we can now assume we've got a good match.  Let's */
               /* go on to take our 'best case' and try to match some   */
               /* more stars to it.                                     */

   for( i = 3; i < N_FIND_COEFFS; i++)   /* zero out higher-order coeffs */
      best_xform[i + 5] = best_xform[i + 5 + N_FIND_COEFFS] = 0.;

                /* We'll start by improving our _linear_ fit.  If asked,  */
                /* we'll also improve the quad and/or cubic fit.          */
   for( order = 1; order <= (int)( params[0] + .5); order++)
      for( i = 0; i < N_PASSES; i++)
         {
         int n_coeffs = (order + 1) * (order + 2) / 2;
         void *lsquare_x = lsquare_init( n_coeffs);
         void *lsquare_y = lsquare_init( n_coeffs);
         double diff[N_FIND_COEFFS], sum_x2 = 0., sum_y2 = 0.;
         double magnitude_top = 0., magnitude_bottom = 0.;
         int n_within_tolerance = 0, err;

         *best_xform = (double)order;
                     /* We're now looking at _all_ image stars... not just */
                     /* the n0 brightest ones.                             */
         for( j = 0; j < n_stars; j++)
            {
            REF_STAR loc, tloc;
            double d, dx, dy, mag = 99.99;
            int closest;
            char tbuff[120];

                           /* Convert the image star to catalog space. */
            xform_point( best_xform, found + j, &loc);
            tloc = loc;
            tloc.dec = -tloc.dec;   /* deal with the silly y convention again */
            deproject_star( &tloc, search_loc[0], search_loc[1]);
            closest = 0;
                           /* As before,  find the nearest catalog star... */
            d = find_closest( n_source, istars, &loc, &closest,
                                          max_allowed_residual * 5.);
            dx = (double)( istars[closest].ra - loc.ra);
            dy = (double)( istars[closest].dec - loc.dec);
                      /* ..and if it's within tolerance,  add it to the */
                      /* least squares fit.  (This _is_ one difference: */
                      /* in the initial pattern matching,  all we had to */
                      /* do was create a transform matching stars A, B to */
                      /* catalog stars C, D,  getting a perfect match for */
                      /* them.  Now,  we have to do a full-fledged least  */
                      /* squares step.)                                   */
            if( d < max_allowed_residual)
               {
               double slopes[N_FIND_COEFFS];

               compute_powers( slopes, found[j].x, found[j].y, order);
               lsquare_add_observation( lsquare_x, dx, 1., slopes);
               lsquare_add_observation( lsquare_y, dy, 1., slopes);
                           /* We may also have to figure out the magnitude */
                           /* parameters,  as follows.                     */
                           /* We should only include those from the  */
                           /* desired band in the magnitude fit:     */
               if( !use_input_mags && found[j].bright > 0)
                  if( !istars[closest].photo_band || !photo_band ||
                               photo_band == istars[closest].photo_band)
                     {
                     double pseudo_mag, sqrt_counts;

                     mag = (double)istars[closest].mag / 100.;
                     sqrt_counts = sqrt( (double)found[j].bright);
                     pseudo_mag = -2.5 * log10( (double)found[j].bright);
                     magnitude_top += (mag - pseudo_mag) * sqrt_counts;
                     magnitude_bottom += sqrt_counts;
                     }
               n_within_tolerance++;
               }
            dx *= 3600.;
            dy *= 3600.;
            if( d < max_allowed_residual)
               {
               sum_x2 += dx * dx;
               sum_y2 += dy * dy;
               }
                     /* Tell the user what's going on. */
            if( d < max_allowed_residual * 3. && i == N_PASSES - 1)
               {                                /* only print on last pass */
               char ra_dec_str[100];
               extern int suppress_display_of_residuals;

               if( found[j].bright > 0)
                  {
                  double pseudo_mag = -2.5 * log10( (double)found[j].bright);

                  mag = best_xform[1] + pseudo_mag;
                  }
               else
                  mag = 0.;
               put_ra_dec_in_str( ra_dec_str, tloc.ra, tloc.dec, 2);
               sprintf( tbuff, "%2d: %4d%5d   %s %4d %5.2lf dist %6.3lf  ",
                  j, istars[closest].zone, istars[closest].number,
                  ra_dec_str, istars[closest].mag, mag, d * 3600.);
               sprintf( tbuff + strlen( tbuff), "dx%6.3lf dy%6.3lf", dx, dy);
               if( istars[closest].ppm_num > 0L && istars[closest].ppm_num < 1000000L)
                  sprintf( tbuff + strlen( tbuff), " PPM %ld\n",
                                               istars[closest].ppm_num);
               else
                  strcat( tbuff, "\n");
               if( !suppress_display_of_residuals)
                  log_printf( log_file, "%s", tbuff);
               }
            if( overlay_file && i == N_PASSES - 1 && found[j].bright)
               {
               if( use_input_mags)        /* for xy lists */
                  log_printf( overlay_file, "s 0%4ld\n", found[j].bright);
               else
                  log_printf( overlay_file, "s 0%4d\n", (int)( mag * 100. + .5));
               stuff_ang_into_buff( tbuff, -tloc.ra);
               stuff_ang_into_buff( tbuff + 3, tloc.dec);
               fwrite( tbuff, 6, 1, overlay_file);
               }
            }
         if( i == N_PASSES - 1)
            {
            best_n_within_tolerance = n_within_tolerance;
            log_printf( log_file, "%d %s\n", n_within_tolerance, strings[109]);
            log_printf( log_file, "rms dx = %.2lf   rms dy = %.2lf  rms err = %.2lf\n",
                                   sqrt( sum_x2 / (double)n_within_tolerance),
                                   sqrt( sum_y2 / (double)n_within_tolerance),
                       sqrt( (sum_x2 + sum_y2) / (double)n_within_tolerance));
            }

                  /* Now do the least squares math,  and adjust the xform. */
         err = lsquare_solve( lsquare_x, diff);
         if( !err)
            for( j = 0; j < n_coeffs; j++)
               best_xform[j + 5] += diff[j];

         err = lsquare_solve( lsquare_y, diff);
         if( !err)
            for( j = 0; j < n_coeffs; j++)
               best_xform[j + N_FIND_COEFFS + 5] += diff[j];

         if( !use_input_mags)
            best_xform[1] = magnitude_top / magnitude_bottom;

         lsquare_free( lsquare_x);
         lsquare_free( lsquare_y);
         }

   free( pairs);
   free( istars);
               /* Return the 'best xform' we've found. */
   memcpy( params, best_xform, N_FIND_PARAMS * sizeof( double));
   params[3] = search_loc[0];
   params[4] = search_loc[1];
   return( best_n_within_tolerance);
}
