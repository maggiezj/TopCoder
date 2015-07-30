#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "gsc_act.h"

/*
   For this sort of work,  it's convenient to swap from the polar (ra, dec)
coordinates to "plate coordinates" (xi, eta).  The following functions
provide transformations in both directions. */

void ra_dec_to_xi_eta( double d_ra, double dec, double dec0, double *xi,
                           double *eta)
{
   double x, y, z, cos_dec0, sin_dec0, cos_dec, sin_dec;

   cos_dec0 = cos( dec0);
   sin_dec0 = sin( dec0);
   cos_dec = cos( dec);
   sin_dec = sin( dec);
   x = sin( d_ra) * cos_dec;
   y = cos_dec0 * sin_dec - cos_dec *  cos( d_ra) * sin_dec0;
   z = sqrt( 1. - x * x - y * y);
   *xi = x / z;
   *eta = y / z;
}

void xi_eta_to_ra_dec( double xi, double eta, double dec0, double *d_ra,
                           double *dec)
{
   double z = sqrt( 1. + xi * xi + eta * eta), x, y, y1, z1;

   x = xi / z;       /* (x, y, z) will be a unit-length vector */
   y = eta / z;
   z = 1. / z;
               /* now rotate by 'dec' radians around the x-axis: */
   y1 = y * cos( dec0) + z * sin( dec0);
   z1 = z * cos( dec0) - y * sin( dec0);
               /* and convert back to polar coords: */
   *dec = asin( y1);
   *d_ra = atan2( x, z1);
}

/*
   As mentioned in the documentation,  GSC-ACT uses a fifth-order fit,
with 21 coefficients on each axis (xi and eta).   */

void compute_gsc_act_powers( double *arr, const double x, const double y,
                     const double mag)
{
   double r2;

   arr[0] = 1.;

   arr[1] = x;
   arr[2] = y;

   arr[3] = x * x;
   arr[4] = x * y;
   arr[5] = y * y;

   arr[6] = x * x * x;
   arr[7] = x * x * y;
   arr[8] = x * y * y;
   arr[9] = y * y * y;

   arr[10] = x * x * x * x;
   arr[11] = x * x * x * y;
   arr[12] = x * x * y * y;
   arr[13] = x * y * y * y;
   arr[14] = y * y * y * y;

   arr[15] = x * x * x * x * x;
   arr[16] = x * x * x * x * y;
   arr[17] = x * x * x * y * y;
   arr[18] = x * x * y * y * y;
   arr[19] = x * y * y * y * y;
   arr[20] = y * y * y * y * y;

   r2 = arr[3] + arr[5];

#ifdef INCLUDE_MAG_PARAMS
   arr[21] = x * r2 * mag;
   arr[22] = y * r2 * mag;
   arr[23] = x * r2 * r2 * mag;
   arr[24] = y * r2 * r2 * mag;
   for( int i = 25; i < 29; i++)
      arr[i] = arr[i - 4] * mag;
#endif
}

/*
   The following function provides a simple way to load up plate
transformation data from the GSC_ACT.DAT file.  You open that file using
a standard fopen( ) call,  and then request data for each plate by name.
You have to supply an array of PLATE_STRUCT structures,  and the length
of that array.  The code is then bright enough to keep the array in order
of use.  This is necessary because you may well have four (sometimes more)
plates covering a given area,  and will be swapping among them briskly. */

int gsc_act_cached_load( FILE *data_file, char *plate_name,
                                    PLATE_STRUCT *plate_array, int n_cached)
{
   int i;
   PLATE_STRUCT temp;
   long header[2];
   char tname[4];

   for( i = 0; i < n_cached; i++)
      if( !strcmp( plate_array[i].plate_name, plate_name))  /* got it */
         {
         if( i)     /* gotta move this plate to the top of the array */
            {
            memcpy( &temp, plate_array + i, sizeof( PLATE_STRUCT));
            memmove( plate_array + 1, plate_array, i * sizeof( PLATE_STRUCT));
            memcpy( plate_array, &temp, sizeof( PLATE_STRUCT));
            }
         return( plate_array[0].n_stars ? 0 : -2);
         }
                     /* gotta load plate data from the file: */
   fseek( data_file, 0L, SEEK_SET);
   fread( header, 2, sizeof( long), data_file);
   if( header[0] != 2076665750)       /* magic/version # */
      return( -1);
   memmove( plate_array + 1, plate_array,
                              (n_cached - 1) * sizeof( PLATE_STRUCT));
   for( i = 0; i < header[1] && fread( tname, 4, 1, data_file); i++)
      if( !strcmp( tname, plate_name))    /* we have a winner */
         {
         fseek( data_file, 2 * sizeof( long) + header[1] * 4 +
                     (long)i * (long)sizeof( PLATE_STRUCT), SEEK_SET);
         fread( plate_array, 1, sizeof( PLATE_STRUCT), data_file);
         return( 0);
         }
                  /* If we can't find the plate,  put a suitable 'dummy' */
                  /* record at the top of the array:                     */
   memset( plate_array, 0, sizeof( PLATE_STRUCT));
   strcpy( plate_array[0].plate_name, plate_name);
   return( -2);            /* plate not found */
}

int apply_transformation( double ra_in, double dec_in, double *ra_out,
                  double *dec_out, PLATE_STRUCT *plate_data)
{
   double xi, eta, power_arr[N_COEFFS];
   int i;

   ra_dec_to_xi_eta( ra_in - plate_data->ra0, dec_in, plate_data->dec0,
                           &xi, &eta);
   compute_gsc_act_powers( power_arr, xi, eta, 0.);
                                   /* magnitude is ignored now */

   for( i = 0; i < N_COEFFS; i++)
      {
      xi  += power_arr[i] * plate_data->coeffs[i];
      eta += power_arr[i] * plate_data->coeffs[i + N_COEFFS];
      }
   xi_eta_to_ra_dec( xi, eta, plate_data->dec0, ra_out, dec_out);
   *ra_out += plate_data->ra0;
   return( 0);
}
