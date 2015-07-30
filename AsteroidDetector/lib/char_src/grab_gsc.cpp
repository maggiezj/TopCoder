#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "watdefs.h"
#include "gsc.h"
#include "crunch.h"
#include "findstar.h"
#include "ppm.h"
#include "tycho.h"
#include "gsc_act.h"
#include "colors.h"
#include "usno_a.h"
#include "ucac2.h"
#include "ucac3.h"
#include "ucac4.h"

#ifdef DEBUG_MEM
#include "checkmem.h"
#endif

#define PI 3.141592653589793
#define MAX_SAO 258997L
#define N_PLATES (12144 / 8 + 1)

int gsc_mag_limit = 9999, gsc_photo_band = 0;
int gsc_bright_mag_limit = 0;
extern int debug_level;
extern FILE *log_file;
extern char *strings[];

int log_printf( FILE *log_file, const char *format, ...);   /* charon.cpp */
int get_sao_info( const char *cd_drive_path, long sao_number, char *buff);
long grab_ppm_data( const char *cd_drive_path, PPM_STAR *ppm_star, long ppm_num);
int find_gsc_star_loc( const char *cd_drive_path, int zone, int number,
                               double *ra_dec);          /* grab_gsc.c */
FILE *local_then_cd_fopen( const char *cd_drive_path, const char *filename);
REF_STAR *grab_ucac2_data( const double x, const double y, const double width,
                const double height, const int cd_drive_letter, int *n_found,
                const int photo_band, const double jd);     /* grab_gsc.c */
REF_STAR *grab_ucac3_data( const double x, const double y, const double width,
                const double height, const int cd_drive_letter, int *n_found,
                const int photo_band, const double jd);     /* grab_gsc.c */
REF_STAR *grab_ucac4_data( const double x, const double y, const double width,
                const double height, const int cd_drive_letter, int *n_found,
                const int photo_band, const double jd);     /* grab_gsc.c */
void exit_2( const int exit_code);                    /* charon.cpp */

int get_sao_info( const char *cd_drive_path, long sao_number, char *buff)
{
   static char zbuff[36];
   static long curr_no = -1;

   if( sao_number < 1L || sao_number > MAX_SAO)
      return( -1);
   if( curr_no != sao_number)
      {
      FILE *ifile;
      long offset = 0L;

      ifile = local_then_cd_fopen( cd_drive_path, "sao\\index.sao");
      if( !ifile)
         return( -1);
      fseek( ifile, sao_number * 3L, SEEK_SET);
      fread( (char *)&offset, 3, 1, ifile);
      fclose( ifile);

      ifile = local_then_cd_fopen( cd_drive_path, "sao\\sao.lmp");
      if( !ifile)
         return( -1);
      fseek( ifile, 16106L + 36L * offset, SEEK_SET);
      fread( zbuff, 36, 1, ifile);
      fclose( ifile);
      }
   memcpy( buff, zbuff, 36);
   curr_no = sao_number;
   return( 0);
}

long grab_ppm_data( const char *cd_drive_path, PPM_STAR *ppm_star, long ppm_num)
{
   long rec_num = ppm_num - 1L;
   FILE *ifile;

   if( ppm_num < 1L)
      return( -1L);
   ifile = local_then_cd_fopen( cd_drive_path, "ppm\\ppm.dat");
   if( !ifile)
      return( -1L);
   if( rec_num > 181720L)
      rec_num--;
   if( rec_num > 378900L)
      rec_num = 378900L;
   fseek( ifile, rec_num * (long)sizeof( PPM_STAR), SEEK_SET);
   do
      {
      if( !(fread( ppm_star, 1, sizeof( PPM_STAR), ifile)))
         {
         fclose( ifile);
         return( -1L);
         }
      rec_num++;
      }
      while( ppm_star->ppm_num < ppm_num);
   fclose( ifile);
   return( (ppm_star->ppm_num == ppm_num) ? rec_num - 1L : -1L);
}

FILE *open_gsc_index( const char *cd_drive_path)
{
   FILE *rval = local_then_cd_fopen( cd_drive_path, "text\\gscdata2.idx");
   if( !rval)
      log_printf( log_file, "Couldn't open GSC index file");
   return( rval);
}

int find_gsc_star_loc( const char *cd_drive_path, int zone, int number,
                               double *ra_dec)
{
   FILE *data_file, *index_file = open_gsc_index( cd_drive_path);
   char *buff = NULL, *tptr;
   long loc[2], bytes_left, buffsize = 32768;
   GSC_STAR star;
   GSC_HEADER hdr;
   char filename[30];
   int dec_zone;
   static int dec_zone_start[25] = {    1,  594, 1178, 1729, 2259, 2781,
                                     3246, 3652, 4014, 4294, 4492, 4615,
                                     4663, 5260, 5838, 6412, 6989, 7523,
                                     8022, 8464, 8840, 9134, 9346, 9490, 9538 };

   if( !index_file)
      return( -1);
   fseek( index_file, 8L * ((long)zone + 731L), SEEK_SET);
   fread( loc, 2, sizeof( long), index_file);
   fclose( index_file);

   while( !buff)
      {
      buffsize >>= 1;
      buff = (char *)calloc( (unsigned)buffsize, sizeof( char));
      }

   for( dec_zone = 0; zone >= dec_zone_start[dec_zone]; dec_zone++)
      ;
   dec_zone--;
   sprintf( filename, "gsc\\%c%04d.lmp",
                  (dec_zone > 11 ? 's' : 'n'), 750 * (dec_zone % 12)
                  - 20 * (dec_zone % 2));
   log_printf( log_file, "GSC zone: %s\n", filename);
   data_file = local_then_cd_fopen( cd_drive_path, filename);
   if( !data_file)
      return( -2);
   fseek( data_file, loc[0], SEEK_SET);
   fread( &hdr, 1, sizeof( GSC_HEADER), data_file);
   hdr.prev_n = 0;
   bytes_left = hdr.filesize - (long)sizeof( GSC_HEADER);
   fread( buff, (long)buffsize > bytes_left ?
                              (unsigned)bytes_left : buffsize, 1, data_file);
   tptr = buff;
   ra_dec[0] = ra_dec[1] = 0.;
   while( bytes_left && ra_dec[0] == 0.)
      {
      unsigned n_valid, rval;

      if( tptr > buff + buffsize - 50)
         {
         unsigned readsize = tptr - buff;

         n_valid = buffsize - readsize;
         memmove( buff, tptr, n_valid);
         if( (long)readsize > bytes_left)
            readsize = (unsigned)bytes_left;
         fread( buff + n_valid, 1, readsize, data_file);
         tptr = buff;
         }
      rval = gsc_star_uncrunch( &star, &hdr, (unsigned char *)tptr);
      tptr += rval;
      bytes_left -= (long)rval;

      if( star.n == number)
         {
         ra_dec[0] = (double)star.x / 100000.;
         ra_dec[1] = (double)star.y / 100000.;
         }
      }
   fclose( data_file);
   free( buff);
   return( ra_dec[0] == 0.);
}

SMALL_AREA *find_small_areas( double x, double y, double width, double height,
                int *n_found, double *minmax)
{
   SMALL_AREA *areas;
   int pass, xoffset;

   minmax[0] = x - width;
   minmax[1] = y - height;
   minmax[2] = x + width;
   minmax[3] = y + height;
   if( !n_found)              /* just after bounds */
      return( NULL);
   for( pass = 0; pass < 2; pass++)
      {
      if( pass)
         {
         areas = (SMALL_AREA *)calloc( *n_found, sizeof( SMALL_AREA));
         if( !areas)
            return( NULL);
         }
      *n_found = 0;
      for( xoffset = -1; xoffset < 2; xoffset++)
         {
         LPT lmin, lmax;

         lmin.x = (long)( minmax[0] * 100000.) + (long)xoffset * THREESIXTY_DEG;
         lmin.y = (long)( minmax[1] * 100000.);
         lmax.x = (long)( minmax[2] * 100000.) + (long)xoffset * THREESIXTY_DEG;
         lmax.y = (long)( minmax[3] * 100000.);
         *n_found += find_regions( &lmin, &lmax,
                                        (pass ? areas + *n_found : NULL),
                                        65000U / sizeof(SMALL_AREA), 1);
         }
      }
   return( areas);
}

static int in_bounds( double *ra, const double dec, const double *minmax)
{
   int rval = 0;

   if( dec > minmax[1] && dec < minmax[3])
      {
      while( *ra > minmax[2])
         *ra -= 360.;
      while( *ra < minmax[0])
         *ra += 360.;
      rval = (*ra > minmax[0] && *ra < minmax[2]);
     }
   return( rval);
}

REF_STAR *grab_gsc_12_data( const double x, const double y,
                             const double width, const double height,
                             const char *gsc_12_filename, int *n_found,
                             const double jd)
{
   FILE *ifile = fopen( gsc_12_filename, "rb");
   int n_stars, count, pass, catalog_found = 0;
   REF_STAR *ret_arr;
   double minmax[4];

   if( !ifile)
      {
      *n_found = -4;
      return( NULL);         /* Couldn't get GSC data */
      }
   find_small_areas( x, y, width, height, NULL, minmax);
   for( pass = 0; pass < 2; pass++)
      {
      char buff[400];

      fseek( ifile, 0L, SEEK_SET);
      memset( buff, 0, sizeof( buff));
      count = 0;
      while( fgets( buff, sizeof( buff), ifile))
         {
         int is_valid = 0;
         REF_STAR star;
         extern char *ref_net_name;

         memset( &star, 0, sizeof( REF_STAR));
                        /* Check for possible GSC 1.2 data: */
         if( buff[25] == '.' && buff[38] == '.' && buff[43] == '.')
            {
            long hr, min, sec, fract, ra, dec;

            sscanf( buff + 17, "%ld %ld %ld.%ld", &hr, &min, &sec, &fract);
            ra = (hr * 360000L + min * 6000L + sec * 100L + fract);
            star.ra = (double)ra / 24000.;
            sscanf( buff + 30, "%ld %ld %ld.%ld", &hr, &min, &sec, &fract);
            dec = (hr * 36000L + min * 600L + sec * 10L + fract);
            star.dec = (double)dec / 36000.;
            if( buff[29] == '-')
               star.dec = -star.dec;
            star.mag = (int)( atof( buff + 47) * 100. + .5);
            star.zone = atoi( buff);
            star.number = atoi( buff + 8);
            if( !catalog_found)
               log_printf( log_file, "GSC 1.2 data detected\n");
            is_valid = catalog_found = 1;
            ref_net_name = "GSC-1.2";
            }
                        /* 13 Jun 2003:  soon,  I hope,  we'll all switch */
                        /* to "final" UCAC2,  and the following lines will */
                        /* be obsolete... but in the interim:              */
         else if( strlen( buff) == 151 && buff[13] == '.' && buff[123] == '.'
                  && atoi( buff + 41) > 1950 && atoi( buff + 41) < 2006
                  && sscanf( buff, "%ld%lf%lf", &star.ppm_num, &star.ra,
                                                   &star.dec) == 3)
            {
            star.mag = (int)( atof( buff + 36) * 100. + .5);
            star.ppm_num = atol( buff);
            if( !catalog_found)
               log_printf( log_file, "UCAC-2 beta data detected\n");
            is_valid = catalog_found = 3;
            ref_net_name = "UCAC2-BETA";
            }
         else if( strlen( buff) == 156 && buff[13] == '.' && buff[128] == '.'
                  && atoi( buff + 41) > 1950 && atoi( buff + 41) < 2006
                  && sscanf( buff, "%ld%lf%lf", &star.ppm_num, &star.ra,
                                                   &star.dec) == 3)
            {
            star.mag = (int)( atof( buff + 36) * 100. + .5);
            if( !catalog_found)
               log_printf( log_file, "UCAC-2 data detected\n");
            is_valid = catalog_found = 3;
            ref_net_name = "UCAC2";
            }
         else if( strlen( buff) == 239 && buff[14] == '.' && buff[46] == '.'
                  && atoi( buff + 60) > 1950 && atoi( buff + 60) < 2009
                  && buff[3] == '-'
                  && sscanf( buff + 4, "%ld%lf%lf", &star.ppm_num, &star.ra,
                                                   &star.dec) == 3)
            {
            const double J2000 = 2451545.;
                     /* Prop motions in .0001 arcsec/yr: */
            const double multiplier =
                        (jd - J2000) / (365.25 * 3.6e+7);

            star.ra  += atof( buff + 96)  * multiplier
                           / cos( star.dec * PI / 180.);
            star.dec += atof( buff + 103) * multiplier;
            star.mag = (int)( atof( buff + 37) * 100. + .5);
            if( !catalog_found)
               log_printf( log_file, "UCAC-3 data detected\n");
            is_valid = catalog_found = 3;
            star.ppm_num += 1000000L * atol( buff);
            ref_net_name = "UCAC3";
            }
         else if( strlen( buff) == 277 && buff[14] == '.' && buff[46] == '.'
                  && atoi( buff + 60) > 1950 && atoi( buff + 60) < 2009
                  && buff[3] == '-'
                  && buff[79] != ' ' && atoi( buff + 77) < 100  /* check pmRA */
                  && buff[83] != ' ' && atoi( buff + 81) < 100  /* check pmDec*/
                  && !memcmp( buff + 57, "  0", 3)
                  && (buff[56] == '0' || buff[56] == '1')
                  && sscanf( buff + 4, "%ld%lf%lf", &star.ppm_num, &star.ra,
                                                   &star.dec) == 3)
            {
            const double J2000 = 2451545.;
                     /* Prop motions in .0001 arcsec/yr: */
            const double multiplier =
                        (jd - J2000) / (365.25 * 3.6e+7);

            star.ra  += atof( buff + 93)  * multiplier
                           / cos( star.dec * PI / 180.);
            star.dec += atof( buff + 100) * multiplier;
            star.mag = (int)( atof( buff + 37) * 100. + .5);
            if( !catalog_found)
               log_printf( log_file, "UCAC-4 data detected\n");
            is_valid = catalog_found = 3;
            star.ppm_num += 1000000L * atol( buff);
            ref_net_name = "UCAC4";
            }
         else if( strlen( buff) > 200 && buff[13] == '.' && buff[139] == '.'
                  && atoi( buff + 76) > 1950 && atoi( buff + 76) < 2006
                  && sscanf( buff, "%ld%lf%lf", &star.ppm_num, &star.ra,
                                                   &star.dec) == 3)
            {
            star.mag = (int)( atof( buff + 57) * 100. + .5);
            if( !catalog_found)
               log_printf( log_file, "UCAC-2 data detected from downloaded data\n");
            is_valid = catalog_found = 3;
            ref_net_name = "UCAC2";
            }
         else if( strlen( buff) > 80 && buff[5] == 'U' && buff[10] == '_'
                  && buff[4] == ' ' && buff[42] == '.')
            {
            star.mag = (int)( atof( buff + 46) * 100. + .5);
            star.ra = atof( buff + 21) * 15. + atof( buff + 24) / 4.
                                             + atof( buff + 27) / 240.;
            star.dec = atof( buff + 34)      + atof( buff + 37) / 60.
                                             + atof( buff + 40) / 3600.;
            star.ppm_num = (atol( buff + 6) / 75 + 24) * 50000000
                                 + atol( buff + 11);
            if( buff[33] == '-')
               star.dec *= -1.;
            if( !catalog_found)
               log_printf( log_file, "USNO-A2.0 data detected from downloaded data\n");
            is_valid = catalog_found = 3;
            ref_net_name = "USNO-A2.0";
            }
                        /* Check for possible Sloan ACR data: */
         else if( buff[24] == '.' && buff[36] == '.' && buff[119] == '.')
            {
            star.ra = atof( buff + 22) * 15.;
            star.dec = atof( buff + 33);
            star.mag = (int)( atof( buff + 45) * 100. + .5);
            star.ppm_num = atol( buff);
            if( !catalog_found)
               log_printf( log_file, "Sloan ACR data detected\n");
            is_valid = catalog_found = 2;
            ref_net_name = "SLOAN-ACR";
            }

         if( is_valid && star.mag && star.mag <= gsc_mag_limit
                                  && star.mag > gsc_bright_mag_limit)
            if( in_bounds( &star.ra, star.dec, minmax))
               {
               if( pass)
                  ret_arr[count] = star;
               count++;
               }

         memset( buff, 0, sizeof( buff));
         }
      n_stars = count;
      if( !pass && n_stars)
         {
         ret_arr = (REF_STAR *)calloc( n_stars, sizeof( REF_STAR));
         if( !ret_arr)
            {
            fclose( ifile);
            *n_found = -5;
            return( NULL);
            }
         }
      }
   fclose( ifile);
   *n_found = n_stars;
   return( ret_arr);
}

REF_STAR *grab_gsc_data( const double x, const double y,
               const double width, const double height,
               const char *cd_drive_path, int *n_found, const int use_gsc_act)
{
   int n_areas, curr_zone = -1, i, curr_plate = 0, j;
   double minmax[4];
   SMALL_AREA *areas;
   FILE *index_file, *data_file = NULL, *gsc_act_file;
   unsigned buffsize;
   long *plate_info = (long *)calloc( N_PLATES, 2 * sizeof( long));
   char *buff, *tptr;
   int n_alloced = 16384;
   int n_per_band[3];
   REF_STAR *ret_arr;
   const int plate_cache_size = 7;
   PLATE_STRUCT *plate_cache = NULL;
   char filename[80];

   if( debug_level)
      log_printf( log_file, "Collecting GSC data: ");
   *n_found = 0;
   ret_arr = (REF_STAR *)calloc( n_alloced, sizeof( REF_STAR));
   if( !plate_info)
      {
      *n_found = -9;
      return( NULL);
      }
   else
      {
      FILE *ifile = fopen( "plates.bin", "rb");

      if( !ifile)
         {
         *n_found = -9;
         return( NULL);
         }
      fread( plate_info, N_PLATES - 1, 2 * sizeof( long), ifile);
      fclose( ifile);
      plate_info[N_PLATES * 2 - 2] = *(long *)"+056";
      plate_info[N_PLATES * 2 - 1] = 2451545L;
      }

   areas = find_small_areas( x, y, width, height, &n_areas, minmax);
   if( debug_level)
      log_printf( log_file, "%d areas\n", n_areas);
   if( !areas)
      {
      *n_found = -1;
      return( NULL);
      }

   index_file = open_gsc_index( cd_drive_path);
   if( !index_file)
      {
      *n_found = -2;
      return( NULL);
      }

   if( use_gsc_act)
      {
      gsc_act_file = fopen( "gsc_act.dat", "rb");
      if( !gsc_act_file)         /* look on CD-ROM: */
         {
         sprintf( filename, "gsc\\gsc_act.dat", cd_drive_path);
         gsc_act_file = local_then_cd_fopen( cd_drive_path, filename);
         }

      if( !gsc_act_file)
         {
         log_printf( log_file,
             "GSC_ACT.DAT not found!  Can't correct to GSC-ACT positions\n");
         *n_found = -1;
         return( NULL);
         }
      plate_cache = (PLATE_STRUCT *)calloc( plate_cache_size,
                                                sizeof( PLATE_STRUCT));
      }


   buffsize = 32000U;
   while( (buff = (char *)malloc( buffsize)) == NULL)
      buffsize /= 2U;
   n_per_band[0] = n_per_band[1] = n_per_band[2] = 0;

   for( i = 0; i < n_areas; i++)
      {
      SMALL_AREA *area_ptr = areas + i;
      long loc[2], bytes_left;
      GSC_STAR star;
      GSC_HEADER hdr;

      if( debug_level > 1)
         log_printf( log_file,
                  "Area %d: small zone %d, ", i, area_ptr->small_zone);
      fseek( index_file, 8L * ((long)area_ptr->small_zone + 731L), SEEK_SET);
      fread( loc, 2, sizeof( long), index_file);
      if( area_ptr->zone != curr_zone)
         {
         if( data_file)
            fclose( data_file);
         curr_zone = area_ptr->zone;
         sprintf( filename, "gsc\\%c%04d.lmp",
                  (curr_zone > 11 ? 's' : 'n'), 750 * (curr_zone % 12)
                  - 20 * (curr_zone % 2));
         if( debug_level > 1)
            log_printf( log_file, "opening '%s' ", filename);
         data_file = local_then_cd_fopen( cd_drive_path, filename);
         if( !data_file)
            {
            *n_found = -3;
            return( NULL);
            }
         }
      fseek( data_file, loc[0], SEEK_SET);
      fread( &hdr, 1, sizeof( GSC_HEADER), data_file);
      hdr.prev_n = 0;
      bytes_left = hdr.filesize - (long)sizeof( GSC_HEADER);
      if( debug_level > 1)
         log_printf( log_file, "file size %ld (compare to %ld, %ld)\n",
                       hdr.filesize, loc[0], loc[1]);
      fread( buff, (long)buffsize > bytes_left ?
                              (unsigned)bytes_left : buffsize, 1, data_file);
      tptr = buff;
      while( bytes_left)
         {
         unsigned n_valid, rval;
         double ra, dec;

         if( tptr > buff + buffsize - 50)
            {
            unsigned readsize = tptr - buff;

            n_valid = buffsize - readsize;
            memmove( buff, tptr, n_valid);
            if( (long)readsize > bytes_left)
               readsize = (unsigned)bytes_left;
            fread( buff + n_valid, 1, readsize, data_file);
            tptr = buff;
            }
         rval = gsc_star_uncrunch( &star, &hdr, (unsigned char *)tptr);
         if( (long)rval > bytes_left)
            {
            *n_found = -4;
            return( NULL);
            }
         tptr += rval;
         bytes_left -= (long)rval;
         ra = (double)star.x / 100000.;
         dec = (double)star.y / 100000.;

         if( star.mag && star.mag <= gsc_mag_limit
                                  && star.mag > gsc_bright_mag_limit)
            if( in_bounds( &ra, dec, minmax))
               {
               if( !star.is_a_dup)
                  {
                  REF_STAR *curr_ptr = ret_arr + (*n_found);

                  if( use_gsc_act)
                     if( !gsc_act_cached_load( gsc_act_file, star.plate + 1,
                                         plate_cache, plate_cache_size))
                        {
                        double ra_out, dec_out;

                        ra *= PI / 180.;
                        dec *= PI / 180.;
                        apply_transformation( ra, dec,
                                      &ra_out, &dec_out, plate_cache);
                        while( ra_out > ra + PI)
                           ra_out -= PI * 2.;
                        while( ra_out < ra - PI)
                           ra_out += PI * 2.;
                        if( use_gsc_act == -1)     /* reversed direction! */
                           {
                           ra_out = ra + ra - ra_out;
                           dec_out = dec + dec - dec_out;
                           }
                        ra = ra_out * 180. / PI;
                        dec = dec_out * 180. / PI;
                        }
                  curr_ptr->ra = ra;
                  curr_ptr->dec = dec;
                  curr_ptr->mag = star.mag;
                  curr_ptr->zone = area_ptr->small_zone;
                  curr_ptr->number = star.n;
                  if( plate_info[curr_plate] != *(long *)( star.plate))
                     {
                     curr_plate = -1;
                     for( j = 0; j < 2 * N_PLATES; j += 2)
                        if( plate_info[j] == *(long *)( star.plate))
                           curr_plate = j;
                     }
                  if( curr_plate >= 0)
                     curr_ptr->ppm_num = plate_info[curr_plate + 1];
                  else
                     curr_ptr->ppm_num = 0;
                  curr_ptr->photo_band = 0;
                  if( star.mag_band <= 18)
                     {
//                   const char *gsc_to_photo_band = "BR BVRVBR VB    V B";
//                   const char *gsc_to_photo_band = "JV BVRVOER      J V";
                     const char *gsc_to_photo_band = "BV?BVRVBRR??????B?V";
                /* 11 Jan 2002:  changed so that 'unknown' bands don't     */
                /* lead to a crash,  and are simply marked as 'unknown'    */

                     curr_ptr->photo_band = gsc_to_photo_band[star.mag_band];
                     if( curr_ptr->photo_band == ' ')
                        curr_ptr->photo_band = 0;
                     else
                        for( j = 0; j < 3; j++)
                           if( "BRV"[j] == curr_ptr->photo_band)
                              n_per_band[j]++;
                     }
                  if( !curr_ptr->photo_band)
                     {
                     log_printf( log_file,
                           "??? Internal error on GSC %d %d: mag band %d\n",
                           area_ptr->small_zone, star.n, star.mag_band);
                     log_printf( log_file, "Please report this to Project Pluto");
                     exit_2( -22);
                     }
                  (*n_found)++;
                  if( *n_found >= n_alloced)
                     {
                     REF_STAR *new_arr;

                     n_alloced += n_alloced;
                     if( debug_level > 2)
                        {
                        log_printf( log_file, "\nReallocing to %ld x %d\n",
                                     n_alloced, sizeof( REF_STAR));
                        if( !debug_check_all_memory( ))
                           log_printf( log_file, "Memory checks out OK\n");
                        else
                           debug_dump_memory( );
                        }
                     new_arr = (REF_STAR *)malloc( n_alloced * sizeof( REF_STAR));
                     if( debug_level > 2)
                        log_printf( log_file, ".");
                     if( new_arr)
                        {
                        memcpy( new_arr, ret_arr, *n_found * sizeof( REF_STAR));
                        if( debug_level > 2)
                           log_printf( log_file, ":");
                        free( ret_arr);
                        if( debug_level > 2)
                           log_printf( log_file, ";");
                        ret_arr = new_arr;
                        }
                     else
                        {
                        if( debug_level)
                          log_printf( log_file, "\nReallocing ");
                        ret_arr = (REF_STAR *)realloc( ret_arr,
                                     n_alloced * sizeof( REF_STAR));
                        if( debug_level)
                           log_printf( log_file, ret_arr ? "successfully\n" : "failed\n");
                        }
                     if( debug_level > 2)
                        log_printf( log_file, ",");
                     if( !ret_arr)
                        {
                        *n_found = -5;
                        return( NULL);
                        }
                     }
                  }
               }
         }
      if( hdr.n_ppm_stars)
         {
         SMALL_PPM_STAR tstar;
         int k, gsc_number;

         fseek( data_file, loc[0] + hdr.filesize, SEEK_SET);
         for( j = 0; j < hdr.n_ppm_stars; j++)
            {
            fread( &tstar, 1, sizeof( SMALL_PPM_STAR), data_file);
            gsc_number = (int)tstar.loc[1];
            if( tstar.ppm_num > 1000000000L)
               tstar.ppm_num -= 1000000000L;
            if( tstar.ppm_num > 410000L)
               tstar.ppm_num %= 410000L;
            for( k = 0; k < (*n_found); k++)
               if( gsc_number == ret_arr[k].number &&
                                       ret_arr[k].zone == area_ptr->small_zone)
                  ret_arr[k].ppm_num = tstar.ppm_num;
            }
         }
      }
   if( debug_level)
      log_printf( log_file, "Freeing, ");
   free( plate_info);
   free( buff);
   free( areas);
   fclose( data_file);
   fclose( index_file);
   if( use_gsc_act)
      {
      fclose( gsc_act_file);
      free( plate_cache);
      }
           /* From which band (B, R, V) do the majority of magnitudes come? */
   for( i = j = 0; j < 3; j++)
      if( n_per_band[j] > n_per_band[i])
         i = j;
            /* OK,  that's the band we want to use: */
   gsc_photo_band = "BRV"[i];
   log_printf( log_file, strings[113], gsc_photo_band);
               /* "GSC magnitudes are mostly in photo band %c\n" */


   if( debug_level)
      log_printf( log_file, "done\n");
   return( ret_arr);
}

double get_tycho_mag( const double vt, const double bt, const char photo_band,
                        char *iband)
{
   double b_v_j = johnson_b_minus_v_from_tycho_b_minus_v( bt - vt);
   double v_j = johnson_v_from_tycho_b_minus_v( bt - vt, vt);
   double mag;

   if( photo_band == 'V')
      mag = v_j;
   else if( photo_band == 'B')
      mag = v_j + b_v_j;
   else        /* compute R from B-V: */
      mag = v_j - b_minus_v_to_v_minus_r( b_v_j);
   if( mag > 0. && mag < 20.)
      *iband = photo_band;
   else
      {
      if( v_j < 20.)
         {
         mag = v_j;
         *iband = 'V';
         }
      else if( v_j < 20. && b_v_j < 20.)
         {
         mag = v_j + b_v_j;
         *iband = 'B';
         }
      else if( vt > 0. && vt < 20.)
         {
         mag = vt;
         *iband = 'v';
         }
      else
         {
         mag = bt;
         *iband = 'b';
         }
      }
   return( mag);
}

char *load_up_act_data( FILE *ifile, int *n_found, int gsc_zone)
{
   long loc[2];
   char *rval;

   if( !ifile)
      return( NULL);
   fseek( ifile, (long)( gsc_zone - 1) * 4L, SEEK_SET);
   fread( loc, 2, sizeof( long), ifile);
   fseek( ifile, 9538L * 4L + loc[0] * 10L, SEEK_SET);
   *n_found = (int)( loc[1] - loc[0]);
   rval = (char *)calloc( *n_found, 10);
   fread( rval, 10, *n_found, ifile);
   return( rval);
}

REF_STAR *grab_tycho_data( const double x, const double y,
                   const double width, const double height,
                   const char *cd_drive_path, int *n_found, double jd,
                   const int photo_band)
{
   int n_areas, i;
   double minmax[4];
   SMALL_AREA *areas;
   FILE *data_file = NULL, *act_file;
   REF_STAR *ret_arr = (REF_STAR *)calloc( 10, sizeof( REF_STAR));
   int n_alloced = 10, version = 6;
   double tycho_pm_factor, tycho_epoch = TYCHO_EPOCH;
   long filesize;

   if( jd < TYCHO_EPOCH)         /* unreasonably small date */
      jd = 2451000.;
   *n_found = 0;

   areas = find_small_areas( x, y, width, height, &n_areas, minmax);
   if( !areas)
      {
      *n_found = -1;
      return( NULL);
      }

   data_file = local_then_cd_fopen( cd_drive_path, "hipp\\lg_tycho.lmp");
   if( !data_file)
      {
      *n_found = -3;
      return( NULL);
      }

   fseek( data_file, 0L, SEEK_END);
   filesize = ftell( data_file);
   if( filesize > 40000000L)
      version = 7;
   if( filesize > 60000000L)
      {
      tycho_epoch = 2451545.;       /* Tycho-2 data is J2000,  not J1991.25 */
      version = 8;
      }
   fseek( data_file, 0L, SEEK_SET);
   tycho_pm_factor = (jd - tycho_epoch) / (3600 * 1000. * 100.);
                     /* cvt from .01 milliarcsec units to degrees */
   tycho_pm_factor /= 365.25;       /* days to years */

   act_file = fopen( "act.dat", "rb");
   if( act_file)
      log_printf( log_file, "Getting proper motions from ACT\n");
   for( i = 0; i < n_areas; i++)
      {
      SMALL_AREA *area_ptr = areas + i;
      int n_act, j;
      char *buff, *tptr;
      char *act_data = load_up_act_data( act_file,
                              &n_act, area_ptr->small_zone);
      long loc[2], bytes_left;

      fseek( data_file, 4L * ((long)area_ptr->small_zone - 1L), SEEK_SET);
      fread( loc, 2, sizeof( long), data_file);
      fseek( data_file, loc[0], SEEK_SET);
      bytes_left = loc[1] - loc[0];
      tptr = buff = (char *)malloc( bytes_left);
      fread( buff, bytes_left, 1, data_file);
      while( bytes_left)
         {
         TYCHO_STAR star;
         int amt_read = parse_tycho_star( tptr, &star, version);
         REF_STAR *curr_ptr = ret_arr + (*n_found);
         double ra, spd;
         char iband;

         ra = tycho_pm_factor * (double)star.pm_vals[0];
         spd = tycho_pm_factor * (double)star.pm_vals[1];
         if( act_data)
            for( j = 0; j < n_act; j++)
               if( star.gsc_num == *(short *)( act_data + j * 10 + 8))
                  {
                  long pm_packed = *(long *)( act_data + j * 10);
                  double pm_ra  = (double)(( pm_packed % 101000L) - 44764L);
                  double pm_dec = (double)( pm_packed / 101000L);
                  short pm_dec_sig = *(short *)( act_data + j * 10 + 6);

                  if( pm_dec_sig & 0x4000)     /* flag for negative pm dec */
                     pm_dec = -pm_dec;
                  pm_ra *= 15. * sin( (double)star.spd * PI / 180.e+7);
                  pm_dec *= 10.;
                              /* both are now in milliarcsecs/century */
                  ra = tycho_pm_factor * pm_ra;
                  spd = tycho_pm_factor * pm_dec;
                  j = n_act;
                  }
                              /* RA/decs are stored in 1e-7 degree units: */
         ra += (double)star.ra * 1.e-7;
         spd += (double)star.spd * 1.e-7;
         curr_ptr->ra = ra;
         curr_ptr->dec = spd - 90.;
         curr_ptr->mag = star.mag;
         curr_ptr->mag = (int)( 100. * get_tycho_mag(
                            (double)star.vt_mag / 1000.,
                            (double)star.bt_mag / 1000., photo_band, &iband));
         curr_ptr->photo_band = iband;
         curr_ptr->zone = star.gsc_zone;
         curr_ptr->number = star.gsc_num;
         curr_ptr->ppm_num = star.ppm_num;
         if( in_bounds( &curr_ptr->ra, curr_ptr->dec, minmax) &&
               curr_ptr->mag < gsc_mag_limit &&
               curr_ptr->mag > gsc_bright_mag_limit)
            (*n_found)++;
         if( *n_found >= n_alloced)
            {
            REF_STAR *new_arr;

            n_alloced += n_alloced / 2;
            new_arr = (REF_STAR *)calloc( n_alloced, sizeof( REF_STAR));
            if( new_arr)
               {
               memcpy( new_arr, ret_arr, *n_found * sizeof( REF_STAR));
               free( ret_arr);
               ret_arr = new_arr;
               }
            else
               ret_arr = (REF_STAR *)realloc( ret_arr,
                                      n_alloced * sizeof( REF_STAR));
            if( !ret_arr)
               {
               *n_found = -5;
               return( NULL);
               }
            }
         tptr += amt_read;
         bytes_left -= amt_read;
         }
      free( buff);
      if( act_data)
         free( act_data);
      }
   if( act_file)
      fclose( act_file);
   free( areas);
   fclose( data_file);
   return( ret_arr);
}

static int parse_tycho2_rec( const char *buff, double *tycho2_data)
{
   int rval = (buff[18] != '.');

   if( !rval)
      {
      *tycho2_data++ = atof( buff + 15);     /* RA */
      *tycho2_data++ = atof( buff + 28);     /* dec */
      *tycho2_data++ = atof( buff + 41);     /* PM in RA */
      *tycho2_data++ = atof( buff + 49);     /* PM in dec */
      *tycho2_data++ = atof( buff + 110);     /* BT */
      *tycho2_data++ = atof( buff + 123);     /* VT */
//    *tycho2_data++ = atof( buff + 75);     /* epoch year in RA */
//    *tycho2_data++ = atof( buff + 83);     /* epoch year in dec */
      }
   return( rval);
}

REF_STAR *grab_tycho2_data( const double x, const double y,
                   const double width, const double height,
                   const char *tycho2_filename, int *n_found, double jd,
                   const int photo_band)
{
   int n_areas, i;
   double minmax[4];
   SMALL_AREA *areas;
   FILE *data_file = NULL, *index_file;
   REF_STAR *ret_arr;
   int n_alloced = 10, tycho_line_size;
   char buff[210];

   if( tycho2_filename[1])
      data_file = fopen( tycho2_filename, "rb");

   if( !data_file)
      {
      sprintf( buff, "%c:\\data\\catalog.dat", *tycho2_filename);
      data_file = fopen( buff, "rb");
      }

   if( !data_file)
      {
      *n_found = -4;
      return( NULL);
      }

        /* The catalogue on the CD-ROM is in DOS (CR/LF) format,  for 208 */
        /* bytes per line.  The data on the CDS site is in Unix (LF only) */
        /* format,  for 207 bytes/line.  The following code figures out   */
        /* with which type we're dealing,  and also verifies that we're   */
        /* really using a Tycho-2 file.                                   */

   fgets( buff, sizeof( buff), data_file);
   tycho_line_size = strlen( buff);

   if( tycho_line_size < 207 || tycho_line_size > 208 ||
            memcmp( buff, "0001 00008 1| |  2.31750494|  2.23184345|", 40))
      {
      *n_found = -2;
      fclose( data_file);
      return( NULL);
      }

   index_file = fopen( "tycho2.idx", "rb");
   if( !index_file)
      {
      *n_found = -6;
      fclose( data_file);
      return( NULL);
      }

   if( jd < TYCHO_EPOCH)         /* unreasonably small date */
      jd = TYCHO_EPOCH;
   *n_found = 0;

   areas = find_small_areas( x, y, width, height, &n_areas, minmax);
   if( !areas)
      {
      *n_found = -1;
      return( NULL);
      }

   ret_arr = (REF_STAR *)calloc( 10, sizeof( REF_STAR));

   for( i = 0; i < n_areas; i++)
      {
      long loc[2], recs_left;
      double tycho2_data[6];
      const double J2000 = 2451545.;
      const double pm_factor = (jd - J2000) / (365.25 * 3600000.);

      fseek( index_file, 4L * ((long)areas[i].small_zone - 1L), SEEK_SET);
      fread( loc, 2, sizeof( long), index_file);
      fseek( data_file, loc[0] * (long)tycho_line_size, SEEK_SET);
      recs_left = loc[1] - loc[0];
      while( recs_left-- && fread( buff, tycho_line_size, 1, data_file))
         if( !parse_tycho2_rec( buff, tycho2_data))
            {
            REF_STAR *curr_ptr = ret_arr + (*n_found);
            double ra, dec, aspect = cos( tycho2_data[1] * PI / 180.), mag;
            char iband;

            ra = tycho2_data[0] + tycho2_data[2] * pm_factor / aspect;
            dec = tycho2_data[1] + tycho2_data[3] * pm_factor;

            curr_ptr->ra = ra;
            curr_ptr->dec = dec;
            mag = get_tycho_mag( tycho2_data[5], tycho2_data[4], photo_band,
                                                     &iband);
            curr_ptr->mag = (int)( mag * 100.);
            curr_ptr->zone = atoi( buff);
            curr_ptr->number = atoi( buff + 4);
            curr_ptr->ppm_num = 0;
            curr_ptr->photo_band = iband;
            if( in_bounds( &curr_ptr->ra, curr_ptr->dec, minmax) &&
                  curr_ptr->mag < gsc_mag_limit &&
                  curr_ptr->mag > gsc_bright_mag_limit)
               (*n_found)++;
            buff[70] = '\0';
            if( *n_found >= n_alloced)
               {
               REF_STAR *new_arr;

               n_alloced += n_alloced / 2;
               new_arr = (REF_STAR *)calloc( n_alloced, sizeof( REF_STAR));
               if( new_arr)
                  {
                  memcpy( new_arr, ret_arr, *n_found * sizeof( REF_STAR));
                  free( ret_arr);
                  ret_arr = new_arr;
                  }
               else
                  ret_arr = (REF_STAR *)realloc( ret_arr,
                                      n_alloced * sizeof( REF_STAR));
               if( !ret_arr)
                  {
                  *n_found = -5;
                  return( NULL);
                  }
               }
            }
      }

   free( areas);
   fclose( data_file);
   fclose( index_file);
   return( ret_arr);
}

static void parse_ucac_data( const char *buff, long *ucac_data, int compressed)
{
   if( !compressed)
      sscanf( buff, "%ld %ld %ld", ucac_data, ucac_data + 1, ucac_data + 2);
   else
      {
      unsigned long tval = ((unsigned long *)buff)[2];

      ucac_data[0] = ((long *)buff)[0];
      ucac_data[1] = ((long *)buff)[1];
      ucac_data[2] = (long)( tval % (unsigned long)1146) + 570L;
      }
}

REF_STAR *grab_ucac1_data( double x, const double y,
                   const double width, const double height,
                   const char *ucac_path, int *n_found, const double jd)
{
   const long ymin = (long)( (y - height + 90.) * 3600000.);  /* use milliarcsecs */
   const long ymax = (long)( (y + height + 90.) * 3600000.);
   REF_STAR *ret_arr = NULL;
   int n_alloced = 10;

   *n_found = 0;
   for( x -= 360.; x - width < 360.; x += 360.)
      if( x + width > 0.)
         {
         const long xmin = (long)( (x - width) * 3600000.);
         const long xmax = (long)( (x + width) * 3600000.);
         long yzone;

         for( yzone = ymin / 1800000L; yzone <= ymax / 1800000L; yzone++)
            {
            char buff[100];
            FILE *ifile;
            long recsize = 0L;

            sprintf( buff, "%sz%03ld", ucac_path, yzone + 1L);
            ifile = fopen( buff, "rb");
            if( !ifile)
               {
               strcat( buff, ".uca");   /* look for compressed version */
               ifile = fopen( buff, "rb");
               recsize = 20L;
               }
            if( !ifile)
               log_printf( log_file, "UCAC file '%s' not found\n", buff);
            else
               {
               long loc = 0L, step, n_recs, ucac_data[3];

               log_printf( log_file, "UCAC file '%s' opened\n", buff);
               if( !recsize)
                  {
                  fgets( buff, sizeof( buff), ifile);
                  recsize = ftell( ifile);
                  }
               fseek( ifile, 0L, SEEK_END);
               n_recs = ftell( ifile) / recsize;
               for( step = 0x40000000; step; step >>= 1)
                  if( loc + step < n_recs)
                     {
                     fseek( ifile, recsize * (loc + step), SEEK_SET);
                     fread( buff, recsize, 1, ifile);
                     parse_ucac_data( buff, ucac_data, recsize == 20L);
                     if( ucac_data[0] < xmin)
                        loc += step;
                     }
               memset( buff, 0, (size_t)recsize);
               *buff = '\0';
               fseek( ifile, loc * recsize, SEEK_SET);
               while( ucac_data[0] < xmax && fread( buff, recsize, 1, ifile))
                  {
                  long ra, spd, mag;

                  parse_ucac_data( buff, ucac_data, recsize == 20L);
                  ra = ucac_data[0];
                  spd = ucac_data[1];
                  mag = ucac_data[2];
                  loc++;
                  if( ra > xmin && ra < xmax && spd > ymin && spd < ymax &&
                            mag < gsc_mag_limit &&
                            mag > gsc_bright_mag_limit)
                     {
                     if( !ret_arr)
                        ret_arr = (REF_STAR *)calloc( n_alloced,
                                                         sizeof( REF_STAR));
                     if( *n_found == n_alloced)
                        {
                        REF_STAR *new_arr;

                        n_alloced += n_alloced / 2;
                        new_arr = (REF_STAR *)calloc( n_alloced,
                                                         sizeof( REF_STAR));
                        if( new_arr)
                           {
                           memcpy( new_arr, ret_arr,
                                              *n_found * sizeof( REF_STAR));
                           free( ret_arr);
                           ret_arr = new_arr;
                           }
                        else
                           ret_arr = (REF_STAR *)realloc( ret_arr,
                                               n_alloced * sizeof( REF_STAR));
                        if( !ret_arr)
                           {
                           *n_found = -5;
                           return( NULL);
                           }
                        }
                     ret_arr[*n_found].mag = mag;
                     ret_arr[*n_found].ra = (double)ra / 3600000.;
                     ret_arr[*n_found].dec = (double)spd / 3600000. - 90.;
                     ret_arr[*n_found].ppm_num = loc;
                     ret_arr[*n_found].zone = (short)( yzone + 1);
                     (*n_found)++;
                     }
                  }
               fclose( ifile);
               }
            }
         }
   return( ret_arr);
}

void remove_trailing_goo( char *buff)
{
   while( *buff && *buff != 10 && *buff != 13)
      buff++;
   *buff = '\0';
}

void get_guide_environment_ptr( char *buff, const char *env_var)
{
   FILE *ifile = fopen( "guide.dat", "rb");

   *buff = '\0';
   if( ifile)
      {
      char tbuff[200];
      const int len = strlen( env_var);

      while( fgets( tbuff, sizeof( tbuff), ifile))
         if( tbuff[len] == '=' && !memcmp( tbuff, env_var, len))
            {
            strcpy( buff, tbuff + len + 1);
            remove_trailing_goo( buff);
            }
      fclose( ifile);
      }
}

REF_STAR *grab_a10_data( const double x, const double y,
               const double width, const double height,
               const int cd_drive_letter, int *n_found, const int photo_band)
{
   const char *a10_name = "a10.dat";
   FILE *ifile;
   long hdr[8], hdr_out[8];
   int data_already_extracted = 0;
   int pass;
   REF_STAR *rval = NULL;

   log_printf( log_file, "CD drive letter '%c' (%d)\n", cd_drive_letter, cd_drive_letter);
   hdr_out[0] = hdr_out[4] = -1L;
   hdr_out[1] = (long)( (x - width) * 360000.);
   hdr_out[2] = (long)( (x + width) * 360000.);
   hdr_out[5] = (long)( (y - height + 90.) * 360000.);
   hdr_out[6] = (long)( (y + height + 90.) * 360000.);
   hdr_out[3] =  1007775886L;       /* magic numbers */
   hdr_out[7] =  2076665750L;       /* magic numbers */
   ifile = fopen( a10_name, "rb");
   if( ifile)
      if( fread( hdr, 8, sizeof( long), ifile))
         if( hdr[0] == -1 && hdr[3] == hdr_out[3] &&
                      hdr[1] < hdr_out[1] && hdr[2] > hdr_out[2] &&
                      hdr[5] < hdr_out[5] && hdr[6] > hdr_out[6])
            data_already_extracted = 1;
                           /* The above means "if a10.dat already exists, */
                           /* has the right magic numbers, and it covers  */
                           /* the area of interest,  then don't bother    */
                           /* re-extracting data;  we're all set." But if */
                           /* we're _not_ ok,  we extract the data,  plus */
                           /* .1 degree of margin on all edges.           */
   if( !data_already_extracted)
      {
      char ax0_path[90], tbuff[90];
      FILE *ofile;
      extern int using_a20;

      if( ifile)
         fclose( ifile);
      sprintf( tbuff, "A%d_PATH", 1 + using_a20);
      get_guide_environment_ptr( ax0_path, tbuff);
      if( *ax0_path)
         log_printf( log_file, "Hunting for Ax.0 data in path '%s'\n", ax0_path);

      memcpy( hdr, hdr_out, 8 * sizeof( long));
      hdr[1] -= 36000L;
      hdr[2] += 36000L;
      hdr[5] -= 36000L;
      hdr[6] += 36000L;
      ofile = fopen( a10_name, "wb");
      fwrite( hdr, 8, sizeof( long), ofile);
      sprintf( tbuff, "%c:\\", cd_drive_letter);
      dos_prompt_extract_a10_data( tbuff, ax0_path, x, y,
                     (width + .1) * 2., (height + .1) * 2., ofile);
      fclose( ofile);
      ifile = fopen( a10_name, "rb");        /* let's try again... */
      }

   for( pass = 0; pass < 2; pass++)
      {
      long ivals[4];
      int count = 0, mag, r_mag, b_mag;
      char iband;

      fseek( ifile, 32L, SEEK_SET);
      while( fread( ivals, 4, sizeof( long), ifile))
         {
         r_mag = (int)( labs( ivals[3]        ) % 1000L) * 10;
         b_mag = (int)( labs( ivals[3] / 1000L) % 1000L) * 10;
         mag = 0;                /* default value */
         if( photo_band == 'B')
            mag = b_mag;
         else if( photo_band != 'V')
            mag = r_mag;
         else
            if( r_mag > 1 && r_mag < 9999 && b_mag > 1 && b_mag < 9999)
               mag = (r_mag * 5 + b_mag * 3 + 4) / 8;
         if( mag)
            iband = photo_band;
         else
            {
            if( b_mag)
               {
               mag = b_mag;
               iband = 'B';
               }
            else
               {
               mag = r_mag;
               iband = 'R';
               }
            }
         if( ivals[1] >= hdr_out[1] && ivals[1] <= hdr_out[2]
                                    && ivals[2] >= hdr_out[5]
                                    && ivals[2] <= hdr_out[6]
                                    && mag > gsc_bright_mag_limit
                                    && mag < gsc_mag_limit && mag)
            {
            if( rval)
               {
               rval[count].ra  = (double)ivals[1] / 360000.;
               rval[count].dec = (double)ivals[2] / 360000. - 90.;
               rval[count].mag = mag;
               rval[count].photo_band = iband;
               rval[count].ppm_num = ivals[0];
               }
            count++;
            }
         }
      *n_found = count;
      if( count && !pass)
         {
         rval = (REF_STAR *)calloc( count, sizeof( REF_STAR));
         if( !rval)
            {
            log_printf( log_file, "RAN OUT OF MEMORY:  Too many A1.0 stars in this area\n");
            log_printf( log_file, "Choose a smaller area or select a brighter mag limit,\n");
            log_printf( log_file, "and try again.\n");
            exit_2( -21);
            }
         }
      }
   fclose( ifile);
   if( !rval || *n_found <= 0)
      rval = grab_gsc_12_data( x, y, width, height, "eso_a2.dat", n_found, 0.);
                  /* Send a zero JD because Ax.0 lacks proper motions,  so */
                  /* the date doesn't matter anyway */
   return( rval);
}

REF_STAR *grab_ucac2_data( const double x, const double y,
               const double width, const double height,
               const int cd_drive_letter, int *n_found, const int photo_band,
               const double jd)
{
   char ucac2_path[80];
   REF_STAR *rval = NULL;

   get_guide_environment_ptr( ucac2_path, "UCAC2_PATH");
   if( *ucac2_path)
      {
      const char *ucac2_filename = "ucac2.txt";
      int retval;
      FILE *ofile = fopen( ucac2_filename, "w");

      retval = extract_ucac2_stars( ofile, x, y,
                               width * 2., height * 2., ucac2_path, 0);
      extract_ucac2_stars( ofile, x, y,         /* add in supplement stars */
                               width * 2., height * 2., ucac2_path, 1);
      fclose( ofile);
      if( retval >= 0)
         rval = grab_gsc_12_data( x, y, width, height, ucac2_filename,
                            n_found, jd);
      else
         {
         log_printf( log_file, "UCAC-2 rval: %d\n", retval);
         *n_found = -1;
         }
      }
   else
      *n_found = -2;
   if( *n_found <= 0)
      {
      rval = grab_gsc_12_data( x, y, width, height, "u2.dat", n_found, jd);
      log_printf( log_file, "Looking for 'u2.dat': rval %p, n_found %d\n", rval, *n_found);
      }
   return( rval);
}

REF_STAR *grab_ucac3_data( const double x, const double y,
               const double width, const double height,
               const int cd_drive_letter, int *n_found,
               const int photo_band, const double jd)
{
   char ucac3_path[80];
   REF_STAR *rval = NULL;

   get_guide_environment_ptr( ucac3_path, "UCAC3_PATH");
   if( *ucac3_path)
      {
      const char *ucac3_filename = "ucac3.txt";
      int retval;
      FILE *ofile = fopen( ucac3_filename, "w");

      retval = extract_ucac3_stars( ofile, x, y,
                               width * 2., height * 2., ucac3_path, 0, 0);
      fclose( ofile);
      if( retval >= 0)
         rval = grab_gsc_12_data( x, y, width, height, ucac3_filename,
                            n_found, jd);
      else
         *n_found = -1;
      }
   else
      *n_found = -2;
#ifdef NO_UCAC3_VIA_VIZIER_YET
   if( *n_found <= 0)
      {
      rval = grab_gsc_12_data( x, y, width, height, "u2.dat", n_found, jd);
      log_printf( log_file, "Looking for 'u2.dat': rval %p, n_found %d\n", rval, *n_found);
      }
#endif
   return( rval);
}

REF_STAR *grab_ucac4_data( const double x, const double y,
               const double width, const double height,
               const int cd_drive_letter, int *n_found,
               const int photo_band, const double jd)
{
   char ucac4_path[180];
   REF_STAR *rval = NULL;

   get_guide_environment_ptr( ucac4_path, "UCAC4_PATH");
   if( *ucac4_path)
      {
      const char *ucac4_filename = "ucac4.txt";
      int retval;
      FILE *ofile = fopen( ucac4_filename, "w");

      retval = extract_ucac4_stars( ofile, x, y,
                               width * 2., height * 2., ucac4_path, 0);
      fclose( ofile);
      if( retval >= 0)
         rval = grab_gsc_12_data( x, y, width, height, ucac4_filename,
                            n_found, jd);
      else
         {
         log_printf( log_file, "UCAC4 error %d\n", retval);
         *n_found = -1;
         }
      log_printf( log_file, "rval for UCAC4 = %p; %d found\n", rval, *n_found);
      }
   else
      *n_found = -2;
#ifdef NO_UCAC4_VIA_VIZIER_YET
   if( *n_found <= 0)
      {
      rval = grab_gsc_12_data( x, y, width, height, "u2.dat", n_found, jd);
      log_printf( log_file, "Looking for 'u2.dat': rval %p, n_found %d\n", rval, *n_found);
      }
#endif
   return( rval);
}
