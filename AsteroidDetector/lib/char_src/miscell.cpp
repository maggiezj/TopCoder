#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include "watdefs.h"
#include "charon.h"
#include "findstar.h"
#include "miscell.h"
#include "date.h"

#ifdef DEBUG_MEM
#include "checkmem.h"
#endif

/* MISCELL.CPP contains a whole slew of functions that were formerly part
   of CHARON.CPP,  but which really ought to work on any OS ('generic
   ANSI C' code.)  It's a bunch of basic utilities to handle grunt work
   such as loading and storing configuration files,  formatting RA/dec
   values,  and so on.  */

#define PI 3.14159265358979323846264338327

int psf_fit_star( FOUND_STAR *found, int cell_size);       /* findstar.c */
int image_pixel_separation = 1, psf_fitting_on = 1;
extern int gsc_mag_limit, gsc_bright_mag_limit, debug_level;
extern FILE *log_file;
int suppress_display_of_residuals = 0;
int language = 'E';
int ignored_north = 2, ignored_south = 2, ignored_east = 2, ignored_west = 2;
char *strings[200];
char *string_file_name = "cstrings.dat";

int reset_language( const int new_language);          /* miscell.cpp */
FILE *local_then_cd_fopen( const char *cd_drive_path, const char *filename);
void exit_2( const int exit_code);                    /* charon.cpp */

static void find_clipped_rect( const IMAGE *image, const int x, const int y,
                                       const int cell_size, int *rect)
{
   const int half_size = cell_size / 2;

   rect[0] = (x > half_size ? x - half_size : 0);
   rect[1] = x + half_size + 1;
   if( rect[1] > image->xsize)
      rect[1] = image->xsize;

   rect[2] = (y > half_size ? y - half_size : 0);
   rect[3] = y + half_size + 1;
   if( rect[3] > image->ysize)
      rect[3] = image->ysize;
}

static PIXEL find_nth_pixel( const IMAGE *image, int x, int y,
                   const int xsize, const int ysize, const int n)
{
   int i, j;
   PIXEL curr_level = 0, bit, *pix_ptr;
   PIXEL max_pixel, min_pixel;

   x -= xsize / 2;         /* move to upper left corner... */
   y -= ysize / 2;
   if( x < 0)              /* then make sure cell is entirely on-screen: */
      x = 0;
   if( y < 0)
      y = 0;
   if( x > image->xsize - xsize)
      x = image->xsize - xsize;
   if( y > image->ysize - ysize)
      y = image->ysize - ysize;

            /* find max and min pixel value within the area: */
   pix_ptr = image->img + x + y * image->xsize;
   max_pixel = min_pixel = *pix_ptr;
   for( j = 0; j < ysize; j++, pix_ptr += image->xsize)
      for( i = 0; i < xsize; i++)
         if( max_pixel < pix_ptr[i])
            max_pixel = pix_ptr[i];
         else if( min_pixel > pix_ptr[i])
            min_pixel = pix_ptr[i];

   for( bit = ((PIXEL)1 << 30); bit; bit >>= 1)
      {
      curr_level |= bit;
      if( curr_level > max_pixel)      /* gone too far... */
         curr_level ^= bit;            /* ...so remove that bit  */
      else if( curr_level >= min_pixel)         /* gotta check in this case: */
         {
         int n_above_thresh = 0;

         pix_ptr = image->img + x + y * image->xsize;
         for( j = 0; j < ysize && n_above_thresh < n;
                                        j++, pix_ptr += image->xsize)
            for( i = 0; i < xsize && n_above_thresh < n; i++)
               if( pix_ptr[i] > curr_level)
                  n_above_thresh++;

         if( n_above_thresh < n)     /* oops!  went too far... */
            curr_level ^= bit;            /* ...so remove that bit  */
         }
      }
   return( curr_level);
}

#ifdef OBSOLETE_NOW
static PIXEL find_threshhold( const IMAGE *image, const int *rect,
                              const double cutoff_point, const PIXEL max_pixel)
{
   int i, j;
   PIXEL curr_level = 0, bit = ((PIXEL)1 << 30);
   const int thresh = (int)( (double)( rect[1] - rect[0]) * (1. - cutoff_point)
                           * (double)( rect[3] - rect[2]) + .5);

   while( bit > max_pixel)
      bit >>= 1;
   for( ; bit; bit >>= 1)
      {
      int n_above_thresh = 0;

      curr_level |= bit;
      for( j = rect[2]; j < rect[3] && n_above_thresh < thresh; j++)
         {
         PIXEL *pix = image->img + j * image->xsize + rect[0];

         for( i = rect[1] - rect[0]; i && n_above_thresh < thresh; i--)
            if( *pix++ > curr_level)
               n_above_thresh++;
         }

      if( n_above_thresh < thresh)     /* oops!  went too far... */
         curr_level ^= bit;            /* ...so remove that bit  */
      }
   return( curr_level);
}
#endif

#define BACK_CELL_SIZE 13

int add_a_star( const IMAGE *image, int x, int y, FOUND_STAR *f,
                      const double cutoff_point, const int cell_size)
{
   int i, j, n_pixels = 0;
   double xcenter = 0., ycenter = 0.;
   long total = 0L;
   PIXEL max_pix = 0, pix, thresh;
   int rect[4];
   FILE *debug_file = fopen( "debug2.dat", "ab");

   find_clipped_rect( image, x, y, abs( cell_size), rect);
   fprintf( debug_file, "Input loc: %d %d; cell size %d\n", x, y, cell_size);

                        /* first,  find the maximum pixel: */
   for( j = rect[2]; j < rect[3]; j++)
      for( i = rect[0]; i < rect[1]; i++)
         if( (pix = image->img[i + j * image->xsize]) > max_pix)
            {
            x = i;
            y = j;
            max_pix = pix;
            }
   fprintf( debug_file, "Moved to: %d %d\n", x, y);


                         /* next,  find threshhold value */
   thresh = find_nth_pixel( image, x, y, BACK_CELL_SIZE, BACK_CELL_SIZE,
                        (int)( (1. - cutoff_point) *
                        (double)( BACK_CELL_SIZE * BACK_CELL_SIZE) + .5));

   fprintf( debug_file, "Output loc: %d %d;  thresh %ld\n", x, y, thresh);
// find_clipped_rect( image, x, y, abs( cell_size), rect);
   for( j = rect[2]; j < rect[3]; j++)
      for( i = rect[0]; i < rect[1]; i++)
         fprintf( debug_file, "%5d%s", image->img[i + j * image->xsize],
                        (i == rect[1] - 1) ? "\n" : " ");

   for( j = rect[2]; j < rect[3]; j++)
      {
      const PIXEL *pix = image->img + j * image->xsize + rect[0];

      for( i = rect[0]; i < rect[1]; i++, pix++)
         if( *pix > thresh)
            {
            const PIXEL diff = *pix - thresh;

            xcenter += (double)( diff * i);
            ycenter += (double)( diff * j);
            total += (long)diff;
            n_pixels++;
            }
      }

   if( n_pixels)
      {
      if( cell_size > 0)
         {
         f->x = xcenter / (double)total;
         f->y = ycenter / (double)total;
         }
      else        /* use Bill Owen's algorithm for comet centroiding: */
         {
         const PIXEL *center_row = image->img + y * image->xsize + x;
         const PIXEL *upper_row = center_row - image->xsize;
         const PIXEL *lower_row = center_row + image->xsize;
         PIXEL s = lower_row[-1], upper, lower, left, right, center;
         double denominator;

         if( s < lower_row[1])
            s = lower_row[1];
         if( s < upper_row[1])
            s = upper_row[1];
         if( s < upper_row[-1])
            s = upper_row[-1];
         left   = (center_row[-1] > s ? center_row[-1] - s : 0);
         center = (center_row[ 0] > s ? center_row[ 0] - s : 0);
         right  = (center_row[ 1] > s ? center_row[ 1] - s : 0);
         upper = (*upper_row > s ? *upper_row - s : 0);
         lower = (*lower_row > s ? *lower_row - s : 0);
         denominator = (double)left + (double)center + (double)right
                     + (double)upper + (double)lower;
         f->x = ((double)right - (double)left ) / denominator;
         f->y = ((double)lower - (double)upper) / denominator;
#ifdef A_WAY_THAT_DIDNT_WORK_WELL
         double a = (double)upper_row[-1] + lower_row[-1] + center_row[-1];
         double b = (double)upper_row[ 0] + lower_row[ 0] + center_row[ 0];
         double c = (double)upper_row[ 1] + lower_row[ 1] + center_row[ 1];

         f->x = .5 * (c - a) / (2 * b - a - c);

         a = (double)upper_row[-1]  + upper_row[0]  + upper_row[1];
         b = (double)center_row[-1] + center_row[0] + center_row[1];
         c = (double)lower_row[-1]  + lower_row[0]  + lower_row[1];
         f->y = .5 * (c - a) / (2 * b - a - c);
#endif

#ifdef NOT_VERY_GOOD_WAY
         PIXEL s = lower_row[-1];
         double denominator;

         if( s < lower_row[1])
            s = lower_row[1];
         if( s < upper_row[1])
            s = upper_row[1];
         if( s < upper_row[-1])
            s = upper_row[-1];
         denominator = (double)*lower_row + (double)*upper_row
                     + (double)center_row[-1] + (double)*center_row
                     + (double)center_row[1] - 5 * (double)s;
         f->x = ((double)center_row[1] - (double)center_row[-1]) / denominator;
         f->y = ((double)lower_row[0] - (double)upper_row[0]) / denominator;
#endif
         fprintf( debug_file, "Initial shifts: %lf %lf\n", f->x, f->y);
         f->x += (double)x;
         f->y += (double)y;
         }
      f->bright = total;
      if( cell_size > 0 && psf_fitting_on)
         {
         FOUND_STAR ftemp = *f;

         if( !psf_fit_star( &ftemp, cell_size))    /* yes,  PSFfing worked */
            *f = ftemp;
         }
      f->x += .5;       /* the '+.5' moves stars to _center_ of pixel */
      f->y += .5;
      }
   fprintf( debug_file, "%d pixels; %lf %lf\n", n_pixels, f->x, f->y);
   fclose( debug_file);
   return( n_pixels);
}

int print_to_screen = 1;

int log_printf( FILE *log_file, const char *format, ...)
{
   va_list marker;
   char tbuff[200];

   va_start( marker, format);
   vsprintf( tbuff, format, marker);
   if( log_file)
      fwrite( tbuff, strlen( tbuff), 1, log_file);
   if( print_to_screen)
      fwrite( tbuff, strlen( tbuff), 1, stdout);
   return( 0);
}

void sprintf_ra_dec( char *str, double ival, const int format, const int is_dec)
{
   long lval;
   int n_places = (format >> 3);

   if( is_dec)
      {
      *str++ = ((ival < 0.) ? '-' : '+');
      if( ival < 0.) ival = -ival;
      }
   else
      {
      if( ival < 0.)
         ival += 360.;
      if( ival > 360.)
         ival -= 360.;
      if( format)
         ival /= 15.;
      }

   switch( format & 7)
      {
      case 0:
         sprintf( str, "%.5lf", ival);
         break;
      case 1:           /* minutes to four places */
         lval = (long)( ival * 60. * 10000. + .5);
         sprintf( str, "%02ldh%02ld.%04ldm", lval / 600000L,
                  (lval / 10000L) % 60L, lval % 10000L);
         if( is_dec)
            {
            str[2] = ' ';
            str[10] = '\'';
            }
         break;
      case 2:        /* seconds to three places */
         ival *= 60. * 60. * 1000.;
         if( !n_places)
            lval = (long)( ival + .5);
         if( n_places == 1)
            lval = (long)( ival + 5.);
         if( n_places == 2)
            lval = (long)( ival + 50.);
         sprintf( str, "%02ldh%02ldm%02ld.%03lds", lval / 3600000L,
                  (lval / 60000L) % 60L, (lval / 1000) % 60L, lval % 1000L);
         if( is_dec)
            {
            str[2] = ' ';
            str[5] = '\'';
            str[11] = '"';
            str[12] = '\0';
            }
         break;
      }
}

void put_ra_dec_in_str( char *str, const double ra, const double dec,
                                                      const int format)
{
   sprintf_ra_dec( str, ra, format, 0);
   strcat( str, "  ");
   sprintf_ra_dec( str + strlen( str), dec, format, 1);
}

int write_out_startup( const char *filename,
                     const CHARON_CONFIG *cconfig, const IMAGE *image)
{
   FILE *startup_file = fopen( filename, "wb");
   char tbuff[160], name_buff[80];
   int i;

   sprintf( tbuff, "%lf %lf %s %lf %lf %d\n", cconfig->max_tilt,
             cconfig->cutoff_point, cconfig->cd_drive_path,
              cconfig->latlon[0], cconfig->latlon[1], cconfig->fit_order);
   fputs( tbuff, startup_file);

   sprintf( tbuff, "%d %d %lf %lf %d %lf %d %lf\n", cconfig->max_stars,
            cconfig->n_stars_used, cconfig->search_dist, cconfig->scale_tolerance,
            cconfig->video_mode, cconfig->time_offset,
            cconfig->saturation_point, cconfig->tilt_angle);
   fputs( tbuff, startup_file);

   strcpy( name_buff, cconfig->target_name);
   for( i = 0; name_buff[i]; i++)
      if( name_buff[i] == ' ')
         name_buff[i] = '_';

   sprintf( tbuff, "%d %d %d %lf %s %d %d %s %d\n", cconfig->ra_dec_format,
             cconfig->color_table_number, cconfig->targeting_on,
             cconfig->max_residual, "dummy",
             image->image_offset, cconfig->cell_size,
             name_buff, psf_fitting_on);

   fputs( tbuff, startup_file);

   sprintf( tbuff, "%s %s %s %d %d %d %d\n", cconfig->report_file_name,
             cconfig->obs_code, cconfig->path_name,
             cconfig->catalog_used, cconfig->file_time_point,
             cconfig->display_flipped, cconfig->pixel_origin_at_top_left);
   fputs( tbuff, startup_file);

   sprintf( tbuff, "%d %lf %.2lf %.2lf %d %d %d %d%c\n", image->is_inverted,
           image->focal_len / 25.4, image->pixel_xsize, image->pixel_ysize,
           cconfig->millisec_per_video_reset, gsc_mag_limit,
            cconfig->ccd_id, gsc_bright_mag_limit, cconfig->photometric_band);
   fputs( tbuff, startup_file);

   sprintf( tbuff, "%d %d %d %d\n", ignored_north, ignored_south,
                     ignored_east, ignored_west);
   fputs( tbuff, startup_file);

   fclose( startup_file);
   return( 0);
}

      /* This function parses angles given in base 60 as six digits plus */
      /* an optional sign and/or fraction,  such as "123456.7" for       */
      /* "12h 34m 56.7s".                                                */

static int grab_base_sixty( const char *istr, double *oval)
{
   long seconds;
   int rval = 0, bytes_scanned;

   if( sscanf( istr, "%ld%n", &seconds, &bytes_scanned) == 1)
      if( bytes_scanned == 6)
         {
         seconds = (seconds / 10000L) * 3600L +
                   ((seconds / 100L) % 100L) * 60L + (seconds % 100L);
         *oval = (double)seconds;
         rval = 6;
         if( istr[rval] == '.')
            {
            *oval += atof( istr + rval);
            rval++;
            while( istr[rval] >= '0' && istr[rval] <= '9')
               rval++;
            }
         }
   return( rval);
}

int parse_ra_dec( const char *istr, double *ra, double *dec)
{
   int comma_loc, bytes_scanned;
   double hr, min, sec = 0., is_negative = 0, deg;

   for( comma_loc = 0; istr[comma_loc] && istr[comma_loc] != ','; comma_loc++)
      ;
   if( !istr[comma_loc])
      return( -1);            /* gotta have a comma... */
   if( strchr( istr + comma_loc + 1, ','))
      return( -2);            /* but gotta have only _one_ comma */

                  /* OK,  let's go after the RA first... */
                  /* try for the six-digit or six-digit-plus-extras first: */
   bytes_scanned = grab_base_sixty( istr, &hr);
   if( bytes_scanned >= 6 && istr[bytes_scanned] == ',' &&
                                           hr >= 0. && hr < 24. * 3600.)
      {
      *ra = hr * 15. / 3600.;    /* cvt RA seconds to decimal degrees */
      istr += bytes_scanned;
      }
   else if( sscanf( istr, "%lf%n", &hr, &bytes_scanned) == 1)
      {
      istr += bytes_scanned;
      if( *istr == 'h')
         {
         *ra = hr * 15.;
         istr++;
         if( sscanf( istr, "%lf%n", &min, &bytes_scanned) == 1 &&
                  istr[bytes_scanned] == 'm')
            {
            *ra += min / 4.;
            istr += bytes_scanned + 1;
            if( sscanf( istr, "%lf%n", &sec, &bytes_scanned) == 1 &&
                     istr[bytes_scanned] == 's')
               {
               *ra += sec / 240.;
               istr += bytes_scanned + 1;
               }
            }
         }
      else              /* assume decimal degrees */
         *ra = hr;
      }
   else
      return( -3);

   if( *istr != ',')
      return( -4);

   istr++;
   if( *istr == '+')
      istr++;
   else if( *istr == '-')
      {
      istr++;
      is_negative = 1;
      }

                  /* try for the six-digit or six-digit-plus-extras first: */
   bytes_scanned = grab_base_sixty( istr, &deg);
   if( bytes_scanned >= 6 && !istr[bytes_scanned] &&
                              deg >= 0. && deg < 90. * 3600.)
      {
      *dec = deg / 3600.;    /* cvt arcseconds to decimal degrees */
      istr += bytes_scanned;
      }
   else if( sscanf( istr, "%lf%n", &deg, &bytes_scanned) == 1)
      {
      *dec = deg;
      istr += bytes_scanned;
      if( sscanf( istr, "%lf%n", &min, &bytes_scanned) == 1 &&
               istr[bytes_scanned] == '\'')
         {
         *dec += min / 60.;
         istr += bytes_scanned + 1;
         if( sscanf( istr, "%lf%n", &sec, &bytes_scanned) == 1 &&
                  istr[bytes_scanned] == '"')
            {
            *dec += sec / 3600.;
            istr += bytes_scanned + 1;
            }
         }
      }
   else
      return( -5);

   if( *istr)
      return( -6);

   if( is_negative)
      *dec = -*dec;
   return( 0);
}

void init_config( CHARON_CONFIG *cconfig)
{
   memset( cconfig, 0, sizeof( CHARON_CONFIG));
   cconfig->cutoff_point = .99;
   cconfig->targeting_on = 1;
   cconfig->max_stars = 500;
   cconfig->n_stars_used = 7,
   strcpy( cconfig->report_file_name, "mpc.log");
            /* It's not part of CHARON_CONFIG,  but gsc_mag_limit */
            /* _is_ stored in the charon.dat file: */
   gsc_mag_limit = 1800;
}


int parse_charon_command( CHARON_CONFIG *cconfig, IMAGE *image,
                                         const char *argv)
{
   int rval = 0;

   switch( argv[1])
      {
      case 'b':
         image->time = extract_jd_from_string( argv + 2);
         if( image->time < 0.)
            {
            log_printf( log_file, strings[57], argv + 2);
            log_printf( log_file, strings[58]);
            exit_2( -24);
            }
         image->time += (image->exposure_length * .5 / 86400.) *
                  (double)( 1 - cconfig->file_time_point);
         image->time += cconfig->time_offset / 86400.;
         break;
      case 'c':
         cconfig->cutoff_point = atof( argv + 2);
         break;
      case 'C':            /* already handled */
         break;
      case 'd':
         {
         extern int is_drift_scan;

         is_drift_scan = 1;
         log_printf( log_file, "Drift-scanned image\n");
         }
         break;
      case 'D':
         cconfig->millisec_per_video_reset = atoi( argv + 2);
         log_printf( log_file, strings[92],
             cconfig->millisec_per_video_reset);
         /* "Will add %d milliseconds to the clock for each video reset\n" */
         break;
      case 'e':
         sscanf( argv + 2, "%lf,%lf", &cconfig->latlon[0],
                                         &cconfig->latlon[1]);
         log_printf( log_file, strings[32], cconfig->latlon[1],
                                            cconfig->latlon[0]);
         cconfig->latlon[0] *= PI / 180.;
         cconfig->latlon[1] *= PI / 180.;
         break;
      case 'f':     /* .PIC image has no focal len */
         image->focal_len = atof( argv + 2);
         break;
//    case 'g':
//       cconfig->observer_code = atoi( argv + 2);
//       break;
      case 'h':
         cconfig->fit_order = atoi( argv + 2);
         log_printf( log_file, strings[31], cconfig->fit_order);
         break;
      case 'i':
         image->is_inverted ^= 1;
         break;
      case 'j':
         {
         double new_offset = atof( argv + 2);

         image->time += (new_offset - cconfig->time_offset) / 86400.;
         log_printf( log_file, strings[30], cconfig->time_offset, new_offset);
         cconfig->time_offset = new_offset;
         }
         break;
      case 'k':
         cconfig->saturation_point = atoi( argv + 2);
         break;
      case 'l':
         gsc_bright_mag_limit = (int)( 100. * atof( argv + 2) + .5);
         log_printf( log_file, strings[93],
                  gsc_bright_mag_limit / 100, gsc_bright_mag_limit % 100);
                          /* "GSC/Ax.0 BRIGHT mag limit reset to %d.%02d\n" */
         break;
      case 'N':
         sscanf( argv + 2, "%d,%d,%d,%d", &ignored_north, &ignored_south,
                     &ignored_east, &ignored_west);
         break;
      case 'n':
         sscanf( argv + 2, "%d,%d", &cconfig->max_stars,
                                       &cconfig->n_stars_used);
         break;
//    case 'P':
//       psf_fitting_on = 0;
//       log_printf( log_file,
//                   "PSF fitting shut off;  centroiding only\n");
//       break;
      case 'r':
         suppress_display_of_residuals = 1;
         log_printf( log_file,
                     "Suppressing display of residuals\n");
         break;
      case 's':
         image->focal_len *= atof( argv + 2);
         break;
      case 't':
         debug_level = atoi( argv + 2);
         break;
      case 'w':
         image_pixel_separation = atoi( argv + 2);
         break;
      case 'y':
         sscanf( argv + 2, "%lf,%lf", &image->pixel_xsize,
                                         &image->pixel_ysize);
         break;
      case 'z':
         image->image_offset = atoi( argv + 2);
         log_printf( log_file, strings[26], image->image_offset);
         break;
#ifdef __WATCOMC__
      case ';':
         {
         extern short _SVGAType;

         _SVGAType = 1;
         }
         break;
      case '!':
         {
         extern int use_pure_vesa;

         use_pure_vesa = 0;
         }
         break;
#endif
      case '=':      /* special flag to start up settings menu */
         break;      /* at program start */
      default:
         rval = -1;
      }
   return( rval);
}

int load_configuration( CHARON_CONFIG *cconfig, IMAGE *image,
                                         const char *config_file_name)
{
   FILE *ifile = fopen( config_file_name, "rb");
   char tbuff[160], dummy[20];
   int i;

   if( !ifile)
      ifile = fopen( "defaults.dat", "rb");
   if( !ifile)
      return( -1);

   if( debug_level)
      log_printf( log_file, "Startup file opened\n");

   if( fgets( tbuff, sizeof( tbuff), ifile))
      sscanf( tbuff, "%lf %lf %s %lf %lf %d", &cconfig->max_tilt,
                &cconfig->cutoff_point, cconfig->cd_drive_path,
                &cconfig->latlon[0], &cconfig->latlon[1],
                &cconfig->fit_order);

   if( debug_level)
      log_printf( log_file, "First line read;");

   if( fgets( tbuff, sizeof( tbuff), ifile))
      sscanf( tbuff, "%d %d %lf %lf %d %lf %d %lf", &cconfig->max_stars,
                       &cconfig->n_stars_used, &cconfig->search_dist,
                       &cconfig->scale_tolerance, &cconfig->video_mode,
                       &cconfig->time_offset, &cconfig->saturation_point,
                       &cconfig->tilt_angle);
   if( debug_level)
      log_printf( log_file, "second line read;");

   if( fgets( tbuff, sizeof( tbuff), ifile))
      sscanf( tbuff, "%d %d %d %lf %s %d %d %s %d", &cconfig->ra_dec_format,
                     &cconfig->color_table_number, &cconfig->targeting_on,
                     &cconfig->max_residual, dummy,
                     &image->image_offset, &cconfig->cell_size,
                     cconfig->target_name, &psf_fitting_on);

   if( fgets( tbuff, sizeof( tbuff), ifile))
      sscanf( tbuff, "%s %s %s %d %d %d %d", cconfig->report_file_name,
                cconfig->obs_code, cconfig->path_name,
                &cconfig->catalog_used, &cconfig->file_time_point,
                &cconfig->display_flipped, &cconfig->pixel_origin_at_top_left);

               /* 6 Feb 2004:  there is no longer an "Ax.0 (data on hard */
               /* drive)" setting.  These,  corresponding to 2 and 5,    */
               /* are now remapped to 1 and 4:  plain old "Ax.0".        */
   if( cconfig->catalog_used == 2 || cconfig->catalog_used == 5)
      cconfig->catalog_used--;
   *dummy = ' ';
   if( fgets( tbuff, sizeof( tbuff), ifile))
      sscanf( tbuff, "%d %lf %lf %lf %d %d %d %d%c", &image->is_inverted,
               &image->focal_len,
               &image->pixel_xsize, &image->pixel_ysize,
               &cconfig->millisec_per_video_reset, &gsc_mag_limit,
               &cconfig->ccd_id, &gsc_bright_mag_limit, dummy);

   cconfig->photometric_band = *dummy;
   if( cconfig->photometric_band < ' ')
      cconfig->photometric_band = ' ';

   if( fgets( tbuff, sizeof( tbuff), ifile))
      sscanf( tbuff, "%d %d %d %d", &ignored_north, &ignored_south,
                     &ignored_east, &ignored_west);

   for( i = 0; cconfig->target_name[i] >= ' '; i++)
      if( cconfig->target_name[i] == '_')
         cconfig->target_name[i] = ' ';
   cconfig->target_name[i] = '\0';

   if( debug_level)
      log_printf( log_file, "Full startup read\n");

   fclose( ifile);
   return( 0);
}

int reset_language( const int new_language)
{
   if( new_language != 'e')            /* make sure new language exists */
      {
      FILE *ifile;

      string_file_name[7] = (char)new_language;
      ifile = fopen( string_file_name, "rb");
      if( ifile)
         fclose( ifile);
      else
         string_file_name[7] = 's';
      }
   language = new_language;
   return( language);
}

int correct_cconfig_with_startup( CHARON_CONFIG *cconfig)
{
   FILE *startup_file = fopen( "startup.mar", "rb");
   char tbuff[150];

   if( !startup_file)
      return( -1);
   if( debug_level)
      log_printf( log_file, "Loading startup file...");
   while( fgets( tbuff, sizeof( tbuff), startup_file))
      switch( atoi( tbuff))
         {
         case 11:
            sscanf( tbuff + 12, "%lf %lf", &cconfig->latlon[0],
                                           &cconfig->latlon[1]);
            cconfig->latlon[0] *= PI / 180.;
            cconfig->latlon[1] *= PI / 180.;
            break;
         case 18:
            {
            int i;

            strcpy( cconfig->cd_drive_path, tbuff + 12);
            for( i = 0; cconfig->cd_drive_path[i] > ' '; i++)
               ;
            cconfig->cd_drive_path[i] = '\0';
            if( i == 1)
               strcat( cconfig->cd_drive_path, ":\\");
            else if( cconfig->cd_drive_path[i - 1] != '\\')
               strcat( cconfig->cd_drive_path, "\\");
            break;
            }
         case 28:
            cconfig->altitude = atof( tbuff + 12);
            break;
         case 51:
            if( language == 'E')       /* i.e.,  not overridden on cmd line */
               reset_language( tbuff[12]);
            break;
         case 63:
            cconfig->a10_drive_letter = tbuff[12];
            break;
         default:
            break;
         }
   fclose( startup_file);
   if( debug_level)
      log_printf( log_file, "done\n");
   return( 0);
}

   /* Some years ago,  I made a Bad Decision.   I downloaded the MPC's  */
   /* list of observer codes,                                           */
   /* http://cfa-www.harvard.edu/iau/lists/ObsCodes.html                */
   /* and renamed it to 'stations.txt'.  I should have left the name    */
   /* unchanged (except to allow for the possibility that the extension */
   /* would be truncated to '.htm').  So now,  the code looks for       */
   /* 'ObsCodes.html',  then 'ObsCodes.htm',  then 'stations.txt'.      */

int get_observer_data( const char *obs_code, char *buff)
{
   FILE *ifile = fopen( "ObsCodes.html", "rb");

   if( !ifile)
      ifile = fopen( "ObsCodes.htm", "rb");
   if( !ifile)
      ifile = fopen( "stations.txt", "rb");
   if( !ifile)
      return( -1);
   while( fgets( buff, 100, ifile))
      if( !memcmp( buff, obs_code, 3))
         {
         fclose( ifile);
         return( 0);
         }
   fclose( ifile);
   return( -2);
}

#define EARTH_MAJOR_AXIS 6378140.
#define EARTH_MINOR_AXIS 6356755.
#define EARTH_AXIS_RATIO (EARTH_MINOR_AXIS / EARTH_MAJOR_AXIS)

int compute_latlon_from_mpc_data( const char *buff, double *latlon)
{
   if( buff[7] == '.' && buff[14] == '.' && buff[23] == '.')
      {                    /* "Standard" MPC parallax constants */
      double lon, x, y;

      log_printf( log_file, strings[34], buff + 30);
      sscanf( buff + 4, "%lf%lf%lf", &lon, &x, &y);
      y /= EARTH_AXIS_RATIO * EARTH_AXIS_RATIO;
      if( lon > 180.)
         lon -= 360.;
      latlon[0] = lon * PI / 180.;
      latlon[1] = atan( y / x);
      }
   else           /* longitude/latitude, maybe alt, given */
      {
      sscanf( buff + 4, "%lf%lf", latlon, latlon + 1);
      latlon[0] *= PI / 180.;
      latlon[1] *= PI / 180.;
      log_printf( log_file, "Unofficial MPC code found\n");
      }
   log_printf( log_file, strings[32], latlon[1] * 180. / PI,
                                      latlon[0] * 180. / PI);
   return( 0);
}

int check_for_guide_cd( const char *cd_drive_path)
{
   static int already_checked = 0, rval;

   if( !already_checked)
      {
      FILE *ifile = local_then_cd_fopen( cd_drive_path, "asteroid\\astnum");

      already_checked = 1;
      if( ifile)
         {
         char tbuff[180];

         fgets( tbuff, sizeof( tbuff), ifile);
         fclose( ifile);
         rval = atoi( tbuff);
         }
      else
         {
         ifile = local_then_cd_fopen( cd_drive_path, "vsop\\vsop.bin");
         if( ifile)
            fclose( ifile);
         rval = (ifile ? 0 : -1);
         }
      }
   return( rval);
}

double extract_jd_from_string( const char *buff)
{
   long year;
   int hour, minute, month, day, got_it_right = 0;
   double second, jd;
   extern char *month_names[12];    /* from DAT.CPP */

   if( sscanf( buff, "%d/%d/%ld %d:%d:%lf",
         &day, &month, &year,  &hour, &minute, &second) == 6)
      if( hour >= 0 && hour < 24 && minute >= 0 && minute < 60 &&
          second >= 0 && second < 60. && day > 0 && day < 32 &&
          month > 0 && month < 13)
         got_it_right = 1;
   if( !got_it_right)
      return( -1.);
   if( year > 0L && year < 50L)
      year += 2000;
   else if( year >= 50L && year < 100L)
      year += 1900;
   jd = (double)dmy_to_day( 0, month, (long)year, 0) +
                  (double)hour / 24. + (double)minute / 1440. +
                  second / 86400. - .5;
   jd += (double)day;
   return( jd);
}

char *write_report_data( const CHARON_CONFIG *cconfig,
         const double jd, const double *target_ra_dec,
         const double target_mag, int report_type,
         const char *ref_net_name)
{
   FILE *rep_file = fopen( cconfig->report_file_name, "a+");
   char *message = NULL;

   if( rep_file)
      {
      long year;
      int month, day, i;
      double d;
      extern char obj_id_for_report[];
      extern int gsc_photo_band;                           /* grab_gsc.cpp */
      char tbuff[100], net_used[20];

      *net_used = '\0';    /* Assume no NET command is given in the file */
      fseek( rep_file, 0L, SEEK_END);
      if( !ftell( rep_file))   /* brand new report file; copy in HEADER.MPC */
         {
         FILE *ifile = fopen( "header.mpc", "r");
         int keep_going = 1;
         time_t t0 = time( NULL);

         while( keep_going && fgets( tbuff, sizeof( tbuff), ifile))
            if( !memcmp( tbuff, "No MPC header", 13))
               keep_going = 0;
            else if( *tbuff != ';' && *tbuff > ' ') /* skip comments, blanks */
               {           /* Charon is supposed to 'fix up' the MPC code: */
               if( !memcmp( tbuff, "COD ", 4))
                  sprintf( tbuff + 4, "%s\n", cconfig->obs_code);
               fwrite( tbuff, strlen( tbuff), 1, rep_file);
               }
         fprintf( rep_file, "COM Charon version %s %s\n", __DATE__, __TIME__);
         fprintf( rep_file, "COM Report file created %s", ctime( &t0));
         fclose( ifile);
         }
      else        /* it's an existing report file;  what NET is being used? */
         {
         fseek( rep_file, 0L, SEEK_SET);
         while( fgets( tbuff, sizeof( tbuff), rep_file))
            if( !memcmp( tbuff, "NET ", 4))
               strcpy( net_used, tbuff + 4);
         fseek( rep_file, 0L, SEEK_END);
         }

             /* If there is no 'NET' statement,  or if we've changed nets, */
             /* the following will detect that and write out a 'NET' line. */
      if( memcmp( ref_net_name, net_used, strlen( ref_net_name)))
         fprintf( rep_file, "NET %s\n", ref_net_name);

             /* Now we can actually worry about writing an MPC/IOTA report: */
      memset( tbuff, ' ', 80);
      d = jd + .5;
      day_to_dmy( (long)d, &day, &month, &year, 0);
      d -= dmy_to_day( 0, month, year, 0);
      sprintf( tbuff + 15, "%4ld %02d %8.5lf ", year, month, d);
      if( report_type == 'o')
         {
                  /* '8' = 'two decimal places' */
         sprintf_ra_dec( tbuff + 32, target_ra_dec[0], 2 | 8, 0);
                  /* '16' = 'single decimal place' */
         sprintf_ra_dec( tbuff + 44, target_ra_dec[1], 2 | 16, 1);
         message = strings[41];
         tbuff[43] = tbuff[55] = ' ';
         }
      if( report_type == 'i')
         {
         sprintf( tbuff + 23, "%9.6lf", d);
                  /* '0' = 'three decimal places' */
         sprintf_ra_dec( tbuff + 32, target_ra_dec[0], 2, 0);
                  /* '8' = 'two decimal place' */
         sprintf_ra_dec( tbuff + 44, target_ra_dec[1], 2 | 8, 1);
         message = strings[40];
         }
      if( report_type == 'p')
         message = strings[98];   /* "Magnitude report written"; */
      for( i = 0; i < 80; i++)
         if( !tbuff[i])
            tbuff[i] = ' ';
      tbuff[34] = tbuff[37] = tbuff[47] = tbuff[50] = ' ';
      memset( tbuff + 56, ' ', 24);
      if( strlen( cconfig->obs_code) == 3)
         {
         tbuff[77] = cconfig->obs_code[0];
         tbuff[78] = cconfig->obs_code[1];
         tbuff[79] = cconfig->obs_code[2];
         }
      tbuff[80] = '\0';
      if( *obj_id_for_report)
         memcpy( tbuff, obj_id_for_report, 12);
      if( tbuff[23] == ' ')    /* leading 0 in day of month */
         tbuff[23] = '0';
      tbuff[14] = 'C';       /* indicate CCD image */
      if( target_mag > 2. && target_mag < 25.)
         if( report_type != 'p')
            {
            sprintf( tbuff + 65, "%4.1lf", target_mag + .05);
            tbuff[69] = ' ';
            }
         else
            {
            sprintf( tbuff + 65, "%4.2lf", target_mag + .005);
            tbuff[70] = ' ';
            }
      if( report_type == 'p')     /* 14 Sep 98:  Added JD */
         {
         sprintf( tbuff + 37, "%14.6lf", jd);
         tbuff[51] = ' ';
         }
      if( cconfig->photometric_band > ' ' && cconfig->photometric_band < 'Z')
         tbuff[70] = cconfig->photometric_band;
               /* If gsc_photo_band is non-zero,  it means that the data */
               /* came from GSC,  and the user doesn't really get to choose */
               /* the photo band... it's "fixed" by what GSC fed him. */
      if( gsc_photo_band)
         tbuff[70] = (char)gsc_photo_band;
      if( !strcmp( cconfig->obs_code, "247"))    /* roving observer: */
         tbuff[14] = 'V';                        /* only change in 1st line */
      fprintf( rep_file, "%s\n", tbuff);
      if( !strcmp( cconfig->obs_code, "247"))    /* roving observer: */
         {                                       /* make a 2nd line  */
         const double latitude = cconfig->latlon[1] * 180. / PI, east_lon =
                   fmod( 3600. + cconfig->latlon[0] * 180. / PI, 360.);

         tbuff[14] = 'v';
         tbuff[32] = '1';
         tbuff[33] = ' ';
         sprintf( tbuff + 34, "%10.6lf %c%9.6lf %5u                247",
                             east_lon, (latitude > 0. ? '+' : '-'),
                             fabs( latitude), (unsigned)cconfig->altitude);
         fprintf( rep_file, "%s\n", tbuff);
         }
      fclose( rep_file);
      }
   return( message);
}

void flip_image( IMAGE *image)
{
   int i, j, k;

   for( j = 0, i = image->ysize - 1; i > j; i--, j++)
      {
      PIXEL *iptr = image->img + i * image->xsize, tval;
      PIXEL *jptr = image->img + j * image->xsize;

      for( k = image->xsize; k; k--)
         {
         tval = *iptr;
         *iptr++ = *jptr;
         *jptr++ = tval;
         }
      }
}

void offset_image( IMAGE *image, const int signed_ints)
{
   if( image->image_offset || signed_ints)
      {
      const int total_pixels = image->xsize * image->ysize;
      int i;

      if( debug_level)
         log_printf( log_file, "Offsetting image...\n");

      for( i = 0; i < total_pixels; i++)
         {
         int new_val = image->img[i] + image->image_offset;

         if( signed_ints && image->img[i])
            {
#ifdef PIXEL_32
            unsigned long tval = (unsigned long)( image->img[i] + 0x80000000);
#else
            unsigned short tval = (unsigned short)( image->img[i] + 0x8000);
#endif

            new_val = (int)tval + image->image_offset;
            }

         if( new_val < image->image_offset || new_val < 0)
            new_val = (PIXEL)0;
         image->img[i] = new_val;
         }
      if( debug_level)
         log_printf( log_file, "Image offset\n");
      }
}

int grab_strings( void)
{
   FILE *ifile = fopen( string_file_name, "rb");
   char buff[100];
   int i, j;

   *buff = '\0';
   for( i = 0; memcmp( buff, "[END]", 5); i++)
      {
      fgets( buff, sizeof( buff), ifile);
      for( j = 0; (unsigned char)buff[j] >= ' '; j++)
         if( buff[j] == '\\' && buff[j + 1] == 'n')
            {
            buff[j] = '\n';
            memmove( buff + j + 1, buff + j + 2, strlen( buff + j));
            }
      buff[j] = '\0';
      strings[i] = (char *)malloc( j + 2);
      strcpy( strings[i], buff);
      }
   fclose( ifile);
   return( 0);
}

int get_palette_name( const int palette_no, char *name)
{
   FILE *ifile = fopen( "palette.dat", "rb");
   int i, j;
   extern int n_color_tables;
   char languages[80];

   fgets( languages, 80, ifile);
   n_color_tables = atoi( languages);
   i = palette_no * (strlen( languages) + 9) + 1;
   for( j = 2; languages[j]; j++)
      if( language == languages[j])
         i += j - 1;
   while( i--)
      fgets( name, 80, ifile);
   fclose( ifile);

   for( i = 0; name[i] != 10 && name[i] != 13; i++)
      ;
   name[i] = '\0';
   return( 0);
}

#define N_KEYWORDS 9

int check_mpc_report_header( const char *header_filename)
{
   FILE *ifile = fopen( header_filename, "rb");
   char buff[100], err_msg[100];
   int line_no = 1, rval = 0, n_lines_found[N_KEYWORDS];

   memset( n_lines_found, 0, N_KEYWORDS * sizeof( int));
   *err_msg = '\0';
   if( !ifile)
      {
      sprintf( err_msg, "MPC report header '%s' was not found!\n",
                                        header_filename);
      rval = -1;
      }
   while( !rval && fgets( buff, sizeof( buff), ifile))
      {
      if( !memcmp( buff, "No MPC header", 13))
         {
         fclose( ifile);
         return( 0);
         }
      if( strlen( buff) > 80)
         {
         sprintf( err_msg,
                "Error at line %d of '%s': line is more than 80 characters\n",
                line_no, header_filename);
         rval = -2;
         }
      else if( *buff != ';' && *buff >= ' ')    /* ignore blanks,  comments */
         {
         static char *keywords[N_KEYWORDS] = {
               "COD ", "CON ", "OBS ", "MEA ", "TEL ", "NET ", "COM ",
               "ACK ", "AC2 " };
         int key_val = -1, i;

         for( i = 0; i < N_KEYWORDS; i++)
            if( !memcmp( buff, keywords[i], 4))
               key_val = i;
         if( key_val == -1)
            {
            sprintf( err_msg,
               "Error at line %d of '%s': line doesn't begin with a keyword\n",
               line_no, header_filename);
            rval = -3;
            }
         else
            {
            n_lines_found[key_val]++;
            }
         }
      line_no++;
      }

   fclose( ifile);

   if( !rval)
      {
      if( n_lines_found[7] != 1)
         {
         sprintf( err_msg,
                   "There should be exactly one ACK line in the report file!");
         rval = -5;
         }
      if( n_lines_found[8] != 1)
         {
         sprintf( err_msg,
                   "There should be exactly one AC2 line in the report file!");
         rval = -6;
         }
      }
   if( *err_msg)
      printf( "%s", err_msg);
   return( rval);
}
