#include <graph.h>
#include <math.h>
#include <conio.h>
#include <stdio.h>
#include <io.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <dos.h>
#include "watdefs.h"
#include "afuncs.h"
#include "findstar.h"
#include "alt_defs.h"
#include "bmouse.h"
#include "mouse.h"
#include "charon.h"
#include "miscell.h"

#ifdef DEBUG_MEM
#include "checkmem.h"

void memory_error( struct _memory_error *err);
#endif

#define MOUSE_SIZE 20

#define PI 3.14159265358979323846264338327
#define TEXT_COLUMNS 27
#define RESIDUAL_THRESH 30
#define RED_INDEX 15
#define YELLOW_INDEX 14
#define GREEN_INDEX 13
#define BLUE_INDEX 12

void find_isophotes( PIXEL *pixels, const int xsize, const int ysize,
            const PIXEL isophote_level, FILE *ofile);       /* isophote.c */
IMAGE *get_blinkable_images( const IMAGE *ref_image, int *n_found,
           const double xfactor, const double yfactor);     /* blink.h */
void show_help_data( void);
void show_image( const PIXEL *buff, const int width, const int height,
      const PIXEL low, const PIXEL high, const double *image_loc,
      const int *disp_loc, const int numcolors);         /* dispimg.cpp */
void mouse_update_buttons( BMOUSE *b, int button);       /* bmouse.cpp */
int add_wcs_header_data( const IMAGE *img);              /* blink.cpp */
int find_obj( double jd, char *obj_id, const char *cd_drive_path, double *ra_dec,
                              double *latlon);
long setup_menu( CHARON_CONFIG *cconfig, IMAGE *img);
int dark_subtract_and_flat_field( IMAGE *image, const IMAGE *dark,
                           const IMAGE *flat);           /* blink.cpp */
int load_kernel_file( const char *filename, double *kernel);
int apply_convolution( IMAGE *image, const int kernel_size, const double *kernel);
int apply_polygon_mask( IMAGE *image, const char *mask_filename);
int reset_language( const int new_language);          /* miscell.cpp */
void mouse_update_buttons( BMOUSE *b, int button);
void exit_2( const int exit_code);                    /* charon.cpp */

FILE *log_file;
FILE *overlay_file;
int use_input_mags, max_background = 0, is_drift_scan = 0;
int debug_level = 0, saturation_point = 65535;
static int started_from_guide = 0;
double get_pixel_dimension( double *params, int true_if_height);
static struct videoconfig vconfig;
extern int using_a20;
static int astnum_value;
extern char *strings[];


static int found_compare( const void *a, const void *b)
{
   if( ((FOUND_STAR *)a)->bright > ((FOUND_STAR *)b)->bright)
      return( 1);
   if( ((FOUND_STAR *)a)->bright < ((FOUND_STAR *)b)->bright)
      return( -1);
   return( 0);
}

static void show_cursor( int x, int y, int hiding)
{
   int dy, dx, count = 0, swap;
   static char restored_pixels[(MOUSE_SIZE * 2 + 1) * 6];

   _setcolor( YELLOW_INDEX);
   for( dy = -1; dy < 2; dy++)
      for( dx = -MOUSE_SIZE; dx <= MOUSE_SIZE; dx++)
         if( abs( dx) > 1)
            for( swap = 0; swap < 2; swap++, count++)
               {
               short x1 = (short)( x + (swap ? dx : dy));
               short y1 = (short)( y + (swap ? dy : dx));
               if( hiding)
                  _setcolor( restored_pixels[count]);
               else
                  {
                  restored_pixels[count] = _getpixel( x1, y1);
                  _setcolor( dy ? 0 : YELLOW_INDEX);
                  }
               _setpixel( x1, y1);
               }
}

static void show_xorred_line( int x1, int y1, int x2, int y2)
{
   _setwritemode( _GXOR);
   _setcolor( (vconfig.numcolors == 16) ? RED_INDEX : 128);
   _moveto( (short)x1, (short)y1);
   _lineto( (short)x2, (short)y2);
   _setwritemode( _GPSET);
}

static double show_zoom_box( int x1, int y1, int x2, int y2, int *zoom_box,
                  const int zooming_out)
{
   int xsize = abs( x2 - x1), i;
   int ysize = abs( y2 - y1);
   int xscr = vconfig.numxpixels, yscr = vconfig.numypixels;
   double rval;

   if( xsize < 10 && ysize < 10)
      {
      zoom_box[0] = zoom_box[2] = x2;
      zoom_box[1] = zoom_box[3] = y2;
      return( zooming_out ? .5 : 1.);
      }
   if( xsize * yscr > ysize * xscr)    /* wider than is high */
      ysize = xsize * yscr / xscr;
   else                                /* tall,  narrow box */
      xsize = ysize * xscr / yscr;
   zoom_box[0] = x1 - xsize;
   zoom_box[1] = y1 - ysize;
   zoom_box[2] = x1 + xsize;
   zoom_box[3] = y1 + ysize;
   for( i = 0; i < 4; i++)
      show_xorred_line( zoom_box[(i & 1) * 2],
                        zoom_box[(i & 2) + 1],
                        zoom_box[i > 1 ? 0 : 2],
                        zoom_box[(i & 1) * 2 + 1]);
   rval = .5 * (double)xscr / (double)xsize;
   if( zooming_out)
      {
      for( i = 0; i < 4; i += 2)
         show_xorred_line( zoom_box[0], zoom_box[i + 1],
                           zoom_box[2], zoom_box[3 - i]);
      rval = 1. / rval;
      }
   return( rval);
}

int directly_read_mouse( BMOUSE *b)
{
   int button = (global_mouse_button >> 1) | (global_mouse_button << 1);

   button = ((button & 3) | (global_mouse_button & 4));
   if( global_mouse_x > (int)b->xmax) global_mouse_x = b->xmin;
   if( global_mouse_y > (int)b->ymax) global_mouse_y = b->ymin;
   if( global_mouse_x < (int)b->xmin) global_mouse_x = b->xmax;
   if( global_mouse_y < (int)b->ymin) global_mouse_y = b->ymax;
   b->dx = (int)( global_mouse_x - b->x);
   b->dy = (int)( global_mouse_y - b->y);
   b->prev_x = b->x;
   b->prev_y = b->y;
   b->x = global_mouse_x;
   b->y = global_mouse_y;
   mouse_update_buttons( b, button);
   return( 0);
}

void exit_2( const int exit_code)
{
   if( started_from_guide)
      {
      log_printf( log_file,  "Exit code %d\n", exit_code);
      log_printf( log_file, strings[88]);
      getch( );
      }
   exit( exit_code);
}

static void place_x_on_screen( char *image_star_buff,
                const short x, const short y, const short color_index)
{
   int i;

   _setcolor( color_index);
   ((short *)image_star_buff)[0] = x;
   ((short *)image_star_buff)[1] = y;
   image_star_buff += 4;
   for( i = 0; i < 40; i++)
      {
      const int delta = 5 + i / 4;
      const int x1 = x + ((i & 1) ? delta : -delta);
      const int y1 = y + ((i & 2) ? delta : -delta);

      *image_star_buff++ = (char)_getpixel( x1, y1);
      _setpixel( x1, y1);
      }
}

static void remove_x_from_screen( const char *image_star_buff)
{
   int i;
   const short x = ((const short *)image_star_buff)[0];
   const short y = ((const short *)image_star_buff)[1];

   image_star_buff += 4;
   for( i = 0; i < 40; i++)
      {
      const int delta = 5 + i / 4;
      const int x1 = x + ((i & 1) ? delta : -delta);
      const int y1 = y + ((i & 2) ? delta : -delta);

      _setcolor( (short)( (unsigned char)*image_star_buff++));
      _setpixel( x1, y1);
      }
}

int n_color_tables;

static void set_palette( const int color_table_number)
{
   int i, len;
   FILE *ifile = fopen( "palette.dat", "rb");
   int red[12], green[12], blue[12];
   char buff[80];
   long *palette = (long *)calloc( vconfig.numcolors, sizeof( long));

   fgets( buff, 80, ifile);
   n_color_tables = atoi( buff);
   len = strlen( buff);
   for( i = color_table_number * (len + 9) + len - 3; i; i--)
      fgets( buff, 80, ifile);
   for( i = 0; i < 11; i++)
      {
      fgets( buff, 80, ifile);
      sscanf( buff, "%d %d %d", blue + i, green + i, red + i);
      }
   fclose( ifile);

   palette[RED_INDEX] = 0x00003f;
   palette[YELLOW_INDEX] = 0x003f3f;
   palette[GREEN_INDEX] = 0x003f00;
   palette[BLUE_INDEX] = 0x3f0000;
   if( vconfig.numcolors == 16)
      for( i = 0; i < 11; i++)
         {
         long rgb = (((long)blue[i] << 16) | ((long)green[i] << 8) | (long)red[i]);

         palette[i] = rgb;
         }
   else
      {
      for( i = 0; i < 64; i++)
         {
         int idx = i * 10 / 63;
         int fraction = (i * 10 * 256 / 63) - idx * 256;
         int r =   red[idx] + (  red[idx + 1] -   red[idx]) * fraction / 256;
         int g = green[idx] + (green[idx + 1] - green[idx]) * fraction / 256;
         int b =  blue[idx] + ( blue[idx + 1] -  blue[idx]) * fraction / 256;

         palette[i + 16] = (((long)b << 16) | ((long)g << 8) | (long)r);
         }
      for( i = 128; i < 256; i++)
         palette[i] = 0x00003f;
      }
   _remapallpalette( palette);
   free( palette);
}


static int short_dist( double dx, double dy)
{
   dx = fabs( dx);
   dy = fabs( dy);
   if( dx > dy)
      dx += dy / 3.;
   else
      dx = dx / 3. + dy;
   return( (int)dx);
}

static void fix_arguments( char **argv, int argc)
{
   int i, loc = 0;
   static char revised_args[256];

   while( argc)
      {
      for( i = 1; i < argc && argv[i][0] != '-'; i++)
         {
         char *new_arg = revised_args + loc;

         if( i == 1)
            strcpy( new_arg, argv[0]);
         strcat( new_arg, " ");
         strcat( new_arg, argv[i]);
         argv[0] = new_arg;
         argv[i][0] = '\0';
         }
      if( i > 1)
         loc += strlen( revised_args + loc) + 1;
      argv += i;
      argc -= i;
      }
}

#ifdef DEBUG_MEM
void memory_error( struct _memory_error *err)
{
   static char *status[6] = { "OK", "Re-freed", "Header overwritten",
         "End overwritten", "Null pointer", "Not allocated" };
   char buff[100];

   if( err->error_no == POINTER_NOT_ALLOCATED)
      return;
   sprintf( buff, "%s: Line %u: %u bytes %s", err->file_allocated,
            err->line_allocated, err->allocation_size, status[err->error_no]);
   log_printf( log_file, "\n%s", buff);

   if( err->error_no && err->error_no != POINTER_NOT_ALLOCATED)
      {
      sprintf( buff, "Found at %s: Line %u", err->file_found, err->line_found);
      log_printf( log_file, "\n%s", buff);
      }
}
#endif

static void show_angle_and_fl( IMAGE *image)
{
   double computed_focal_length, angle_between_axes;

   computed_focal_length =
            image->focal_len / get_pixel_dimension( image->xform, 0);
   computed_focal_length /= (180. / PI) * 3600.;
   log_printf( log_file, strings[42], computed_focal_length / 25.4,
                                                computed_focal_length);

   computed_focal_length =
            image->focal_len / get_pixel_dimension( image->xform, 1);
   computed_focal_length /= (180. / PI) * 3600.;
   log_printf( log_file, strings[43], computed_focal_length / 25.4,
                                                computed_focal_length);

   angle_between_axes = (180. / PI) * fabs(
           atan2( image->xform[6 + N_FIND_COEFFS], image->xform[6]) -
           atan2( image->xform[7 + N_FIND_COEFFS], image->xform[7]));
   if( angle_between_axes > 180.)
      angle_between_axes -= 180.;
   log_printf( log_file, strings[44], angle_between_axes);
}

static int advance_clock( int n_milliseconds)
{
   struct _dostime_t t;
   int hour, minute, second, hsecond;

   _dos_gettime( &t);
   hsecond = t.hsecond + n_milliseconds / 10;
   second = t.second;
   minute = t.minute;
   hour = t.hour;

   second += hsecond / 100;
   hsecond %= 100;
   if( hsecond < 0)
      {
      hsecond += 100;
      second--;
      }

   minute += second / 60;
   second %= 60;
   if( second < 0)
      {
      second += 60;
      minute--;
      }

   hour += minute / 60;
   minute %= 60;
   if( minute < 0)
      {
      minute += 60;
      hour--;
      }

   t.hsecond = hsecond;
   t.second = second;
   t.minute = minute;
   t.hour = hour;
   _dos_settime( &t);
   log_printf( log_file, strings[91], n_milliseconds);
                          /* "Clock advanced %d milliseconds\n" */
   return( 0);
}

char *ref_net_name = NULL;

static void create_ax0_desig( char *buff, unsigned long ax0_num)
{
   unsigned long zone, running_num;

   if( ax0_num < 50000000u * 24u * 2u)
      {
      zone = ax0_num / 50000000u;
      running_num = ax0_num % 50000000u;
      }
   else
      {
      zone = (ax0_num - 2400000000u) / 4000000u;
      running_num = ax0_num % 4000000u;
      *buff++ = 'S';
      }
   sprintf( buff, "A%d_%04d_%08d", zone / 24 + 1, (zone % 24) * 75,
                     running_num);
}

void make_star_desig( const REF_STAR *star, char *desig)
{
   if( !memcmp( ref_net_name, "UCAC2", 5))
      sprintf( desig, "2UCAC%7ld", star->ppm_num);
   else if( !memcmp( ref_net_name, "UCAC3", 5))
      sprintf( desig, "3UCAC%03ld-%06ld", star->ppm_num / 1000000L,
                                          star->ppm_num % 1000000L);
   else if( !memcmp( ref_net_name, "UCAC4", 5))
      sprintf( desig, "4UCAC%03ld-%06ld", star->ppm_num / 1000000L,
                                          star->ppm_num % 1000000L);
   else if( *ref_net_name == 'U' && ref_net_name[1] == 'C')     /* UCAC1 */
      sprintf( desig, "UCAC%4d%7ld", star->zone, star->ppm_num);
   else if( *ref_net_name == 'G' || *ref_net_name == 'A'
                                 || *ref_net_name == 'T')    /* GSC */
      sprintf( desig, "GSC %4d%5d", star->zone, star->number);
   else if( *ref_net_name == 'U' && ref_net_name[1] == 'S')    /* USNO */
      create_ax0_desig( desig, star->ppm_num);
   else                 /* should never get here */
      sprintf( desig, "??? %s", ref_net_name);
}

REF_STAR *grab_ucac2_data( const double x, const double y, const double width,
                const double height, const int cd_drive_letter, int *n_found,
                const int photo_band, const double jd);     /* grab_gsc.c */
REF_STAR *grab_ucac3_data( const double x, const double y, const double width,
                const double height, const int cd_drive_letter, int *n_found,
                const int photo_band, const double jd);     /* grab_gsc.c */
REF_STAR *grab_ucac4_data( const double x, const double y, const double width,
                const double height, const int cd_drive_letter, int *n_found,
                const int photo_band, const double jd);     /* grab_gsc.c */
int extract_ucac2_stars( FILE *ofile, const double ra, const double dec,
                  const double width, const double height, const char *path,
                  const int is_supplement);

static REF_STAR *grab_catalog_stars( const double xc, const double yc,
                 const CHARON_CONFIG *cconfig,
                 const char *gsc_12_filename, const double jd,
                 int *n_stars_found)
{
   REF_STAR *ref_stars;
   int n_stars;
   const double width = cconfig->search_dist / cos( yc * PI / 180.);

   if( gsc_12_filename)
      {
      if( gsc_12_filename[strlen( gsc_12_filename) - 1] == '\\')
         {
         ref_stars = grab_ucac1_data( xc, yc, width, cconfig->search_dist,
                              gsc_12_filename, &n_stars, jd);
         ref_net_name = "UCAC1";
         }
      else
         {                                        /* tycho/hipparcos */
         ref_stars = grab_tycho2_data( xc, yc, width, cconfig->search_dist,
                                 gsc_12_filename, &n_stars, jd,
                                 cconfig->photometric_band);
         if( n_stars >= 0)
            ref_net_name = "TYCHO-2";
         else     /* see if it's GSC 1.2 or Sloan ACR */
            {
            ref_stars = grab_gsc_12_data( xc, yc, width, cconfig->search_dist,
                              gsc_12_filename, &n_stars, jd);
            if( n_stars <= 0)
               {
               const char *ucac2_filename = "ucac2.txt";
               FILE *ofile = fopen( ucac2_filename, "w");

               n_stars = extract_ucac2_stars( ofile, xc, yc,
                            width, cconfig->search_dist, gsc_12_filename, 0);
               fclose( ofile);
               if( n_stars > 0)
                  ref_stars = grab_gsc_12_data( xc, yc,
                            width, cconfig->search_dist, ucac2_filename,
                            &n_stars, jd);
               }
            }
         }
      }
   else switch( cconfig->catalog_used)
      {
      case CATALOG_GSC_11:   /* "standard" GSC 1.1 */
      case CATALOG_GSC_ACT:  /* "standard" GSC-ACT */
         {
         const int dir_used = (cconfig->catalog_used == CATALOG_GSC_ACT)
                                           - (astnum_value >= 53122);

         ref_net_name = ((cconfig->catalog_used == CATALOG_GSC_ACT)
                             ? "GSC-ACT" : "GSC-1.1");
         if( *cconfig->cd_drive_path == '*')
            {
            log_printf( log_file,
                   "Sorry, you need the Guide CD to use the %s catalog.\n",
                   ref_net_name);
            exit_2( -1);
            }
         ref_stars = grab_gsc_data( xc, yc, width, cconfig->search_dist,
                 cconfig->cd_drive_path, &n_stars, dir_used);
         }
         break;
      case CATALOG_A2:        /* look on CD for A2.0 */
      case CATALOG_A1:        /* look on CD for A1.0 */
         ref_stars = grab_a10_data( xc, yc, width, cconfig->search_dist,
           cconfig->a10_drive_letter, &n_stars, cconfig->photometric_band);
         if( ref_stars && n_stars > 0)
            {
            static char *ax0_name = "USNO-xx         ";
            char tbuff[40];                  /* set up ref net name using    */
            int i;                           /* the desig of the first star, */
                                             /* which always starts w/ "A1", */
            ref_net_name = ax0_name;         /* "A2", "SA1", or "SA2"        */
            create_ax0_desig( tbuff, ref_stars->ppm_num);
            for( i = 0; tbuff[i] != '_'; i++)
               ref_net_name[i + 5] = tbuff[i];
            ref_net_name[i + 5] = '\0';
            }
         break;
      case CATALOG_UCAC2:
         ref_stars = grab_ucac2_data( xc, yc, width, cconfig->search_dist,
              cconfig->a10_drive_letter, &n_stars, cconfig->photometric_band, jd);
         ref_net_name = "UCAC2";
         break;
      case CATALOG_TYCHO2:
         ref_net_name = ((astnum_value >= 56691) ? "TYCHO-2" : "ACT");
         if( *cconfig->cd_drive_path == '*')
            {
            log_printf( log_file,
                  "Sorry, you need the Guide CD to get Tycho data.\n");
            exit_2( -2);
            }

         ref_stars = grab_tycho_data( xc, yc, width, cconfig->search_dist,
                                 cconfig->cd_drive_path, &n_stars, jd,
                                 cconfig->photometric_band);
         break;
      case CATALOG_UCAC3:
         ref_stars = grab_ucac3_data( xc, yc, width, cconfig->search_dist,
                              cconfig->a10_drive_letter, &n_stars,
                              cconfig->photometric_band, jd);
         ref_net_name = "UCAC3";
         break;
      case CATALOG_UCAC4:
         ref_stars = grab_ucac4_data( xc, yc, width, cconfig->search_dist,
                              cconfig->a10_drive_letter, &n_stars,
                              cconfig->photometric_band, jd);
         ref_net_name = "UCAC4";
         break;
      default:
         ref_net_name = NULL;       /* just to dodge a compiler warning */
         log_printf( log_file, "??? use_a10 = %d\n", cconfig->catalog_used);
         log_printf( log_file, "Internal error;  please contact Project Pluto");
         exit_2( -3);
         break;
      }

   if( n_stars <= 0)
      {
      int msg;

      switch( n_stars)
         {
         case -9:          /* couldn't open PLATES.BIN */
            msg = 50;
            break;
         case -1:          /* ran out of memory */
         case -5:
            msg = 51;
            break;
         case -2:          /* couldn't open GSCDATA2.IDX */
            msg = 52;
            break;
         case -3:          /* couldn't open GSC data on CD-ROM drive */
            msg = 53;
            break;
         case -4:          /* Error in reading the GSC data */
            msg = 54;
            break;
         case -12:          /* Couldn't identify GSC plate */
            msg = 55;
            break;
         default:
            msg = 36;
            if( n_stars)
               log_printf( log_file, "Error %d\n", n_stars);
            break;
         }
      log_printf( log_file, strings[msg]);
      exit_2( -4);
      }
   *n_stars_found = n_stars;
   log_printf( log_file, strings[95],  n_stars, ref_net_name);
                           /* "%d %s stars found\n" */
   return( ref_stars);
}

static void check_alt_azzes( const double *solar_ra_dec, const double *ra_dec,
                 const double jd, const double *latlon)
{
   int i;

   for( i = 0; i < 2; i++)
      {
      int is_troublesome_alt = 0;
      double alt_az[2], loc_epoch[2], hr_ang;

      full_ra_dec_to_alt_az( (DPT *)( i ? ra_dec : solar_ra_dec),
                (DPT *)alt_az, (DPT *)loc_epoch,
                (DPT *)latlon, jd, &hr_ang);
      log_printf( log_file, strings[111 + i], alt_az[1] * 180. / PI,
                                             alt_az[0] * 180. / PI + 180.);
      if( alt_az[1] * 180. / PI > -12. && !i)
         is_troublesome_alt = 1;
      if( alt_az[1] * 180. / PI < 10. && i)
         is_troublesome_alt = 1;
      if( is_troublesome_alt)
         {
         log_printf( log_file, strings[80 + i]);
         getch( );
         }
      }
}

int check_default_mpc_report_header( void)
{
   FILE *ifile = fopen( "header.mpc", "rb"), *ofile;

   if( ifile)
      fclose( ifile);
   else
      {
      char buff[64];
      int bytes_read;

      ifile = fopen( "header.def", "rb");
      ofile = fopen( "header.mpc", "wb");
      while( bytes_read = fread( buff, 1, 64, ifile))
         fwrite( buff, 1, bytes_read, ofile);
      fclose( ifile);
      fclose( ofile);
      printf( "\nWARNING:  You will have to set up your HEADER.MPC file to\n");
      printf( "reflect the header you want used when sending reports to the\n");
      printf( "Minor Planet Center.  If you check HEADER.MPC with a text\n");
      printf( "editor,  you will see full instructions on how to set this up.\n");
      printf( "In the meantime,  you can proceed using an example header...\n");
      printf( "just be sure to fix it before sending in MPC reports.\n\n");
      printf( "Hit any key to continue:\n");
      getch( );
      }
   return( 0);
}

#ifdef NOT_USED_ANYMORE
static void isophote_draw( const int is_move, const double x1,
               const double y1, const double step, const double hair_len)
{
   static double x0 = 0., y0 = 0., remains = 0.;

   if( !is_move)
      {
      const double dx = x1 - x0, dy = y1 - y0;
      const double len = sqrt( dx * dx + dy * dy);

      _moveto( (short)x0, (short)y0);
      _lineto( (short)x1, (short)y1);
      while( remains < len)
         {
         double x2 = x0 + remains * dx / len, y2 = y0 + remains * dy / len;

         _moveto( (short)x2, (short)y2);
         x2 -= hair_len * dy / len;
         y2 += hair_len * dx / len;
         _lineto( (short)x2, (short)y2);
         remains += step;
         }
      remains -= len;
      }
   x0 = x1;
   y0 = y1;
}
#endif


static void draw_circle( const short x, const short y, const short radius)
{
   _ellipse( _GBORDER, x - radius, y - radius, x + radius, y + radius);
}

static int reset_video_mode( const int new_mode,
                      CHARON_CONFIG *cconfig,
                      BMOUSE *bmouse, int *disp_loc, int *mouse_port)
{
   int rval = _setvideomode( (short)new_mode);
   static int first_time = 1;

   if( !rval)                    /* If video mode setting fails,  maybe */
      {                          /* trying VESA (or not trying it,  if  */
      extern short _SVGAType;    /* we were already VESA) will work     */

      _SVGAType ^= 1;
      rval = _setvideomode( (short)new_mode);
      }

   if( rval)
      {
      if( first_time)
         log_printf( log_file, "Video mode set, ");
      cconfig->video_mode = new_mode;
      _getvideoconfig( &vconfig);
      if( first_time)
         log_printf( log_file, "%d x %d pixels, ", vconfig.numxpixels,
                                                   vconfig.numypixels);
      set_palette( cconfig->color_table_number);
      if( first_time)
         log_printf( log_file, "palette set, ");
      disp_loc[2] = vconfig.numxpixels;
      disp_loc[3] = vconfig.numypixels;

      memset( bmouse, 0, sizeof( BMOUSE));
      bmouse->xmax = (unsigned)vconfig.numxpixels - 1;
      bmouse->ymax = (unsigned)vconfig.numypixels - 1;
      bmouse->x = bmouse->xmax >> 1;
      bmouse->y = bmouse->ymax >> 1;
      if( !*mouse_port)
         init_mouse( bmouse);
      if( first_time)
         log_printf( log_file, "mouse set\n");
      }
   first_time = 0;
   return( rval);
}

void main( int argc, char **argv)
{
   int i, c = 0, show_gsc_stars = 1, signed_ints = 0;
   int n_stars, n_gsc_stars = 0;
   int ra_dec_from_command_line = 0;
   int n_matched, do_initial_setup_menu = (argc == 1);
   PIXEL min_used, max_used, range;
   PIXEL ten_percent_pt, pixel_under_cursor = (PIXEL)0;
   int fixed_target_pixel = 0;
   double xc = 0, yc = 0, target_pixel[2];
   double xfactor, yfactor;
   double start_loc[7];
   FOUND_STAR *found = NULL;
   int disp_loc[4] = {0, 0, 319, 199};
   double img_loc[4] = {0, 0, 0, 0};
   int resid_thresh = 1;
   int text_on = 1;
   int mouse_port = 0, use_cdelta_values = 1;
   int n_video_sets = 0, err_code, kernel_size = 0;
   double *kernel = NULL;
   char *mask_filename = NULL, *isophote_file_name = NULL;
   IMAGE *image, *blink_images = NULL;
   int n_blinkable = 0, curr_blink = -1;
   int xpixel_col = 0, ypixel_col = 0, mag_col = 0;
   BMOUSE bmouse;
   REF_STAR *gsc_stars = NULL;
   double ra_dec[4], solar_ra_dec[2], original_pixel_sizes[2];
   double original_focal_len = 0.;
   FILE *ifile;
   char tbuff[180];
   char *message = NULL, message_buff[80];
   char *gsc_12_filename = NULL;
   char *config_file_name = "charon.dat";
   char *dark_filename = NULL, *flat_filename = NULL;
   long total_pixels;
   CHARON_CONFIG cconfig;
#ifdef __WATCOMC__
   extern short _SVGAType;
   char *env_ptr = getenv( "MOUSE");

   _SVGAType = 0;
   setvbuf( stdout, NULL, _IONBF, 0);
   log_file = fopen( "charon.log", "w");
   setvbuf( log_file, NULL, _IONBF, 0);

   log_printf( log_file, "CHARON version %s %s\n", __DATE__, __TIME__);

   if( env_ptr)
      mouse_port = atoi( env_ptr);
#endif
#ifdef DEBUG_MEM
   set_memory_error_function( memory_error);
#endif
   for( i = 0; i < argc; i++)
      if( argv[i][0] == '-')
         {
         int drop_this_argument = 0, j;

         if( argv[i][1] == 'd' && argv[i][2])
            {
            dark_filename = argv[i] + 2;
            if( *dark_filename == ' ')
               dark_filename++;
            }
         else if( argv[i][1] == 'F' && argv[i][2])
            {
            flat_filename = argv[i] + 2;
            if( *flat_filename == ' ')
               flat_filename++;
            }
         else if( argv[i][1] == 'D' && argv[i][2])
            {
            debug_level = atoi( argv[i] + 2);
            log_printf( log_file, "Debug level set to %d\n", debug_level);
            drop_this_argument = 1;
            }
         else if( argv[i][1] == 'C' && argv[i][2])
            {
            config_file_name = argv[i] + 2;
            log_printf( log_file, "Using configuration file %s\n",
                                    config_file_name);
            drop_this_argument = 1;
            }
         else if( argv[i][1] == 'L' && argv[i][2])
            {
            reset_language( argv[i][2]);
            log_printf( log_file, "Reset to language %c\n", argv[i][2]);
            drop_this_argument = 1;
            }
         else if( !stricmp( argv[i], "-sxc"))
            {
            extern int colour_starlight_express;

            colour_starlight_express = 1;
            log_printf( log_file,
                     "Starlight XPress images are assumed to be in color\n");
            drop_this_argument = 1;
            }
         else if( !stricmp( argv[i], "-guide"))
            {
            started_from_guide = 1;
            do_initial_setup_menu = 1;
            drop_this_argument = 1;
            }
#ifdef __WATCOMC__
         else if( !stricmp( argv[i], "-vesa"))
            {
            _SVGAType = 1;
            drop_this_argument = 1;
            }
#endif
         if( drop_this_argument)
            {
            argc--;
            for( j = i; j < argc; j++)
               argv[j] = argv[j + 1];
            i--;
            }
         }

   check_default_mpc_report_header( );

   if( check_mpc_report_header( "header.mpc"))
      exit_2( -5);

   if( debug_level)
      log_printf( log_file, "Initializing configuration...");
   init_config( &cconfig);
   cconfig.video_mode = _MRES256COLOR;
   image = (IMAGE *)calloc( 1, sizeof( IMAGE));

   if( debug_level)
      log_printf( log_file, "Loading %s...", config_file_name);

   err_code = load_configuration( &cconfig, image, config_file_name);
   if( err_code)
      {
      log_printf( log_file, "ERROR %d in loading configuration!\n", err_code);
      exit_2( -6);
      }

   if( debug_level)
      log_printf( log_file, "done\nLoading STARTUP.MAR...");

   correct_cconfig_with_startup( &cconfig);

   if( debug_level)
      log_printf( log_file, "Grabbing strings...");
   grab_strings( );
   if( debug_level)
      log_printf( log_file, "done\n");

   if( strlen( cconfig.obs_code) == 3)
      if( !get_observer_data( cconfig.obs_code, tbuff))
         if( strcmp( cconfig.obs_code, "247"))     /* roving observer */
            compute_latlon_from_mpc_data( tbuff, cconfig.latlon);

#ifdef NOT_NEEDED_NOW
   astnum_value = check_for_guide_cd( cconfig.cd_drive_path);
   if( astnum_value < 0)
      {
      log_printf( log_file, "%s\n%s\n", strings[23], strings[24]);
      log_printf( log_file, "Path used: %s\n",  cconfig.cd_drive_path);
      exit_2( -7);
      }
#endif

   if( argc > 1 && argv[1][0] != '-')     /* file name on command line */
      strcpy( cconfig.path_name, argv[1]);


   if( debug_level)
      log_printf( log_file, "Fixing arguments...\n");
   if( argc > 2)
      fix_arguments( argv + 2, argc - 2);
   if( debug_level)
      log_printf( log_file, "Arguments fixed\n");

   ra_dec[2] = ra_dec[3] = 0.;      /* assume non-moving object at first */

   for( i = 0; i < argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'G':
               if( argv[i][2])
                  {
                  gsc_12_filename = argv[i] + 2;
                  log_printf( log_file, strings[99], gsc_12_filename);
                            /* "Using GSC 1.2 file %s\n" */
                  }
               break;
            case 'K':
               kernel_size = load_kernel_file( argv[i] + 2, NULL);
               if( kernel_size > 1)
                  {
                  kernel = (double *)calloc( kernel_size * kernel_size,
                                                      sizeof( double));
                  load_kernel_file( argv[i] + 2, kernel);
                  }
               break;
            case 'l':
               image->is_xy_list = 1;
               sscanf( argv[i] + 2, "%d,%d", &xpixel_col, &ypixel_col, &mag_col);
               log_printf( log_file, "Expecting text file\n");
               break;
            case 'm':
               log_printf( log_file, "Ignoring pixel sizes in FITS header\n");
               use_cdelta_values = 0;
               break;
            case 'M':
               if( !argv[i][3])        /* single character */
                  mouse_port = atoi( argv[i] + 2);
               else
                  mask_filename = argv[i] + 2;
               break;
            case 'o':
               strcpy( cconfig.target_name, argv[i] + 2);
               ra_dec_from_command_line = 1;
               break;
            case 'p':
               if( err_code = parse_ra_dec( argv[i] + 2, &xc, &yc))
                  {
                  log_printf( log_file, "Couldn't parse '%s': %d\n",
                                         argv[i] + 2, err_code);
                  exit_2( -8);
                  }
               ra_dec_from_command_line = 1;
               break;
            case 'Q':
               sscanf( argv[i] + 2, "%lf,%lf", target_pixel, target_pixel + 1);
               fixed_target_pixel = 1;
               break;
            case 'u':
               overlay_file = fopen( argv[i] + 2, "ab");
               break;
            case '@':
               signed_ints = 1;
               break;
            case '=':
               do_initial_setup_menu = 1;
               break;
            default:
               if( parse_charon_command( &cconfig, image, argv[i]))
                  {
                  log_printf( log_file, strings[25], argv[i]);
                  exit_2( -9);
                  }
               break;
            }

   if( debug_level)
      log_printf( log_file, "Arguments read\n");

   if( do_initial_setup_menu)
      {
      if( debug_level)
         log_printf( log_file, "Entering setup menu...");
      setup_menu( &cconfig, image);
      if( debug_level)
         log_printf( log_file, "done\n");
      image->focal_len *= 25.4;  /* cvt to mm */
      write_out_startup( config_file_name, &cconfig, image);
      image->focal_len /= 25.4;  /* cvt back to inches */
      }

   original_pixel_sizes[0] = image->pixel_xsize;
   original_pixel_sizes[1] = image->pixel_ysize;
   original_focal_len = image->focal_len;

   log_printf( log_file, "is_xy_list = %d\n", image->is_xy_list);
   if( !image->is_xy_list)
      if( load_image( cconfig.path_name, image))
         {
         log_printf( log_file, strings[33], cconfig.path_name);
         exit_2( -10);
         }

   if( dark_filename)
      {
      IMAGE dark;

      memcpy( &dark, image, sizeof( IMAGE));
      if( load_image( dark_filename, &dark))
         {
         log_printf( log_file, strings[33], dark_filename);
         exit_2( -11);
         }
      dark_subtract_and_flat_field( image, &dark, NULL);
      free( dark.filename);
      free( dark.img);
      }


   if( kernel_size > 1)
      {
      log_printf( log_file, "Applying %d-square kernel...", kernel_size);
      apply_convolution( image, kernel_size, kernel);
      log_printf( log_file, "done\n");
      }

   if( !use_cdelta_values)
      {
      image->pixel_xsize = original_pixel_sizes[0];
      image->pixel_ysize = original_pixel_sizes[1];
      image->focal_len = original_focal_len;
      }
              /* Bev Ewen-Smith's Starlight XPress images,  and certain */
              /* FITS images, contain an RA/dec.  Unless that RA/dec is */
              /* overridden by the 'p' or 'o' command line switches,    */
              /* we'll make use of it:                                  */
   if( image->ra_dec_from_header && !ra_dec_from_command_line)
      {
      xc = image->xform[0];
      yc = image->xform[1];
      }

   if( image->time)
      image->time += (image->exposure_length * .5 / 86400.) *
                                    (double)( 1 - cconfig.file_time_point);
   image->time += cconfig.time_offset / 86400.;
   image->xform[0] = 1.;                     /* linear xformation */
   target_pixel[0] = target_pixel[1] = 0.;

   if( !xc && !yc)
      {
      int rval = find_obj( image->time, cconfig.target_name,
                    cconfig.cd_drive_path, ra_dec, cconfig.latlon);

      if( !rval)
         {
         xc = ra_dec[0];
         yc = ra_dec[1];
         log_printf( log_file, strings[29], xc, yc);
         }
      else
         {
         if( rval == -99)
            log_printf( log_file,
                  "Charon needs a Guide CD to find the object '%s'.\n",
                                                   cconfig.target_name);
         else
            {
            log_printf( log_file, strings[28], cconfig.target_name);
            log_printf( log_file, "rval %d\n", rval);
            }
         exit_2( -12);
         }
      }

   if( !xc && !yc)
      {
      log_printf( log_file, "%s\n%s\n", strings[21], strings[22]);
      exit_2( -13);
      }

   ra_dec[0] = -xc * PI / 180.;
   ra_dec[1] = yc * PI / 180.;
   ra_dec[2] *= -PI / 180.;
   ra_dec[3] *= PI / 180.;
   find_obj( image->time, "sun", cconfig.cd_drive_path,
                                 solar_ra_dec, cconfig.latlon);
   solar_ra_dec[0] *= -PI / 180.;
   solar_ra_dec[1] *=  PI / 180.;

   if( image->time)
      check_alt_azzes( solar_ra_dec, ra_dec, image->time, cconfig.latlon);

   if( !image->pixel_xsize || !image->pixel_ysize)
      {
      log_printf( log_file, "%s\n", strings[60]);
      exit_2( -14);
      }
   else
      {
      double xsize_in_arcseconds = 3600. * (180. / PI) *
                        image->pixel_xsize / (image->focal_len * 25400.);
      double ysize_in_arcseconds = 3600. * (180. / PI) *
                        image->pixel_ysize / (image->focal_len * 25400.);
      log_printf( log_file, strings[100], xsize_in_arcseconds,
                                          ysize_in_arcseconds);
                         /* "Pixels are %.2lf x %.2lf arcseconds\n" */
      log_printf( log_file, strings[101],
                          xsize_in_arcseconds * (double)image->xsize / 60.,
                          ysize_in_arcseconds * (double)image->ysize / 60.);
                 /* "Total field is %.2lf x %.2lf arcminutes\n" */
      }

   if( !image->focal_len)
      {
      log_printf( log_file, "%s\n", strings[59]);
      exit_2( -15);
      }

   if( !image->is_xy_list)
      offset_image( image, signed_ints);

   if( mask_filename)
      {
      log_printf( log_file, "Applying mask file '%s'...", mask_filename);
      i = apply_polygon_mask( image, mask_filename);
      log_printf( log_file, i ? "MASK FILE FAILED\n" : "done\n");
      }

   using_a20 = (cconfig.catalog_used == CATALOG_A2);

   found = calloc( sizeof( FOUND_STAR), cconfig.max_stars + 2);
   if( image->time)
      log_printf( log_file, "%s %.5lf\n", strings[123], image->time);
   gsc_stars = grab_catalog_stars( xc, yc, &cconfig, gsc_12_filename,
                                            image->time, &n_gsc_stars);

   start_loc[0] = xc;
   start_loc[1] = yc;
   start_loc[2] = cconfig.search_dist;
   start_loc[3] = cconfig.scale_tolerance * cconfig.scale_tolerance;
   start_loc[4] = cconfig.max_residual;
   start_loc[5] = cconfig.max_tilt;   /* default is any tilt is acceptable */
   start_loc[6] = cconfig.tilt_angle; /* default is 0 (aligned N/S) */
   image->xform[0] = (double)cconfig.fit_order;

   if( image->is_xy_list)
      {
      ifile = fopen( cconfig.path_name, "rb");
      i = 0;
      if( image->image_type == 'p')        /* .MAP file */
         {
         while( fgets( tbuff, sizeof( tbuff), ifile))
            if( !memcmp( tbuff, "date ", 5))
               while( fgets( tbuff, sizeof( tbuff), ifile))
                  if( sscanf( tbuff, "%lf %lf %ld", &found[i].x,
                            &found[i].y, &found[i].bright) == 3)
                     {
                     if( debug_level > 2)
                        log_printf( log_file, "%3d: %lf %lf %ld \n",
                                  i, found[i].x, found[i].y, found[i].bright);
                     if( signed_ints)
                        found[i].bright = -found[i].bright;
                     found[i].x *= image->focal_len;
                     found[i].y *= image->focal_len;
                     i++;
                     }
         }
      else if( image->image_type == 'n')        /* .FIN file */
         {
         double xpixel, ypixel, mag;

         image->focal_len *= 25.4;    /* convert focal len from inches to mm */
         xfactor = image->pixel_xsize / (image->focal_len * 1000.);
         yfactor = image->pixel_ysize / (image->focal_len * 1000.);
         xfactor *= (180. / PI) * 3600.;
         yfactor *= (180. / PI) * 3600.;
         while( fgets( tbuff, sizeof( tbuff), ifile))
            if( sscanf( tbuff + 5, "%lf%lf%lf", &xpixel, &ypixel, &mag) == 3)
               {
               found[i].x = xpixel * xfactor;
               found[i].y = ypixel * yfactor;
               found[i].bright = -(long)( atof( tbuff + 25) * 100. + .5);
               i++;
               }
         }
      else if( xpixel_col)
         {
         image->focal_len *= 25.4;    /* convert focal len from inches to mm */
         xfactor = image->pixel_xsize / (image->focal_len * 1000.);
         yfactor = image->pixel_ysize / (image->focal_len * 1000.);
         xfactor *= (180. / PI) * 3600.;
         yfactor *= (180. / PI) * 3600.;
         log_printf( log_file, "Getting data...\n");
         while( fgets( tbuff, sizeof( tbuff), ifile) && i < cconfig.max_stars)
            {
            double x = atof( tbuff + xpixel_col - 1);
            double y = atof( tbuff + ypixel_col - 1);
            double mag = atof( tbuff + mag_col - 1);

            if( x && y && strlen( tbuff) > ypixel_col)
               {
                  /* "focal length" becomes "arcsec/pixel" */
               found[i].x = x * xfactor;
               found[i].y = y * yfactor;
               found[i].bright = (long)( mag * 100. + .5);
               i++;
               log_printf( log_file, "%lf %lf\n", x, y);
               }
            }
         }
      fclose( ifile);
      n_stars = i;
      if( debug_level)
         log_printf( log_file, "%d map stars;  rescaling by %.5lf\n",
                              n_stars, image->focal_len);
      if( n_stars < cconfig.n_stars_used)
         cconfig.n_stars_used = n_stars;

      if( image->is_inverted)
         for( i = 0; i < n_stars; i++)
            found[i].y = -found[i].y;

      if( debug_level)
         log_printf( log_file, "Sorting...");
      qsort( found, n_stars, sizeof( FOUND_STAR), found_compare);
      if( debug_level)
         {
         log_printf( log_file, "Done\n");
         for( i = 0; i < cconfig.n_stars_used; i++)
            log_printf( log_file, "%3d: %lf %lf %ld \n", i, found[i].x,
                                    found[i].y, found[i].bright);
         }

      use_input_mags = 1;
      n_matched = get_initial_match( found, n_stars, image->xform,
                start_loc, cconfig.n_stars_used, gsc_stars, n_gsc_stars, 0);
//    image->focal_len = 100.;
      if( n_matched > 1)
         show_angle_and_fl( image);
      if( n_matched > 1 && image->image_type == 'n')        /* .FIN file */
         {
         double xpixel, ypixel;
         FILE *ofile;

         ifile = fopen( image->filename, "rb");
         strcpy( strchr( image->filename, '.'), ".fi2");
         ofile = fopen( image->filename, "wb");
         i = 0;
         while( fgets( tbuff, sizeof( tbuff), ifile))
            {
            if( sscanf( tbuff, "%lf%lf", &xpixel, &ypixel) == 2)
               {
               double ra_dec[2];
               char tstr[80];

               pixel_to_ra_dec( ra_dec, (double *)&found[i], image->xform);
               put_ra_dec_in_str( tstr,
                       ra_dec[0], ra_dec[1], cconfig.ra_dec_format);
               strcat( tstr, "  ");
               fwrite( tstr, strlen( tstr), 1, ofile);
               i++;
               }
            fwrite( tbuff, strlen( tbuff), 1, ofile);
            }
         fclose( ifile);
         fclose( ofile);
         }

      exit_2( 0);
      }
                                /* now we _know_ we're handling an image */
   image->focal_len *= 25.4;    /* convert focal len from inches to mm */
   xfactor = image->pixel_xsize / (image->focal_len * 1000.);
   yfactor = image->pixel_ysize / (image->focal_len * 1000.);
   total_pixels = (long)image->xsize * (long)image->ysize;
   ten_percent_pt = find_histo( image, cconfig.cutoff_point);
   saturation_point = cconfig.saturation_point;
                     /* 22 Dec 99:  these two lines should replace the */
                     /* #ifdef OLD_FLIP_CODE commented out below       */
   if( image->is_inverted)
      flip_image( image);
   n_stars = find_stars( image->img, image->xsize, image->ysize,
                   ten_percent_pt, found, cconfig.max_stars, cconfig.cell_size);

   log_printf( log_file, strings[35], n_stars);

   xfactor *= (180. / PI) * 3600.;
   yfactor *= (180. / PI) * 3600.;

   for( i = 0; i < n_stars; i++)
      {
      found[i].x *= xfactor;
      found[i].y *= yfactor;
      }
   if( n_stars > 3)
      {
      extern int gsc_photo_band;

      n_matched = get_initial_match( found, n_stars, image->xform,
                  start_loc, cconfig.n_stars_used, gsc_stars, n_gsc_stars,
                  gsc_photo_band ? gsc_photo_band : cconfig.photometric_band);
      }
   else
      n_matched = 0;

   for( i = 0; i < n_stars; i++)
      {
      found[i].x /= xfactor;
      found[i].y /= yfactor;
      }

   if( n_matched > 1)
      show_angle_and_fl( image);
   else
      {
      printf( "\nCharon failed to get a good pattern match.  Check your\n");
      printf( "settings and image,  and try again.\n");
      printf( "%s\n", strings[39]);
                        /* For debugging,  it may be helpful to allow */
                        /* one to continue by hitting 'c': */
      if( getch( ) != 'c')
         exit_2( -20);
      }

   if( n_matched < 5)
      {
      printf( strings[114]);
      printf( strings[115]);
      printf( strings[116]);
      printf( strings[117]);
/*             WARNING:  Charon probably did not accurately\n
               pattern-match catalog stars to your image. Check the\n
               display carefully to confirm that the match is an\n
               accurate one.\n                                          */
      printf( "%s\n", strings[39]);
      getch( );
      }

   min_used = max_used = image->img[0];
   for( i = image->xsize * image->ysize - 1; i; i--)
      {
      if( image->img[i] > max_used) max_used = image->img[i];
      if( image->img[i] < min_used) min_used = image->img[i];
      }
   image->low_end = find_histo( image, .6);
   image->high_end = find_histo( image, .99);
   if( image->low_end == image->high_end)
      image->high_end++;
   log_printf( log_file, strings[37], min_used, max_used);
   log_printf( log_file, strings[38], image->low_end, image->high_end);
   image->tilt_angle = atan2( image->xform[7], image->xform[6]);
   log_printf( log_file, strings[110],   /* "Tilt angle: %.2lf degrees\n" */
                  image->tilt_angle * 180. / PI);
   log_printf( log_file, strings[39]);
   if( getch( ) == 27)
      exit_2( -16);
   if( !reset_video_mode( cconfig.video_mode, &cconfig, &bmouse,
                                         disp_loc, &mouse_port))
      {
      log_printf( log_file, strings[74]);
      exit_2( -17);
      }
   n_video_sets++;

   if( mouse_port)
      {
      mouse_install( mouse_port);
      log_printf( log_file, "Mouse installed\n");
      bmouse.mouse_installed = 2;
      global_mouse_x = bmouse.x;
      global_mouse_y = bmouse.y;
      }

   img_loc[2] = image->ysize * disp_loc[2] / disp_loc[3];
   img_loc[3] = image->ysize;

   while( c != 27)
      {
      double img_width = img_loc[2] - img_loc[0];
      double img_height = img_loc[3] - img_loc[1];
      int redrawn = 1, nearest_target, nearest_gsc, prev_button = 0;
      int prev_image_star = -1;
      char image_star_buff[44];
      double pixel[2],  zoom_factor = 1.;
      double image_to_display_x =
                  (double)disp_loc[2] / img_width;
      double image_to_display_y =
                  (double)disp_loc[3] / img_height;
      double target_mag = 0., target_ra_dec[2], clicked_ra_dec[2];

      log_printf( log_file, "Showing image; ");
      show_image( image->img, image->xsize, image->ysize,
                   image->low_end, image->high_end, img_loc, disp_loc,
                   vconfig.numcolors);
//    log_printf( log_file, "Image shown ");

      if( isophote_file_name &&
                       (ifile = fopen( isophote_file_name, "rb")) != NULL)
         {
         int is_move = 0, n_crossings = 0, is_edge_obj = 0, n_inside = 0.;
         double x1 = 0., y1 = 0., radius = 0.;
         double prev_x = 0., prev_y = 0.;
         double avg_x = 0., avg_y = 0., avg_r = 0.;
         double avg_r2 = 0., avg_xr = 0., avg_yr = 0.;

         while( fgets( tbuff, sizeof( tbuff), ifile))
            if( *tbuff == 'N')
               {
               if( !is_edge_obj)
                  {
                  _setcolor( (n_crossings & 1) ? BLUE_INDEX : GREEN_INDEX);
                  draw_circle( (short)x1, (short)y1, (short)radius);
                  _moveto( (short)( x1 - radius / 3.), (short)y1);
                  _lineto( (short)( x1 + radius / 3.), (short)y1);
                  _moveto( (short)x1, (short)( y1 - radius / 3.));
                  _lineto( (short)x1, (short)( y1 + radius / 3.));
                  if( n_crossings & 1)
                     {
                     avg_x += x1;
                     avg_y += y1;
                     avg_r += radius;
                     avg_r2 += radius * radius;
                     avg_xr += x1 * radius;
                     avg_yr += y1 * radius;
                     n_inside++;
                     }
                  }
               is_edge_obj = (tbuff[3] == '0');
               if( !is_edge_obj)
                  {
                  double area;

                  sscanf( tbuff + 10, "%lf%lf%lf", &area, &x1, &y1);
                  radius = sqrt( fabs( area) / PI) * image_to_display_x;
                  x1 = (x1 - img_loc[0]) * image_to_display_x + .5;
                  y1 = (y1 - img_loc[1]) * image_to_display_y + .5;
                  fgets( tbuff, sizeof( tbuff), ifile);
                  }
               _setcolor( is_edge_obj ? BLUE_INDEX : RED_INDEX);
               is_move = 1;
               n_crossings = 0;
               }
            else
               {
               double x, y;
               const double xcenter = (disp_loc[0] + disp_loc[2]) * .5;
               const double ycenter = (disp_loc[1] + disp_loc[3]) * .5;

               sscanf( tbuff, "%lf%lf", &x, &y);
               x = (x / 1000. - img_loc[0]) * image_to_display_x + .5;
               y = (y / 1000. - img_loc[1]) * image_to_display_y + .5;
#ifdef NOT_USED_ANYMORE
               isophote_draw( is_move, x, y, 5., 3.);
#else
               if( is_move)
                  _moveto( (short)x, (short)y);
               else
                  _lineto( (short)x, (short)y);
#endif
               if( !is_move)
                  if( y >= ycenter && prev_y < ycenter ||
                                         y < ycenter && prev_y >= ycenter)
                     if( x + (prev_x - x) * (ycenter - y) / (prev_y - y) < xcenter)
                        n_crossings++;
               is_move = 0;
               prev_x = x;
               prev_y = y;
               }
         fclose( ifile);
         if( n_inside > 1)
            {
            double determ;

            avg_x /= (double)n_inside;
            avg_y /= (double)n_inside;
            avg_r /= (double)n_inside;
            avg_r2 /= (double)n_inside;
            avg_xr /= (double)n_inside;
            avg_yr /= (double)n_inside;
            determ = avg_r2 - avg_r * avg_r;
            if( determ)
               {
               FOUND_STAR *f = found + n_stars - 1;

               x1 = (avg_r2 * avg_x - avg_xr * avg_r) / determ;
               y1 = (avg_r2 * avg_y - avg_yr * avg_r) / determ;
               _setcolor( BLUE_INDEX);
               _moveto( (short)(x1 - 15.), (short)(y1 - 15.));
               _lineto( (short)(x1 + 15.), (short)(y1 + 15.));
               _moveto( (short)(x1 - 15.), (short)(y1 + 15.));
               _lineto( (short)(x1 + 15.), (short)(y1 - 15.));
               draw_circle( (short)x1, (short)y1, 5);
               f->x = x1 / image_to_display_x + img_loc[0];
               f->y = y1 / image_to_display_y + img_loc[1];
               f->bright = -3;
               }
            }
         }
      if( cconfig.targeting_on && curr_blink == -1)
         {
         _setcolor( RED_INDEX);
         for( i = 0; i < n_stars; i++)
            {
            double r = 5. * sqrt(
                       (double)found[i].bright / (double)found->bright);
            double x = found[i].x - img_loc[0];
            double y = found[i].y - img_loc[1];

            x *= image_to_display_x;
            y *= image_to_display_y;
            x += .5;
            y += .5;
            r += .5;
            draw_circle( (short)x, (short)y, (short)r + 1);
            }
         if( !fixed_target_pixel)
            {
            ra_dec_to_pixel( target_pixel, xc, yc, image->xform);
            target_pixel[0] /= xfactor;
            target_pixel[1] /= yfactor;
            }
         pixel[0] = (target_pixel[0] - img_loc[0]) *
                                      image_to_display_x;
         pixel[1] = (target_pixel[1] - img_loc[1]) *
                                      image_to_display_y;
         _moveto( 0, (short)pixel[1]);
         _lineto( vconfig.numxpixels, (short)pixel[1]);
         _moveto( (short)pixel[0], 0);
         _lineto( (short)pixel[0], vconfig.numypixels);
         }
      _setcolor( GREEN_INDEX);
      if( show_gsc_stars)
         for( i = 0; i < n_gsc_stars; i++)
            {
            double x, y, r = 5.5, gsc_pixel[2];

            ra_dec_to_pixel( gsc_pixel, gsc_stars[i].ra, gsc_stars[i].dec,
                                   image->xform);
            gsc_pixel[0] /= xfactor;
            gsc_pixel[1] /= yfactor;

            x = gsc_pixel[0] - img_loc[0];
            y = gsc_pixel[1] - img_loc[1];
            x *= image_to_display_x;
            y *= image_to_display_y;
            _moveto( (short)( x - r), (short)y);
            _lineto( (short)( x + r), (short)y);
            _moveto( (short)x, (short)( y - r));
            _lineto( (short)x, (short)( y + r));
            }

      if( message)
         {
         _settextposition( vconfig.numtextcols, 1);
         _settextcolor( YELLOW_INDEX);
         _outtext( message);
         message = NULL;
         }
      range = image->high_end - image->low_end;

      _setcolor( YELLOW_INDEX);
      if( mouse_port)
         directly_read_mouse( &bmouse);
      else
         mouse_read( &bmouse);
      show_cursor( bmouse.x, bmouse.y, 0);
      while( !(bmouse.accum_released & 5) && !kbhit( ))
         {
         if( mouse_port)
            directly_read_mouse( &bmouse);
         else
            mouse_read( &bmouse);
         pixel[0] = bmouse.x;
         pixel[1] = bmouse.y;
         if( (bmouse.button | prev_button) & 5)
            {
            pixel[0] = bmouse.pressed_x;
            pixel[1] = bmouse.pressed_y;
            }
         pixel[0] = img_loc[0] + pixel[0] / image_to_display_x;
         pixel[1] = img_loc[1] + pixel[1] / image_to_display_y;

         if( bmouse.dx || bmouse.dy || redrawn || bmouse.button != prev_button)
            {
            int dist, best_dist = 32767;
            int zoom_box[4];
            double x_value_to_display = pixel[0];
            double y_value_to_display = pixel[1];

            if( cconfig.pixel_origin_at_top_left & 1)
               y_value_to_display = (double)image->ysize - y_value_to_display;
            if( cconfig.pixel_origin_at_top_left & 2)
               x_value_to_display = (double)image->xsize - x_value_to_display;

            _settextcolor( YELLOW_INDEX);
            _settextposition( 1, vconfig.numtextcols - TEXT_COLUMNS);
            sprintf( tbuff, "%7.2lf %7.2lf ", x_value_to_display,
                                              y_value_to_display);
            if( pixel[0] >= 0. && pixel[1] >= 0. &&
                    pixel[0] < (double)image->xsize &&
                    pixel[1] < (double)image->ysize)
               {
               pixel_under_cursor = image->img[
                           (int)pixel[0] + (int)pixel[1] * image->xsize];
               sprintf( tbuff + 16, "%10u", pixel_under_cursor);
               }
            else
               strcpy( tbuff + 16, "          ");
            if( text_on)
               _outtext( tbuff);

            pixel[0] *= xfactor;
            pixel[1] *= yfactor;
            pixel_to_ra_dec( ra_dec, pixel, image->xform);
            put_ra_dec_in_str( tbuff, ra_dec[0], ra_dec[1], cconfig.ra_dec_format);
            if( bmouse.button != 2)
               memcpy( clicked_ra_dec, ra_dec, 2 * sizeof( double));
            _settextposition( 2, vconfig.numtextcols - TEXT_COLUMNS);
            if( text_on)
               _outtext( tbuff);
            pixel[0] /= xfactor;
            pixel[1] /= yfactor;

            nearest_gsc = -1;
            for( i = 0; i < n_gsc_stars; i++)
               {
               double gsc_pixel[2];

               ra_dec_to_pixel( gsc_pixel, gsc_stars[i].ra, gsc_stars[i].dec,
                                               image->xform);
               gsc_pixel[0] /= xfactor;
               gsc_pixel[1] /= yfactor;
               if( (dist = short_dist( gsc_pixel[0] - pixel[0],
                               gsc_pixel[1] - pixel[1])) < best_dist)
                  {
                  best_dist = dist;
                  nearest_gsc = i;
                  }
               }
            if( nearest_gsc != -1)
               {
               REF_STAR *tstar = gsc_stars + nearest_gsc;

               make_star_desig( tstar, tbuff);
               strcat( tbuff, "         ");

               sprintf( tbuff + 18, "mag %2d.%02d%c",
                      tstar->mag / 100, tstar->mag % 100, tstar->photo_band);
               if( cconfig.catalog_used == CATALOG_A1
                        || cconfig.catalog_used == CATALOG_A2)
                  {                                /* for Ax.0,  show mag to */
                  tbuff[26] = tstar->photo_band;   /* .1 precision,  not .01 */
                  tbuff[27] = '\0';
                  }
               _settextcolor( GREEN_INDEX);
               _settextposition( 3, vconfig.numtextcols - TEXT_COLUMNS);
               if( text_on)
                  _outtext( tbuff);
               put_ra_dec_in_str( tbuff, tstar->ra, tstar->dec,
                                 cconfig.ra_dec_format);
               _settextposition( 4, vconfig.numtextcols - TEXT_COLUMNS);
               if( text_on)
                  _outtext( tbuff);
               }

            nearest_target = -1;
            best_dist = 32767;
            for( i = 0; i < n_stars; i++)
               if( (dist = short_dist( found[i].x - pixel[0],
                               found[i].y - pixel[1])) < best_dist)
                  {
                  best_dist = dist;
                  nearest_target = i;
                  }
            if( nearest_target != -1)
               {
               double star_pixel[2], x_resid, y_resid;

               star_pixel[0] = found[nearest_target].x * xfactor;
               star_pixel[1] = found[nearest_target].y * yfactor;
               pixel_to_ra_dec( target_ra_dec, star_pixel, image->xform);
               put_ra_dec_in_str( tbuff,
                       target_ra_dec[0], target_ra_dec[1], cconfig.ra_dec_format);
               _settextcolor( RED_INDEX);
               _settextposition( 5, vconfig.numtextcols - TEXT_COLUMNS);
               if( text_on)
                  _outtext( tbuff);
#ifndef TEMPORARY_CHECK_FOR_RMN_PROBLEM
               y_resid = 3600. * (target_ra_dec[1] - yc);
               x_resid = 3600. * (target_ra_dec[0] - xc) * cos( yc * PI / 180.);
#else
               x_resid = found[nearest_target].x - pixel[0];
               y_resid = found[nearest_target].y - pixel[1];
#endif
               if( fabs( y_resid) < RESIDUAL_THRESH &&
                            fabs( x_resid) < RESIDUAL_THRESH || !resid_thresh)
                  {
                  char *format_text = "%6.2lf %6.2lf";

                  if( fabs( x_resid) > 99. || fabs( y_resid) > 99.)
                     format_text = "%6.1lf %6.1lf";
                  if( fabs( x_resid) > 999. || fabs( y_resid) > 999.)
                     format_text = "%6.0lf %6.0lf";
                  sprintf( tbuff, format_text, x_resid, y_resid);
                  }
               else
                  strcpy( tbuff, "             ");
               _settextposition( 6, vconfig.numtextcols - TEXT_COLUMNS);
               if( text_on)
                  _outtext( tbuff);

               _settextposition( 6, vconfig.numtextcols - 14);
               if( found[nearest_target].bright > 0)
                  {
                  target_mag = image->xform[1] - 2.5 *
                                log10( (double)found[nearest_target].bright);

                  sprintf( tbuff, " mag %5.2lf     ", target_mag);
                  }
               else
                  {
                  sprintf( tbuff, " no mag      ");
                  target_mag = 0.;
                  }

               if( bmouse.button == 2)
                  {
                  double dist, posn_ang, radians[4];

                  radians[0] = clicked_ra_dec[0];
                  radians[1] = clicked_ra_dec[1];
                  radians[2] = ra_dec[0];
                  radians[3] = ra_dec[1];
                  for( i = 0; i < 4; i++)
                     radians[i] *= PI / 180.;
                  calc_dist_and_posn_ang( radians, radians + 2,
                                                  &dist, &posn_ang);

                  dist *= 3600. * (180. / PI);
                  posn_ang = -posn_ang; /* 27 Apr 98: GJG pointed out error */
                  if( posn_ang < 0.)
                     posn_ang += PI + PI;
                  sprintf( tbuff, "%8.2lf\" %5.1lf", dist, posn_ang * 180. / PI);
                  }

               if( text_on)
                  _outtext( tbuff);

               }

            show_cursor( bmouse.prev_x, bmouse.prev_y, 1);
            if( prev_image_star != -1 && cconfig.targeting_on)
               remove_x_from_screen( image_star_buff);
            if( prev_button == 2)
               show_xorred_line( bmouse.pressed_x, bmouse.pressed_y,
                           bmouse.prev_x, bmouse.prev_y);
            if( bmouse.button == 2)
               show_xorred_line( bmouse.pressed_x, bmouse.pressed_y,
                           bmouse.x, bmouse.y);
            zoom_factor = 1.;
            if( prev_button & 5)
               zoom_factor = show_zoom_box( bmouse.pressed_x, bmouse.pressed_y,
                            bmouse.prev_x, bmouse.prev_y, zoom_box,
                            prev_button == 4);
            if( bmouse.button & 5)
               zoom_factor = show_zoom_box( bmouse.pressed_x, bmouse.pressed_y,
                           bmouse.x, bmouse.y, zoom_box, bmouse.button == 4);
            if( nearest_target != -1 && cconfig.targeting_on)
               {
               double x = found[nearest_target].x - img_loc[0];
               double y = found[nearest_target].y - img_loc[1];

               x = x * image_to_display_x + .5;
               y = y * image_to_display_y + .5;
               prev_image_star = nearest_target;
               place_x_on_screen( image_star_buff, (short)x, (short)y,
                                                   RED_INDEX);
               }
            show_cursor( bmouse.x, bmouse.y, 0);
            prev_button = bmouse.button;
            }
         redrawn = 0;
         }

      if( kbhit( ))
         {
         c = tolower( getch( ));
         if( !c)
            c = getch( ) + 256;
         }
      else           /* mouse click */
         {
         img_width  /= zoom_factor;
         img_height /= zoom_factor;
         img_loc[0] = pixel[0] - img_width / 2.;
         img_loc[1] = pixel[1] - img_height / 2.;
         img_loc[2] = img_loc[0] + img_width;
         img_loc[3] = img_loc[1] + img_height;
         c = 0;
         }

      switch( c)
         {
         case '+':
         case '-':
         case '*':
         case '/':
            {
            PIXEL delta = range / 5 + 1;
            double new_high, new_low;

            if( c == '+' || c == '*')
               new_high = (double)image->high_end - (double)delta;
            else
               new_high = (double)image->high_end + (double)delta;
            if( c == '+' || c == '/')
               new_low = (double)image->low_end - (double)delta;
            else
               new_low = (double)image->low_end + (double)delta;
            if( new_low > new_high - 2.)
               new_low = new_high - 2.;

            if( new_low < 0)
               image->low_end = 0;
            else if( new_low > (double)MAX_PIXEL)
               image->low_end = MAX_PIXEL;
            else
               image->low_end = (PIXEL)new_low;
            if( new_high < 0)
               image->high_end = 0;
            else if( new_high > (double)MAX_PIXEL)
               image->high_end = MAX_PIXEL;
            else
               image->high_end = (PIXEL)new_high;
            }
            break;
         case 'b':
            if( !blink_images)
               {
               blink_images = get_blinkable_images( image, &n_blinkable,
                                          xfactor, yfactor);
               if( n_blinkable)
                  {
                  for( i = 0; i < n_blinkable; i++)
                     {
                     offset_image( blink_images + i, signed_ints);
                     match_contrasts( blink_images + i, image);
                     }
                  sprintf( message_buff, strings[102], n_blinkable);
                                       /* "%d blink image(s) loaded" */
                  message = message_buff;
                  }
               else
                  message = strings[103];  /* "No blinkable image found" */
               }
            else
               {
               for( i = 0; i < 2; i++)
                  {
                  if( curr_blink != -1)
                     {
                     PIXEL *tpixel = blink_images[curr_blink].img;
                     char *tname = blink_images[curr_blink].filename;
                     int tlow = blink_images[curr_blink].low_end;
                     int thigh = blink_images[curr_blink].high_end;
                     int j;
                     double new_pixel[2], tdouble;

                     pixel[0] = (img_loc[0] + img_loc[2]) * xfactor / 2.;
                     pixel[1] = (img_loc[1] + img_loc[3]) * yfactor / 2.;
                     pixel_to_ra_dec( ra_dec, pixel, image->xform);

                     blink_images[curr_blink].img      = image->img;
                     blink_images[curr_blink].filename = image->filename;
                     blink_images[curr_blink].low_end  = image->low_end;
                     blink_images[curr_blink].high_end  = image->high_end;
                     image->img = tpixel;
                     image->filename = tname;
                     image->low_end = tlow;
                     image->high_end = thigh;
                     for( j = 0; j < N_FIND_PARAMS; j++)
                        {
                        tdouble = blink_images[curr_blink].xform[j];

                        blink_images[curr_blink].xform[j] = image->xform[j];
                        image->xform[j] = tdouble;
                        }

                     ra_dec_to_pixel( new_pixel, ra_dec[0], ra_dec[1],
                                                image->xform);
                     new_pixel[0] /= xfactor;
                     new_pixel[1] /= yfactor;
                     pixel[0] /= xfactor;
                     pixel[1] /= yfactor;
                     for( j = 0; j < 4; j++)
                        img_loc[j] += new_pixel[j & 1] - pixel[j & 1];
                     xc = image->xform[3];
                     yc = image->xform[4];

                                          /* And now,  to be bold... */
                     ten_percent_pt = find_histo( image, cconfig.cutoff_point);
                     saturation_point = cconfig.saturation_point;
//                   n_stars = find_stars( image->img, image->xsize,
//                               image->ysize, ten_percent_pt, found,
//                               cconfig.max_stars, cconfig.cell_size);
                     }
                  if( !i)     /* only on first pass! */
                     {
                     curr_blink++;
                     if( curr_blink == n_blinkable)
                        curr_blink = -1;
                     }
                  }
               message = image->filename;
               }
            break;
         case 'c':
         case ZKEY_F1:
            if( n_stars < cconfig.max_stars)
               n_stars++;
            add_a_star( image, (int)pixel[0], (int)pixel[1],
                                found + n_stars - 1, cconfig.cutoff_point,
                                (c == ZKEY_F1 ? -3 : cconfig.cell_size));
            message = strings[104];   /* "Star added at cursor" */
            if( c == ZKEY_F1)
               message = "Centroided using Bill Owen's comet method";
            break;
#ifdef NOT_YET
         case ZKEY_F2:
            message = "Click on other end of trail:";
            break;
#endif
         case ZKEY_F3:
            if( isophote_file_name)
               {
               unlink( isophote_file_name);
               isophote_file_name = NULL;
               }
            break;
         case ZKEY_F4:
            {
            FILE *ofile;

            isophote_file_name = "isophote.txt";
            ofile = fopen( isophote_file_name, "r+b");
            if( !ofile)
               ofile = fopen( isophote_file_name, "wb");
            else
               fseek( ofile, 0L, SEEK_END);
            find_isophotes( image->img, image->xsize, image->ysize,
                        pixel_under_cursor, ofile);
            fclose( ofile);
            }
            break;
         case 'e':
            resid_thresh ^= 1;
            message = (resid_thresh ? "30 arcsec limit on residuals" :
                                 "No limit on residuals");
            break;
         case 'f':
            cconfig.ra_dec_format = (cconfig.ra_dec_format + 1) % 3;
            message = strings[14 + cconfig.ra_dec_format];
            break;
         case 'g':
            show_gsc_stars ^= 1;
            message = strings[show_gsc_stars ? 122 : 121];
                   /* "GSC stars turned on" : "GSC stars turned off"; */
            break;
         case 'h':
            {
            extern int display_method;

            message = (display_method ? "Switched to linear contrast"
                                      : "Switched to logarithmic contrast");
            display_method ^= 1;
            }
            break;
         case 'l':
            text_on ^= 1;
            message = strings[118];      /* "Legend toggled" */
            break;
         case 'r':
            cconfig.color_table_number = (cconfig.color_table_number + 1)
                              % n_color_tables;
            set_palette( cconfig.color_table_number);
            get_palette_name( cconfig.color_table_number, message_buff);
            message = message_buff;
            break;
         case 't':
            cconfig.targeting_on ^= 1;
            message = strings[cconfig.targeting_on ? 120 : 119];
                         /* "Image stars shown" :  "Image stars shut off" */
            break;
         case ZKEY_NUM_DEL:
         case '.':
            if( nearest_target >= 0)
               {
               memmove( found + nearest_target, found + nearest_target + 1,
                         (n_stars - nearest_target) * sizeof( FOUND_STAR));
               n_stars--;
               message = "Image star deleted";
               }
            break;
         case ZKEY_NUM_5:
         case '5':
            img_loc[0] += img_width / 4;
            img_loc[2] -= img_width / 4;
            img_loc[1] += img_height / 4;
            img_loc[3] -= img_height / 4;
            break;
         case '0':
         case ZKEY_NUM_0:
            img_loc[0] -= img_width / 2;
            img_loc[2] += img_width / 2;
            img_loc[1] -= img_height / 2;
            img_loc[3] += img_height / 2;
            break;
         case ZKEY_NUM_1: case ZKEY_NUM_2: case ZKEY_NUM_3:
         case ZKEY_NUM_4:                  case ZKEY_NUM_6:
         case ZKEY_NUM_7: case ZKEY_NUM_8: case ZKEY_NUM_9:
            img_loc[0] += ZKEY_CURSOR_DX( c) * img_width / 4;
            img_loc[2] += ZKEY_CURSOR_DX( c) * img_width / 4;
            img_loc[1] += ZKEY_CURSOR_DY( c) * img_height / 4;
            img_loc[3] += ZKEY_CURSOR_DY( c) * img_height / 4;
            break;
         case ZKEY_NUM_CTRL_1: case ZKEY_NUM_CTRL_2: case ZKEY_NUM_CTRL_3:
         case ZKEY_NUM_CTRL_4:                       case ZKEY_NUM_CTRL_6:
         case ZKEY_NUM_CTRL_7: case ZKEY_NUM_CTRL_8: case ZKEY_NUM_CTRL_9:
            {
            static short keys[9] = { ZKEY_NUM_CTRL_1, ZKEY_NUM_CTRL_2,
                    ZKEY_NUM_CTRL_3, ZKEY_NUM_CTRL_4, ZKEY_NUM_CTRL_5,
                    ZKEY_NUM_CTRL_6, ZKEY_NUM_CTRL_7, ZKEY_NUM_CTRL_8,
                    ZKEY_NUM_CTRL_9 };
            for( i = 0; i < 9; i++)
               if( c == keys[i])
                  {
                  img_loc[0] += i % 3 - 1;
                  img_loc[2] += i % 3 - 1;
                  img_loc[1] += 1 - i / 3;
                  img_loc[3] += 1 - i / 3;
                  }
            }
            break;
         case ' ':
            img_loc[0] = img_loc[1] = 0;
            img_loc[2] = image->ysize * disp_loc[2] / disp_loc[3];
            img_loc[3] = image->ysize;
            break;
         case ZKEY_TAB:
            img_loc[0] = target_pixel[0] - 20 *  disp_loc[2] / disp_loc[3];
            img_loc[1] = target_pixel[1] - 20;
            img_loc[2] = target_pixel[0] + 20 *  disp_loc[2] / disp_loc[3];
            img_loc[3] = target_pixel[1] + 20;
            break;
         case 'p':      /* photometry:  two-decimal mags,  no posn */
         case 'o':
         case 'i':
            if( curr_blink > -1)
               message = "Can't do astrometry on 'blink' images!";
            else
               {
               double jd = image->time, drift;

               if( is_drift_scan)
                  {
                  if( fabs( image->tilt_angle) > 3. * PI / 4.)
                     drift = (image->xsize - found[nearest_target].x)
                                       * xfactor;
                  else if( image->tilt_angle > PI / 4.)  /* 45-135 degrees */
                     drift = found[nearest_target].y * yfactor;
                  else if( image->tilt_angle < -PI / 4.)  /* -45--135 degrees */
                     drift = (image->ysize - found[nearest_target].y)
                                       * yfactor;
                  else                                /* -45 to 45 degrees */
                     drift = found[nearest_target].x * xfactor;

                  drift /= 15 * cos( target_ra_dec[1] * PI / 180.);
                  jd += fabs( drift) / 86400.;
                  }
               message = write_report_data( &cconfig, jd, target_ra_dec,
                               target_mag, c, ref_net_name);
               }
            break;
         case 'x':
            cconfig.pixel_origin_at_top_left ^= 2;
            break;
         case 'y':
            cconfig.pixel_origin_at_top_left ^= 1;
            break;
         case 'z':
            {
            FILE *ofile = fopen( "images.dat", "a");
            int low_end, high_end, j;
            const int linear_terms[9] = { 3, 4, 5, 6, 7,
                    5 + N_FIND_COEFFS, 6 + N_FIND_COEFFS, 7 + N_FIND_COEFFS, 1 };

            fprintf( ofile, "L");
            for( i = 0; i < 9; i++)
               {
               double oval = image->xform[linear_terms[i]];

               if( i >= 2)
                  oval *= 1000000.;
               fprintf( ofile, " %9.5lf", oval);
               }
            fprintf( ofile, "\n");
            if( image->xform[1] > 1.5)    /* include quad terms */
               {                          /* (8,9,10, 18,19,20) */
               fprintf( ofile, "Q");
               for( j = 0; j <= N_FIND_COEFFS; j += N_FIND_COEFFS)
                  for( i = 8; i <= 10; i++)
                     fprintf( ofile, " %9.5lf", image->xform[i + j] * 1.e+9);
               fprintf( ofile, "\n");
               }
            if( image->xform[1] > 2.5)    /* include cubic terms */
               {                          /* (11-14, 21-24) */
               fprintf( ofile, "C");
               for( j = 0; j <= N_FIND_COEFFS; j += N_FIND_COEFFS)
                  for( i = 11; i <= 14; i++)
                     fprintf( ofile, " %9.5lf", image->xform[i + j] * 1.e+12);
               fprintf( ofile, "\n");
               }

            fprintf( ofile, " ");
            for( i = 0; i < 4; i++)
               {
               double corner_pixel[2], corner_ra_dec[2];

               corner_pixel[0] = ((i & 1) ? image->xsize : 0);
               corner_pixel[1] = ((i & 2) ? image->ysize : 0);
                     /* SBIG images are 'flipped'... */
               if( isdigit( image->image_type))
                  corner_pixel[1] = image->ysize - corner_pixel[1];
               if( image->is_inverted)
                  corner_pixel[1] = image->ysize - corner_pixel[1];
               corner_pixel[0] *= xfactor;
               corner_pixel[1] *= yfactor;
               pixel_to_ra_dec( corner_ra_dec, corner_pixel, image->xform);
               fprintf( ofile, "%8.4lf %8.4lf  ",
                                     corner_ra_dec[0], corner_ra_dec[1]);
               }
            low_end = (image->low_end - image->image_offset) & 0xffff;
            high_end = (image->high_end - image->image_offset) & 0xffff;
            if( signed_ints)
               {
               low_end ^= 0x8000;
               high_end ^= 0x8000;
               }
            fprintf( ofile, "%5u %5u  %5d %5d %9ld %d %s\n",
               low_end, high_end,
               image->xsize, image->ysize, image->data_offset,
               image->byte_order, image->filename);
            fclose( ofile);
            message = strings[61];
            }
            break;
         case 13:
            _setvideomode( _DEFAULTMODE);
            image->focal_len /= 25.4;
            setup_menu( &cconfig, image);
            image->focal_len *= 25.4;
            if( !reset_video_mode( cconfig.video_mode, &cconfig, &bmouse,
                                         disp_loc, &mouse_port))
               reset_video_mode( _MRES256COLOR, &cconfig, &bmouse,
                                         disp_loc, &mouse_port);
            n_video_sets++;
            img_loc[2] = image->ysize * disp_loc[2] / disp_loc[3];
            img_loc[3] = image->ysize;
            img_loc[0] = img_loc[1] = 0;
            break;
         case 'v': case 'V':
            {
            extern const int video_modes[];
            int n_modes, curr_idx, mode_set;

            for( i = 0; video_modes[i]; i++)
               if( video_modes[i] == cconfig.video_mode)
                  curr_idx = i;
            n_modes = i;
            do
               {
               if( c == 'v')
                  curr_idx = (curr_idx + 1) % n_modes;
               else
                  curr_idx = (curr_idx + n_modes - 1) % n_modes;
               cconfig.video_mode = video_modes[curr_idx];
               mode_set = reset_video_mode( video_modes[curr_idx], &cconfig,
                                  &bmouse, disp_loc, &mouse_port);
               } while( !mode_set);
            message = message_buff;
            sprintf( message, "%d x %d pixels, %d colors set",
                        vconfig.numxpixels, vconfig.numypixels,
                        vconfig.numcolors);
            }
            break;
         case 'w': case 'W':
            if( add_wcs_header_data( image))
               {
               sprintf( message_buff, "Couldn't write WCS data! '%c'",
                                 image->image_type);
               message = message_buff;
               }
            break;
         case '?':
            _setvideomode( _DEFAULTMODE);
            show_help_data( );
            getch( );
            _setvideomode( cconfig.video_mode);
            set_palette( cconfig.color_table_number);
            n_video_sets++;
            break;
         default:
            break;
         }
      }
   if( mouse_port)
      mouse_deinstall( );
   _setvideomode( _DEFAULTMODE);

   write_out_startup( config_file_name, &cconfig, image);

   if( gsc_stars)
      free( gsc_stars);
   if( found)
      free( found);
   if( kernel)
      free( kernel);
   if( image)
      {
      if( image->filename)
         free( image->filename);
      if( image->img)
         free( image->img);
      free( image);
      }
   for( i = 0; strings[i]; i++)
      free( strings[i]);

   if( blink_images)
      {
      for( i = 0; i < n_blinkable; i++)
         {
         free( blink_images[i].filename);
         free( blink_images[i].img);
         }
      free( blink_images);
      }

   if( n_video_sets && cconfig.millisec_per_video_reset)
      advance_clock( n_video_sets * cconfig.millisec_per_video_reset);

   if( isophote_file_name)
      unlink( isophote_file_name);
#ifdef DEBUG_MEM
   if( n_allocs || n_bytes_allocated)
      log_printf( log_file, "\nUnfreed  memory: %lu allocations %lu bytes",
         n_allocs, n_bytes_allocated);
   debug_dump_memory( );
#endif
   exit_2( 0);
}
