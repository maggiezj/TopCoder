#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <string.h>
#include <graph.h>
#include <ctype.h>
// #include "io.h"
#include "watdefs.h"
#include "alt_defs.h"
#include "findstar.h"
#include "charon.h"
#include "miscell.h"

#ifdef DEBUG_MEM
#include "checkmem.h"
#endif

#define DISPMEM ((char *)0xb8000)
#define MONOMEM ((char *)0xb0000)
#define N_ENTRIES 28
#define N_LINES_USED 14
#define N_CATALOGS 10
#define SCREENLOC( I) (screen_buffer + (I) * 80)
#define PI 3.141592653589793238462643383279502884197169399375

extern int debug_level, psf_fitting_on;
extern int gsc_mag_limit;
extern char *string_file_name, *strings[];
extern FILE *log_file;

void blast_string( char *str, int posn, int size, int attr)
{
   char *tptr = DISPMEM + posn * 2;

   while( *str && size--)
      {
      *tptr++ = *str;
      *tptr++ = (char)attr;
      if( *str)
         str++;
      }
   while( size--)
      {
      *tptr++ = ' ';
      *tptr++ = (char)attr;
      }
}

static const char *grab_from_end( const char *tptr)
{
   int i;

   for( i = 0; tptr[i] && i < 79; i++)
      ;
   while( i && (unsigned char)tptr[i] <= ' ')
      i--;         /*  back up over trailing space */
   while( i && (unsigned char)tptr[i] != ' ')
      i--;         /* next,  back up over real data */
   return( tptr + i + 1);
}


#include <process.h>

static void file_dialogue( char *filename)
{
   char *argv[20], tbuff[180];
   int n_args = 3, i, loc = 0;
   extern int language;
   FILE *ifile;

   argv[0] = "di";
   argv[1] = "-ipath.dat";
   argv[2] = "-lz";
   argv[2][2] = language;
   sprintf( tbuff, "-s%s", filename);
   for( i = 0; tbuff[i]; i++)
      if( tbuff[i] == '\\')       /* yes,  there is a path */
         loc = i;
   if( loc)
      {
      tbuff[loc] = '\0';
      argv[n_args++] = tbuff;
      }
   argv[n_args++] = NULL;
   spawnv( P_WAIT, argv[0], (const char * const *)argv);

   ifile = fopen( "path.dat", "rb");
   if( ifile)
      {
      if( fgets( tbuff, sizeof( tbuff), ifile) && *tbuff > ' ')
         {
         for( i = 0; tbuff[i] > ' '; i++)
            ;
         tbuff[i] = '\0';
         strcpy( filename, tbuff);
         }
      fclose( ifile);
      }
}

#define LOC_REPORT_FILE_NAME       0
#define LOC_CUTOFF_POINT           1
#define LOC_IMAGE_OFFSET           2
#define LOC_SCALE_TOLERANCE        3
#define LOC_SEARCH_DIST            4
#define LOC_OBSERVER_CODE          5
#define LOC_TIME_OFFSET            6
#define LOC_MAX_RESIDUAL           7
#define LOC_N_STARS_USED           8
#define LOC_COLOR_TABLE            9
#define LOC_RA_DEC_FORMAT         10
#define LOC_VIDEO_MODE            11
#define LOC_FIT_ORDER             12
#define LOC_CELL_SIZE             13
#define LOC_SATURATION_POINT      14
#define LOC_FOCAL_LENGTH          15
#define LOC_PIXEL_XSIZE           16
#define LOC_PIXEL_YSIZE           17
#define LOC_IS_INVERTED           18
#define LOC_CATALOG_USED          19
#define LOC_FILE_TIME_POINT       20
#define LOC_MAG_LIMIT             21
#define LOC_SELECT_FILE           22
#define LOC_TARGET_NAME           23
#define LOC_PSF_FITTING           24
#define LOC_MAX_TILT              25
#define LOC_PHOTO_BAND            26
#define LOC_SAVE_SETTINGS         27

const int video_modes[11] = {
             _MRES256COLOR,
             _VRES256COLOR, _SVRES256COLOR, _XRES256COLOR, _YRES256COLOR,
             _ZRES256COLOR,
             _VRES16COLOR,  _SVRES16COLOR,  _XRES16COLOR,  _YRES16COLOR,
              0 };

long setup_menu( CHARON_CONFIG *cconfig, IMAGE *img)
{
   char buff[150];
   char *screen_buffer = (char *)calloc( 80, N_ENTRIES);
   FILE *ifile = fopen( string_file_name, "rb"), *ofile;
   int i, c = 0, curr_entry = 0, help_shown = -1;
   int help_offset_line = 0, new_help_offset_line = 0, paged_help = 0;
   int top_shown = 0, video_mode_index = 0;
   long rval = 0;
   extern int n_color_tables;
   static const int catalog_remaps[N_CATALOGS] =
                    { 70, 71, -1, 73, 78, -1, 82, 124, 128, 129 };

   if( debug_level)
      log_printf( log_file,  "Entering settings menu...");
   for( i = 0; video_modes[i]; i++)
      if( cconfig->video_mode == video_modes[i])
         video_mode_index = i;

   sprintf( SCREENLOC( LOC_REPORT_FILE_NAME), strings[0], cconfig->report_file_name);
   sprintf( SCREENLOC( LOC_CUTOFF_POINT), strings[1], cconfig->cutoff_point);
   sprintf( SCREENLOC( LOC_IMAGE_OFFSET), strings[2], img->image_offset);
   sprintf( SCREENLOC( LOC_SCALE_TOLERANCE), strings[3], cconfig->scale_tolerance - 1.);
   sprintf( SCREENLOC( LOC_SEARCH_DIST), strings[4], cconfig->search_dist);
   sprintf( SCREENLOC( LOC_OBSERVER_CODE), strings[5], cconfig->obs_code);
   sprintf( SCREENLOC( LOC_TIME_OFFSET), strings[6], cconfig->time_offset);
   sprintf( SCREENLOC( LOC_MAX_RESIDUAL), strings[7], cconfig->max_residual);
   sprintf( SCREENLOC( LOC_N_STARS_USED), strings[8], cconfig->n_stars_used);
   sprintf( SCREENLOC( LOC_CELL_SIZE), strings[46], cconfig->cell_size);
   sprintf( SCREENLOC( LOC_SATURATION_POINT), strings[56], cconfig->saturation_point);
   sprintf( SCREENLOC( LOC_FOCAL_LENGTH), strings[9], img->focal_len);
   sprintf( SCREENLOC( LOC_PIXEL_XSIZE), strings[10], img->pixel_xsize);
   sprintf( SCREENLOC( LOC_PIXEL_YSIZE), strings[11], img->pixel_ysize);
   strcpy(  SCREENLOC( LOC_SAVE_SETTINGS), strings[64]);
   sprintf( SCREENLOC( LOC_MAG_LIMIT), strings[84], gsc_mag_limit / 100,
                                    abs( gsc_mag_limit % 100));
   sprintf( SCREENLOC( LOC_TARGET_NAME), strings[87], cconfig->target_name);
   sprintf( SCREENLOC( LOC_MAX_TILT), strings[89],
                                   cconfig->max_tilt * 180. / PI + .001);
   if( cconfig->tilt_angle)
      {
      char tbuff[90], *tptr = SCREENLOC( LOC_MAX_TILT);

      sprintf( tbuff, "%.2lf +/- ", cconfig->tilt_angle * 180. / PI + .001);
      tptr += strlen( tptr) - 1;
      while( *tptr != ' ')
         tptr--;
      memcpy( tptr - strlen( tbuff), tbuff, strlen( tbuff));
      }

   if( debug_level)
      log_printf( log_file,  "strings set...");

   while( c != 27)
      {
      static const int video_strings[10] =
            { 17, 18, 19, 20, 125, 127,
              75, 76, 77, 126 };
      int video_string_no = video_strings[video_mode_index];

      get_palette_name( cconfig->color_table_number, buff);
      strcpy( SCREENLOC( LOC_COLOR_TABLE    ), buff);
      strcpy( SCREENLOC( LOC_RA_DEC_FORMAT  ), strings[14 + cconfig->ra_dec_format]);
      strcpy( SCREENLOC( LOC_VIDEO_MODE     ), strings[video_string_no]);
      strcpy( SCREENLOC( LOC_FIT_ORDER      ), strings[47 + cconfig->fit_order - 1]);
      strcpy( SCREENLOC( LOC_IS_INVERTED    ), strings[62 + img->is_inverted]);
      strcpy( SCREENLOC( LOC_CATALOG_USED   ), strings[catalog_remaps[cconfig->catalog_used]]);
      strcpy( SCREENLOC( LOC_FILE_TIME_POINT), strings[67 + cconfig->file_time_point]);
      i = strlen( cconfig->path_name) - 78;
      if( i < 0)
         i = 0;
      strcpy( SCREENLOC( LOC_SELECT_FILE), cconfig->path_name + i);
      strcpy( SCREENLOC( LOC_PSF_FITTING), strings[psf_fitting_on ? 85 : 86]);
      sprintf( SCREENLOC( LOC_PHOTO_BAND), strings[90],
                                               cconfig->photometric_band);

      if( debug_level)
         log_printf( log_file,  "\nReady to put on screen...\n");
                     /* Make sure the range shown is within bounds: */
      if( top_shown > curr_entry)
         top_shown = curr_entry;
      if( top_shown < curr_entry - N_LINES_USED + 1)
         top_shown = curr_entry - N_LINES_USED + 1;
      if( top_shown < 0)
         top_shown = 0;
      if( top_shown > N_ENTRIES - N_LINES_USED)
         top_shown = N_ENTRIES - N_LINES_USED;
      for( i = 0; i < N_LINES_USED; i++)
         {
         int j;
         char *tptr = SCREENLOC( i + top_shown);

         for( j = strlen( tptr); j < 80; j++)
            tptr[j] = '\0';
         blast_string( SCREENLOC( i + top_shown), i * 80, 80,
                       (i + top_shown == curr_entry) ? 31 : 15);
         }

      if( help_shown != curr_entry || new_help_offset_line != help_offset_line)
         {
         int line, search_char = 'A' + curr_entry;

         if( curr_entry > 25)
            search_char = 'a' + (curr_entry - 26);
         if( help_shown != curr_entry)    /* new subject means start at top: */
             new_help_offset_line = 0;
         help_shown = curr_entry;
         help_offset_line = new_help_offset_line;
         blast_string( "", N_LINES_USED * 80, (25 - N_LINES_USED) * 80, 15);
         fseek( ifile, 0L, SEEK_SET);
         line = (N_LINES_USED + 1) - new_help_offset_line;
         paged_help = 0;
         while( fgets( buff, sizeof( buff), ifile))
            if( *buff == '~' && buff[1] == (char)search_char)
               while( fgets( buff, sizeof( buff), ifile) && *buff != '~')
                  {
                  if( line < 25 && line > N_LINES_USED)
                     {
                     buff[strlen( buff) - 2] = '\0';
                     buff[80] = '\0';
                     blast_string( buff, line * 80, 80, 15);
                     }
                  if( line == 25)
                     {
                     blast_string( strings[83], 24 * 80, 80, 63);
                     paged_help = 1;
                     }
                  line++;
                  }
         }

      c = getch( );
      if( !c)
         c = getch( ) + 256;
      switch( c)
         {
         case ZKEY_NUM_4: case ZKEY_NUM_6:
            if( curr_entry >= N_LINES_USED)
               curr_entry -= N_LINES_USED;
            else if( curr_entry + N_LINES_USED < N_ENTRIES)
               curr_entry += N_LINES_USED;
            break;
         case ZKEY_NUM_8:
            curr_entry--;
            if( curr_entry == -1)
               curr_entry += N_ENTRIES;
            break;
         case ZKEY_NUM_2:
            curr_entry++;
            if( curr_entry == N_ENTRIES)
               curr_entry = 0;
            break;
         case ZKEY_NUM_3:
            if( paged_help)
               new_help_offset_line++;
            break;
         case ZKEY_NUM_9:
            if( new_help_offset_line > 0)
               new_help_offset_line--;
            break;
         case ZKEY_CTRL_TAB:
            memcpy( MONOMEM, DISPMEM, 4000);
            break;
         case 8:           /* backspace */
            {
            char *tptr = SCREENLOC( curr_entry) + 20;

            if( curr_entry == LOC_TARGET_NAME)
               while( tptr[-1] != ':')
                  tptr--;
            if( *tptr)
               {
               memmove( tptr + 1, tptr, strlen( tptr) - 1);
               *tptr = ' ';
               }
            }
            break;
         case 13:
            if( curr_entry == LOC_COLOR_TABLE)
               cconfig->color_table_number = (cconfig->color_table_number + 1)
                           % n_color_tables;
            else if( curr_entry == LOC_RA_DEC_FORMAT)
               cconfig->ra_dec_format = (cconfig->ra_dec_format + 1) % 3;
            else if( curr_entry == LOC_VIDEO_MODE)
               {
               video_mode_index++;
               if( !video_modes[video_mode_index])
                  video_mode_index = 0;
               }
            else if( curr_entry == LOC_FIT_ORDER)
               cconfig->fit_order = (cconfig->fit_order % 3) + 1;     /* 1, 2, or 3 */
            else if( curr_entry == LOC_IS_INVERTED)
               img->is_inverted ^= 1;
            else if( curr_entry == LOC_CATALOG_USED)
               {
               int i;
               static const char catalog_loop[9] = {
                           CATALOG_TYCHO2,
                           CATALOG_A1,
                           CATALOG_A2,
                           CATALOG_GSC_11,
                           CATALOG_GSC_ACT,
                           CATALOG_UCAC2,
                           CATALOG_UCAC3,
                           CATALOG_UCAC4,
                           CATALOG_TYCHO2 };

               for( i = 0; catalog_loop[i] != cconfig->catalog_used; i++)
                  ;
               cconfig->catalog_used = catalog_loop[i + 1];
               }
            else if( curr_entry == LOC_FILE_TIME_POINT)
               cconfig->file_time_point =
                          (cconfig->file_time_point + 1) % 3;
            else if( curr_entry == LOC_SELECT_FILE)
               {
               file_dialogue( cconfig->path_name);
               help_shown = -1;
               }
            else if( curr_entry == LOC_PSF_FITTING)
               psf_fitting_on ^= 1;
            else if( curr_entry == LOC_SAVE_SETTINGS)
               c = 27;
            rval ^= 1L << curr_entry;
            break;
         default:
            if( c >= ' ' && c <= 'z')
               if( curr_entry == LOC_PHOTO_BAND)
                  cconfig->photometric_band = toupper( c);
               else
                  {
                  char *tptr = SCREENLOC( curr_entry) + 20;
                  int len;

                  if( curr_entry == LOC_TARGET_NAME)
                     while( tptr[-1] != ':')
                        tptr--;
                  len = strlen( tptr) - 1;
                  if(( rval >> curr_entry) & 1L)
                     memmove( tptr, tptr + 1, len + 1);
                  else
                     memset( tptr, ' ', len + 1);
                  tptr[len] = (char)c;
                  rval |= 1L << curr_entry;
                  }
            break;
         }
      }

   if( debug_level)
      log_printf( log_file,  "Now out of setting menu:");

   cconfig->video_mode = video_modes[video_mode_index];
   if( debug_level)
      log_printf( log_file,  "mode %d, ", cconfig->video_mode);
   strcpy( cconfig->report_file_name,
                         grab_from_end( SCREENLOC( LOC_REPORT_FILE_NAME)));
   if( debug_level)
      log_printf( log_file,  "report file %s\n", cconfig->report_file_name);
   cconfig->cutoff_point = atof( grab_from_end( SCREENLOC( LOC_CUTOFF_POINT)));
   if( debug_level)
      log_printf( log_file,  "cutoff %d, ", cconfig->cutoff_point);
   img->image_offset = atoi( grab_from_end( SCREENLOC( LOC_IMAGE_OFFSET)));
   cconfig->scale_tolerance = atof( grab_from_end( SCREENLOC( LOC_SCALE_TOLERANCE))) + 1.;
   cconfig->search_dist = atof( grab_from_end( SCREENLOC( LOC_SEARCH_DIST)));
   strcpy( cconfig->obs_code,
                        grab_from_end( SCREENLOC( LOC_OBSERVER_CODE)));
   if( debug_level)
      log_printf( log_file,  "obs code %s, ", cconfig->obs_code);
   cconfig->time_offset = atof( grab_from_end( SCREENLOC( LOC_TIME_OFFSET)));
   memcpy( buff, SCREENLOC( LOC_MAX_TILT), 80);
   buff[80] = '\0';
   for( i = 0; buff[i] != ':'; i++)
      ;
   if( strstr( buff, "+/-"))
      sscanf( buff + i + 1, "%lf +/- %lf", &cconfig->tilt_angle,
                                         &cconfig->max_tilt);
   else
      {
      cconfig->tilt_angle = 0.;
      cconfig->max_tilt = atof( buff + i + 1);
      }
   cconfig->max_tilt *= PI / 180.;
   cconfig->tilt_angle *= PI / 180.;
   if( debug_level)
      log_printf( log_file,  "tilt %lf (%lf) ", cconfig->max_tilt, cconfig->tilt_angle);
   cconfig->max_residual = atof( grab_from_end( SCREENLOC( LOC_MAX_RESIDUAL)));
   cconfig->n_stars_used = atoi( grab_from_end( SCREENLOC( LOC_N_STARS_USED)));
   cconfig->cell_size = atoi( grab_from_end( SCREENLOC( LOC_CELL_SIZE)));
   cconfig->saturation_point = atoi( grab_from_end( SCREENLOC( LOC_SATURATION_POINT)));
   if( debug_level)
      log_printf( log_file,  "sat point %d, ", cconfig->saturation_point);
   img->focal_len =    atof( grab_from_end( SCREENLOC( LOC_FOCAL_LENGTH)));
   img->pixel_xsize =  atof( grab_from_end( SCREENLOC( LOC_PIXEL_XSIZE)));
   img->pixel_ysize =  atof( grab_from_end( SCREENLOC( LOC_PIXEL_YSIZE)));
   gsc_mag_limit = (int)
                 ( atof( grab_from_end( SCREENLOC( LOC_MAG_LIMIT))) * 100.);
   strcpy( buff, SCREENLOC( LOC_TARGET_NAME));
   for( i = 0; buff[i] && buff[i] != ':'; i++)
      ;
   for( i++; buff[i] == ' '; i++)
      ;
   if( buff[i])
      strcpy( cconfig->target_name, buff + i);
   if( debug_level)
      log_printf( log_file,  "target %s\n", cconfig->target_name);

   fclose( ifile);

   for( i = 0; i < 80 * N_ENTRIES; i++)
      if( (unsigned char)screen_buffer[i] < ' ')
         screen_buffer[i] = ' ';

   ofile = fopen( "settings.dat", "wb");
   if( debug_level)
      log_printf( log_file,  "Writing settings... ");
   for( i = 0; i < N_ENTRIES; i++)
      {
      fwrite( SCREENLOC( i), 80, 1, ofile);
      fprintf( ofile, "\n");
      }
   fclose( ofile);

   free( screen_buffer);
   blast_string( "", 0, 25 * 80, 15);        /* a.k.a "clear screen" */
   if( debug_level)
      log_printf( log_file,  "Returning from setup menu: %d\n", rval);
   return( rval);
}

void show_help_data( void)
{
   FILE *ifile = fopen( string_file_name, "rb");
   char tbuff[100];
   int i;

   while( fgets( tbuff, sizeof( tbuff), ifile))
      if( tbuff[0] == '~' && tbuff[1] == '?')
         {
         tbuff[0] = ' ';
         for( i = 0; i < 25 && tbuff[0] != '~'; i++)
            if( fgets( tbuff, sizeof( tbuff), ifile))
               {
               tbuff[strlen( tbuff) - 2] = '\0';
               tbuff[80] = '\0';
               blast_string( tbuff, i * 80, 80, 15);
               }
         }
   fclose( ifile);
}
