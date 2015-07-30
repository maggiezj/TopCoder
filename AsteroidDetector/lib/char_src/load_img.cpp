#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <conio.h>
#include "watdefs.h"
#include "findstar.h"
#include "date.h"

#ifdef DEBUG_MEM
#include "checkmem.h"
#endif

int log_printf( FILE *log_file, const char *format, ...);   /* charon.cpp */
void exit_2( const int exit_code);                    /* charon.cpp */

extern FILE *log_file;
extern int debug_level;

#define PI 3.141592653589793238462643383279502884
#define ST4_WIDTH  192
#define ST4_HEIGHT 165
#define ST5_WIDTH  323
#define ST5_HEIGHT 242
#define ST6_WIDTH  375
#define ST6_HEIGHT 242
#define ST7_WIDTH  765
#define ST7_HEIGHT 510
#define ST7_PIXEL_SIZE_X       9.
#define ST7_PIXEL_SIZE_Y       9.
#define ST6_PIXEL_SIZE_X      23.
#define ST6_PIXEL_SIZE_Y      27.
#define ST5_PIXEL_SIZE_X      10.
#define ST5_PIXEL_SIZE_Y      10.
#define ST4_PIXEL_SIZE_X      13.75
#define ST4_PIXEL_SIZE_Y      16.
#define ST9_PIXEL_SIZE_X      20.
#define ST9_PIXEL_SIZE_Y      20.
#define ST9_WIDTH            512
#define ST9_HEIGHT           512
#define PIC_OFFSET 290
#define AST_PIXEL_SIZE_X       9.
#define AST_PIXEL_SIZE_Y       9.
#define STARLIGHT_PIXEL_WIDTH    12.7
#define STARLIGHT_PIXEL_HEIGHT   16.6
// #define STARLIGHT_PIXEL_WIDTH    9.8
// #define STARLIGHT_PIXEL_HEIGHT   12.6
#define STARLIGHT_HX516_SIZE      7.4
#define CB_FILESIZE 63007L

/* find_histo( ) is used to determine the pixel value below which a certain */
/* fraction of pixels in an image lie.  It's used,  for example,  to set    */
/* the initial contrast for an image;  60% of the pixels are supposed to be */
/* black,  1% white.  find_histo( ) determines where those points lie.      */

#ifdef PIXEL_32
#define N_HISTO_PASSES 3
#else
#define N_HISTO_PASSES 2
#endif

PIXEL find_histo( const IMAGE *img, const double fraction)
{
   extern int ignored_north, ignored_south, ignored_east, ignored_west;
   long i, j;
   int pass;
   const int real_xsize = img->xsize - (ignored_east + ignored_west);
   const int real_ysize = img->ysize - (ignored_south + ignored_north);
   const long target = (long)
                  ( fraction * (double)real_xsize * (double)real_ysize + .5);
   long *histo = (long *)calloc( 256, sizeof( long));
   PIXEL level1, level2, rval = 0;

   if( debug_level)
      log_printf( log_file, "find_histo: ");
   for( pass = N_HISTO_PASSES; pass; pass--)
      {
      int bitshift = pass * 8;
      long base_count = 0L;

      if( debug_level)
         log_printf( log_file, " %d", pass);
      level1 = (rval >> bitshift) << bitshift;
      level2 = level1 + ((PIXEL)1 << bitshift) - (PIXEL)1;
      bitshift -= 8;
      memset( histo, 0, 256 * sizeof( long));
      for( j = ignored_north; j < img->ysize - ignored_south; j++)
         {
         const PIXEL *tptr = img->img + ignored_west + img->xsize * j;

         for( i = real_xsize; i; i--, tptr++)
            if( *tptr < level1)
               base_count++;
            else if( *tptr <= level2)
               histo[(*tptr >> bitshift) & 0xff]++;
         }
      for( j = 0; base_count + histo[j] < target; j++)
         base_count += histo[j];
      rval |= ((PIXEL)j << bitshift);
      }
   free( histo);
   if( debug_level)
      log_printf( log_file, " done\n");
   return( rval);
}

/* When blinking images,  you ideally want the 'blink' images to be loaded */
/* with a contrast that matches the original image.  The following function */
/* accomplishes this by looking at the 'low_end' and 'high_end' of the      */
/* source_img,  finding out what fraction of the pixels in source_img are   */
/* under those points,  and then making sure that the same fractions apply  */
/* to 'img'.                                                                */

int match_contrasts( IMAGE *img, const IMAGE *source_img)
{
   int i;

   for( i = 0; i < 2; i++)       /* first low,  then high end */
      {
      long fraction = 0L, j;
      long total_pixels = (long)source_img->xsize * (long)source_img->ysize;
      PIXEL reference = (i ? source_img->high_end : source_img->low_end);

      for( j = 0; j < total_pixels; j++)
         if( source_img->img[j] < reference)
            fraction++;
      reference = find_histo( img,
                     (double)fraction / (double)total_pixels);
      if( i)
         img->high_end = reference;
      else
         img->low_end = reference;
      }
   return( 0);
}

int colour_starlight_express = 0;

int load_image( const char *filename, IMAGE *img)
{
   FILE *ifile = fopen( filename, "rb");
   int month, year = -1, hour = 0, minute = 0;
   int is_compressed = 0, i, j, signed_ints = 0;
   double day, seconds = 0.;
   IMAGE rval;
   const char *extension = strchr( filename, '.');
   char ibuff[180];
   long filesize;
   int bits_per_pixel = 0;
   int real_line_size = 0;
   extern char *strings[];

   if( debug_level)
      log_printf( log_file, "Loading image %s\n", filename);
   if( !ifile || !extension)
      return( -1);
   fseek( ifile, 0L, SEEK_END);
   filesize = ftell( ifile);
   fseek( ifile, 0L, SEEK_SET);
   memcpy( &rval, img, sizeof( IMAGE));
   rval.image_type = tolower( filename[ strlen( filename) - 1]);
   rval.byte_order = 0;
   rval.ra_dec_from_header = 0;
   fread( ibuff, 80, 1, ifile);

   if( !memcmp( ibuff, "SIMPLE  =", 9))   /* a FITS file w/o correct exten */
      {
      log_printf( log_file, strings[107], "FITS");
                          /* "Image identified as %s file\n" */
      extension = ".fit";
      rval.image_type = 't';
      rval.byte_order = 1;
      }
   if( !memcmp( ibuff, "MIRA", 4))
      {
      extension = ".mim";
      rval.image_type = 'm';
      }

                        /* It so happens that the PixCel 255 is basically */
                        /* a modified SBIG one,  so we can do this:       */
   if( !memcmp( ibuff, "PixCel255 Compressed Image", 26))
      strcpy( ibuff, "ST-6 Compressed Image");

   if( !memcmp( ibuff, "ST-", 3) && !memcmp( ibuff + 4, " Compressed Image",17))
      {
      char *ext_buff = ".stx";

      extension = ext_buff;
      rval.image_type = ext_buff[3] = ibuff[3];
      log_printf( log_file, strings[107], extension);
                                  /* "Image identified as %s file\n" */
      }

                           /* now check for B E-S format: */
   fseek( ifile, 291000L, SEEK_SET);
   if( fgets( ibuff, 80, ifile))
      if( !memcmp( ibuff, "Exposure length  :", 18))
         {
         extension = ".sx";
         rval.image_type = 'x';
         rval.xsize = 500;
         rval.ysize = 291;
         rval.pixel_xsize = STARLIGHT_PIXEL_WIDTH;
         rval.pixel_ysize = STARLIGHT_PIXEL_HEIGHT;
         rval.exposure_length = atof( ibuff + 18);
         fgets( ibuff, 80, ifile);
         sscanf( ibuff + 18, "%d/%d/%lf %d:%d:%lf",
               &year, &month, &day, &hour, &minute, &seconds);
         fgets( ibuff, 80, ifile);     /* image subject */
         fgets( ibuff, 80, ifile);
         rval.focal_len = atof( ibuff + 18) / (25.4 * 1000.);
         rval.is_inverted = 1;
         fgets( ibuff, 80, ifile);        /* rotation... skip it */
         fgets( ibuff, 80, ifile);
         rval.xform[0] = atof( ibuff + 18) * 180. / PI;
         fgets( ibuff, 80, ifile);
         rval.xform[1] = atof( ibuff + 18) * 180. / PI;
         rval.ra_dec_from_header = 1;
         }
   fseek( ifile, 0L, SEEK_SET);

   if( filesize == 295800 || filesize == 261120 || filesize == 660000L)
      {                                /* assume Starlight XPress */
      extension = ".sx";
      rval.image_type = 'x';
      if( filesize == 660000L)            /* HX516 */
         {
         rval.xsize = 660;
         rval.ysize = 500;
         rval.pixel_xsize = rval.pixel_ysize = STARLIGHT_HX516_SIZE;
         }
      else                                /* MX5 or SX */
         {
         rval.xsize = 510;
                 /* It looks as if,  realistically,  only about 248 lines */
                 /* contain clean data: */
         rval.ysize = 248;
//       rval.ysize = 256;
                /* In a "290-line" image,  the first 34 lines are garbage: */
         if( filesize == 295800L)
            fseek( ifile, (290L - 256L + 8L) * 510L * 2L, SEEK_SET);
                        /* Then skip eight lines: */
         else
            fseek( ifile, 8L * 510L * 2L, SEEK_SET);

         rval.pixel_xsize = STARLIGHT_PIXEL_WIDTH;
         rval.pixel_ysize = STARLIGHT_PIXEL_HEIGHT;
         }
      log_printf( log_file, strings[107], "Starlight XPress SX");
      }

   if( !memicmp( extension, ".st", 3))
      switch( extension[3])
         {
         case '4':
            rval.ysize = ST4_HEIGHT;
            rval.xsize = ST4_WIDTH;
            rval.pixel_xsize = ST4_PIXEL_SIZE_X;
            rval.pixel_ysize = ST4_PIXEL_SIZE_Y;
            break;
         case '5':
         case '6':
         case '7':
         case '8':
         case '9':
            {
            char *tptr, *tbuff;

            if( extension[3] == '4')
               {
               rval.ysize = ST4_HEIGHT;
               rval.xsize  = ST4_WIDTH;
               rval.pixel_xsize = ST4_PIXEL_SIZE_X;
               rval.pixel_ysize = ST4_PIXEL_SIZE_Y;
               }
            if( extension[3] == '5')
               {
               rval.ysize = ST5_HEIGHT;
               rval.xsize  = ST5_WIDTH;
               rval.pixel_xsize = ST5_PIXEL_SIZE_X;
               rval.pixel_ysize = ST5_PIXEL_SIZE_Y;
               }
            if( extension[3] == '6')
               {
               rval.ysize = ST6_HEIGHT;
               rval.xsize  = ST6_WIDTH;
               rval.pixel_xsize = ST6_PIXEL_SIZE_X;
               rval.pixel_ysize = ST6_PIXEL_SIZE_Y;
               }
            if( extension[3] == '7' || extension[3] == '8')
               {
               rval.ysize = ST7_HEIGHT;
               rval.xsize  = ST7_WIDTH;
               rval.pixel_xsize = ST7_PIXEL_SIZE_X;
               rval.pixel_ysize = ST7_PIXEL_SIZE_Y;
               }
            if( extension[3] == '9')
               {
               rval.ysize = ST9_HEIGHT;
               rval.xsize  = ST9_WIDTH;
               rval.pixel_xsize = ST9_PIXEL_SIZE_X;
               rval.pixel_ysize = ST9_PIXEL_SIZE_Y;
               }

            tbuff = (char *)calloc( 2049, 1);
            fread( tbuff, 2048, 1, ifile);
            if( tptr = strstr( tbuff, "Focal_length = "))
               rval.focal_len = atof( tptr + 15);
            if( tptr = strstr( tbuff, "Focal_Length = "))
               rval.focal_len = atof( tptr + 15);
            if( tptr = strstr( tbuff, "Height = "))
               rval.ysize = atoi( tptr + 9);
            if( tptr = strstr( tbuff, "Width = "))
               rval.xsize = atoi( tptr + 8);
            if( tptr = strstr( tbuff, "Date = "))
               sscanf( tptr + 7, "%d/%lf/%d", &month, &day, &year);
            if( tptr = strstr( tbuff, "Time = "))
               sscanf( tptr + 7, "%d:%d:%lf", &hour, &minute, &seconds);
            if( tptr = strstr( tbuff, "X_pixel_size = "))
               rval.pixel_xsize = 1000. * atof( tptr + 15);
            if( tptr = strstr( tbuff, "Y_pixel_size = "))
               rval.pixel_ysize = 1000. * atof( tptr + 15);
            if( tptr = strstr( tbuff, "Note = Sekunden je Reihe"))
               {
               extern int is_drift_scan;

               is_drift_scan = 1;
               log_printf( log_file, "Drift-scanned image\n");
               }

            rval.exposure_length = .01;   /* was .001; fixed 17 Dec 97 */
            if( tptr = strstr( tbuff, "Number_exposures = "))
               rval.exposure_length *= atof( tptr + 18);
            if( tptr = strstr( tbuff, "Each_exposure = "))
               rval.exposure_length *= atof( tptr + 15);
            is_compressed = (tbuff[5] == 'C');
            if( !memcmp( tbuff, "PixCel255 Compressed Image", 26))
               is_compressed = 1;
            if( tbuff[3] == '4')
               {
               rval.ysize = ST4_HEIGHT - 1;           /* dunno why,  though! */
               is_compressed = (tbuff[6] == 'C');      /* allow for the 'X' */
               }
            free( tbuff);
            }
            break;
         default:
            fclose( ifile);
            return( -2);
            break;
         }
   else if( !stricmp( extension, ".pic"))
      {
      short size[2];

      fseek( ifile, 68L, SEEK_SET);
      fread( size, 2, sizeof( short), ifile);
      rval.xsize = size[0];
      rval.ysize = size[1];
      }
   else if( !stricmp( extension, ".ast"))    /* Mike Kretlow file */
      {
      long size[2];
      float jd;

      rval.image_type = 'a';        /* 't' is taken by .LST */
      fread( size, 2, sizeof( long), ifile);
      rval.xsize = (int)size[0];
      rval.ysize = (int)size[1];
      rval.pixel_xsize = AST_PIXEL_SIZE_X;
      rval.pixel_ysize = AST_PIXEL_SIZE_Y;
      fseek( ifile, 36L, SEEK_SET);
      fread( &jd, 1, sizeof( float), ifile);
      rval.time = (double)jd;
      fseek( ifile, 512L, SEEK_SET);      /* data after the 512-byte hdr */
      }
   else if( !stricmp( extension, ".fit") || !stricmp( extension, ".mim"))
      {
#ifdef MORE_TROUBLE_THAN_ITS_WORTH
      double tval;
#endif
      int pixel_size_from_cdelt = 0, time_found = 0;
      int origin_flag = 0;

      rval.image_type = extension[1];
      for( i = 0; fread( ibuff, 80, 1, ifile) &&
                                           memcmp( ibuff, "END ", 4); i++)
         {
         j = atoi( ibuff + 11);
         if( !memcmp( ibuff, "NAXIS", 5))
            {
            if( ibuff[5] == '1')
               rval.xsize = j;
            else if( ibuff[5] == '2')
               rval.ysize = j;
                        /* Cookbook 252 pixels are 25.5x19.7 microns */
            if( rval.xsize == 252 && rval.ysize == 242)
               {
               rval.pixel_xsize = 25.5;
               rval.pixel_ysize = 19.7;
               }
            if( rval.xsize == 242 && rval.ysize == 252)
               {
               rval.pixel_xsize = 19.7;
               rval.pixel_ysize = 25.5;
               }
            }

#ifdef MORE_TROUBLE_THAN_ITS_WORTH
         if( !memcmp( ibuff, "CDELT", 5) || !memcmp( ibuff, "CRDELT", 6))
            {
            int axis;

            tval = fabs( atof( ibuff + 11));
            if( ibuff[1] == 'R')
               {
               pixel_size_from_cdelt |= 1;
               axis = ibuff[6];
               }
            else
               {
               pixel_size_from_cdelt |= 2;
               axis = ibuff[5];
               tval *= PI / 180.;     /* CDELT=in degrees, cvt to radians */
               }
            tval *= 25.4 * 1000. * rval.focal_len;     /* in microns */
            if( axis == '1')
               rval.pixel_xsize = tval;
            else if( axis == '2')
               rval.pixel_ysize = tval;
            }
#endif

         if( !memcmp( ibuff, "BITPIX ", 7))
            bits_per_pixel = j;

         if( !memcmp( ibuff, "DATASEC ", 8))
            {
            int tval[4];

            if( sscanf( ibuff + 12, "%d:%d,%d:%d", tval, tval + 1,
                                                tval + 2, tval + 3) == 4)
               {
               extern int ignored_north, ignored_south;
               extern int ignored_east, ignored_west;

               ignored_east = tval[0] + 2;
               ignored_west = rval.xsize - tval[1] + 2;
               ignored_north = tval[2] + 2;
               ignored_south = rval.ysize - tval[3] + 2;
               }
            }

                     /* Modified 10 Aug 98,  to account for a weird date */
                     /* format change.  Some are '10/08/1998'; others are*/
                     /* '1998-08-10' */
         if( !memcmp( ibuff, "DATE-OBS", 8))
            {
                     /* Some SOHO images have 'yyyy/mm/dd': */
            if( ibuff[15] == '/' && ibuff[18] == '/' &&
                  sscanf( ibuff + 11, "%d/%d/%lf", &year, &month, &day) == 3
                  && year > 1990 && year < 2100 && month >= 1 && month <= 12)
               ;
            else
               if( sscanf( ibuff + 11, "%lf/%d/%d", &day, &month, &year) != 3)
                  if( sscanf( ibuff + 11, "%d-%d-%lf", &year, &month, &day) != 3)
                     day = 0.;
            if( day)
               if( day > 31 || month > 12)
                  {
 log_printf( log_file,  "WARNING:  the DATE-OBS field in the FITS header is not in a\n");
 log_printf( log_file,  "correct format.  It should be either 'DD/MM/YYYY' or 'YY-MM-DD'.\n");
 log_printf( log_file,  "(In either case,  Charon will accept a two or four-digit year;\n");
 log_printf( log_file,  "but not all programs will be so forgiving.)\n");
 log_printf( log_file,  "Hit a key: \n");
                  getch( );
                  day = 0.;
                  }
                              /* Vince Gardiner (Pictor 416) case: */
                              /* DATE-OBS= '1999-06-25T19:38:45'   */
                              /* Added 16 July 1999                */
            if( ibuff[15] == '-' && ibuff[18] == '-' && ibuff[21] == 'T'
                     && ibuff[24] == ':' && ibuff[27] == ':')
                if( sscanf( ibuff + 22, "%d:%d:%lf",
                          &hour, &minute, &seconds) == 3)
                   time_found = 1;
            }

         if( !memcmp( ibuff, "TIME-OBS", 8) || !memcmp( ibuff, "UT      ", 8))
            {
            int n_found = sscanf( ibuff + 11, "%d:%d:%lf",
                                                 &hour, &minute, &seconds);

            if( origin_flag == 1 && n_found != 3)      /* R Berry odd time format */
               n_found = sscanf( ibuff + 11, "%d/%d/%lf",
                                                 &hour, &minute, &seconds);
            if( n_found != 3)
               {
log_printf( log_file, "ERROR:  the observation time is in an invalid format.  This _may_\n");
log_printf( log_file, "keep Charon from processing the image correctly.\n");
log_printf( log_file, "Hit any key: \n");
               getch( );
               }
            time_found = 1;
            }

                        /* R A Kowalski,  for one,  uses "TIME    " */
                        /* Nicolas Minary uses "UT-START"           */
         if( !time_found)
            if( !memcmp( ibuff, "TIME    ",8) || !memcmp( ibuff, "UT-START",8))
               sscanf( ibuff + 11, "%d:%d:%lf", &hour, &minute, &seconds);


         if( !memcmp( ibuff, "DATETIME", 8))
            {
            int m1;

            for( m1 = 0; m1 < 12 && memicmp( set_month_name( m1 + 1, NULL),
                                             ibuff + 14, 3);
                        m1++)
               ;
            if( m1 < 12 && ibuff[13] == ' ' && ibuff[17] == ' ')
               {
               month = m1 + 1;
               sscanf( ibuff + 18, "%lf %d:%d:%lf %d", &day,
                         &hour, &minute, &seconds, &year);
               }
            }

         if( !memcmp( ibuff, "EXPTIME ", 8))
            {
            char *tptr = ibuff + 9;

            if( *tptr == '\'')
               tptr++;
            rval.exposure_length = atof( tptr);
            }
         if( !memcmp( ibuff, "EXPOSURE", 8))
            {
            char tbuff[80];

            if( sscanf( ibuff + 11, "%lf%s", &rval.exposure_length, tbuff) == 2)
               if( !strcmp( tbuff, "ms"))
                  {
                  log_printf( log_file,
                           "Non-standard exposure length of %.3lf ms\n",
                           rval.exposure_length);
                  rval.exposure_length /= 1000.;
                  }
            }

                              /* Charles Juels sent in an image where */
                              /* the target RA/dec was in the header: */
         if( !memcmp( ibuff, "OBJCTRA = '", 11))
            {
            double hr, min, sec;
            char end_apostrophe;

            if( sscanf( ibuff + 11, "%lf %lf %lf%c", &hr, &min, &sec,
                    &end_apostrophe) == 4 && end_apostrophe == '\'')
               {
               rval.xform[0] = (hr + min / 60. + sec / 3600.) * 15.;
               fread( ibuff, 80, 1, ifile);
               sscanf( ibuff + 12, "%lf %lf %lf", &hr, &min, &sec);
               rval.xform[1] = (hr + min / 60. + sec / 3600.);
               if( ibuff[11] == '-')
                  rval.xform[1] *= -1.;
               rval.ra_dec_from_header = 1;
               }
            }
#ifdef CODE_REMOVED_TOO_MUCH_TROUBLE
                               /* Paul Boltwood,  and a few others,  have   */
                               /* images with the target RA/dec also in the */
                               /* header,  except as 'CRVAL1' & '2'         */
                  /* Then I removed this 'cuz some images,  notably ones    */
                  /* processed w/PinPoint,  have CRVAL1=RA,  not dec;  and  */
                  /* do _not_ have CRVAL1 and CRVAL2 right next to each     */
                  /* other.                                                 */
         if( !memcmp( ibuff, "CRVAL1  = ", 10))
            {
            rval.xform[1] = atof( ibuff + 10);     /* dec comes first (?) */
            fread( ibuff, 80, 1, ifile);
            rval.xform[0] = atof( ibuff + 10);
            rval.ra_dec_from_header = 1;
            }
#endif

         if( !memcmp( ibuff, "BZERO   =", 9))
            {
            double bzero = atof( ibuff + 10);

            if( bits_per_pixel == 8 &&
                               bzero > 127.9 && bzero < 128.1)
               signed_ints = 1;
            else if( bits_per_pixel == 16 &&
                               bzero > 32767.9 && bzero < 32768.1)
               signed_ints = 1;
            else if( bits_per_pixel == 32 &&
                               bzero > 2147483647.9 && bzero < 2147483648.1)
               signed_ints = 1;
            else if( bzero > 0.1 || bzero < -0.1)
               {
log_printf( log_file,
          "WARNING:  This image has a non-standard value for BZERO.  This\n");
log_printf( log_file,
          "will _probably_ not be a problem,  but the image _may_ be shown\n");
log_printf( log_file,
          "or processed incorrectly.  Contact Project Pluto if that happens.\n");
log_printf( log_file,
          "Hit any key:\n");
               getch( );
               }
            }
         if( !memcmp( ibuff, "ORIGIN  = 'AIP245 1.00 by Richard Berry'", 40))
            {
            log_printf( log_file,  "Recognized AIP245 1.00,  by Richard Berry;\n");
            log_printf( log_file,  "assuming non-standard observation time format\n");
            origin_flag = 1;
            }
         }
      if( extension[1] == 'f')         /* FITS file */
         fseek( ifile, (long)(i / 36 + 1) * 36L * 80L, SEEK_SET);
      else                       /* MIRA:  file is in subdir */
         {
         char img_filename[200];
         int name_loc = 0;

         fclose( ifile);
         for( i = 0; filename[i]; i++)
            if( filename[i] == ':' || filename[i] == '\\')
               name_loc = i;
         memcpy( img_filename, filename, name_loc);
         strcpy( img_filename + name_loc, "\\pix");
         strcat( img_filename, filename + name_loc);
         ifile = fopen( img_filename, "rb");
         log_printf( log_file, "Opening '%s'\n", img_filename);
         if( !ifile)
            {
            log_printf( log_file, "Couldn't open image file '%s'\n",
                                            img_filename);
            exit_2( -23);
            }
         }

      if( pixel_size_from_cdelt == 3)
         {
         rval.focal_len = 9. / rval.pixel_xsize;
         rval.pixel_xsize *= rval.focal_len;
         rval.pixel_ysize *= rval.focal_len;
         }
      }
   else if( !stricmp( extension, ".lst"))
      {
      rval.is_xy_list = 1;
//    rval.focal_len = .21115;
      }
   else if( !stricmp( extension, ".map"))
      {
      rval.is_xy_list = 1;
      rval.focal_len = .21115 / 10.637;
      }
   else if( !stricmp( extension, ".fin"))
      {
      rval.is_xy_list = 1;
      }
   else if( !stricmp( extension, ".sx"))
      {
      }
   else if( !stricmp( extension, ".bmp"))
      {

      fseek( ifile, 0L, SEEK_SET);
      fread( ibuff, 80, 1, ifile);
      rval.data_offset = *(long *)( ibuff + 10);
      fseek( ifile, rval.data_offset, SEEK_SET);
      bits_per_pixel = ibuff[28];
      rval.xsize = *(long *)( ibuff + 18);
      rval.ysize = *(long *)( ibuff + 22);
               /* BMPs are weird:  the line size is always rounded up to */
               /* the next multiple of four: */
      real_line_size = rval.xsize * bits_per_pixel / 8;
      real_line_size = (real_line_size + 3) & (~3);
      rval.image_type = 'b';
      }
   else if( filesize == CB_FILESIZE &&
                             !stricmp( extension, ".pa"))
      {                        /* Cookbook cameras are two-part */
      rval.xsize = 252;
      rval.ysize = 242;
      rval.pixel_xsize = 25.5;
      rval.pixel_ysize = 19.7;
      rval.image_type = 'k';        /* 'Kookbook' */
      }
   else if( filesize == CB_FILESIZE &&
                             !stricmp( extension, ".p1"))
      {                        /* "Three-parter" */
      rval.xsize = 378;
      rval.ysize = 242;
      rval.pixel_xsize = 25.5 * 2. / 3.;  /* not _totally_ sure of sizes! */
      rval.pixel_ysize = 19.7;
      rval.image_type = 'k';        /* 'Kookbook' */
      }
   else
      {
      fclose( ifile);
      return( -3);
      }

   if( !rval.is_xy_list)
      {
      if( debug_level)
         {
         log_printf( log_file, "Image is %ld x %ld bytes\n",
                                       rval.xsize, rval.ysize);
#ifdef DEBUG_MEM
         debug_dump_memory( );
#endif
         log_printf( log_file, "Allocating memory for image...");
         }
      rval.img = (PIXEL *)calloc( rval.xsize * rval.ysize, sizeof( PIXEL));
      if( debug_level)
         log_printf( log_file, "Pointer returned is %p\n", rval.img);
      if( !rval.img)
         {
         fclose( ifile);
         return( -4);
         }
      if( debug_level)
         log_printf( log_file, "Memory allocated successfully\n");
      }

   switch( rval.image_type)
      {
      case '4':
         {
         unsigned char tbuff[ST4_WIDTH];
         PIXEL *tptr = rval.img;

         rval.data_offset = 0;
         fseek( ifile, 0L, SEEK_SET);
         for( i = 0; i < ST4_HEIGHT; i++)
            {
            fread( tbuff, ST4_WIDTH, 1, ifile);
            for( j = 0; j < ST4_WIDTH; j++)
               *tptr++ = tbuff[j];
            }
         fread( tbuff, ST4_WIDTH, 1, ifile);
         sscanf( (char *)tbuff + 89, "%lf", &rval.focal_len);
         sscanf( (char *)tbuff + 21, "%d %d %lf", &year, &month, &day);
         rval.byte_order = 2;       /* just a byte per pixel... */
         }
         break;
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
         {
         rval.data_offset = 2048;
         fseek( ifile, 2048L, SEEK_SET);
         if( !is_compressed)
#ifdef PIXEL_32
            {
            unsigned short *uptr = (unsigned short *)rval.img;

            fread( rval.img, rval.xsize, rval.ysize * sizeof( unsigned short),
                                             ifile);
            for( i = rval.xsize * rval.ysize - 1; i >= 0; i--)
               rval.img[i] = (PIXEL)uptr[i];
            }
#else
            fread( rval.img, rval.xsize, rval.ysize * sizeof( unsigned short),
                                             ifile);
#endif
         else
            {
            PIXEL *tptr = rval.img;
            short len;
            char *tbuff = (char *)calloc( rval.xsize, sizeof( unsigned short));
            char *cptr;

            rval.byte_order = 3;       /* compressed SBIG form: */
            for( i = 0; i < rval.ysize; i++)
               {
               fread( &len, sizeof( short), 1, ifile);
               fread( tbuff, len, 1, ifile);
               if( len == 2 * rval.xsize)
                  for( j = 0; j < rval.xsize; j++)
                     *tptr++ = ((unsigned short *)tbuff)[j];
               else
                  {
                  *tptr++ = *(unsigned short *)tbuff;
                  cptr = tbuff + 2;
                  for( j = (PIXEL)( rval.xsize - 1); j; j--, tptr++)
                     if( *cptr == -128)
                        {
                        *tptr = *(unsigned short *)( cptr + 1);
                        cptr += 3;
                        }
                     else
                        *tptr = tptr[-1] + (*cptr++);
                  }
               }
            free( tbuff);
            }
         }
         break;
      case 'c':
         {
         char *iptr;

         fseek( ifile, 80L, SEEK_SET);
         fread( ibuff, 90, 1, ifile);
         for( iptr = ibuff, i = 0; i < 70; i++, iptr++)
            if( iptr[2] == '/' && iptr[5] == '/')
               {
               if( atoi( iptr + 6) > 80 && atoi( iptr + 6) < 2010)
                  if( !memcmp( ibuff + 15, "DD/MM/YY", 7))
                     sscanf( iptr, "%d/%lf/%d", &month, &day, &year);
                  else
                     sscanf( iptr, "%lf/%d/%d", &day, &month, &year);
               }
            else if( iptr[2] == ':' && iptr[5] == ':' && iptr[8] == '.')
               sscanf( iptr, "%d:%d:%lf", &hour, &minute, &seconds);
         if( debug_level)
            log_printf( log_file, "Date: %d %d %d %d:%02d:%02d\n",
                            month, (int)day,
                            year, hour, minute, (int)seconds);
         rval.exposure_length = atof( ibuff + 75);
         if( debug_level)
            log_printf( log_file, "Exposure length:  %.2lf seconds\n",
                          rval.exposure_length);
         rval.data_offset = (long)PIC_OFFSET;
         fseek( ifile, (long)PIC_OFFSET, SEEK_SET);
         if( !is_compressed)
            for( i = rval.ysize - 1; i >= 0; i--)
               {
               PIXEL *tptr = rval.img + i * rval.xsize;

               fread( tptr, rval.xsize, sizeof( unsigned short), ifile);
#ifdef PIXEL_32
               for( j = rval.xsize - 1; j >= 0; j--)
                  tptr[j] = *((unsigned short *)tptr + j);
#endif
               }
         if( debug_level)
            log_printf( log_file, "Image data read\n");
         }
         break;
      case 'f':            /* FITS img */
      case 'm':            /* MIRA img */
      case 'x':            /* Starlight XPress img */
      case 'a':            /* Mike Kretlow .AST file */
      case 'b':            /* BMP file */
         {
         rval.data_offset = ftell( ifile);
         for( i = rval.ysize - 1; i >= 0; i--)
            {
            unsigned char *tptr = (unsigned char *)( rval.img + i * rval.xsize);
#ifdef PIXEL_32
            PIXEL *tpixel_ptr = (PIXEL *)tptr;
#endif

            if( rval.image_type == 'b')
               {
               fread( tptr, real_line_size, 1, ifile);
               if( bits_per_pixel == 8)
                  for( j = rval.xsize - 1; j >= 0; j--)
                     ((PIXEL *)tptr)[j] = (PIXEL)tptr[j];
               else
                  for( j = rval.xsize - 1; j >= 0; j--)
                     ((PIXEL *)tptr)[j] =
                          (PIXEL)(tptr[j * 3] + tptr[j * 3 + 1] + tptr[j * 3 + 2]);
               }
            else if( bits_per_pixel == 8)    /* J Greaves fix,  I hope */
               {
               fread( tptr, rval.xsize, 1, ifile);
               for( j = rval.xsize - 1; j >= 0; j--)
                  *(PIXEL *)( tptr + j * 4) = (PIXEL)tptr[j];
               }
#ifdef PIXEL_32
            else if( bits_per_pixel == 32 || bits_per_pixel == -32)
               {
               fread( tptr, rval.xsize, sizeof( unsigned long), ifile);
               if( rval.byte_order)
                  {
                  for( j = rval.xsize; j; j--)
                     {
                     unsigned char temp = *tptr;    /* byte swapping time */

                     tptr[0] = tptr[3];
                     tptr[3] = temp;
                     temp = tptr[1];
                     tptr[1] = tptr[1];
                     tptr[2] = temp;
                     if( bits_per_pixel == -32)
                        *(long *)tptr = (long)( *(float *)tptr * 100.);
                     tptr += 4;
                     }
                  tptr -= 4 * rval.xsize;
                  }
               }
#endif
            else              /* 16-bit image */
               {
               fread( tptr, rval.xsize, sizeof( unsigned short), ifile);
               if( rval.byte_order)
                  {
                  for( j = rval.xsize; j; j--)
                     {
                     unsigned char temp = *tptr;          /* byte swapping time */

                     tptr[0] = tptr[1];
                     tptr[1] = temp;
                     tptr += 2;
                     }
                  tptr -= 2 * rval.xsize;
                  }
#ifdef PIXEL_32
               for( j = rval.xsize - 1; j >= 0; j--)
                  tpixel_ptr[j] = *((unsigned short *)tptr + j);
#endif
               }
            }
#ifdef WHY_IS_THIS_HERE
         if( rval.image_type == 'f')
            for( i = rval.ysize * rval.xsize; i; i--, iptr++)
               if( *iptr == (PIXEL)32768)
                  *iptr = 0;
#endif
         if( colour_starlight_express)
            {
                        /* In colour SX images,  only every other line and */
                        /* every other column is to be taken seriously. */
            rval.xsize /= 2;
            rval.ysize /= 2;
            for( j = 0; j < rval.ysize; j++)
               for( i = 0; i < rval.xsize; i++)
                  rval.img[i + j * rval.xsize] =
                              rval.img[i + i + (4 * j + 2) * rval.xsize];
            }
         }
         break;
      case 'k':         /* Kookbook camera, in two or three parts */
         {
         char part_ii_name[60], *end_ptr;
         FILE *part_ii;

         fseek( ifile, 7L, SEEK_SET);
         strcpy( part_ii_name, filename);
         end_ptr = part_ii_name + strlen( filename) - 1;
         if( *end_ptr == 'a')       /* CB245 image,  252x242 pixels */
            {
            fread( rval.img, 252u * 125u, sizeof( short), ifile);
            *end_ptr = 'b';
            part_ii = fopen( part_ii_name, "rb");
            if( part_ii)
               {                               /* 242-125=117 lines in part ii */
               fseek( part_ii, 7L, SEEK_SET);
               fread( rval.img + 252u*125u, 252u*117u, sizeof(short), part_ii);
               fclose( part_ii);
               }
            }
         else if( *end_ptr == '1')    /* 378x242 image,  three parts */
            {
            fread( rval.img, 378u * 81u, sizeof( short), ifile);
            *end_ptr = '2';
            part_ii = fopen( part_ii_name, "rb");
            if( part_ii)
               {                               /* 242-125=117 lines in part ii */
               fseek( part_ii, 7L, SEEK_SET);
               fread( rval.img + 378u*81u, 378u*81u, sizeof(short), part_ii);
               fclose( part_ii);
               }
            *end_ptr = '3';
            part_ii = fopen( part_ii_name, "rb");
            if( part_ii)
               {                               /* 242-125=117 lines in part ii */
               fseek( part_ii, 7L, SEEK_SET);
               fread( rval.img + 378u*162u, 378u*80u, sizeof(short), part_ii);
               fclose( part_ii);
               }
            }
         }
         break;
      case 'p':         /* .MAP file */
      case 't':         /* .LST file */
      case 'n':         /* .FIN file */
         break;
      default:
         return( -5);
      }
   fclose( ifile);

   if( signed_ints)
      {
      PIXEL *tptr = rval.img;

      if( bits_per_pixel == 32)
         for( i = rval.xsize * rval.ysize; i; i--)
            *tptr++ ^= 0x80000000;
      else
         for( i = rval.xsize * rval.ysize; i; i--)
            *tptr++ ^= 0x8000;
      }

   if( !rval.is_xy_list && year != -1)
      {
      char buff[80];

      if( year > 50 && year < 100)
         year += 1900;
      if( year >= 0 && year < 50)
         year += 2000;
      day += (double)dmy_to_day( 0, month, (long)year, 0) + (double)hour / 24.
               + (double)minute / 1440. + seconds / 86400. - .5;
      full_ctime( buff, day, 0);
      log_printf( log_file, strings[65], buff);   /* "Date extracted from header: %s\n" */
      rval.time = day;
      }
   else
      log_printf( log_file, strings[66]);         /* "No date given in header\n" */
   rval.filename = (char *)malloc( strlen( filename) + 1);
   strcpy( rval.filename, filename);
   memcpy( img, &rval, sizeof( IMAGE));
   if( debug_level)
      log_printf( log_file, "load_image( ) completed\n");
   return( 0);
}
