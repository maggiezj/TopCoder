#include <stdlib.h>
#include <string.h>
#include <graph.h>
#include <math.h>
#include <dos.h>
#include "watdefs.h"
#include "findstar.h"

#ifdef DEBUG_MEM
#include "checkmem.h"
#endif

#ifdef __WATCOMC__
#define HDR_OFFSET 6
#define DISPMEM ((char *)0xA0000)
#else
#define HDR_OFFSET 4
#define DISPMEM ((char far *)0xA0000000)
#endif

void show_image( const PIXEL *buff, const int width, const int height,
      const PIXEL low, const PIXEL high, const double *image_loc,
      const int *disp_loc, const int numcolors);         /* dispimg.cpp */
static void dump_16_color_line( int x, int y, char *ibuff);

int use_pure_vesa = 1;

static void set_vesa_page( short page)
{
   union REGS r;

   r.h.ah = 0x4f;
   r.h.al = 5;
   r.w.bx = 0;
   r.w.dx = page;
   int386( 0x10, &r, &r);
}

static void dump_vesa_line( const int line, int width, char *pixels)
{
   long offset = (long)line * (long)width;
   long offset1 = offset + (long)width;
   static int curr_page = -1;
   int page = (int)( offset >> 16);
   int page1 = (int)( offset1 >> 16);

   if( page != curr_page || line == 0)
      {
      set_vesa_page( page);
      curr_page = page;
      }
   if( page == page1)
      memcpy( DISPMEM + (offset & 0xffff), pixels, width);
   else
      {
      int width1 = ((long)page1 << 16) - offset;

      memcpy( DISPMEM + (offset & 0xffff), pixels, width1);
      set_vesa_page( page1);
      memcpy( DISPMEM, pixels + width1, width - width1);
      curr_page = page1;
      }
}

int display_method = 0;

#define DISPLAY_LINEAR           0
#define DISPLAY_LOGARITHMIC      1
#define DISPLAY_HISTOGRAM        2

void show_image( const PIXEL *buff, const int width, const int height,
      const PIXEL low, const PIXEL high, const double *image_loc,
      const int *disp_loc, const int numcolors)
{
   int i, j;
   const PIXEL *tptr;
   const int disp_height = disp_loc[3] - disp_loc[1];
   const int disp_width  = disp_loc[2] - disp_loc[0];
   const double image_height = image_loc[3] - image_loc[1];
   const double image_width  = image_loc[2] - image_loc[0];
   int image_y, prev_image_y = -2;
   int *x = calloc( disp_width, sizeof( int));
   unsigned char *obuff = malloc( 16 + disp_width);
   unsigned char *remap_tbl = NULL;
   long multiplier;
   const int sixteen = (numcolors == 16);
   const int zero_offset = (sixteen ? 0 : 16);
   const char white_pixel_value = (sixteen ? 10 : 79);
   const int max_tbl_size = 70000;
   static const unsigned char dither_matrix[64] = {
        0 << 2, 32 << 2,  8 << 2, 40 << 2,  2 << 2, 34 << 2, 10 << 2, 42 << 2,
       48 << 2, 16 << 2, 56 << 2, 24 << 2, 50 << 2, 18 << 2, 58 << 2, 26 << 2,
       12 << 2, 44 << 2,  4 << 2, 36 << 2, 14 << 2, 46 << 2,  6 << 2, 38 << 2,
       60 << 2, 28 << 2, 52 << 2, 20 << 2, 62 << 2, 30 << 2, 54 << 2, 22 << 2,
        3 << 2, 35 << 2, 11 << 2, 43 << 2,  1 << 2, 33 << 2,  9 << 2, 41 << 2,
       51 << 2, 19 << 2, 59 << 2, 27 << 2, 49 << 2, 17 << 2, 57 << 2, 25 << 2,
       15 << 2, 47 << 2,  7 << 2, 39 << 2, 13 << 2, 45 << 2,  5 << 2, 37 << 2,
       63 << 2, 31 << 2, 55 << 2, 23 << 2, 61 << 2, 29 << 2, 53 << 2, 21 << 2 };

   *(short *)obuff = (short)disp_width;
   *(short *)(obuff + 2) = 1;               /* 1 pixel high... */
   *(short *)(obuff + 4) = 8;               /* ...eight bits deep */
   obuff += HDR_OFFSET;
   for( i = 0; i < disp_width; i++)
      {
      x[i] = (int)( image_loc[0] + (double)i * image_width / disp_width);
      if( x[i] < 0 || x[i] >= width)
         x[i] = -1;
      }

   if( sixteen)
      multiplier = 10L * 65536L / (long)( high - low);
   else
      multiplier = 63L * 65536L / (long)( high - low);
   if( high - low < max_tbl_size && !sixteen)
      {
      remap_tbl = (unsigned char *)calloc( high - low, sizeof( char));
      if( display_method == DISPLAY_LINEAR)
         for( i = 0; i < high - low; i++)
            remap_tbl[i] = (char)( 16 + 63L * (long)i / (long)( high - low));
      else if( display_method == DISPLAY_LOGARITHMIC)
         {
         const double log_low = log( (double)low);
         const double log_mult = 63. / (log( (double)high) - log_low);

         for( i = 0; i < high - low; i++)
            {
            const double tval = (log( (double)( i + low)) - log_low) * log_mult;

            remap_tbl[i] = (char)( 16 + (short)tval);
            }
         }
      remap_tbl -= low;
      }

   for( i = 0; i < disp_height; i++)
      {
      const unsigned char *dither_ptr = dither_matrix + 8 * (i & 7);

      image_y = (int)( image_loc[1] + image_height * (double)i / disp_height);
      if( image_y < 0 || image_y >= height)
         image_y = -1;
      if( image_y != prev_image_y || sixteen)
         {
         prev_image_y = image_y;
         memset( obuff, zero_offset, disp_width);
         if( image_y != -1)
            {
            tptr = buff + width * image_y;
            for( j = disp_width; j; j--, obuff++, x++)
               if( *x >= 0)
                  {
                  PIXEL ival = tptr[*x];

                  if( ival >= high)
                     *obuff = white_pixel_value;
                  else if( ival > low)
                     if( remap_tbl)
                        *obuff = remap_tbl[ival];
                     else
                        {
                        long tval = multiplier * (long)( ival - low);

                        *obuff = (unsigned char)( zero_offset +
                                    ((unsigned short *)&tval)[1]);
                        if( sixteen)
                           if( ((unsigned char *)&tval)[1] > dither_ptr[j & 7])
                              (*obuff)++;
                        }
                  }
            obuff -= disp_width;
            x -= disp_width;
            }
         }
      if( sixteen)
         for( j = 0; j < disp_width; j += 32)
            dump_16_color_line( j, i, (char *)obuff + j);
      else if( use_pure_vesa)
         dump_vesa_line( (short)( disp_loc[1] + i), disp_width, (char *)obuff);
      else
         _putimage( (short)disp_loc[0], (short)( disp_loc[1] + i),
                        (char *)obuff - HDR_OFFSET, _GPSET);
      }
   free( obuff - HDR_OFFSET);
   free( x);
   if( remap_tbl)
      free( remap_tbl + low);
}

static void dump_16_color_line( int x, int y, char *ibuff)
{
   char image[16 + HDR_OFFSET], *tptr = image + HDR_OFFSET;
   static const unsigned char bits[8] = { 128, 64, 32, 16, 8, 4, 2, 1 };
// static const unsigned char bits[8] = { 1, 2, 4, 8, 16, 32, 64, 128 };
   int i, j, plane, bit;
   const int n_pixels = 32;

   *(short *)image = (short)( n_pixels);
   *(short *)(image + 2) = 1;               /* 1 pixel high... */
   *(short *)(image + 4) = 4;               /* ...four bits deep */
   memset( tptr, 0, n_pixels / 2);
   for( i = 0; !ibuff[i] && i < n_pixels; i++)
      ;
   if( i != n_pixels)         /* can skip 'blanks' */
      for( bit = 8, plane = 0; plane < 4; plane++, bit >>= 1)
         {
         for( i = 0; i < n_pixels / 8; i++)
            for( j = 0; j < 8; j++, ibuff++)
               if( *ibuff & bit)
                  tptr[i] |= bits[j];
         ibuff -= n_pixels;
         tptr += n_pixels / 8;
         }
   _putimage( (short)x, (short)y, image, _GPSET);
}
