#include "gsc.h"

#define N_LARGE_REGIONS   732

static char zone_size[24] = {48, 47, 45, 43, 40, 36,
                            32, 27, 21, 15, 9, 3,
                            48, 47, 45, 43, 40, 36,
                            32, 27, 21, 15, 9, 3};

static unsigned char n_subdivs[N_LARGE_REGIONS / 4] = {
          0x55, 0x55, 0xa9, 0xaa, 0x5a, 0x55, 0x55, 0xa5, 0xaa, 0xaa,
          0xaa, 0x55, 0x55, 0x55, 0xaa, 0xaa, 0x56, 0x55, 0x55, 0xa5,
          0xaa, 0xaa, 0xaa, 0x55, 0x55, 0xa9, 0xaa, 0x6a, 0x55, 0x40,
          0x55, 0xa9, 0xaa, 0xaa, 0x5a, 0x55, 0xaa, 0xaa, 0x6a, 0x15,
          0x00, 0x55, 0xa9, 0xaa, 0xaa, 0x9a, 0xaa, 0xaa, 0xaa, 0x56,
          0x01, 0x54, 0x95, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0x6a, 0x15,
          0x40, 0x55, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0x56, 0x05, 0x54,
          0xa5, 0xaa, 0xaa, 0xaa, 0xaa, 0x55, 0x55, 0xa5, 0xaa, 0xaa,
          0xaa, 0x56, 0x55, 0xa9, 0xaa, 0xaa, 0x55, 0xa5, 0xaa, 0x56,
          0xaa, 0x5a, 0x50, 0x95, 0xaa, 0xaa, 0x56, 0x55, 0x95, 0xaa,
          0xaa, 0xaa, 0x5a, 0x15, 0x40, 0x55, 0xaa, 0xaa, 0x5a, 0x55,
          0xa5, 0xaa, 0xaa, 0xaa, 0x56, 0x05, 0x50, 0xa5, 0xaa, 0xaa,
          0x56, 0xa5, 0xaa, 0xaa, 0xaa, 0x5a, 0x15, 0x40, 0x95, 0xaa,
          0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0x56, 0x05, 0x54, 0xa5,
          0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0x56, 0x55, 0x55, 0xa9,
          0xaa, 0xaa, 0xaa, 0xaa, 0xaa, 0x5a, 0x55, 0x55, 0xaa, 0xaa,
          0xaa, 0xaa, 0xaa, 0x5a, 0x55, 0x95, 0xaa, 0xaa, 0xaa, 0xaa,
          0x5a, 0x55, 0xaa, 0xaa, 0xaa, 0x6a, 0x55, 0xaa, 0xaa, 0x6a,
          0xaa, 0xaa, 0xaa };

#define LARGE_SIZE 750000L

/*

   find_gsc_regions( ) takes as input an RA/dec rectangle,  running from
'min' to 'max',  and determines the GSC tiles that would cover that
rectangle.  It returns the tile numbers in the array 'ret_arr',  up
to a maximum of 'max_ret'.

   'min' and 'max' are in units of 1/100000 degree,  because those are the
source units of the GSC.  'level' should be set to 1 to get GSC tile
numbers,  or 2 to get large area numbers.  The concept of a GSC large area
is discussed below.

   The GSC divides the sky at three levels.  The first level is into
24 zones of declination,  each 7.5 degrees high.  Logically,  these would
run from the north celestial pole to the south,  or vice versa.  Instead,
the first zone covers dec 0 to +7.5;  the next,  +7.5 to +15,  and so
on northward to the pole.  Zone 13 covers dec 0 to -7.5;  zone 14,
dec -7.5 to -15.0;  and so on,  to the south pole.

   Each zone is subdivided in RA into three to 48 "large areas".  The
two zones bordering the celestial equator are each subdivided into 48
large areas,  each of which is 7.5 degrees wide... so they are basically
"square".  The number of large areas in each zone was chosen to maintain
that property,  as best as was possible;  thus,  the zones bordering the
poles are only subdivided into three zones,  each covering a 120-degree
slice of RA.  The number of subdivisions is stored in the 'zone size'
array,  but can be figured by computing 48. * cos( center_dec),  where
'center_dec' is the declination of the center of the zone,  and rounding
(not truncating!) the result.  I was lazy and stored the result in a
24-element array instead.

   And now we come to third,  oddest level of GSC subdivision.  Each
'large area' is subdivided into a 2x2,  3x3,  or 4x4 grid of actual
GSC tiles,  or 'small areas'.  The subdivision into 'large areas' was
algorithmically based,  but that into small areas depends on where you
are in the sky.  In sparse regions,  such as Virgo,  the subdivision is
usually 2x2.  In dense regions,  such as Scorpius,  it's usually 4x4.
The subdivision scale is encoded in the 'n_subdivs' array,  above.  This
code was originally written in 1992,  at a time when devoting 792 bytes
to _anything_ would be considered excessive.  So instead,  the data was
compressed into two bits per subdivision,  or 4:1,  so the array really
only requires 792/4 = 198 bytes.  (I mention this should anyone wonder
why I went to such extreme measures.  Back in 1992,  expending 792 bytes
was being a little profligate with memory.)

*/

int find_regions( const LPT *min, const LPT *max, SMALL_AREA FAR *ret_arr,
                                    const int max_ret, const int level)
{
   int zone, rval = 0, loc = 1, large_loc = 0, i, loc1;
   long south;

   if( min->x > THREESIXTY_DEG || max->x < 0L)
      return( 0);
   for( zone = 0; zone < 24; zone++, south += LARGE_SIZE)
      {
      if( zone < 12)
         south = (long)zone * LARGE_SIZE;
      else
         south = -(long)(zone - 11) * LARGE_SIZE;
      if( south < max->y && south + LARGE_SIZE > min->y)
         {
         loc1 = loc;
         for( i = 0; i < zone_size[zone]; i++)
            {
            int subdivs, i1, j1;
            long west, east, step;

            i1 = large_loc + i;
            j1 = (i1 & 3) * 2;
            subdivs = ((n_subdivs[i1 / 4] >> j1) & 3) + 2;
            east = THREESIXTY_DEG * (long) i      / (long)zone_size[zone];
            west = THREESIXTY_DEG * (long)(i + 1) / (long)zone_size[zone];
            step = west - east;
            if( west > min->x && east < max->x && rval < max_ret)
               {           /* first store large area */
               if( level & 2)
                  {
                  if( ret_arr)
                     {
                     ret_arr->zone = zone;
                     ret_arr->x1 = east;
                     ret_arr->x2 = west;
                     ret_arr->y1 = south;
                     ret_arr->y2 = south + LARGE_SIZE;
                     ret_arr->small_zone = -(large_loc + i + 1);
                     ret_arr++;
                     }
                  rval++;
                  }
               if( level & 1)
                  for( i1 = 0; i1 < subdivs; i1++)
                     for( j1 = 0; j1 < subdivs && rval < max_ret; j1++)
                        {
                        long x1, y1, x2, y2;

                        x1 = east + (long) i1      * step / (long)subdivs;
                        x2 = east + (long)(i1 + 1) * step / (long)subdivs;
                        y1 = south + (long)j1 * LARGE_SIZE / (long)subdivs;
                        y2 = y1 + LARGE_SIZE / subdivs;
                        if( y1<max->y && y2>min->y && x1 < max->x && x2 > min->x)
                           {
                           if( ret_arr)
                              {
                              ret_arr->zone = zone;
                              ret_arr->x1 = x1;
                              ret_arr->x2 = x2;
                              ret_arr->y1 = y1;
                              ret_arr->y2 = y2;
                              if( zone < 12)
                                 ret_arr->small_zone = loc1 + i1 + j1 * subdivs;
                              else
                                 ret_arr->small_zone = loc1 + i1 +
                                                   (subdivs - j1 - 1) * subdivs;
                              ret_arr++;
                              }
                           rval++;
                           }
                        }     /* end of search for small zones */
               }
            loc1 += subdivs * subdivs;
            }
         }
      for( i = 0; i < zone_size[zone]; i++, large_loc++)
         {
         int subdiv = n_subdivs[large_loc / 4] >> (2 * (large_loc & 3));

         subdiv = (subdiv & 3) + 2;
         loc += subdiv * subdiv;
         }
      }     /* end of for( zone...) loop */
   return( rval);
}

/*
   small_zone_to_large_data( ) won't be of much use to most of the world.
Guide uses it in situations where it knows of a GSC tile number,  and
wants to know in which large area and zone it would fit,  and where it
would fit within the large area.  This helps for some of the bizarre
indexing schemes it uses.  */

int small_zone_to_large_data( const int small_zone, int *ret_data)
{
   int zone, i, loc, sm_start, subdivs;
   static short zone_start[25] = {    1,  594, 1178, 1729, 2259, 2781,
                                   3246, 3652, 4014, 4294, 4492, 4615,
                                   4663, 5260, 5838, 6412, 6989, 7523,
                                   8022, 8464, 8840, 9134, 9346, 9490, 9538 };

   if( small_zone < 1 || small_zone > 9537)
      return( -1);
   for( zone = 24; zone_start[zone] > small_zone; zone--)
      ;
   sm_start = zone_start[zone];
   for( i = loc = 0; i < zone; i++)
      loc += zone_size[i];
   for( i = loc; sm_start <= small_zone; i++)
      {
      int shift = (i & 3) * 2;

      subdivs = ((n_subdivs[i / 4] >> shift) & 3) + 2;
      sm_start += subdivs * subdivs;
      }
   sm_start -= subdivs * subdivs;
   if( ret_data)
      {
      *ret_data++ = zone;
      *ret_data++ = i;
      *ret_data++ = small_zone - sm_start;
      *ret_data++ = subdivs;
      }
   return( zone);
}

/*
   find_gsc_rect( ) takes as input the number of a large area (assumed
to run from -1 to -732) or of a GSC tile (assumed to run from 1 to 9537),
and computes its bounding rectangle (stored as rect[0...3]).  One is
more commonly interested in a reverse procedure,  where you need to know
which tiles cover a given area... but a few parts of Guide use the following
handy function. */

int find_gsc_rect( long *rect, const short region)
{
   int i, rval = -1;

   if( region < 0 && region >= -732)  /* looking for a _large_ region */
      {
      int count = -region - 1, zone = 0, subdivs;

      for( zone = 0; count >= zone_size[zone]; zone++)
         count -= zone_size[zone];
      subdivs = zone_size[zone];
      rval = zone;
      if( zone >= 12)            /* allow for southern hemisphere flip */
         zone = 11 - zone;
      for( i = 0; i < 2; i++)
         {
         *rect++ = (long)( count + i) * THREESIXTY_DEG / (long)subdivs;
         *rect++ = (long)( zone + i) * LARGE_SIZE;
         }
      }
   if( region > 0 && region <= 9537)
      {
      long lg_rect[4];
      int zone_data[4], zone = small_zone_to_large_data( region, zone_data);
      int x, y;

      find_gsc_rect( lg_rect, (short)-zone_data[1]);
      x = zone_data[2] % zone_data[3];
      y = zone_data[2] / zone_data[3];
      if( zone >= 12)                  /* again,  s hemisphere flip */
         y = zone_data[3] - 1 - y;
      lg_rect[2] -= lg_rect[0];
      lg_rect[3] -= lg_rect[1];
      for( i = 0; i < 2; i++)
         {
         *rect++ = lg_rect[0] + lg_rect[2] * (long)(x+i) / (long)zone_data[3];
         *rect++ = lg_rect[1] + lg_rect[3] * (long)(y+i) / (long)zone_data[3];
         }
      rval = *zone_data;
      }
   return( rval);
}

#ifdef TEST_CODE

/*
   A little 'example' program for exercising two of the above functions.
I found it handy in debugging,  but it also gives an indication of how
these functions are used... assuming you care,  which is unlikely;  the
find_gsc_regions( ) function is _much_ more likely to get actual use.
*/

#include <stdlib.h>
#include <stdio.h>

void main( int argc, char **argv)
{
   int i, start = atoi( argv[1]), end = atoi( argv[2]);

   for( i = start; i <= end; i++)
      {
      int lg_data[4];
      long rect[4];

      if( i > 0)
         small_zone_to_large_data( i, lg_data);
      else
         lg_data[0] = lg_data[1] = lg_data[2] = lg_data[3] = 0;
      find_gsc_rect( rect, (short)i);
      printf( "%5d: %2d %3d %2d %1d: %9ld %9ld  %9ld %9ld\n",
                              i, lg_data[0], lg_data[1],
                              lg_data[2], lg_data[3],
                              rect[0], rect[2], rect[1], rect[3]);
      }
}
#endif
