#include <string.h>
#include "crunch.h"

#define THREESIXTY_DEG 36000000L

int gsc_star_uncrunch( GSC_STAR *star, GSC_HEADER *hdr,
                                         const unsigned char *buff)
{
   int rval;
   unsigned short cmag, temp;
   long dx, dy;

   memset( star, 0, sizeof( GSC_STAR));
   star->plate_id_no = (short)(*buff & 0x7);
   if( *buff & 0x10)
      {
      star->n = hdr->prev_n;
      star->is_a_dup = 1;
      rval = 5;
      }
   else
      {
      star->n = hdr->prev_n + 1;
      rval = 8;
      }
   if( star->plate_id_no == 7)
      star->plate_id_no = (unsigned short)buff[rval++];
   if( *buff & 0x08)
      star->obj_type = (unsigned short)buff[rval++];
   memcpy( star->plate, hdr->plates + 5 * star->plate_id_no, 5);
   temp = (unsigned short)*buff >> 5;
   if( temp == 7)    /* mag error didn't fit */
      {
      star->mag_err = buff[rval++];
      if( star->mag_err == 255)     /* didn't fit a byte either */
         {
         star->mag_err = *(unsigned short *)(buff + rval);
         rval += 2;
         }
      }
   else
      star->mag_err = hdr->magerr_0 + temp;
   cmag = ((unsigned short)buff[1] | ((unsigned short)buff[2] << 8)) & 1023;
   if( cmag == 1023)    /* mag maxed out & didn't fit in ten bits */
      {
      cmag = (unsigned short)buff[rval++] << 8;
      cmag |= (unsigned short)buff[rval++];
      }
   else
      cmag += 600;
   star->mag = cmag;
   if( buff[2] & 4)
      star->mag_band = (unsigned short)buff[rval++];
   else
      star->mag_band = hdr->magband_0;
   temp = ((unsigned short)buff[2] >> 3) & 7;
   if( temp == 7)       /* posn error didn't fit */
      {
      star->posn_err = (unsigned short)buff[rval++];
      if( star->posn_err == 255)    /* posn error didn't fit unsigned char */
         {
         star->posn_err = (unsigned short)buff[rval]
                       | ((unsigned short)buff[rval + 1] << 8);
         rval += 2;
         }
      }
   else
      star->posn_err = temp;
   if( star->n == hdr->prev_n)
      {           /* difference encoding time */
      dx = (long)buff[3] | (((long)buff[2] << 2) & 0xf00);
      dy = (long)buff[4];
      if( dx == 1023L || dy == 255L)      /* need a bigger difference */
         {
         dx = ((long)buff[rval  ] << 8) + (long)buff[rval + 1] - 32768L;
         if( dx != 32767L)
            {
            dy = ((long)buff[rval+2] << 8) + (long)buff[rval + 3] - 32768L;
            rval += 4;
            }
         else
            {
            dx = *(long *)(buff + rval + 2);
            dy = *(long *)(buff + rval + 6);
            rval += 10;
            }
         }
      else
         {
         dx -= 512L;
         dy -= 128L;
         }
      star->x = hdr->prev_x + dx;
      star->y = hdr->prev_y + dy;
      }
   else           /* star position isn't difference encoded */
      {
      dx = ((long)buff[5] << 16) | ((long)buff[6] << 8) | (long)buff[7];
      dx &= 0x7fffff;
      dy = (((long)buff[2] >> 6) << 17) | ((long)buff[3] << 8) | (long)buff[4];
      if( buff[5] & 0x80)
         dy |= 0x10000;
      star->x = hdr->ra0  + dx;
      star->y = hdr->dec0 + dy;
      }
   if( star->x < 0L)
      star->x += THREESIXTY_DEG;
   if( star->x >= THREESIXTY_DEG)
      star->x -= THREESIXTY_DEG;
   hdr->prev_x = star->x;
   hdr->prev_y = star->y;
   hdr->prev_n = star->n;
   return( rval);
}
