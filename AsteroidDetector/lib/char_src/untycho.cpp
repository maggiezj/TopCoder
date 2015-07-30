#include <string.h>
#include "watdefs.h"
#include "tycho.h"

static int extract_mag_err( const int packed_err)
{
   int rval;

   if( packed_err < 15)
      rval = (short)packed_err;
   else
      {
      int count = (packed_err >> 3) - 1;

      rval = (packed_err & 7) + 8;
      rval <<= count;
      }
   return( rval);
}

/*
  Some changes were made for version 8.  I wanted to include the error
measurements for VT and BT.  Tycho-1 (used in Guides <= 7) included
Johnson V magnitude data;  for Tycho-2 (Guides >= 8),  this has to
be computed from VT and BT.  Also,  Tycho-2 includes epoch data for
the RA and dec.  (In Tycho-1,  all data was epoch 1991.25.)

   So,  if the 'version number' indicates Guide 8 data,  the code takes
the "Johnson V" data and extracts sigmas for VT and BT,  and then computes
Johnson V from BT and VT data.  It grabs epoch data for the RA and dec
(a byte apiece),  filling in '91' for both values for Tycho-1 (Guide <= 7)
data.

   Of older interest:  you'll see a couple of "if( version_number >= 7)"
lines,  because a few changes were made between Guides 6 and 7.  (The
dataset was still Tycho-1,  but I'd realized that I'd omitted some items
that I really should have included.)  Thus,  for a Guide 7 or later disk,
the code reads in parallax,  BT,  VT,  and some sigmas.  (With Guide 6,
this data was stored in a separate file... a really bad move on my part.
For details on accessing that file,  see TYC_EXT.CPP.)
*/

int parse_tycho_star( const char FAR *data_ptr, TYCHO_STAR *star,
                      const int version_number)
{
   char flags = data_ptr[14];
   const char FAR *tptr;
   int i;

   memset( star, 0, sizeof( TYCHO_STAR));
   FMEMCPY( star, data_ptr, 14);
   tptr = data_ptr + 15;
   if( version_number >= 7)
      {
      const short FAR *extras = (const short FAR *)tptr;

      star->parallax =     extras[0];
      star->parallax_err = extras[1];
      star->posn_err[0] =  extras[2];
      star->posn_err[1] =  extras[3];
      star->bt_mag =       extras[4];
      star->vt_mag =       extras[5];
               /* recompute Johnson V magnitude,  using method from p. 57, */
               /* _Intro & Guide to the Data_,  and extract BT & VT errors */
      if( version_number >= 8)
         {
         const short bt_minus_vt = star->bt_mag - star->vt_mag;

         star->vt_mag_err = (short)extract_mag_err( star->mag & 63);
         star->bt_mag_err = (short)extract_mag_err( (star->mag >> 6) & 63);
         star->mag = star->vt_mag - (short)( bt_minus_vt * 9 / 100L);
         if( bt_minus_vt < 500)
            {
            if( bt_minus_vt < 100)
               star->mag -= (short)( 7 * (bt_minus_vt - 100) / 150);
            else
               star->mag -= (short)( (bt_minus_vt - 100) / 80);
            }
         else        /* 'redder' end */
            {
            if( bt_minus_vt < 1400)
               star->mag -= 5;
            else
               star->mag -= (short)( 5 + (bt_minus_vt - 1400) / 40);
            }
         if( !star->bt_mag)         /* got nothing but vt to go on */
            star->mag = star->vt_mag;
         else if( !star->vt_mag)    /* got nothing but bt to use */
            star->mag = star->bt_mag;
//       star->mag /= 10;
         }
      else                          /* mags now in .001 units,  not .01 */
         star->mag *= 10;
      tptr += 6 * sizeof( short);
      }

   if( version_number >= 8)
      {
      star->epoch_ra = *tptr++;
      star->epoch_dec = *tptr++;
      }
   else
      star->epoch_ra = star->epoch_dec = 91;

   if( flags & 1)
      {
      FMEMCPY( star->gsc_id, tptr, 3);
      tptr += 3;
      }
   if( flags & 2)
      {
      FMEMCPY( &star->ppm_num, tptr, 3);
      tptr += 3;
      }
   if( flags & 4)
      {
      FMEMCPY( &star->hipp_num, tptr, 3);
      tptr += 3;
      }

   if( (flags & 128) || ( star->hipp_num & 0x800000))   /* variable flag */
      {
      FSTRCPY( star->var_desig, tptr);
      tptr += FSTRLEN( tptr) + 1;
      star->hipp_num &= 0xfffff;
      }

   for( i = 0; i < ((version_number >= 7) ? 4 : 2); i++)
      {
      long oval;

      if( star->hipp_num)
         {
         if( *(short FAR *)tptr == 32767)
            {
            oval = *(long FAR *)( tptr + 2);
            tptr += 6;
            }
         else
            {
            oval = *(short FAR *)tptr;
            tptr += 2;
            }
         }
      else
         {
         if( *tptr < 127 && *tptr > -127)
            oval = (long)*tptr++;
         else if( *tptr == 127)
            {
            oval = *(short FAR *)( tptr + 1);
            tptr += 3;
            }
         else
            {
            oval = *(long FAR *)( tptr + 1);
            tptr += 5;
            }
         if( i < 2)
            oval *= 100;      /* cvt to .01 milliarcsec/yr */
         }
      if( i < 2)
         star->pm_vals[i] = oval;
      else
         star->pm_err[i - 2] = (short)oval;
      }

   if( flags & 8)
      {
      FMEMCPY( &star->hd_num, tptr, 3);
      tptr += 3;
      }
   if( flags & 16)
      {
      star->yale = *(short FAR *)tptr;
      star->flamsteed = tptr[2];
      star->bayer = tptr[3];
      star->bayer_suffix = tptr[4];
      star->constell = tptr[5];
      tptr += 6;
      }
   if( flags & 32)
      {
      FMEMCPY( &star->sao_num, tptr, 3);
      tptr += 3;
      }
   if( flags & 64)
      {
      star->spectrum[0] = (char)( *tptr++ & 0x7f);
      if( tptr[-1] & 0x80)
         star->spectrum[1] = *tptr++;
      }
                  /* Workaround for %!@$ Xi Sco: */
   if( star->gsc_zone == 5619 && star->gsc_num == 1257)
      {
      star->hipp_num = 78727L;
      star->bayer = 14;
      star->constell = 73;
      }
                  /* Workaround for %!@$ Xi UMa: */
   if( star->gsc_zone == 2520 && star->gsc_num == 2634)
      {
      star->hipp_num = 55203L;
      star->bayer = 14;
      star->constell = 83;
      star->flamsteed = 53;
//    star->sao_num = 62484L;
      if( star->gsc_id[0] == '1')
         {
         star->bayer_suffix = 1;
         star->yale = 4375;
         star->hd_num = 98231L;
         }
      else
         {
         star->bayer_suffix = 2;
         star->yale = 4374;
         star->hd_num = 98230L;
         }
      }
   return( tptr - data_ptr);
}

#ifdef TEST_PROGRAM

#include <stdio.h>
#include <stdlib.h>

void main( int argc, char **argv)
{
   FILE *ifile;
   char filename[80], FAR *buff;
   unsigned loc, bytes_read, buff_size;
   unsigned gsc_zone = atoi( argv[2]);
   unsigned gsc_num  = atoi( argv[3]);
   long offsets[2];

   sprintf( filename, "%s:\\hipp\\lg_tycho.lmp", argv[1]);
   ifile = fopen( filename, "rb");
   if( !ifile)
      {
      printf( "%s not opened\n", filename);
      exit( 0);
      }

   fseek( ifile, 4L * ((long)gsc_zone - 1L), SEEK_SET);
   fread( offsets, 2, sizeof( long), ifile);
   fseek( ifile, offsets[0], SEEK_SET);
   buff_size = (unsigned)( offsets[1] - offsets[0]);
   buff = (char FAR *)malloc( buff_size);
   fread( buff, buff_size, 1, ifile);
   fclose( ifile);
   for( loc = 0; loc < buff_size; loc += bytes_read)
      {
      TYCHO_STAR star;

      bytes_read = parse_tycho_star( buff + loc, &star, 8);
      if( !gsc_num || star.gsc_num == (int)gsc_num)
         {
         short packed_mag_err = ((short *)( buff + loc))[2];

         printf( "Star GSC %u %u%s found\n", gsc_zone, star.gsc_num, star.gsc_id);
         printf( "RA: %.7lf degrees\n", (double)star.ra / 1.e+7);
         printf( "dec: %.7lf degrees\n", (double)star.spd / 1.e+7 - 90.);
         printf( "Johnson V: %d.%03d\n", star.mag / 1000, star.mag % 1000);
         if( star.ppm_num)
            printf( "PPM %ld\n", star.ppm_num);
         if( star.hipp_num)
            printf( "HIP %ld\n", star.hipp_num);
         if( star.sao_num)
            printf( "SAO %ld\n", star.sao_num);
         if( star.hd_num)
            printf( "HD %ld\n", star.hd_num);
         if( star.yale)
            printf( "Yale (HR) %d\n", star.yale);
         printf( "Proper motion in RA: %.2lf +/- %d milliarcseconds/year\n",
               star.pm_vals[0] / 100., star.pm_err[0]);
         printf( "Proper motion in dec: %.2lf +/- %d milliarcseconds/year\n",
               star.pm_vals[1] / 100., star.pm_err[1]);
         printf( "Parallax: %d +/- %d milliarcseconds\n",
               star.parallax, star.parallax_err);
         if( star.spectrum[0])
            printf( "Spectral type %s\n", star.spectrum);
         if( star.var_desig[0])
            printf( "Variable desig %s\n", star.var_desig);
         if( star.bayer || star.bayer_suffix || star.constell)
            printf( "Bayer %d %d %d\n",
                             star.bayer, star.bayer_suffix, star.constell);
         printf( "Epoch RA: %d        Epoch dec: %d\n",
                     star.epoch_ra + 1900, star.epoch_dec + 1900);
         printf( "BT mag: %d.%03d +/- %d.%03d     VT mag: %d.%03d +/- %d.%03d \n",
                     star.bt_mag / 1000, star.bt_mag % 1000,
                     star.bt_mag_err / 1000, star.bt_mag_err % 1000,
                     star.vt_mag / 1000, star.vt_mag % 1000,
                     star.vt_mag_err / 1000, star.vt_mag_err % 1000);
         printf( "Packed mag err: %d (%d %d)\n", packed_mag_err,
                        packed_mag_err >> 6, packed_mag_err & 63);
         }
      }
   free( buff);
}
#endif
