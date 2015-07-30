/* A10_EXTR.CPP         Copyright 1997,  Project Pluto                  */
/* All rights reserved.  This source code may be used freely for        */
/* non-commercial uses.                                                 */

/* Under MSC,     compile as:  cl -W4 -Ox a10_extr.cpp            */
/* Under Watcom,  compile as:  wcl386 -W4 -Ox a10_extr.cpp        */
/* If all you want is the example DOS program,  add -c -DMAKE_EXE */

/* This includes some routines to figure out which CDs would be needed */
/* to cover a given area,  plus the actual code to do the extraction.  */
/* The choice between A1.0 and A2.0 is controlled via the global int   */
/* using_a20.  If you stick in an SA1.0 or SA2.0 disk,  the code will  */
/* automatically recognize that fact and extract data from it.         */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "usno_a.h"

#define MAX_RA (360L * 3600L * 100L)
#define SPACING (MAX_RA / n_entries)
#define ZONE_SIZE (MAX_RA / 48L)
#define BUFFSIZE 1000
#define EXTRACT struct extract
#define SWAP( A, B, TEMP)   { TEMP = A;  A = B;  B = TEMP; }
#define SA10_DISK   (1 << 14)

static long n_entries = 1440L;

/* Yes,  I realize globals are evil... but: */

int using_a20 = 0;

EXTRACT
   {
   long ra1, ra2, spd1, spd2;
   long start_rec;
   int zone_no, cd_no;
   };

   /* The 24 declination zones of the A1.0 catalog are spread out in a */
   /* somewhat random manner across ten CDs,  as shown in this table:  */

static const char cd_table[24] =
                   { 1, 1, 6, 5, 3, 2,  1, 4, 6, 5, 7, 10,
                     8, 7, 8, 9, 9, 4, 10, 3, 2, 6, 2,  3 };

   /* Those of the A2.0 are spread out as follows,  over eleven CDs: */

static const char cd_table_a20[24] =
          { 1, 1, 9, 7, 5, 4,  3, 2, 1, 6, 7, 10,
            9, 8, 8,11,10,11,  6, 4, 2, 3, 3,  2 };
       /*   0    15    30     45    60    75     */
       /*  90   105   120    135   150   165     */


/* find_extracts() determines what pieces of A2.0 would be needed to
cover the desired RA/dec rectangle.  It handles issues such as the
rectangle crossing 0h or 24h in RA.  The resulting array of 'extracts'
is sorted by CD and,  within CD,  by zone number.  That makes it a little
easier to determine which CDs would be needed to cover a given area.
Of course,  if you're getting SAx.0 data,  or data from DVD or the
hard drive,  the sort order is irrelevant.         */

static int find_extracts( EXTRACT *ret_arr, const long center_ra,
          const long center_spd, const long width, const long height)
{
   EXTRACT total, *tptr = ret_arr;
   int rval, i, j;
   long zone, delta = 0L;

   total.ra1 = center_ra - width / 2L;
   total.ra2 = center_ra + width / 2L;
   total.spd1 = center_spd - height / 2L;
   total.spd2 = center_spd + height / 2L;

   for( zone = 0; zone < 24L; zone++)
      if( zone * ZONE_SIZE <= total.spd2 && (zone + 1L) * ZONE_SIZE > total.spd1)
         {
         memcpy( tptr, &total, sizeof( EXTRACT));
         if( tptr->spd2 > (zone + 1L) * ZONE_SIZE)
            tptr->spd2 = (zone + 1L) * ZONE_SIZE;
         if( tptr->spd1 < zone * ZONE_SIZE)
            tptr->spd1 = zone * ZONE_SIZE;
         tptr->zone_no = (int)zone;
         tptr->cd_no = (using_a20 ? cd_table_a20[tptr->zone_no] :
                                    cd_table[tptr->zone_no]);
         tptr++;
         }
   rval = tptr - ret_arr;
   if( total.ra1 < 0L)
      delta = MAX_RA;
   if( total.ra2 > MAX_RA)
      delta = -MAX_RA;
   if( delta)              /* wraps around 0h: duplicating time: */
      {
      memcpy( tptr, ret_arr, rval * sizeof( EXTRACT));
      for( i = 0; i < rval; i++, tptr++)
         {
         tptr->ra1 += delta;
         tptr->ra2 += delta;
         }
      rval *= 2;
      tptr = ret_arr;
      for( i = 0; i < rval; i++, tptr++)
         {
         if( tptr->ra1 < 0L)
            tptr->ra1 = 0L;
         if( tptr->ra2 < 0L)
            tptr->ra2 = 0L;
         if( tptr->ra1 > MAX_RA)
            tptr->ra1 = MAX_RA;
         if( tptr->ra2 > MAX_RA)
            tptr->ra2 = MAX_RA;
         }
      }
               /* bubblesort in order of CD # (small array) */
               /* within a CD,  sort by zone # */

   for( i = 0; i < rval; i++)
      for( j = i + 1; j < rval; j++)
         if( ret_arr[j].cd_no < ret_arr[i].cd_no ||
             (ret_arr[j].cd_no == ret_arr[i].cd_no &&
              ret_arr[j].zone_no < ret_arr[i].zone_no))
            {
            total = ret_arr[j];
            ret_arr[j] = ret_arr[i];
            ret_arr[i] = total;
            }

   return( rval);
}

/* Set_file_extents() puzzles out,  for a given EXTRACT,  which stretch
of an Ax.0 zone needs to be examined,  by reading index data from the
sax0.ind or ax0.ind file.  This is the best way to set the extents,
but if the .ind files aren't available,  the code will try setting
extents using the .acc files distributed with Ax.0 (using the function
set_file_extents_using_acc()).  And if there aren't any .acc files,
it simply does a binary search to find out where the data starts. */

static int set_file_extents( int n_extents, EXTRACT *ret_arr)
{
   char *filename = (n_entries == 180L ? "sa10.idx" : "a10.idx");
   FILE *ifile;

   if( using_a20)
      filename = (n_entries == 180L ? "sa20.idx" : "a20.idx");
   ifile = fopen( filename, "rb");
   if( !ifile)
      return( -2);

   while( n_extents--)
      {
      const int start = (int)( ret_arr->ra1 / SPACING);
      const long loc = ((n_entries + 1L) * (long)ret_arr->zone_no + start);

      fseek( ifile, loc * 4L, SEEK_SET);
      fread( &ret_arr->start_rec, sizeof( long), 1, ifile);
      ret_arr++;
      }

   fclose( ifile);
   return( 0);
}

/* For various historical reasons,  this code provides _three_ different
ways to figure out which part of an Ax.0 file covers the desired span of
RA values.  set_file_extents(),  the default,  uses an "accelerator"
index file I put together for use with Guide.  If that file isn't present,
set_file_extents_using_acc() looks for the .ACC file supplied by USNO
for a given zone.  And if _that_ fails,  set_file_extents_by_searching()
does a binary search to narrow down the extents. */

static int set_file_extents_using_acc( EXTRACT *eptr, const char *filename)
{
   char buff[50];
   FILE *ifile;

   strcpy( buff, filename);
   strcpy( buff + strlen( buff) - 3, "acc");
   ifile = fopen( buff, "rb");
   if( !ifile)
      eptr->start_rec = 0L;
   else
      {
      const long n_acc_entries = 96L;
      const long acc_spacing = MAX_RA / n_acc_entries;
      const int start = (int)( eptr->ra1 / acc_spacing);

      fseek( ifile, start * 30L, SEEK_SET);
      fread( buff, 30, 1, ifile);
      eptr->start_rec = atol( buff + 5) - 1L;

      fclose( ifile);
      }
   return( ifile == NULL);
}

static void set_file_extents_by_searching( EXTRACT *eptr, FILE *ifile)
{
   long total_recs;
   long loc = 0, step, loc1;

   fseek( ifile, 0L, SEEK_END);
   total_recs = ftell( ifile) / 12L;

   for( step = 0x40000000; step; step >>= 1)
      if( (loc1 = step + loc) < total_recs)
         {
         char tbuff[4], tchar;

         fseek( ifile, loc1 * 12L, SEEK_SET);
         fread( tbuff, 1, sizeof( long), ifile);
         tchar = tbuff[0];
         tbuff[0] = tbuff[3];
         tbuff[3] = tchar;
         tchar = tbuff[1];
         tbuff[1] = tbuff[2];
         tbuff[2] = tchar;
         if( *(long *)tbuff < eptr->ra1)
            loc = loc1;
         }
   eptr->start_rec = loc;
}

/* When Ax.0 or SAx.0 stars are read in,  they're read as three long
integers.  The catalog number is implicit in the star's location in the
file.  But when written out to 'ofile',  you'll skip over some stars,
and the catalog number must be explicitly stored.  That's why the
following code writes out _four_ longs per star:  the first is the
catalog number,  the rest are the "usual" three long ints,  unchanged
except for the byte order being fixed.

   That first long int 'id' provides the version (A1.0, A2.0, SA1.0,
or SA2.0),  the zone (0-23,  from south to north pole),  and offset
within the zone.  All that data is packed into 'id' as follows:

   For A1.0,  id = offset +  zone       * 50 million
   For A2.0,  id = offset + (zone + 24) * 50 million
   For SA1.0, id = offset +  zone       *  4 million + 1.2 billion
   For SA2.0, id = offset + (zone + 24) *  4 million + 1.2 billion

   All this takes advantage of the fact that,  for Ax.0,  no zone
has more than 50 million stars;  for SAx.0,  no star has more than
four million stars.

   The return value of extract_ax0_data( ) is a bitmask indicating
which CDs cover the desired area.  Use is made of this in the
find_needed_ax0_cds() function below,  to figure out what CD(s) to
prompt the user for.     */

static int extract_ax0_data( const char *cd_path, const double ra_degrees,
            const double dec_degrees, const double width_degrees,
            const double ht_degrees, FILE *ofile)
{
   EXTRACT *ret_arr = (EXTRACT *)calloc( 48, sizeof( EXTRACT));
   int i, n_ret, rval = 0, use_acc_files = 0, catalog_offset;
   long *buff;
   char sa10_file[20];
   FILE *sa_test;

   n_ret = find_extracts( ret_arr, (long)( ra_degrees * 360000.),
                                  (long)( (dec_degrees + 90.) * 360000.),
                                   (long)( width_degrees * 360000.),
                                   (long)( ht_degrees * 360000.));
   if( !cd_path)        /* just figuring out which CDs we need */
      {
      for( i = 0; i < n_ret; i++)
         rval |= (1 << ret_arr[i].cd_no);
      free( ret_arr);
      return( rval);
      }

   n_entries = 1440L;

      /* if ZONE0450.CAT is on the CD,  and is smaller than 100MB, */
      /* it's SA 1.0 or 2.0 */

   sprintf( sa10_file, "%s\\zone0450.cat", cd_path);
   sa_test = fopen( sa10_file, "rb");
   if( sa_test)
      {
      fseek( sa_test, 0L, SEEK_END);
      if( ftell( sa_test) < 100000000L)
         {
         n_entries = 180L;
         rval = SA10_DISK;
         }
      fclose( sa_test);
      }

   if( rval == SA10_DISK)
      catalog_offset = (using_a20 ? 2496000000u : 2400000000u);
   else
      catalog_offset = (using_a20 ? 1200000000u : 0u);

   if( set_file_extents( n_ret, ret_arr))
      use_acc_files = 1;

   buff = (long *)calloc( BUFFSIZE, 3 * sizeof( long));
   for( i = 0; i < n_ret; i++)
      {
      FILE *ifile;
      char filename[40];
      EXTRACT *eptr = ret_arr + i;

      sprintf( filename, "%s\\zone%04d.cat", cd_path, eptr->zone_no * 75);
      ifile = fopen( filename, "rb");
      if( ifile)
         {
         long loc;
         int stop_reading = 0;

         if( use_acc_files)
            if( set_file_extents_using_acc( eptr, filename))
               set_file_extents_by_searching( eptr, ifile);
         fseek( ifile, eptr->start_rec * 12L, SEEK_SET);
         loc = eptr->start_rec + catalog_offset +
                  (rval == SA10_DISK ? 4000000 : 50000000) * eptr->zone_no;
         while( !stop_reading && ofile)
            {
            int j, k, n_read;

            n_read = fread( buff, 12, BUFFSIZE, ifile);
            if( !n_read)
               stop_reading = 1;
            for( j = 0; !stop_reading && j < n_read; j++)
               {
               long swapped[4], ra, spd;

               memcpy( swapped + 1, buff + j * 3, 3 * sizeof( long));
               for( k = 1; k < 4; k++)    /* A1.0 data is "wrong-endian" */
                  {
                  char tchar, *swapptr = (char *)( swapped + k);

                  SWAP( swapptr[0], swapptr[3], tchar);
                  SWAP( swapptr[1], swapptr[2], tchar);
                  }
               ra = swapped[1];
               spd = swapped[2];

               if( ra > eptr->ra2)
                  stop_reading = 1;
               else if( spd > eptr->spd1 && spd < eptr->spd2 &&
                                                        ra > eptr->ra1)
                  {
                  swapped[0] = loc + j + 1L;
                  fwrite( swapped, 4, sizeof( long), ofile);
                  }
               }
            loc += n_read;
            }
         if( rval >= 0)
            rval |= (1 << eptr->cd_no);
         fclose( ifile);
         }
      }
   free( buff);
   free( ret_arr);
   return( rval);
}

/* If you call extract_ax0_data with a NULL path and output file,  all
it does is to figure out which CDs are needed;  that info is returned as
a bitmask.  The following convenience function takes advantage of this.
The list of CDs covering a given area is declination-dependent only. */

int find_needed_ax0_cds( const double dec_degrees, const double ht_degrees)
{
   return( extract_ax0_data( NULL, 1., dec_degrees, 1., ht_degrees, NULL));
}

/* extract_ax0_from_multipath() allows the cd_path string to specify
multiple places where Ax.0 data might be found.  For example,  if it
was 'x:\;d:\a2",  the function would attempt to get Ax.0 data from
the root of drive x:.  If that failed,  it would try to find data in
the folder d:\a2.

   Aside from that,  it's indistinguishable from extract_ax0_data(). */

int extract_ax0_from_multipath( const char *cd_path,
            const double ra_degrees, const double dec_degrees,
            const double width_degrees,
            const double ht_degrees, FILE *ofile)
{
   char curr_path[80];
   int i, rval = 0, curr_rval;
   int disks_left = find_needed_ax0_cds( dec_degrees, ht_degrees);

   while( disks_left && *cd_path && rval >= 0)
      {
      for( i = 0; cd_path[i] && cd_path[i] != ';'; i++)
         curr_path[i] = cd_path[i];
      cd_path += i;
      if( *cd_path == ';')
         cd_path++;
      if( i == 1)       /* drive letter only */
         curr_path[i++] = ':';
      if( curr_path[i - 1] == '\\')       /* remove trailing backslash */
         i--;
      curr_path[i] = '\0';
      curr_rval = extract_ax0_data( curr_path, ra_degrees, dec_degrees,
                        width_degrees, ht_degrees, ofile);
      if( curr_rval < 0)
         rval = curr_rval;
      else
         {
         disks_left &= ~curr_rval;
         rval |= curr_rval;
         }
      }
   return( rval);
}

int build_prompt_string( const int needed_mask,
        const char *prompt_string, const char *cd_path, char *output_prompt)
{
   char tbuff[50], needed_cds[13];
   int i, n_needed = 0;

   for( i = 0; i < 12; i++)            /* count the number of disks needed */
      if( (needed_mask >> i) & 1)
         needed_cds[n_needed++] = (char)i;
   *tbuff = '\0';
   if( n_needed > 2)
      {
            /* If multiple CDs will be needed,  list 'em all up front: */
      for( i = 0; i < n_needed - 1; i++)
         sprintf( tbuff + strlen( tbuff), "%d, ", needed_cds[i]);
      sprintf( tbuff + strlen( tbuff), "or %d", needed_cds[i]);
      }
   else if( n_needed == 2)
      sprintf( tbuff, "%d or %d", needed_cds[0], needed_cds[1]);
   else if( n_needed == 1)
      sprintf( tbuff, "%d", needed_cds[0]);
   for( i = 0; prompt_string[i] != '%'; i++)
      output_prompt[i] = prompt_string[i];
   output_prompt[i] = ' ';
   strcpy( output_prompt + i + 1, tbuff);
   sprintf( output_prompt + strlen( output_prompt), prompt_string + i + 2,
                                     cd_path);
   return( 0);
}

#ifdef OBSOLETE_WINDOWS_CODE
#include <windows.h>          /* For MessageBox( ) prototype */

long windows_extract_a10_data( const char *cd_path, const char *path,
            const double ra_degrees,
            const double dec_degrees, const double width_degrees,
            const double ht_degrees, FILE *ofile, const char *prompt_string)
{
   char path_buff[3];
   int needed_mask = find_needed_ax0_cds( dec_degrees, ht_degrees);

   if( path && *path)
      needed_mask &= ~extract_ax0_from_multipath( path, ra_degrees,
                           dec_degrees, width_degrees, ht_degrees, ofile);

   path_buff[0] = (char)cd_letter;
   path_buff[1] = ':';
   path_buff[2] = '\0';

   while( needed_mask)
      {
      int cd_mask = extract_ax0_data( path_buff, ra_degrees, dec_degrees,
                                       width_degrees, ht_degrees, ofile);

      if( cd_mask == SA10_DISK)              /* SA1.0 CD */
         return( 0);
      if( cd_mask == -2)
         return( -2);
      needed_mask &= ~cd_mask;      /* mask out the disks we don't need */
      if( needed_mask)
         {
         char buff[100], title[10];

         build_prompt_string( needed_mask, prompt_string, cd_path,
                                                            buff);
         sprintf( title, "A%d.0", using_a20 + 1);
         if( MessageBox( NULL, buff, title, MB_OKCANCEL) == IDCANCEL)
            return( -3);
         }
      }
   return( 0);
}
#else

#include <conio.h>
#include <dos.h>

long dos_prompt_extract_a10_data( const char *cd_path, const char *path,
            const double ra_degrees,
            const double dec_degrees, const double width_degrees,
            const double ht_degrees, FILE *ofile)
{
   int needed_mask = find_needed_ax0_cds( dec_degrees, ht_degrees);

   if( path && *path)
      needed_mask &= ~extract_ax0_from_multipath( path, ra_degrees,
                           dec_degrees, width_degrees, ht_degrees, ofile);

            /* The following code allowed for situations where Ax.0 was */
            /* _not_ installed to the hard drive.  If the above call    */
            /* looked for Ax.0 data in the 'multipath' and didn't get   */
            /* all the needed zones,  it prompted the user for the      */
            /* disk(s) required.  Very slick,  but no longer necessary, */
            /* and I ran into maintenance issues... it was simplest to  */
            /* drop support for this now very rarely used feature.      */
#ifdef OBSOLETE_CODE
   while( needed_mask)
      {
      char buff[100];
      char *prompt = "Please insert Ax.0 CD number %d in drive %c:";
      int cd_mask;

      prompt[15] = (char)( '1' + using_a20);
      build_prompt_string( needed_mask, prompt, cd_path, buff);
      printf( "%s\n", buff);
      if( getch( ) == 27)
         return( -1);
      else
         {
         struct _find_t c_file;
         char *path = "z:\\*.*";

         *path = (char)cd_letter;
         _dos_findfirst( path, _A_VOLID, &c_file);
         printf( "Volume name: %s\n", c_file.name);
         }

      cd_mask = extract_ax0_data( path_buff, ra_degrees, dec_degrees,
                                       width_degrees, ht_degrees, ofile);

      if( cd_mask == SA10_DISK)              /* SA1.0 CD */
         {
         printf( "SA%d.0 data extracted\n", using_a20 + 1);
         return( 0);
         }
      if( cd_mask == -2)
         {
         printf( "A%d.0/SA%d.0 index not found!\n", using_a20 + 1,
                                                    using_a20 + 1);
         return( -2);
         }
      needed_mask &= ~cd_mask;
      }
#endif
   return( 0);
}
#endif

#ifdef MAKE_EXE

#include <conio.h>

void main( int argc, char **argv)
{
   int cd_letter = argv[5][0];
   FILE *ofile = NULL;
   double center_ra = atof( argv[1]);
   double center_dec = atof( argv[2]);
   double width = atof( argv[3]);
   double height = atof( argv[4]);

   if( argc > 6)
      ofile = fopen( argv[6], "ab");
   dos_prompt_extract_a10_data( cd_letter, center_ra, center_dec,
                                width, height, ofile);
   if( ofile)
      fclose( ofile);
   printf( "Please reinsert the Guide CD and hit any key:");
   getch( );
}
#endif
