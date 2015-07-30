#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "watdefs.h"
#include "comets.h"
#include "date.h"
#include "ppm.h"
#include "afuncs.h"

#ifdef DEBUG_MEM
#include "checkmem.h"
#endif

#define N_GCVS_RECS 31992
#define PI 3.14159265358979323
#define GAUSS_K .01720209895
#define SOLAR_GM (GAUSS_K * GAUSS_K)
#define EARTH_MAJOR_AXIS 6378140.
#define EARTH_MINOR_AXIS 6356755.
#define EARTH_AXIS_RATIO (EARTH_MINOR_AXIS / EARTH_MAJOR_AXIS)
#define EARTH_MAJOR_AXIS_IN_AU (EARTH_MAJOR_AXIS / (1000. * AU_IN_KM))
#define EARTH_MINOR_AXIS_IN_AU (EARTH_MINOR_AXIS / (1000. * AU_IN_KM))

extern FILE *log_file;
char obj_id_for_report[20];
char *constell_names =
"AndAntApsAqrAqlAraAriAurBooCaeCamCncCVnCMaCMiCapCarCasCenCepCetChaCirColCom\
CrACrBCrvCrtCruCygDelDorDraEquEriForGemGruHerHorHyaHyiIndLacLeoLMiLepLibLup\
LynLyrMenMicMonMusNorOctOphOriPavPegPerPhePicPscPsAPupPyxRetSgeSgrScoSclSct\
SerSexTauTelTriTrATucUMaUMiVelVirVolVul";

static int messier_to_ngc[109] = {
        1952, 7089, 5272, 6121, 5904, 6405, 6475, 6523,
        6333, 6254, 6705, 6218, 6205, 6402, 7078, 6611,
        6618, 6613, 6273, 6514, 6531, 6656, 6494, 6603,
           0, 6694, 6853, 6626, 6913, 7099,  224,  221,
         598, 1039, 2168, 1960, 2099, 1912, 7092,    0,
        2287, 1976, 1982, 2632,    0, 2437, 2422, 2548,
        4472, 2323, 5194, 7654, 5024, 6715, 6809, 6779,
        6720, 4579, 4621, 4649, 4303, 6266, 5055, 4826,
        3623, 3627, 2682, 4590, 6637, 6681, 6838, 6981,
        6994,  628, 6864,  650, 1068, 2068, 1904, 6093,
        3031, 3034, 5236, 4374, 4382, 4406, 4486, 4501,
        4552, 4569, 4571, 6341, 2447, 4736, 3351, 3368,
        3587, 4192, 4254, 4321, 5457, 5866,  581, 4594,
        3379, 4258, 6171, 3556, 3992 };

int get_earth_loc( const double t_millenia, double *results);
static int grab_elements( const char *cd_drive_path, long epoch,
                         uint32_t *elems, int ast);
long grab_ppm_data( const char *cd_drive_path, PPM_STAR *ppm_star, long ppm_num);
int get_sao_info( const char *cd_drive_path, long sao_number, char *buff);
int calc_asteroid_posn( double ut, ELEMENTS *class_elem, double *ra_dec,
                                double *latlon);
int find_gsc_star_loc( const char *cd_drive_path, int zone, int number,
                               double *ra_dec);          /* grab_gsc.c */
int log_printf( FILE *log_file, const char *format, ...);   /* charon.cpp */
int parse_ra_dec( const char *istr, double *ra, double *dec);  /* miscell.c */
FILE *local_then_cd_fopen( const char *cd_drive_path, const char *filename);
void get_guide_environment_ptr( char *buff, const char *env_var);

FILE *local_then_cd_fopen( const char *cd_drive_path, const char *filename)
{
   FILE *rval;
   int i;

   for( i = strlen( filename); i && filename[i - 1] != '\\'; i--)
      ;
   rval = fopen( filename + i, "rb");     /* try in current directory first */
   if( !rval && i)
      rval = fopen( filename, "rb");      /* then in hard-drive install */
   while( !rval && *cd_drive_path)        /* then in path */
      {
      char fullname[100];

      for( i = 0; cd_drive_path[i] && cd_drive_path[i] != ';'; i++)
         ;
      memcpy( fullname, cd_drive_path, i);
      cd_drive_path += i;
      if( *cd_drive_path == ';')
         cd_drive_path++;
      if( fullname[i - 1] != '\\')
         fullname[i++] = '\\';
      strcpy( fullname + i, filename);
      rval = fopen( fullname, "rb");
      }
   return( rval);
}

static int grab_elements( const char *cd_drive_path, long epoch,
                                       uint32_t *elems, int ast)
{
   char filename[80];
   FILE *ifile;
   int rval = 0;

   sprintf( filename, "asteroid\\asteroid\\%ld\\%ld.ast",
                      epoch / 1000L, epoch);
   ifile = local_then_cd_fopen( cd_drive_path, filename);
   if( !ifile)
      return( -1);
   if( fseek( ifile, 6L + (long)(ast - 1) * 24L, SEEK_SET))
      rval = -2;
   if( fread( elems, sizeof( uint32_t), 6, ifile) != 6)
      rval = -3;
   fclose( ifile);
   return( rval);
}

int find_variable( long var_no, const char *cd_drive_path, double *ra_dec)
{
   char tbuff[190];
   FILE *ifile;
   long loc = 0, step = 32768L, n_recs;
   int found_it = 0, reclen, using_gcvs4 = 0;

   sprintf( tbuff, "variable\\gcvs3.dat", cd_drive_path);
   ifile = local_then_cd_fopen( cd_drive_path, "variable\\gcvs3.dat");
   if( !ifile)
      {
      ifile = local_then_cd_fopen( cd_drive_path, "variable\\gcvs4.dat");
      using_gcvs4 = 1;
      }
   if( !ifile)
      return( -1);
   fgets( tbuff, 190, ifile);
   reclen = strlen( tbuff);
               /* 1 Nov 2000:  more recent GCVS dropped the 3-digit */
               /* variable star constellation ID #                  */
   if( reclen == 180)
      using_gcvs4 = 0;
   fseek( ifile, 0L, SEEK_END);
   n_recs = ftell( ifile) / (long)reclen;
   while( step && !found_it)
      {
      if( loc + step < n_recs)
         {
         fseek( ifile, (long)reclen * (long)( loc + step) + (long)using_gcvs4,
                                             SEEK_SET);
         fread( tbuff, 40, 1, ifile);
         tbuff[6] = '\0';
         if( atol( tbuff) <= var_no)
            loc += step;
         if( atol( tbuff) == var_no)
            found_it = 1;
         }
      step >>= 1;
      }
   fclose( ifile);
   tbuff[39] = '\0';
   tbuff[6] = ' ';
   printf( "%s\n", tbuff);
   if( using_gcvs4)
      memmove( tbuff + 16, tbuff + 18, 21);
   if( found_it)
      if( reclen == 180)         /* 1 Nov 2000:  recent GCVS w/full precision */
         {
         long ra = atol( tbuff + 20);
         long dec = atol( tbuff + 29);

         ra_dec[0] = (double)( ra / 10000L) * 15. +
                  (double)( (ra % 10000L) / 100L) / 4. +
                  (double)( ra % 100) / 240.;
         if( tbuff[26] == '.')                     /* tenths of RA sec */
            ra_dec[0] += atof( tbuff + 26) / 240.;
         ra_dec[1] = (double)( dec / 10000L) +
                     (double)( (dec % 10000L) / 100L) / 60. +
                     (double)( dec % 100) / 3600.;
         if( tbuff[28] == '-')
            ra_dec[1] = -ra_dec[1];

         return( 0);
         }
      else        /* older GCVS with lousy precision */
         {
         long ra = atol( tbuff + 16);
         long dec = atol( tbuff + 25);

         ra_dec[0] = (double)( ra / 10000L) * 15. +
                  (double)( (ra % 10000L) / 100L) / 4. +
                  (double)( ra % 100) / 240.;
         if( tbuff[29] == '.')
            ra_dec[1] = (double)( dec / 100L) +
                     (double)( dec % 100L) / 60.;
         else if( tbuff[29] >= '0' && tbuff[29] <= '9')
            ra_dec[1] = (double)( dec / 10000L) +
                     (double)( (dec % 10000L) / 100L) / 60. +
                     (double)( dec % 100) / 3600.;
         printf( "%c: %ld\n", tbuff[29], dec);
         if( tbuff[24] == '-')
            ra_dec[1] = -ra_dec[1];
         return( 0);
         }
   return( -2);      /* variable not found */
}

static int get_ngc_or_ic( const int is_ngc, int cat_no,
                    const char *cd_drive_path, double * const ra_dec)
{
   char buff[80];
   FILE *ifile;

   if( cat_no < 1)
      return( -1);
   if( is_ngc)
      {
      if( cat_no > 7840)
         return( -1);
      cat_no += 5386;
      }
   else        /* is IC */
      if( cat_no > 5386)
         return( -1);
   ifile = local_then_cd_fopen( cd_drive_path, "ngc2000\\ngc2000.dat");
   if( !ifile)
      return( -2);
   fseek( ifile, 102L * (long)( cat_no - 1), SEEK_SET);
   fread( buff, 50, 1, ifile);
   fclose( ifile);
   ra_dec[0] = atof( buff + 13) + atof( buff + 16) / 60.;
   ra_dec[1] = atof( buff + 23) + atof( buff + 26) / 60.;
   ra_dec[0] *= 15;
   if( buff[22] == '-')
      ra_dec[1] *= -1.;
   return( 0);
}

static void precess_1950_to_2000( double *ra_dec)
{
   double matrix[9], old_posn[2], new_posn[2];

   old_posn[0] =  ra_dec[0] * PI / 180.;
   old_posn[1] =  ra_dec[1] * PI / 180.;
   setup_precession( matrix, 1950., 2000.);
   precess_ra_dec( matrix, new_posn, old_posn, 0);
   ra_dec[0] =  new_posn[0] * 180. / PI;
   ra_dec[1] =  new_posn[1] * 180. / PI;
}

int calc_asteroid_posn( double ut, ELEMENTS *class_elem, double *ra_dec,
                                double *latlon)
{
   double td = ut + td_minus_ut( ut) / 86400.;
   double earth_loc[6], target_loc[4], dist = 0., x, y, z;
   const double sin_obliq_2000 = .397777156;
   const double cos_obliq_2000 = .917482062;
   int i, j;

   get_earth_loc( (td - 2451545.0) / 365250., earth_loc);
   for( i = 0; i < 2; i++)
      {
      if( class_elem)
         comet_posn( class_elem, td - dist / AU_PER_DAY, target_loc);
      else
         for( j = 0; j < 4; j++)
            target_loc[j] = 0.;
      dist = 0.;
      for( j = 0; j < 3; j++)
         {
         target_loc[j] -= earth_loc[j];
         dist += target_loc[j] * target_loc[j];
         }
      dist = sqrt( dist);
      }

   x = target_loc[0];
   y = target_loc[1] * cos_obliq_2000 - target_loc[2] * sin_obliq_2000;
   z = target_loc[1] * sin_obliq_2000 + target_loc[2] * cos_obliq_2000;

   if( latlon && latlon[1] < 2.)
      {
      double hr_ang = green_sidereal_time( ut) + latlon[0];
      double u, x0, y0;

      u = atan( EARTH_AXIS_RATIO * sin( latlon[1]) / cos(latlon[1]));
      y0 = EARTH_MINOR_AXIS_IN_AU * sin( u);
      x0 = EARTH_MAJOR_AXIS_IN_AU * cos( u);

      x -= cos( hr_ang) * x0;
      y -= sin( hr_ang) * x0;
      z -=                y0;
      dist = sqrt( x * x + y * y + z * z);
      }

   ra_dec[0] = atan2( y, x) * 180. / PI;
   if( ra_dec[0] < 0.) ra_dec[0] += 360.;
   ra_dec[1] = asin( z / dist) * 180. / PI;
   return( 0);
}

#define HUND_MILLION (1.e+8)

static uint32_t grab_elem( const char *buff)
{
   uint32_t rval = 0;

   while( *buff == ' ')
      buff++;
   while( *buff != ' ')
      {
      if( *buff >= '0' && *buff <= '9')
         rval = rval * 10 + (uint32_t)( *buff - '0');
      else
         if( *buff != '.')
            return( (uint32_t)-1);
      buff++;
      }
   return( rval);
}

static long extract_mpc_epoch( const char *buff)
{
   long year;
   int monthday[2], i;

   year = 1900L + (buff[0] - 'J') * 100L + buff[1] * 10L + buff[2] - 11 * '0';
   for( i = 0; i < 2; i++)
      {
      monthday[i] = buff[3 + i] - '0';
      if( monthday[i] > 10)
         monthday[i] += 10 - 'A' + '0';
      }
   return( dmy_to_day( monthday[1], monthday[0], year, 0) - 1L);
}

static int setup_elems_from_mpcorb( uint32_t *elems, const char *buff)
{
   double semimaj;

   elems[0] = grab_elem( buff +  26) * 100L;  /* mean anom */
   elems[3] = grab_elem( buff +  48) * 100L;  /* asc node */
   elems[4] = grab_elem( buff +  37) * 100L;  /* arg per */
   elems[5] = grab_elem( buff +  59) * 100L;  /* incl */
   elems[2] = grab_elem( buff +  70) * 10L;   /* eccentric */
   semimaj = atof( buff + 93);
   if( semimaj > 21.)
      semimaj = 21. + (semimaj - 21.) / 4.;
   if( semimaj < 42.9)
      elems[1] = (uint32_t)( semimaj * HUND_MILLION);
   else
      elems[1] = 0L;
   return( 0);
}

static int extended_hex_char( const int ival)
{
   int rval;

   if( ival < 10)
      rval = ival + '0';
   else if( ival < 10 + 26)
      rval = ival + 'A' - 10;
   else
      rval = ival + 'a' - 36;
   return( rval);
}

static int setup_obj_id_for_report( char *obj_id_for_report,
                              const char *obj_name)
{
   int i, j, rval = 0;
   char buff[20];

   for( i = j = 0; obj_name[i] && j < 19; i++)  /* create a version of the name */
      if( obj_name[i] != ' ')                   /* with all spaces eliminated */
         buff[j++] = obj_name[i];
   buff[j] = '\0';

   for( i = 0; isdigit( buff[i]); i++)
      ;
   if( buff[i] && i == 4)         /* provisional desig */
      {
      int year = atoi( buff);
      int sub_designator;

      memset( obj_id_for_report, ' ', 12);
      for( i = 0; i < 4; i++)
         {
         const char *surveys[4] = { "P-L", "T-1", "T-2", "T-3" };

         if( !stricmp( buff + 4, surveys[i]))
            {
            memcpy( obj_id_for_report + 8, buff, 4);
            memcpy( obj_id_for_report + 5, "PLST1ST2ST3S" + i * 3, 3);
            return( rval);
            }
         }

      sprintf( obj_id_for_report + 5, "%c%02d",
                  'A' - 10 + year / 100, year % 100);
      obj_id_for_report[6] = buff[2];    /* decade */
      obj_id_for_report[7] = buff[3];    /* year */

      i = 4;
      while( buff[i] == ' ')
         i++;                 /* skip spaces separating year & half-desig */

      obj_id_for_report[8] = toupper( buff[i]);    /* prelim desig is */
      i++;
      if( isalpha( buff[i]))                       /* a little scrambled */
         {
         obj_id_for_report[11] = toupper( buff[i]);
         i++;
         }
      else
         obj_id_for_report[11] = '0';

      sub_designator = atoi( buff + i);
      obj_id_for_report[10] = (char)( '0' + sub_designator % 10);
      if( sub_designator < 100)
         obj_id_for_report[9] = (char)( '0' + sub_designator / 10);
      else if( sub_designator < 360)
         obj_id_for_report[9] = (char)( 'A' + sub_designator / 10 - 10);
      else
         obj_id_for_report[9] = (char)( 'a' + sub_designator / 10 - 36);
      }
   else if( !buff[i])       /* simple numbered asteroid */
      {
      const int asteroid_number = atoi( buff);

      sprintf( obj_id_for_report, " %04d       ", asteroid_number % 10000);
      *obj_id_for_report = extended_hex_char( asteroid_number / 10000);
      }
   else     /* strange ID that isn't decipherable */
      rval = -1;
   return( rval);
}

static int find_object_from_mpcorblike_file( const char *filename,
                                            const char *obj_id, char *buff)
{
   FILE *ifile = fopen( filename, "rb");
   int i, len, got_it = 0;

   printf( "Looking within '%s': %s\n", filename,
                        (ifile ? "Got it" : "didn't find it"));
   if( !ifile)             /* perhaps a path to MPCORB data is  */
      {                    /* specified in 'guide.dat'?  Let's check: */
      get_guide_environment_ptr( buff, "MPCORB_PATH");

      if( *buff)
         {
         strcat( buff, filename);
         ifile = fopen( buff, "rb");
         printf( "Looking within '%s': %s\n", buff,
                        (ifile ? "Got it" : "didn't find it"));
         }
      }
   if( !ifile)
      return( 0);

   len = strlen( obj_id);
   for( i = 0; isdigit( obj_id[i]); i++)
      ;

   while( !got_it && fgets( buff, 256, ifile))
      if( strlen( buff) > 200 && buff[16] == '.' && buff[138] == '.')
         {
         if( !obj_id[i] && atoi( obj_id))         /* numbered object */
            {
            int number_from_mpcorb = atoi( buff + 1);

            if( *buff >= '0' && *buff <= '9')
               number_from_mpcorb += (*buff - '0') * 10000;
            else if( *buff >= 'A' && *buff <= 'Z')
               number_from_mpcorb += (*buff - 'A' + 10) * 10000;
            else if( *buff >= 'a' && *buff <= 'z')
               number_from_mpcorb += (*buff - 'a' + 36) * 10000;
            if( atoi( obj_id) == number_from_mpcorb)
               got_it = 1;
            }
         else
            {
            char *name_text = buff + 174;          /* Name can start in */
                                                   /* column 174 or 175 */
            if( *name_text == ')' || *name_text == ' ')
               name_text++;
            got_it = (!memicmp( name_text, obj_id, len)
                                          && name_text[len] <= ' ');
            }
         }

   if( got_it)
      log_printf( log_file, "Asteroid orbital data found in '%s'\n",
                                    filename);
   fclose( ifile);
   return( got_it);
}

int find_obj( double jd, char *obj_id, const char *cd_drive_path, double *ra_dec,
                              double *latlon)
{
   char buff[256];
   ELEMENTS class_elem;
   int i, j, target_no = -1, obj_type, rval = 0, got_it = 0;
   long obj_num = atol( obj_id + 3);
   FILE *ifile;
   extern char *strings[];

   while( *obj_id == ' ')                 /* skip any leading spaces */
      obj_id++;

   if( !parse_ra_dec( obj_id, ra_dec, ra_dec + 1))
      return( rval);

   if( !strcmp( obj_id, "sun"))
      return( calc_asteroid_posn( jd, NULL, ra_dec, latlon));

   if( *cd_drive_path == '*')    /* if there's no Guide CD... */
      if( *obj_id != 'A' && *obj_id != 'a' && *obj_id != 'C' && *obj_id != 'c')
         return( -99);     /* we can only handle asteroids and comets */

   if( !memicmp( obj_id, "ppm", 3))
      {
      PPM_STAR ppm_star;

      if( grab_ppm_data( cd_drive_path, &ppm_star, obj_num) != -1L)
         {
         ra_dec[0] = (double)ppm_star.ra * 15. / (3600. * 1000.);
         ra_dec[1] = (double)ppm_star.dec / (3600. * 100.);
         }
      else
         rval = -1;
      return( rval);
      }

   if( !memicmp( obj_id, "sao", 3))
      {
      long sao_data[9];

      if( !get_sao_info( cd_drive_path, obj_num, (char *)sao_data))
         for( i = 0; i < 2; i++)
            ra_dec[i] = (180. / PI) * (double)sao_data[i] /
                                         (double)(1L << 28);
      else
         rval = -1;
      return( rval);
      }

   if( !memicmp( obj_id, "gsc", 3))
      {
      int zone, number;

      sscanf( obj_id + 3, "%d %d", &zone, &number);
      return( find_gsc_star_loc( cd_drive_path, zone, number, ra_dec));
      }

   if( !memicmp( obj_id, "sn", 2))
      {
      obj_id += 2;
      ifile = local_then_cd_fopen( cd_drive_path, "text\\supernov.nam");
      if( !ifile)
         return( -1);
      while( fgets( buff, sizeof( buff), ifile))
         if( !memcmp( buff, obj_id, strlen( obj_id)))
            {
            double hr, min, sec = 0., deg;

            fclose( ifile);
            buff[44] = buff[57] = '\0';
            sscanf( buff + 32, "%lf %lf %lf", &hr, &min, &sec);
            ra_dec[0] = (double)hr * 15. + min / 4. + sec / 240.;
            sec = 0.;
            sscanf( buff + 46, "%lf %lf %lf", &deg, &min, &sec);
            ra_dec[1] = (double)deg + min / 60. + sec / 3600.;
            if( buff[45] == '-')
               ra_dec[1] = -ra_dec[1];
            if( buff[66] != '2')
               precess_1950_to_2000( ra_dec);
            return( 0);
            }
      fclose( ifile);
      return( -1);
      }

   if( !memicmp( obj_id, "ngc", 3))
      return( get_ngc_or_ic( 1, atoi( obj_id + 3), cd_drive_path, ra_dec));

   if( !memicmp( obj_id, "ic", 2))
      return( get_ngc_or_ic( 0, atoi( obj_id + 2), cd_drive_path, ra_dec));

   if( *obj_id == 'm')
      if( obj_id[1] == ' ' || isdigit( obj_id[1]))
         {
         int mess_no = atoi( obj_id + 1);

         if( mess_no > 0 && mess_no < 110)
            if( messier_to_ngc[mess_no - 1])
               return( get_ngc_or_ic( 1, messier_to_ngc[mess_no - 1],
                                          cd_drive_path, ra_dec));
         return( -1);
         }
   if( !memicmp( obj_id, "nsv", 3))
      {
      int deg, hr;
      double min, seconds;

      if( obj_num < 1L || obj_num > 14811L)
         return( -1);
      ifile = local_then_cd_fopen( cd_drive_path, "variable\\nsv3.dat");
      if( !ifile)
         return( -1);
      if( obj_num > 10360)
         obj_num++;
      fseek( ifile, 97L * (long)( obj_num - 1), SEEK_SET);
      fread( buff, 70, 1, ifile);
      fclose( ifile);
      buff[70] = '\0';
      printf( "%s\n", buff);
      sscanf( buff + 8, "%d %lf %lf", &hr, &min, &seconds);
      ra_dec[0] = (double)hr * 15. + min / 4. + seconds / 240.;
      sscanf( buff + 19, "%d %lf", &deg, &min);
      ra_dec[1] = (double)deg + min / 60.;
      if( buff[18] == '-')
         ra_dec[1] = -ra_dec[1];
      precess_1950_to_2000( ra_dec);
      return( 0);
      }

   switch( obj_type = *obj_id++)
      {
      case 'A':        /* asteroid from built-in data or MPCORB-like file */
      case 'a':
         {
         int got_it;
         uint32_t l_elems[6];
         long l_epoch;

         while( *obj_id == ' ')              /* skip leading spaces */
            obj_id++;
         got_it = find_object_from_mpcorblike_file( "neatod.dat",
                              obj_id, buff);
         if( !got_it)
            got_it = find_object_from_mpcorblike_file( "mpcorb.dat",
                              obj_id, buff);
         if( !got_it)
            got_it = find_object_from_mpcorblike_file( "mpcorbcr.dat",
                              obj_id, buff);

         if( got_it)       /* asteroid found in an MPCORB-like file: */
            {
            setup_elems_from_mpcorb( l_elems, buff);
            l_epoch = extract_mpc_epoch( buff + 20);
            }
         else       /* didn't get it from MPCORB:  look to hard drive: */
            {
            for( i = 0; isdigit( obj_id[i]); i++)
               ;
            if( !obj_id[i])       /* it's a numbered object */
               target_no = atoi( obj_id);
            else                 /* it's an unnumbered object;  we gotta */
               {                 /* find out where it is in the file     */
               ifile = local_then_cd_fopen( cd_drive_path, "asteroid\\astnames.dat");
               if( ifile)
                  {
                  int len;

                  fgets( buff, sizeof( buff), ifile);
                  len = strlen( buff);
                  for( i = 2; fread( buff, len, 1, ifile) && target_no < 0; i++)
                     {
                     for( j = len; j && buff[j - 1] <= ' '; j--)
                     ;
                     buff[j] = '\0';
                     if( !stricmp( buff + 12, obj_id))
                        target_no = i;
                     }
                  fclose( ifile);
                  }
               }
            if( target_no <= 0)
               return( -2);
            else
               {
               l_epoch = (((long)jd + 100L) / 200L) * 200L;
               if( grab_elements( cd_drive_path, l_epoch, l_elems, target_no))
                  return( -3);
               got_it = 1;
               }
            }

         if( got_it)       /* we got elements,  either from MPCORBlike files */
            {              /* or Guide's built-in asteroid data:             */
            setup_obj_id_for_report( obj_id_for_report, obj_id);
            setup_elems_from_ast_file( &class_elem,
                      l_elems, (double)l_epoch + .5);
            log_printf( log_file, "epoch = %ld.5; '%s'\n", l_epoch,
                                 obj_id_for_report);
            }
         }
         break;
      case 'c': case 'C':
         while( *obj_id == ' ')              /* skip leading spaces */
            obj_id++;
         class_elem.perih_time = 0.;       /* select "worst" case elements */
         for( i = 0; i < 2 && target_no; i++)
            {
            char extracted_name[50];

            extracted_name[0] = '\0';
            if( !i)
               strcpy( buff, "comets.dat");
            else
               strcpy( buff, "asteroid\\cometg.dat");
            if( *cd_drive_path == '*')
               {
               strcpy( buff, "soft02cm.txt");
               if( i)
                  *buff = '\0';
               }
            if( *buff)
               ifile = local_then_cd_fopen( cd_drive_path, buff);
            else
               ifile = NULL;
            while( ifile && fgets( buff, sizeof( buff), ifile))
               {
               char *tptr;

               buff[41] = '\0';
               got_it = 0;
               tptr = strstr( buff, obj_id);
               if( tptr)
                  {
                  int end_char = tptr[strlen( obj_id)];

                  if( tptr > buff && tptr[-1] == '(' && end_char == ')')
                     got_it = 1;
                  if( end_char == ' ')
                     got_it = 1;
                  }
               if( got_it)
                  {
                  int year, month;
                  double day, mean_anom, epoch;

                  sscanf( buff + 42, "%lf %d %d", &day, &month, &year);
                  day += (double)dmy_to_day( 0, month, year, 0) - .5;
                  if( fabs( day - jd) < fabs( class_elem.perih_time - jd))
                     {
                     int j;

                     class_elem.perih_time = day;
                     for( j = 60; buff[j] != '.' && j < 72; j++)
                        ;
                     if( j != 72)
                        mean_anom = atof( buff + 60);
                     else                 /* reference stored in bytes  */
                        mean_anom = 0.;   /* 60-71;  assume MA = 0      */
                        /* normally goes,  MA=0,  i.e.,  comet-style elems */

                     sscanf( buff + 72,
                            "%lf %lf %lf %lf %lf %lf %lf %lf",
                            &class_elem.q, &class_elem.ecc,
                            &class_elem.incl, &class_elem.arg_per,
                            &class_elem.asc_node, &epoch,
                            &class_elem.abs_mag, &class_elem.slope_param);
                     class_elem.is_asteroid = (buff[152] == 'A');
                     class_elem.incl *= PI / 180.;
                     class_elem.arg_per *= PI / 180.;
                     class_elem.asc_node *= PI / 180.;
                     derive_quantities( &class_elem, SOLAR_GM);
                     class_elem.perih_time -= (mean_anom / 360.) *
                                                ( class_elem.t0 * 2. * PI);
                     strcpy( extracted_name, buff);
                     target_no = 0;
                     }
                  }
               }
            if( ifile)
               fclose( ifile);
            strcpy( buff, extracted_name);
            }
         if( target_no)       /* no object found */
            return( -4);

         for( i = got_it = 0; !got_it && buff[i]; i++)
            if( buff[i] == '(')
               {
               int periodic_num = atoi( buff + i + 1), j;

               for( j = i + 1; buff[j] >= '0' && buff[j] <= '9'; j++)
                  ;
               if( periodic_num && buff[j] == 'P' && buff[j + 1] == ')')
                  {
                  got_it = 1;
                  sprintf( obj_id_for_report, "%04dP       ", periodic_num);
                  }
               }

         for( i = 0; !got_it && buff[i]; i++)
            if( !memcmp( buff + i, "19", 2) || !memcmp( buff + i, "20", 2))
               if( !setup_obj_id_for_report( obj_id_for_report, buff + i))
                  {
                  if( i > 1 && buff[i-1] == '/')
                     if( buff[i-2] == 'C' || buff[i-2] == 'P')
                     obj_id_for_report[4] = buff[i - 2];
                  got_it = 1;
                  }
         printf( strings[105], class_elem.perih_time);
                     /* "Time of perihelion %.2lf\n" */
         break;
      case 'v': case 'V':        /* variable */
         {
         long id_num;

         while( *obj_id == ' ')              /* skip leading spaces */
            obj_id++;
         id_num = decipher_var_desig( obj_id);
         printf( "id num %ld\n", id_num);
         if( id_num > 0L)
            {
            while( *obj_id > ' ')
               obj_id++;
            for( i = 0; i < 88; i++)
               if( !memicmp( obj_id + 1, constell_names + i * 3, 3))
                  {
                  id_num += 10000L * (long)(i + 1);
                  printf( "id num %ld\n", id_num);
                  if( !find_variable( id_num, cd_drive_path, ra_dec))
                     {
                     precess_1950_to_2000( ra_dec);
                     return( 0);
                     }
                  }
            }
         return( -1);         /* no such luck */
         }
         break;
      default:
         return( -1);
         break;
      }

   rval = calc_asteroid_posn( jd, &class_elem, ra_dec, latlon);
   if( rval)
      log_printf( log_file, strings[106], rval);
                  /* "ERROR in calc_asteroid_posn: %d\n" */
   else        /* movement over .01 day: */
      calc_asteroid_posn( jd + .01, &class_elem, ra_dec + 2, latlon);
   if( rval)
      log_printf( log_file, strings[106], rval);
                  /* "ERROR in calc_asteroid_posn: %d\n" */
   return( rval);
}
