#include <malloc.h>
#include <string.h>
#include <search.h>
#include <stdlib.h>
#include <conio.h>
#include <ctype.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <dos.h>
#include <time.h>
#include <io.h>
#include <i86.h>
#include <direct.h>
#include <stdint.h>
#include "alt_defs.h"
#include "watdefs.h"

#define CTRL(X) ((X)-64)
#define SECONDS_MASK 31
#define MINUTE_MASK 63
#define HOUR_MASK 31
#define DATE_MASK 31
#define MONTH_MASK 15
#define YEAR_MASK 127
#define SWAP( A, B, TEMP) { TEMP=A;    A=B;    B=TEMP;   }

#define DIRECT_ENT struct direct_ent

DIRECT_ENT
   {
   long size;
   unsigned date, time, attrib;
   char name[9];
   char ext[4];
   };

// int xscr = 80, yscr = 25, display_mode = 0;
int xscr = 80, yscr = 50, display_mode = 8;
static char path_str[132];
static char sole_extension[10];
char *dispmem = (char *)0xb8000;

#define MINIMAL_N_ALLOCED 10000

int main( int argc, char *argv[]);

char **strings;

#define TEXT_MODE_80x25          0
#define TEXT_MODE_100x37         1
#define TEXT_MODE_80x60          2
#define TEXT_MODE_132x44         3
#define TEXT_MODE_132x50         4
#define TEXT_MODE_132x25         5
#define TEXT_MODE_132x28         6
#define TEXT_MODE_80x28          7
#define TEXT_MODE_80x50          8
#define N_TEXT_MODES             9

int set_mode( int mode, int *xsize, int *ysize)
{
   static const unsigned char modes[N_TEXT_MODES] =
                 { 3, 42, 38, 34, 34, 34, 34, 3, 3};
   static const unsigned char xsizes[N_TEXT_MODES] =
                 { 80, 100, 80, 132, 132, 132, 132, 80, 80 };
   static const unsigned char ysizes[N_TEXT_MODES] =
                 { 25, 37, 60, 44, 50, 25, 28, 28, 50 };
            /* for a 16-high font, use 0x14 */
            /* for a 14-high font, use 0x11 */
            /* for an 8-high font, use 0x12 */
   static const unsigned char font_sizes[N_TEXT_MODES] =
                         { 0x14, 0, 0, 0, 0x12, 0x14, 0x11, 0x11, 0x12 };
   union _REGS regset;

   regset.w.ax = modes[mode];
   int386( 0x10, &regset, &regset);
   if( font_sizes[mode])
      {
      regset.h.ah = 0x11;
      regset.h.al = font_sizes[mode];
      regset.h.bl = 0x00;
      int386( 0x10, &regset, &regset);
      }
   if( xsize)
      *xsize = xsizes[mode];
   if( ysize)
      *ysize = ysizes[mode];
   return( 0);
}

static void show_message( const char *message)
{
   char FAR *loc = dispmem + xscr * 2 * (yscr / 2) + xscr;

   loc -= 2 * (strlen( message) / 2);        /* center the message */
   while( *message)
      {
      *loc++ = *message++;
      *loc++ = 47;               /* green backgnd */
      }
}

static void update_legend( const char *path_str, const int n_files,
                                    const int64_t n_bytes)
{
   char FAR *loc = dispmem + xscr * 2 * (yscr - 2);
   char buff[132];
   char tbuff[20];
   int i;

   for( i = 0; i < xscr; i++)
      {
      *loc++ = path_str[i];
      *loc++ = 47;               /* green backgnd */
      }

   sprintf( buff, strings[0], n_files);   /*  "%3d Files" */
               /* Show the total number of bytes in those files... */
   sprintf( tbuff, "%12lld", n_bytes);
               /* ...then spread out digits in groups of 3,  such as */
               /* "314 159 265 358"            */
   memmove( tbuff + 12, tbuff + 9, 4);
   memmove( tbuff + 8, tbuff + 6, 3);
   memmove( tbuff + 4, tbuff + 3, 3);
   tbuff[3] = tbuff[7] = tbuff[11] = ' ';
   strcat( buff, tbuff);
   strcat( buff, strings[1] + 4);         /* " bytes" */
   for( i = 0; buff[i]; i++)
      {
      *loc++ = buff[i];
      *loc++ = 31;
      }

   loc += 6;      /* skip three spaces */

   strcpy( buff, strings[2]);    /* "Hit '?' for a list of commands" */
   for( i = 0; buff[i]; i++)
      {
      *loc++ = buff[i];
      *loc++ = 95;
      }
}

#define MAX_HISTORY 100

char *prev_directories[MAX_HISTORY];
int cursor_locs[MAX_HISTORY];

void save_cursor_loc( const int cursor_loc)
{
   char path[_MAX_PATH];
   int i;

   getcwd( path, _MAX_PATH);
   for( i = 0; i < MAX_HISTORY && prev_directories[i]; i++)
      if( !strcmp( path, prev_directories[i]))
         {
         cursor_locs[i] = cursor_loc;
         return;
         }
   if( i == MAX_HISTORY)
      {
      for( i = 0; i < MAX_HISTORY; i++)
         {
         free( prev_directories[i]);
         prev_directories[i] = NULL;
         }
      i = 0;
      }
   prev_directories[i] = (char *)malloc( strlen( path) + 1);
   strcpy( prev_directories[i], path);
   cursor_locs[i] = cursor_loc;
}

int get_cursor_loc( void)
{
   char path[_MAX_PATH];
   int cursor_loc = 0, i;

   getcwd( path, _MAX_PATH);
   for( i = 0; i < MAX_HISTORY && prev_directories[i]; i++)
      if( !strcmp( path, prev_directories[i]))
         cursor_loc = cursor_locs[i];
   return( cursor_loc);
}

static DIRECT_ENT *get_directory( char *fileid, int *n_files)
{
   int found_one = 0;
   struct find_t c_file;
   int rval, i, n_alloced = MINIMAL_N_ALLOCED;
   DIRECT_ENT *dir = (DIRECT_ENT *)calloc( MINIMAL_N_ALLOCED, sizeof( DIRECT_ENT));

   dir = (DIRECT_ENT *)calloc( MINIMAL_N_ALLOCED, sizeof( DIRECT_ENT));
   *n_files = 0;
   rval = _dos_findfirst(fileid, _A_SUBDIR | _A_NORMAL |
                                   _A_HIDDEN | _A_SYSTEM | _A_ARCH, &c_file);
   while( !rval)
      {
      int use_it = 1;

      if( *n_files > n_alloced - 3)
         {
         n_alloced += n_alloced / 2;
         dir = (DIRECT_ENT *)realloc( dir, n_alloced * sizeof( DIRECT_ENT));
         if( !dir)
            {
            printf( "Failed to alloc %d directory entries\n", n_alloced);
            exit( 0);
            }
         }
      for( i = 0; c_file.name[i] != '.' && c_file.name[i]; i++)
         dir[*n_files].name[i] = c_file.name[i];
      dir[*n_files].name[i] = '\0';
      if( c_file.name[i] == '.')       /* extension */
         strcpy( dir[*n_files].ext, c_file.name + i + 1);
      else
         dir[*n_files].ext[0] = '\0';
      if( c_file.name[0] == '.' && (c_file.attrib & _A_SUBDIR))
         {            /* nameless directory */
         if( found_one)
            strcpy( dir[*n_files].name, "..");
         else
            strcpy( dir[*n_files].name, ".");
         found_one = 1;
         }
      if( sole_extension[0] && stricmp( sole_extension, dir[*n_files].ext))
         use_it = ( c_file.attrib & _A_SUBDIR);
      dir[*n_files].date = c_file.wr_date;
      dir[*n_files].time = c_file.wr_time;
      dir[*n_files].attrib = c_file.attrib;
      dir[*n_files].size = c_file.size;
      rval = _dos_findnext(&c_file);
      if( use_it)
         (*n_files)++;
      }
   getcwd( path_str, xscr);
   return( dir);
}

static char **load_strings( char *string_file_name)
{
   FILE *ifile = fopen( string_file_name, "rb");
   int n_lines, i, loc = 0, filesize;
   char buff[100];
   char *text;
   char **rval;

   if( !ifile)
      {
      *string_file_name = 'e';
      ifile = fopen( string_file_name, "rb");
      }
   if( !ifile)
      ifile = fopen( "c:\\utils\\e_di.dat", "rb");
   if( !ifile)
      return( NULL);

   for( n_lines = 0; fgets( buff, sizeof( buff), ifile); n_lines++)
      ;
   n_lines++;
   filesize = (int)ftell( ifile);
   rval = (char **)malloc( filesize + n_lines * sizeof( char **));
   text = (char *)( rval + n_lines);
   fseek( ifile, 0L, SEEK_SET);
   for( n_lines = 0; fgets( buff, sizeof( buff), ifile); n_lines++)
      {
      char *tptr = strstr( buff, "[]");

      if( *tptr)
         *tptr = '\0';
      for( i = 0; buff[i] != 10 && buff[i] != 13 && buff[i]; i++)
         ;
      buff[i] = '\0';
      rval[n_lines] = text + loc;
      strcpy( rval[n_lines], buff);
      loc += i + 1;
      }
   rval[n_lines] = NULL;
   fclose( ifile);
   return( rval);
}

void make_dirent_string( char *s, DIRECT_ENT *dir)
{
   unsigned day, month, year, min, hour;
   char *month_str = strings[3];
   int trail_loc = 38;
   char mstr[4];
            /* "JAN FEB MAR APR MAY JUN JUL AUG SEP OCT NOV DEC" */

   day = dir->date & DATE_MASK;
   month = (dir->date >> 5) & MONTH_MASK;
   year = (dir->date >> 9) & YEAR_MASK;
   year = (year+80) % 100;
   min = (dir->time >> 5) & MINUTE_MASK;
   hour = (dir->time >> 11) & HOUR_MASK;
   mstr[3] = '\0';
   memcpy( mstr, month_str + (month - 1) * 4, 3);
   if (dir->attrib & _A_SUBDIR)
       sprintf(s, "%-8s %-3s  <dir>  %2d%3s%02d %2d:%02d      ", dir->name, dir->ext,
          day, mstr, year,  hour, min);
   else
       sprintf(s, "%-8s %-3s%8ld %2d%3s%02d %2d:%02d      ", dir->name, dir->ext,
          dir->size, day, mstr, year,  hour, min);
   if( dir->attrib & _A_RDONLY)
      s[trail_loc--] = 'r';
// if( dir->attrib & _A_ARCH)
//    s[trail_loc--] = 'a';
   if( dir->attrib & _A_HIDDEN)
      s[trail_loc--] = 'h';
   if( dir->attrib & _A_SYSTEM)
      s[trail_loc--] = 's';
   if( s[38] == ' ')          /* No flags set;  show the seconds */
      sprintf( s + 34, ":%02d      ", (dir->time & SECONDS_MASK) * 2);
}

void show_cursor( int entry, int color, int width)
{
   char FAR *loc;
   int i, page_ht = yscr - 2;

   loc = dispmem + 2 * xscr * (entry % page_ht) + width * 2 * (entry / page_ht);
   for( i = 0; i < width; i++)
      loc[i + i + 1] = (char)color;
}

void disp_entry( int entry, DIRECT_ENT *dir, int width)
{
   char FAR *loc;
   char buff[50];
   int i, page_ht = yscr - 2;

   loc = dispmem + 2 * xscr * (entry % page_ht) + width * 2 * (entry / page_ht);
   if( dir)
      make_dirent_string( buff, dir);
   else
      memset( buff, ' ', width);
   if( width == 16)
      if( dir->attrib & _A_SUBDIR)
         buff[15] = ' ';
      else
         memset( buff + 13, ' ', 3);

   for( i = 0; i < width; i++)
      {
      *loc++ = buff[i];
      *loc++ = 15;      /* background color */
      }
}

int compare_name( const DIRECT_ENT *a, const DIRECT_ENT *b)
{
   int rval;

   rval = (b->attrib & _A_SUBDIR) - (a->attrib & _A_SUBDIR);
   if( !rval)
      rval = strcmp( a->name, b->name);
   if( !rval)
      rval = strcmp( a->ext, b->ext);
   return( rval);
}

int compare_ext( const DIRECT_ENT *a, const DIRECT_ENT *b)
{
   int rval;

   rval = (b->attrib & _A_SUBDIR) - (a->attrib & _A_SUBDIR);
   if( !rval)
      rval = strcmp( a->ext, b->ext);
   if( !rval)
      rval = strcmp( a->name, b->name);
   return( rval);
}

int compare_size( const DIRECT_ENT *a, const DIRECT_ENT *b)
{
   int rval = 0;

   rval = (b->attrib & _A_SUBDIR) - (a->attrib & _A_SUBDIR);
   if( !rval)
      if( a->size > b->size)
         rval = 1;
      else if( a->size < b->size)
         rval = -1;
   return( rval);
}

int compare_date( const DIRECT_ENT *a, const DIRECT_ENT *b)
{
   if( a->date != b->date)
      return( (a->date < b->date) ? -1 : 1);
   else
      return( a->time - b->time);
}

int compare_time( const DIRECT_ENT *a, const DIRECT_ENT *b)
{
   return( a->time - b->time);
}

void show_help_info( void)
{
   int i, j;

   for( i = 0; i < xscr * yscr * 2; i += 2)
      {
      dispmem[i] = ' ';
      dispmem[i + 1] = 15;
      }

   for( i = 0; strings[i + 7][0]; i++)
      for( j = 0; strings[i + 7][j]; j++)
         dispmem[i * 2 * xscr + j * 2] = strings[i + 7][j];
}

#ifdef __WATCOMC__
void _chdrive( const unsigned drive_letter)
{
   unsigned n_drives;

   _dos_setdrive( drive_letter, &n_drives);
}
#endif

int main( int argc, char *argv[])
{
   char FAR *backup;
   char *fileid = "*.*", *initial_path = NULL;
   char original_path[200];
   char *string_file_name = "e_di.dat";
   FILE *output_file = NULL;
   int i, c, cursor_loc = 0, page = 0, n_files, width = 40;
   int displayed_page = 0, update = 1, resort = 1;
   int return_to_home = 0, is_for_charon = 0;
   int sort_order = 1;
   DIRECT_ENT *dir;
   int( *sort_func)( const DIRECT_ENT *elem1, const DIRECT_ENT *elem2);
   time_t t0 = 0;

   getcwd( original_path, 200);
   for( i = 1; i < argc; i++)
      if( argv[i][0] == '-')
         switch( argv[i][1])
            {
            case 'r':
               return_to_home = 1;
               break;
            case 's':
               initial_path = argv[i] + 2;
               break;
            case 'f':
               output_file = fopen( argv[i] + 2, "wb");
               break;
            case '.':
               strcpy( sole_extension, argv[i] + 2);
               break;
            case 'i':
               output_file = fopen( argv[i] + 2, "wb");
               return_to_home = 1;
               is_for_charon = 1;
               break;
            case 'l':
               *string_file_name = argv[i][2];
               break;
            default:
               break;
            }
      else
         fileid = argv[i];
   strings = load_strings( string_file_name);
   if( initial_path)
      {
      _chdrive( toupper( *initial_path) - 'A' + 1);
      chdir( initial_path);
      }

   n_files = page = cursor_loc = 0;
   backup = (char FAR *)FMALLOC( 2 * xscr * yscr); /* store original screen */
   FMEMCPY( backup, dispmem, 2 * xscr * yscr);
   FMEMSET( dispmem, 0, 2 * xscr * yscr);
   sort_func = compare_ext;
   dir = get_directory( fileid, &n_files);
   c = 0;
   while( c != 27)
      {
      if( resort && sort_func)
         {
         int j, gap = 3280;

         while( gap /= 3)
            for( j = 0; j < gap; j++)
               for( i = j; i + gap < n_files; )
                  if( sort_func( dir + i, dir + i + gap) * sort_order > 0)
                     {
                     DIRECT_ENT temp = dir[i + gap];

                     dir[i + gap] = dir[i];
                     dir[i] = temp;
                     if( i >= gap)
                        i -= gap;
                     }
                  else
                     i += gap;
         }
      if( page != displayed_page || update || resort)
         {
         char FAR *loc;
         int64_t n_bytes = 0;

         for( i = 0; i < (yscr - 2) * (xscr / width); i++)
            if( page + i < n_files)
               disp_entry( i, dir + page + i, width);
            else
               disp_entry( i, NULL, width);
         displayed_page = page;
         loc = dispmem + 2 * xscr * (yscr - 2);
         for( i = 0; i < n_files; i++)
            n_bytes += (int64_t)dir[i].size;
         update_legend( path_str, n_files, n_bytes);
         update = 0;
         resort = 0;
         }
      show_cursor( cursor_loc - page, 112, width);
      if( width == 16)
         {
         char FAR *loc = dispmem + 2 * xscr * yscr - 48 * 2;
         char buff[180];

         make_dirent_string( buff, dir + cursor_loc);
         for( i = 0; buff[i] && i < 40; i++)
            {
            *loc++ = buff[i];
            *loc++ = 15;
            }
         }
      while( !kbhit( ))
         if( time( NULL) != t0)
            {
            char tstr[50];
            char FAR *bottom = dispmem + 2 * xscr * yscr;

            t0 = time( NULL);
            sprintf( tstr, "%02ld:%02ld:%02ld", (t0 / 3600L) % 24L,
                        (t0 / 60L) % 60L, t0 % 60L);
            for( i = 0; i < 8; i++)
               {
               bottom[i + i - 16] = tstr[i];
               bottom[i + i - 15] = 64 + 15;
               }
            }
      c = getch( );
      if( c == 0)
         c = getch( ) + 256;
      show_cursor( cursor_loc - page, 15, width);
      switch( c)
         {
         case ZKEY_NUM_8:
            cursor_loc--;
            break;
         case ZKEY_NUM_2:
            cursor_loc++;
            break;
         case ZKEY_NUM_4:
            cursor_loc -= (yscr - 2);
            break;
         case ZKEY_NUM_0:                /* top of current column */
         case ZKEY_NUM_5:                /* middle of current column */
         case ZKEY_NUM_9:                /* top quarter */
         case ZKEY_NUM_3:                /* last quarter */
            cursor_loc /= (yscr - 2);
            cursor_loc *= (yscr - 2);
            if( c == ZKEY_NUM_5)
               cursor_loc += (yscr - 2) / 2;
            else if( c == ZKEY_NUM_9)
               cursor_loc += (yscr - 2) / 4;
            else if( c == ZKEY_NUM_3)
               cursor_loc += 3 * (yscr - 2) / 4;
            break;
         case ZKEY_NUM_6:
            cursor_loc += (yscr - 2);
            break;
         case ZKEY_NUM_7:
            cursor_loc = 0;
            break;
         case ZKEY_NUM_1:
            cursor_loc = n_files - 1;
            break;
         case 9:        /* tab */
            c = getch( );
            save_cursor_loc( cursor_loc);
            _chdrive( toupper( c) - 'A' + 1);
            free( dir);
            dir = get_directory( fileid, &n_files);
            cursor_loc = get_cursor_loc( );
            resort = 1;
            break;
         case ZKEY_ALT_R:
            save_cursor_loc( cursor_loc);
            _chdrive( toupper( original_path[0]) - 'A' + 1);
            chdir( original_path);
            free( dir);
            dir = get_directory( fileid, &n_files);
            cursor_loc = get_cursor_loc( );
            resort = 1;
            break;
         case ZKEY_ALT_N:
            sort_func = compare_name;
            resort = 1;
            break;
         case ZKEY_ALT_Z:
            sort_func = NULL;
            resort = 1;
            break;
         case ZKEY_ALT_E:
            sort_func = compare_ext;
            resort = 1;
            break;
         case ZKEY_ALT_S:
            sort_func = compare_size;
            resort = 1;
            break;
         case ZKEY_ALT_T:
            sort_func = compare_time;
            resort = 1;
            break;
         case ZKEY_ALT_C:
            display_mode = (display_mode + 1) % N_TEXT_MODES;
            set_mode( display_mode, &xscr, &yscr);
            update = 1;
            break;
         case ZKEY_ALT_D:
            sort_func = compare_date;
            resort = 1;
            break;
         case ZKEY_ALT_W:
            {
            char name[90];
#ifdef __WATCOMC__
            int pmode = ((dir[cursor_loc].attrib & _A_RDONLY) ?
                                    S_IREAD | S_IWRITE : S_IREAD);

            sprintf( name, "%s.%s", dir[cursor_loc].name, dir[cursor_loc].ext);
            dir[cursor_loc].attrib ^= _A_RDONLY;
            chmod( name, pmode);
#else
            int pmode = ((dir[cursor_loc].attrib & _A_RDONLY) ?
                                    _S_IREAD | _S_IWRITE : _S_IREAD);

            sprintf( name, "%s.%s", dir[cursor_loc].name, dir[cursor_loc].ext);
            dir[cursor_loc].attrib ^= _A_RDONLY;
            _chmod( name, pmode);
#endif
            cursor_loc++;
            update = 1;
            }
            break;
         case ZKEY_ALT_V:
            width = 16 + 40 - width;
            update = 1;
            break;
         case ZKEY_ALT_X:
            show_message( "Confirm reset all files as readable,  unhidden,  non-system (Y/N):");
            if( getch( ) == 'y')
               {
               int i;
               char name[90];

               for( i = 0; i < n_files; i++)
                  if( dir[i].attrib & _A_RDONLY)
                     {
                     sprintf( name, "%s.%s", dir[i].name, dir[i].ext);

                     dir[i].attrib ^= _A_RDONLY;
#ifdef __WATCOMC__
                     chmod( name, S_IREAD | S_IWRITE);
#else
                     _chmod( name, _S_IREAD | _S_IWRITE);
#endif
                     update = 1;
                     }
               }
            break;
         case '-':
            {
            int i, j;
            DIRECT_ENT temp;

            for( i = 0, j = n_files - 1; i < j; i++, j--)
               {
               memcpy( &temp, dir + i, sizeof( DIRECT_ENT));
               memcpy( dir + i, dir + j, sizeof( DIRECT_ENT));
               memcpy( dir + j, &temp, sizeof( DIRECT_ENT));
               }
            sort_order *= -1;
            }
            update = 1;
            break;
         case '.': case '\\': case '/':
            save_cursor_loc( cursor_loc);
            if( c == '.')
               chdir( "..");
            else
               chdir( "\\");
            free( dir);
            dir = get_directory( fileid, &n_files);
            cursor_loc = get_cursor_loc( );
            resort = 1;
            break;
         case ' ':         /* edit or change directory */
            if( dir[cursor_loc].attrib & _A_SUBDIR)
               {
               char name[90];

               save_cursor_loc( cursor_loc);
               sprintf( name, "%s.%s", dir[cursor_loc].name, dir[cursor_loc].ext);
               chdir( name);
               free( dir);
               dir = get_directory( fileid, &n_files);
               cursor_loc = get_cursor_loc( );
               resort = 1;
               }
            else
               {
               char command[90];

               sprintf( command, "be %s.%s", dir[cursor_loc].name,
                                             dir[cursor_loc].ext);
               if( !output_file)
                  system( command);
               else
                  {
                  getcwd( command, 180);
                  i = strlen( command);
                  if( command[i - 1] == '\\')
                     command[i - 1] = '\0';
                  fprintf( output_file, "%s\\%s.%s\n", command,
                                    dir[cursor_loc].name, dir[cursor_loc].ext);
                  fprintf( output_file, "%s %s %s\n", command,
                                    dir[cursor_loc].name, dir[cursor_loc].ext);
                  fclose( output_file);
                  c = 27;
                  }
               }
            update = 1;
            break;
         case ZKEY_F1:          /* fixr */
         case ZKEY_F2:          /* fdump */
            if( !is_for_charon)
               {
               char command[90];

               sprintf( command, "%s %s.%s", (c == ZKEY_F1) ?  "fixr" : "fdump",
                           dir[cursor_loc].name, dir[cursor_loc].ext);
               system( command);
               update = 1;
               }
            break;
         case '?':
            show_help_info( );
            getch( );
            update = 1;
            break;
         case ZKEY_F8:
            FMEMCPY( dispmem, backup, 2 * xscr * yscr);
            getch( );
            update = 1;
            break;
         case ZKEY_CTRL_N:      /* delete _everything_ */
            if( !is_for_charon)
               {
               show_message( "Confirm delete ALL FILES IN THIS DIRECTORY");
               if( getch( ) == 'y')
                  {
                  for( i = n_files - 1; i; i--)
                     {
                     char fullname[90];
                     int rval;

                     sprintf( fullname, "%s.%s", dir[i].name, dir[i].ext);
                     if( dir[i].attrib & _A_SUBDIR)
                        rval = rmdir( fullname);
                     else
                        rval = unlink( fullname);
                     if( !rval)
                        {
                        n_files--;
                        memmove( dir + i, dir + i + 1,
                                    (n_files - i) * sizeof( DIRECT_ENT));
                        }
                     }
                  save_cursor_loc( cursor_loc);
                  chdir( "..");
                  free( dir);
                  dir = get_directory( fileid, &n_files);
                  cursor_loc = get_cursor_loc( );
                  resort = 1;
                  }
               update = 1;
               }
            break;
         case ZKEY_NUM_DEL:
            if( !is_for_charon)
               {
               char fullname[90];
               int rval;

               sprintf( fullname, "%s.%s", dir[cursor_loc].name,
                                           dir[cursor_loc].ext);
               if( dir[cursor_loc].attrib & _A_SUBDIR)
                  rval = rmdir( fullname);
               else
                  rval = unlink( fullname);
               if( !rval)
                  {
                  n_files--;
                  memmove( dir + cursor_loc, dir + cursor_loc + 1,
                           (n_files - cursor_loc) * sizeof( DIRECT_ENT));
                  }
               update = 1;
               }
            break;
         default:
            {
            int i, got_it = 0;
            int look_for_name = (sort_func == compare_name);

            if( c >= 'A' && c <= 'Z')
               look_for_name ^= 1;

            c = toupper( c);
            for( i = 1; i < n_files && !got_it; i++)
               {
               int new_cursor_loc = (cursor_loc + i) % n_files;

               if( look_for_name)
                  got_it = (dir[new_cursor_loc].name[0] == c);
               else
                  got_it = (dir[new_cursor_loc].ext[0] == c);
               if( got_it)
                  cursor_loc = new_cursor_loc;
               }
            }
            break;
         }
      if (cursor_loc > n_files - 1) cursor_loc = n_files - 1;
      if (cursor_loc < 0) cursor_loc = 0;
      while( page > cursor_loc)
         page -= (yscr - 2);
      if( page < 0)
         page = 0;
      while( page + (yscr - 2) * (xscr / width) <= cursor_loc)
         page += (yscr - 2);
      }
// if( display_mode)
//    set_mode( TEXT_MODE_80x25, NULL, NULL);
   FMEMCPY( dispmem, backup, 2 * xscr * yscr);
   if( return_to_home)
      {
      _chdrive( toupper( original_path[0]) - 'A' + 1);
      chdir( original_path);
      }
   printf( "\n%s    ", strings[4]);
   printf( "%s", strings[5]);
   return( 0);
}
