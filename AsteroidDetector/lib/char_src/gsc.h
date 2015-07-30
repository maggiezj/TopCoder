#define LPT struct lpt

#ifndef FAR
#define FAR
#endif

LPT
   {
   long x, y;
   };

#define STAR struct star

STAR
   {
   unsigned short x, y;
   unsigned char mag, obj_type;
   };

#define SMALL_AREA struct small_area

SMALL_AREA
   {
   long x1, x2, y1, y2;
   int zone, small_zone, n_stars, n_pieces;
   STAR FAR *stars, FAR *stars2;
   STAR FAR **pieces;
   long *offsets;
   };

#define THREESIXTY_DEG 36000000L
#define ONEEIGHTY_DEG  18000000L
#define NINETY_DEG      9000000L

int find_regions( const LPT *min, const LPT *max, SMALL_AREA FAR *ret_arr,
                                    const int max_ret, const int level);
int find_gsc_regions( const LPT *min, const LPT *max, short *ret_arr,
                                    const int max_ret, const int level);
int find_bsc_regions( LPT *min, LPT *max, int n_divisions, int FAR *ret_arr,
                              SMALL_AREA *sm_arr);
int find_gsc_rect( long *rect, const short region);
