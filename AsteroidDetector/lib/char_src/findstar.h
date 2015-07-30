#define N_FIND_COEFFS 66
#define N_FIND_PARAMS 137
      /* The above allows for a tenth-order fit,  meaning 66 coefficients
         in x and another 66 in y,  plus five more parameters.       */

#define FOUND_STAR struct found_star
#define REF_STAR struct ref_star

FOUND_STAR
   {
   double x, y;
   long bright;
// short zone, number;
   };

REF_STAR
   {
   double ra, dec;
   long ppm_num;
   short zone, number, mag;
   char photo_band;
   };

#ifdef PIXEL_32
#define PIXEL unsigned long
#else
#define PIXEL unsigned short
#endif

#define MAX_PIXEL ((PIXEL)-1)
#define IMAGE struct image

IMAGE
   {
   char *filename;
   PIXEL *img, low_end, high_end;
   int xsize, ysize, image_type, byte_order;
   int is_inverted, is_xy_list, image_offset;
   double pixel_xsize, pixel_ysize;
   double focal_len, time, exposure_length, tilt_angle;
   double xform[N_FIND_PARAMS];
   long data_offset;
   int ra_dec_from_header, is_drift_scan;
   };

int find_stars( PIXEL *image, int xsize, int ysize,
            PIXEL limit, FOUND_STAR *found, int n_max, int cell_size);
int get_initial_match( const FOUND_STAR *found, const int n_stars,
           double *params, const double *search_loc, const int n0,
           const REF_STAR *gsc_stars, const int n_gsc_stars,
           const int photo_band);
REF_STAR *grab_gsc_12_data( const double x, const double y,
                             const double width, const double height,
                             const char *gsc_12_filename, int *n_found,
                             const double jd);
REF_STAR *grab_gsc_data( const double x, const double y,
               const double width, const double height,
               const char *cd_drive_path, int *n_found, const int use_gsc_act);
REF_STAR *grab_a10_data( const double x, const double y,
               const double width, const double height,
               const int cd_drive_letter, int *n_found, const int photo_band);
REF_STAR *grab_tycho_data( const double x, const double y,
                   const double width, const double height,
                   const char *cd_drive_path, int *n_found, double jd,
                   const int photo_band);
REF_STAR *grab_tycho2_data( const double x, const double y,
                   const double width, const double height,
                   const char *tycho2_filename, int *n_found, double jd,
                   const int photo_band);
REF_STAR *grab_ucac1_data( double x, const double y,
                   const double width, const double height,
                   const char *ucac_path, int *n_found, const double jd);
int ra_dec_to_pixel( double *pixel, const double ra, const double dec,
                                       const double *params);
int pixel_to_ra_dec( double *ra_dec, const double *pixel, const double *params);

PIXEL find_histo( const IMAGE *img, const double fraction);
int load_image( const char *filename, IMAGE *img);           /* load_img.c */
int match_contrasts( IMAGE *img, const IMAGE *source_img);   /* load_img.c */

#define CATALOG_A1            1
#define CATALOG_A2            4
#define CATALOG_GSC_11        0
#define CATALOG_GSC_ACT       6
#define CATALOG_UCAC2         7
#define CATALOG_UCAC3         8
#define CATALOG_UCAC4         9
#define CATALOG_TYCHO2        3
