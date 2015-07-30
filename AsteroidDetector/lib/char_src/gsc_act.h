#ifdef INCLUDE_MAG_PARAMS
#define N_COEFFS 29
#else
#define N_COEFFS 21
#endif

#define PLATE_STRUCT struct plate_struct

#pragma pack(4)

PLATE_STRUCT
   {
   char plate_name[4];
   long n_stars, n_rejected_in, n_rejected_out;
   double rms_in, rms_out;
   double ra0, dec0, coeffs[N_COEFFS * 2];
   };
#pragma pack()

void ra_dec_to_xi_eta( double d_ra, double dec, double dec0, double *xi,
                           double *eta);
void xi_eta_to_ra_dec( double xi, double eta, double dec0, double *d_ra,
                           double *dec);
void compute_gsc_act_powers( double *arr, const double x, const double y,
                     const double mag);
int gsc_act_cached_load( FILE *data_file, char *plate_name,
                                    PLATE_STRUCT *plate_array, int n_cached);
int apply_transformation( double ra_in, double dec_in, double *ra_out,
                  double *dec_out, PLATE_STRUCT *plate_data);
