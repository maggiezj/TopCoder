int add_a_star( const IMAGE *image, int x, int y, FOUND_STAR *f,
                      const double cutoff_point, const int cell_size);
int log_printf( FILE *log_file, const char *format, ...);
void sprintf_ra_dec( char *str, double ival, const int format,
                                                      const int is_dec);
void put_ra_dec_in_str( char *str, const double ra, const double dec,
                                                      const int format);
int parse_ra_dec( const char *istr, double *ra, double *dec);
int write_out_startup( const char *filename,
                     const CHARON_CONFIG *cconfig, const IMAGE *image);
void init_config( CHARON_CONFIG *cconfig);
int parse_charon_command( CHARON_CONFIG *cconfig, IMAGE *image,
                                         const char *argv);
int load_configuration( CHARON_CONFIG *cconfig, IMAGE *image,
                                         const char *config_file_name);
int correct_cconfig_with_startup( CHARON_CONFIG *cconfig);
int get_observer_data( const char *obs_code, char *buff);
int compute_latlon_from_mpc_data( const char *buff, double *latlon);
int check_for_guide_cd( const char *cd_drive_path);
double extract_jd_from_string( const char *buff);
char *write_report_data( const CHARON_CONFIG *cconfig,
         const double jd, const double *target_ra_dec,
         const double target_mag, int report_type,
         const char *ref_net_name);
void flip_image( IMAGE *image);
void offset_image( IMAGE *image, const int signed_ints);
int grab_strings( void);
int get_palette_name( const int palette_no, char *name);
int check_mpc_report_header( const char *header_filename);
