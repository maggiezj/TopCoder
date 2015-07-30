#define CHARON_CONFIG struct charon_config

CHARON_CONFIG
   {
   char path_name[120], cd_drive_path[80];
   char report_file_name[80], target_name[40];
   double cutoff_point, time_offset, max_residual;
   double search_dist, scale_tolerance, max_tilt, tilt_angle, latlon[2];
   double altitude;
   int millisec_per_video_reset, a10_drive_letter;
   int n_stars_used, max_stars, video_mode;
   char obs_code[4];
   int ra_dec_format, photometric_band;
   int cell_size, fit_order, saturation_point;
   int catalog_used, file_time_point, pixel_origin_at_top_left;
   int targeting_on, color_table_number, display_flipped, ccd_id;
   };
