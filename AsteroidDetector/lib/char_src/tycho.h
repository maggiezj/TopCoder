#define TYCHO_STAR               struct tycho_star

#pragma pack(2)
TYCHO_STAR
   {
   short gsc_zone, gsc_num, mag;
   unsigned long ra, spd;
   long ppm_num, hipp_num, hd_num, sao_num;
   long pm_vals[2];
   unsigned char flamsteed;
   char gsc_id[4], constell, bayer, bayer_suffix, spectrum[3];
   char var_desig[15];
   short yale;
            /* The following items are only available in Guide 7 CDs: */
   short parallax, parallax_err;
   short pm_err[2], posn_err[2];
   short bt_mag, vt_mag;
   short bt_mag_err, vt_mag_err;
            /* The following items are only available in Guide 8 CDs: */
   short epoch_ra, epoch_dec;
   };
#pragma pack()

   /* All Tycho/Hipparcos data is for the epoch J1991.25,  8.75 Julian years */
   /* before J2000 = 2451545.  A Julian year is exactly 365.25 days.  So:    */

#define TYCHO_EPOCH (2451545.0 - 365.25 * 8.75)

int parse_tycho_star( const char FAR *data_ptr, TYCHO_STAR *star,
                      const int version_number);
