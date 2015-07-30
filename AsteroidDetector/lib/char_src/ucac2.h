      /* Forcing byte-aligned packing is probably not _essential_,  but... */
#pragma pack( 1)

#define UCAC2_STAR struct ucac2_star

UCAC2_STAR
   {
   long ra, dec;                  /* RA/dec at J2000.0,  ICRS,  in milliarcsec  */
   short mag;                     /* internal UCAC mag,  red bandpass,  .01 mag */
   signed char e_ra, e_dec;       /* sigmas in RA and dec at central epoch      */
   signed char nobs;              /* Number of UCAC observations of this star   */
   signed char e_posn;            /* Position error for UCAC observations only  */
   signed char ncat;              /* # of catalog positions used for pmRA, pmDC */
   signed char catalog_flag;      /* ID of major catalogs use in pmRA, pmDE     */
   short epoch_ra;                /* Central epoch for mean RA, minus 1975      */
   short epoch_dec;               /* Central epoch for mean DE, minus 1975      */
   long pm_ra;                    /* Proper motion in RA (no cos DE), mas/year  */
   long pm_dec;                   /* Proper motion in DE, mas/year              */
   signed char e_pm_ra;           /* s.e. of pmRA (*cos DEm)                    */
   signed char e_pm_dec;          /* s.e. of pmDE                               */
   signed char q_pm_ra;           /* Goodness of fit for pmRA (in .05 units)    */
   signed char q_pm_dec;          /* Goodness of fit for pmDE (in .05 units)    */
   long twomass_id;               /* 2MASS pts_key star identifier              */
   short j_mag;                   /* 2MASS J  magnitude (.001 mags)             */
   short h_mag;                   /* 2MASS H  magnitude (.001 mags)             */
   short ks_mag;                  /* 2MASS Ks magnitude (.001 mags)             */
   signed char twomass_ph_qual;   /* 2MASS modified ph_qual flag                */
   signed char twomass_cc_flag;   /* 2MASS modified cc_flg                      */
   };
#pragma pack( )

/* Note: sizeof( UCAC2_STAR) = 44 bytes */

void flip_ucac2_star( UCAC2_STAR *star);
int extract_ucac2_stars( FILE *ofile, const double ra, const double dec,
                  const double width, const double height, const char *path,
                  const int is_supplement);
int extract_ucac2_info( const long ucac2_number, UCAC2_STAR *star,
                     const char *path);
int write_ucac2_star( const long offset, char *obuff,
                                          const UCAC2_STAR *star);
