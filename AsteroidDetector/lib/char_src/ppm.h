#define PPM_STAR struct ppm_star
#define SMALL_PPM_STAR struct small_ppm_star
#define PPM_SIZE 135

#pragma pack(1)

PPM_STAR
   {
   long ra;                         /* in .001 seconds */
   long dec;                        /* in .01 arcseconds */
   long sao_num, hd_num, ppm_num;
   unsigned short dm_number;        /* includes dm_letter? */
   unsigned short agk3_number;
   unsigned short gsc_zone, gsc_number, gsc_dist;     /* NOT ACTUALLY SET */
   short pm_ra, pm_dec;             /* in .0001 sec & .001" per year */
   unsigned short epoch_ra, epoch_dec;  /* in years/100 since 1900 */
   unsigned char mag;                       /* in tenths of a mag */
   char spectrum[2];
   unsigned char dm_cat;
   char dm_zone;
   char number_obs, notes;
   unsigned char sig_ra, sig_dec;            /* in arcsec/100 */
   unsigned char sig_pm_ra, sig_pm_dec;      /* in arcsec/10000 */
   char agk3_zone;
   };

      /* Smaller PPM record for drawing... reduced to 18 bytes/rec */
SMALL_PPM_STAR
   {
   long sao_num;
   short loc[3];
   long ppm_num;
   char pm_ra, pm_dec;                /* in .01 sec & .1" per year */
   unsigned char spectrum;
   unsigned char mag;                 /* in tenths of a mag */
   };
#pragma pack()
