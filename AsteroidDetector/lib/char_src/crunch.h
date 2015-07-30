#define GSC_STAR struct gsc_star
#define GSC_HEADER struct gsc_header

#pragma pack(1)
GSC_STAR
   {
   long x, y;
   unsigned short n, is_a_dup, posn_err, mag_err, mag;
   unsigned short plate_id_no, mag_band, obj_type;
   char plate[5];
   };

GSC_HEADER
   {
   long unused[5], n_ppm_stars;
   long ra_min, ra_max, dec_min, dec_max, ra0, dec0;
   long prev_x, prev_y, filesize;
   unsigned short n_stars, n_plates, prev_n, total_lines, n_entries;
   unsigned short area_no, zone_no, magerr_0, magband_0;
   char plates[60];
   };
#pragma pack()

#define X_RANGE (1L << 23)
#define Y_RANGE (1L << 19)

int gsc_star_parse( GSC_STAR *star, GSC_HEADER *hdr, char *buff);
int gsc_star_crunch( GSC_STAR *star, GSC_HEADER *hdr, char *buff);
int gsc_star_uncrunch( GSC_STAR *star, GSC_HEADER *hdr,
                                            const unsigned char *buff);
