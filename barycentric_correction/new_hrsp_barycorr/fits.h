#ifndef FITS_DEFINED

   #define FITS_DEFINED

   #define CARD_BYTES       80
   #define RECORD_CARDS     36
   #define RECORD_BYTES   2880   /* 36 cards times 80 bytes */

   #define HEAD_RECS        20
   #define HEAD_CARDS      720  /* 20 records times 36 cards   */
   #define HEAD_BYTES    57600  /* 20 records times 2880 bytes */

   #define MAXNAXIS 3
   #define MAXFITSCOLS 200
   #define MAXCOLNAME   32
   #define MAXCOLDISP    8

   #define MAX_FITS_KEYLEN 8
   #define MAX_FITS_VALLEN 70

   #define MAX_REGRESSION_NDIM   2
   #define MAX_REGRESSION_ARG    8
   #define MAX_REGRESSION_COEFF  100

   typedef char fitskey_type[MAX_FITS_KEYLEN+1];
   typedef char fitsval_type[MAX_FITS_VALLEN+1];

   typedef struct
   {
      char dat[HEAD_BYTES];
      int nrec,ncard;
   }
   fits_head;

   typedef struct
   {
      file_type file;

      /* ---- Primary header (image) ------ */
      fits_head head;
      int bitpix;  /* Bits per pixel (negative for floating format) */
      int naxis;   /* Number of axes */
      int npix[MAXNAXIS+1]; /* Number of pixels per axis */
      double crpix[MAXNAXIS+1]; /* Reference pixel */
      double crval[MAXNAXIS+1]; /* Coordinate at reference pixel */
      double cdelt[MAXNAXIS+1]; /* Coordinate increment per pixel */
      int pixsize; /* Pixel size in bytes */
      int totpix;   /* Total number of pixels */

      /* ---- Extension header (binary table) ------ */
      fits_head extend;
      int nrows,ncols,rowlen;
      int colstart[MAXFITSCOLS+1],colsize[MAXFITSCOLS+1];
      char coltype[MAXFITSCOLS+1],coldisp[MAXFITSCOLS+1][MAXCOLDISP+1];
      char colname[MAXFITSCOLS+1][MAXCOLNAME+1];
      int coloffset[MAXFITSCOLS+1];

      int binstart;  /* Offset of the binary data block */
      int totbinrec;  /* Total number of binary data records */
      unsigned char *bin;   /* Binary FITS data */
      float *pix; /* 32-bit real image pixels */
   }
   fits_file;

   typedef struct
   {
      int n;
      double min,max,tot;
      int pixmin[MAXNAXIS+1],pixmax[MAXNAXIS+1];
   }
   imgstat;

   typedef struct
   {
      int ndim;
      int ycol,xcol[MAX_REGRESSION_NDIM+1],deg[MAX_REGRESSION_NDIM+1];
      char arg[MAX_REGRESSION_ARG+1];
      int scol,fcol,rcol,nsel;
      double rms,kappa;
      int ma;
      double a[MAX_REGRESSION_COEFF+1];
   }
   regre_block;

#endif

/* ------------------------ Function Declaration ------------------------ */
void clear_fits_header(fits_head *head);
void clear_primary_header(fits_file *fits);
void clear_extension_header(fits_file *fits);
int fits_key_char_ok(char c);
int extract_fits_keyname(char *s,fitskey_type key);
char *normalize_fits_keyname(fitskey_type key,fitskey_type normkey);
char *collect_keyname(fits_file *fits,fits_head *head,int row,
                      fitskey_type key);
int find_fits_card(fits_head *head,char *text);
int find_fits_comment(fits_head *head,char *com);
int find_fits_keyword(fits_head *head,fitskey_type key);
int get_fits_keyword(fits_file *fits,fits_head *head,fitskey_type key);
void read_keyword_value(fits_file *fits,fits_head *head,int row,
                        fitsval_type val);
void read_keyword_textual(fits_file *fits,fits_head *head,int row,
                          char *dest,int maxlen);
double read_keyword_double(fits_file *fits,fits_head *head,int row);
int read_keyword_float(fits_file *fits,fits_head *head,int row);
int read_keyword_integer(fits_file *fits,fits_head *head,int row);
char read_keyword_logical(fits_file *fits,fits_head *head,int row);
void get_keyword_textual(fits_file *fits,fits_head *head,char *key,
                         char *str,int maxlen);
char get_keyword_logical(fits_file *fits,fits_head *head,char *key);
int get_keyword_integer(fits_file *fits,fits_head *head,char *key);
double get_keyword_double(fits_file *fits,fits_head *head,char *key);
float get_keyword_float(fits_file *fits,fits_head *head,char *key);
void check_keyword_logical(fits_file *fits,fits_head *head,char *key,char a);
void check_keyword_integer(fits_file *fits,fits_head *head,char *key,int a);
void remove_fits_cards(fits_file *fits,fits_head *head,int first,int last);
void overwrite_existing_card(fits_head *head,int seq,char *text);
void insert_new_card(fits_file *fits,fits_head *head,char *text);
void write_fits_card(fits_file *fits,fits_head *head,char *text);
void insert_blank_card(fits_file *fits,fits_head *head);
void insert_comment_card(fits_file *fits,fits_head *head,char *com);
void add_comment(char *text,char *com);
void write_keyword_logical(fits_file *fits,fits_head *head,char *key,
                           char a,char *com);
void write_keyword_integer(fits_file *fits,fits_head *head,char *key,
                           int a,char *com);
void write_keyword_double(fits_file *fits,fits_head *head,char *key,
                          double a,char *com);
void write_keyword_float(fits_file *fits,fits_head *head,char *key,
                         float a,char *com);
void write_keyword_textual(fits_file *fits,fits_head *head,char *key,
                           char *a,char *com);
int find_regression_block(fits_head *head,char *key);
int get_regression_block(fits_file *fits,fits_head *head,char *key);
void remove_regression_block(fits_file *fits,fits_head *head,char *key);
void insert_integer_value(fits_file *fits,fits_head *head,char *var,int val);
void insert_double_value(fits_file *fits,fits_head *head,char *var,double val);
void insert_integer_list(fits_file *fits,fits_head *head,char *var,
                         int *x,int n);
int find_fits_array(fits_head *head,char *var);
void remove_fits_array(fits_file *fits,fits_head *head,char *var);
void insert_double_array(fits_file *fits,fits_head *head,char *var,
                         double *x,int n);
void write_double_array(fits_file *fits,fits_head *head,char *var,
                        double *x,int n);
void write_regression_block(fits_file *fits,fits_head *head,char *key,
                            regre_block *reg);
char *get_regression_var(fits_file *fits,fits_head *head,char *key,
                         int row,char *var);
int extract_integer_value(fits_file *fits,fits_head *head,char *key,
                          int row,char *var);
void extract_integer_list(fits_file *fits,fits_head *head,char *key,
                          int row,char *var,int *a,int n);
double extract_double_value(fits_file *fits,fits_head *head,char *key,
                          int row,char *var);
void extract_double_array(fits_file *fits,fits_head *head,char *key,
                          int row,char *var,double *a,int n);
void get_double_array(fits_file *fits,fits_head *head,char *var,
                      double *a,int n);
void read_regression_block(fits_file *fits,fits_head *head,char *key,
                          regre_block *reg);
void set_totbinrec(fits_file *fits,int size);
void set_img_totbinrec(fits_file *fits);
void set_tbl_totbinrec(fits_file *fits);
void allocate_whole_binary(fits_file *fits);
void free_binary_data(fits_file *fits);
float collect_pixel_value(fits_file *fits,int x,int y);
void collect_pixel_row(fits_file *fits,int row,double *val);
void deposit_pixel_value(fits_file *fits,int x,int y,float val);
void deposit_pixel_row(fits_file *fits,int row,double *val);
void invert_image_pixels(fits_file *fits);
void get_raw_fits_table_data(fits_file *fits);
void get_normal_table_data(fits_file *fits);
void vset_fits_name(fits_file *fits,char *name,va_list ap);
void set_fits_name(fits_file *fits,char *name, ...);
void open_fits_file(fits_file *fits);
void create_fits_file(fits_file *fits);
void seek_fits_file(fits_file *fits,int offset,int whence);
void read_fits_file(fits_file *fits,int size,int count,char *buf);
void write_fits_file(fits_file *fits,int size,int count,char *buf);
void close_fits_file(fits_file *fits);
void load_fits_header(fits_file *fits,fits_head *head,char *signature);
void load_primary_header(fits_file *fits);
void load_extension_header(fits_file *fits);
void save_primary_header(fits_file *fits);
void save_extension_header(fits_file *fits);
void load_whole_binary(fits_file *fits);
void load_image_pixels(fits_file *fits);
void load_table_data(fits_file *fits);
void save_whole_binary(fits_file *fits);
void save_image_pixels(fits_file *fits);
void save_table_data(fits_file *fits);
void load_image_header(fits_file *fits);
void default_display_format(fits_file *fits,int col);
void load_table_header(fits_file *fits);
void save_image_header(fits_file *fits);
void save_table_header(fits_file *fits);
void load_image(fits_file *fits);
void load_table(fits_file *fits);
void load_raw_image(fits_file *fits);
void load_raw_table(fits_file *fits);
void write_image_date(fits_file *fits);
void save_image(fits_file *fits);
void save_raw_image(fits_file *fits);
void save_table(fits_file *fits);
void save_raw_table(fits_file *fits);
int vfits_table_ok(char *name,va_list ap);
int fits_table_ok(char *name, ...);
void vread_image_header(fits_file *fits,char *name,va_list ap);
void read_image_header(fits_file *fits,char *name, ...);
void vread_table_header(fits_file *fits,char *name,va_list ap);
void read_table_header(fits_file *fits,char *name, ...);
fits_head *vread_fits_header(fits_file *fits,char *name,va_list ap);
fits_head *read_fits_header(fits_file *fits,char *name, ...);
void read_binary_data(fits_file *fits);
void read_image_pixels(fits_file *fits);
void vread_image(fits_file *fits,char *name,va_list ap);
void vread_table(fits_file *fits,char *name,va_list ap);
void read_image(fits_file *fits,char *name,...);
void read_table(fits_file *fits,char *name,...);
void vread_raw_image(fits_file *fits,char *name,va_list ap);
void read_raw_image(fits_file *fits,char *name,...);
void vread_raw_table(fits_file *fits,char *name,va_list ap);
void read_raw_table(fits_file *fits,char *name,...);
void vwrite_image(fits_file *fits,char *name,va_list ap);
void vwrite_table(fits_file *fits,char *name,va_list ap);
void write_image(fits_file *fits,char *name, ...);
void write_table(fits_file *fits,char *name, ...);
void vwrite_raw_image(fits_file *fits,char *name,va_list ap);
void write_raw_image(fits_file *fits,char *name, ...);
void vwrite_raw_table(fits_file *fits,char *name,va_list ap);
void write_raw_table(fits_file *fits,char *name, ...);
void check_image_1d(fits_file *fits);
void check_image_2d(fits_file *fits);
void check_same_image_size(fits_file *x,fits_file *y);
void check_same_image_step(fits_file *x,fits_file *y);
void check_fft_image(fits_file *x);
void check_real_image(fits_file *fits);
int pixel_inside(fits_file *fits,int x,int y);
double world_coordinate(fits_file *x,int axis,double pix);
double pixel_coordinate(fits_file *x,int axis,double val);
double get_axis_start(fits_file *x,int axis);
void add_constant(char *xname,double c,char *yname);
void write_image_descriptor_integer(char *img,char *key,int a,char *com);
void write_image_descriptor_double(char *img,char *key,double a,char *com);
void write_image_descriptor_float(char *img,char *key,float a,char *com);
void write_image_descriptor_logical(char *img,char *key,char a,char *com);
void write_image_descriptor_textual(char *img,char *key,char *a,char *com);
void update_image_minmax(char *name);
void normalize_image(char *xname,char *yname);
void locate_pixel(fits_file *fits,int k,int *pix);
void stat_image_fits(fits_file *fits,imgstat *s);
void stat_image_file(char *name,imgstat *s);
void multiply_image_fits(fits_file *xfits,fits_file *yfits,fits_file *zfits);
void divide_image_fits(fits_file *xfits,fits_file *yfits,fits_file *zfits);
void subtract_image_fits(fits_file *xfits,fits_file *yfits,fits_file *zfits);
void subtract_image_file(char *xname,char *yname,char *zname,int pmod);
void divide_image_file(char *xname,char *yname,char *zname,int pmod);
void create_image(char *name,fits_head *head,int bitpix,int naxis,...);
void read_column_format(char *form,char *type,int *size);
void vcreate_fits_table(char **coldef,int nrows,char *name,va_list ap);
int find_table_column(fits_file *t,char *c);
int get_table_column(fits_file *t,char *c);
void get_several_table_columns(fits_file *t,char *list,int *c,int n);
void read_table_column_double(fits_file *t,int col,double *x);
void read_table_column_integer(fits_file *t,int col,int *x);
void write_table_column_double(fits_file *t,int col,double *x);
void write_table_column_integer(fits_file *t,int col,int *x);
