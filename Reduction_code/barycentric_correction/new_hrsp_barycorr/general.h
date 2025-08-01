#ifndef GENERAL_DEFINED

   #define GENERAL_DEFINED

   #define MAX_LOGROW   512
   #define MAX_DIRNAME  60
   #define MAX_FILENAME 20
   #define MAX_PATHNAME 81
   #define MAX_SYSCOM   256

   #define MAX_IMGNO    999

   typedef char logrow_type[MAX_LOGROW+1];
   typedef char dirname_type[MAX_DIRNAME+1];
   typedef char filename_type[MAX_FILENAME+1];
   typedef char pathname_type[MAX_PATHNAME+1];
   typedef char syscom_type[MAX_SYSCOM+1];

   typedef struct
   {
      pathname_type path;
      dirname_type dir;
      filename_type name;
      FILE *dat;
   }
   file_type;

#endif

/* ------------------------ Function Declaration ------------------------ */
void get_abort(char *msg, ...);
void *get_space(int n);
void get_free(void *a);
void check_memory(void);
int vbuilds(char *buf,int bufsize,char *msg,va_list ap);
int builds(char *s,int n,char *msg, ...);
void vlogfile(char *msg,va_list ap);
void logfile(char *msg, ...);
void vreport(char *msg,va_list ap);
void report(char *msg, ...);
void vinform(char *msg,va_list ap);
void inform(char *msg, ...);
void get_error(char *msg, ...);
void vget_system(char *com,va_list ap);
void get_system(char *com, ...);
char *herc_dir(void);
char *keep_tmpfiles(void);
char *today(void);
void vset_file_name(file_type *file,pathname_type path,va_list ap);
void set_file_name(file_type *file,pathname_type path, ...);
void vremove_file(pathname_type path,va_list ap);
void remove_file(pathname_type path, ...);
int vfile_exists(pathname_type path,va_list ap);
int file_exists(pathname_type path, ...);
void get_fopen(file_type *file,char *mode);
void check_file_open(file_type *file);
void  get_fseek(file_type *file,int offset,int whence);
void  get_fread(char *buf,int size,int count,file_type *file);
void  get_fwrite(char *buf,int size,int count,file_type *file);
void get_fclose(file_type *file);
void vopen_file(file_type *file,char *mode,pathname_type path,va_list ap);
void open_file(file_type *file,char *mode,pathname_type path, ...);
int number_of_lines(pathname_type path, ...);
int nint(double x);
int upper_case_ok(char s);
int lower_case_ok(char s);
void get_upper_case(char *s);
/* */
double roundx(double x,int k);
int alpha_ok(char s);
int digit_ok(char s);
int alphanum_ok(char c);
int ascii_ok(char s);
int text_char_ok(char s);
int ascii_string_ok(char *s,int n);
int text_string_ok(char *s,int n);
int blank_string(char *s);
int same_text(char *s,char *w);
int same_case_text(char *s,char *w);
int power_of_two(int n);
char sign_character(double x);
void get_lower_case(char *s);
void remove_trailing_blanks(char *s);
int filename_match_ok(filename_type filename,filename_type template);
char *next_tmp_filename(filename_type name);
void remove_tmp_file(pathname_type path, ...);
int count_items(char *s);
void collect_int(char *s,int n,int *a);
void list_integers(int *a,int n,char *s,int max);
void copy_integer_array(int *a,int n,int *b);
void copy_double_array(double *a,int n,double *b);
void copy_string_array(char *a,int n,int m,char *b);
