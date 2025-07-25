#ifndef GENERAL_DEFINED

   #define GENERAL_DEFINED

   #define MAX_LOG_ROW   512
   #define MAX_DIR_NAME  128
   #define MAX_FILE_NAME 128
   #define MAX_PATH_NAME 512
   #define MAX_SYS_COM   512

   #define MAX_IMGNO    999

   #define INVALID_TYPE   0
   #define STANDARD_TYPE  1
   #define KIWISPEC_TYPE  2

   typedef unsigned char byte;

   typedef char log_row_type[MAX_LOG_ROW+1];
   typedef char dir_name_type[MAX_DIR_NAME+1];
   typedef char file_name_type[MAX_FILE_NAME+1];
   typedef char path_name_type[MAX_PATH_NAME+1];
   typedef char sys_com_type[MAX_SYS_COM+1];

   typedef struct
   {
      path_name_type path;
      dir_name_type dir;
      file_name_type name;
      FILE *dat;
   }
   file_type;

#endif

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
void inform_title(char *msg);
void inform_table(char *msg);
void get_error(char *msg, ...);
int vget_system(int error_status,char *com,va_list ap);
void get_system(char *com, ...);
void get_system_err(int error_status,char *com, ...);
int go_system(char *com, ...);
char *hrsp_dir(void);
char *keep_tmpfiles(void);
char *today(void);
void split_path
  (path_name_type pathname,dir_name_type dirname,file_name_type filename);
void parse_file_name(char *full,char *fnam,char *fext);
void vset_file_name(file_type *file,path_name_type path,va_list ap);
void set_file_name(file_type *file,path_name_type path, ...);
void vremove_file(path_name_type path,va_list ap);
void remove_file(path_name_type path, ...);
int vfile_exists(path_name_type path,va_list ap);
int file_exists(path_name_type path, ...);
void get_fopen(file_type *file,char *mode);
void check_file_open(file_type *file);
void  get_fseek(file_type *file,int offset,int whence);
void  get_fread(byte *buf,int size,int count,file_type *file);
void  get_fwrite(byte *buf,int size,int count,file_type *file);
void get_fclose(file_type *file);
void vopen_file(file_type *file,char *mode,path_name_type path,va_list ap);
void open_file(file_type *file,char *mode,path_name_type path, ...);
int number_of_lines(path_name_type path, ...);
int nint(double x);
int upper_case_ok(char s);
int lower_case_ok(char s);
void get_upper_case(char *s);
double roundx(double x,int k);
int imin(int a,int b);
int imax(int a,int b);
int end_of_line_ok(char s);
int white_spc_ok(char s);
int empty_spc_ok(char s);
int quot_chr_ok(char s);
int alpha_ok(char s);
int digit_ok(char s);
int alphanum_ok(char c);
int ascii_ok(char s);
int extended_ascii_ok(unsigned char s);
int text_char_ok(char s);
int ascii_string_ok(char *s,int n);
int extended_ascii_string_ok(unsigned char *s,int n);
int text_string_ok(char *s,int n);
int blank_string(char *s);
int same_text(char *s,char *w);
int same_case_text(char *s,char *w);
int power_of_two(int n);
char sign_char(double x);
double signed_double(char sgn,double x);
double use_sign(double s,double x);
void get_lower_case(char *s);
void trim_right(char *s);
int file_name_match_ok(file_name_type file_name,file_name_type template);
char *next_tmp_file_name(file_name_type name);
void remove_tmp_file(path_name_type path, ...);
int count_items(char *s);
void collect_int(char *s,int n,int *a);
void list_integers(int *a,int n,char *s,int max);
void copy_integer_array(int *a,int n,int *b);
void copy_double_array(double *a,int n,double *b);
void copy_string_array(char *a,int n,int m,char *b);
int parse_uint(char *intstr,char **endptr,int def);
int str_to_int_def(char *a,int def);
double str_to_dbl_def(char *a,double def);
char str_to_log(char *a);
char *unquote(char *s,char *u,int nmax);
int collect_single_integer(char *s,char **e,int n,char *id);
int collect_integer_list(char *s,int *sel,int n,char *id);
int numeric_string_ok(char *str,int n);
int date_string_ok(char *datestr,int *year,int *month,int *day);
int time_string_ok(char *timestr,int dec,int *hour,int *min,double *sec);
