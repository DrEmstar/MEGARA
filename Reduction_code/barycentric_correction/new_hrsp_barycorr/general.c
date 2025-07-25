#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <math.h>

#include "general.h"

#define MAX_MEMPTR 50

static char date_today[] = "0000-00-00";

static void *memptr[MAX_MEMPTR];
static int memindx=0;

void get_abort(char *msg, ...)
{
   va_list ap;

   va_start(ap,msg);

   fprintf(stderr,"\n");
   vfprintf(stderr,msg,ap);
   fprintf(stderr,"\n\nProgram aborted.\n");

   va_end(ap);

   exit(1);
}

void *get_space(int n)
{
   void *a;

   if (memindx >= MAX_MEMPTR) get_abort("get_space: Too many memory blocks!");
   a = malloc((size_t)n);
   if (a == NULL) get_abort("get_space: Cannot allocate %d bytes!",n);
   memptr[memindx++] = a;
   return(a);
}

void get_free(void *a)
{
   int k,rem;

   if (memindx < 1)
      get_abort("get_free: No more blocks to free!");
   for (k=0;k<memindx && a!=memptr[k];k++);
   if (a != memptr[k])
      get_abort("get_free: Trying to free a non-existing block!");
   free(a);
   memindx--;
   rem = memindx-k;
   if (rem > 0) memmove(&memptr[k],&memptr[k+1],rem*sizeof(void *));
}

void check_memory(void)
{
   if (memindx != 0)
      get_abort("check_memory: Memory allocation table not empty!");
}

int vbuilds(char *buf,int bufsize,char *msg,va_list ap)
{
   char *temp;
   int flag=0;

   temp = (char *)get_space(bufsize+1);
   vsnprintf(temp,(size_t)(bufsize+1),msg,ap);
   memmove(buf,temp,(size_t)bufsize);
   buf[bufsize-1]=0;
   if (strlen(temp) > bufsize-1) flag = 1;
   get_free(temp);
   return(flag);
}

int builds(char *s,int n,char *msg, ...)
{
   va_list ap;
   int flag;

   va_start(ap,msg);
   flag = vbuilds(s,n,msg,ap);
   va_end(ap);

   return(flag);
}

void vlogfile(char *msg,va_list ap)
{
   FILE *dat;
   log_row_type row;
   int retval;
   va_list aq;

   dat = fopen("hrsp.log","a");
   if (dat == NULL)
   {
      fprintf(stderr,"\nvlogfile: Cannot open 'hrsp.log'!\n");
      exit(1);
   }

   va_copy(aq,ap);
   retval = vbuilds(row,sizeof(log_row_type),msg,aq);
   va_end(aq);

   if (retval != 0)
   {
      fprintf(stderr,"\nvlogfile: The following message is too long"
                     "to be written into the log file:\n");
      vfprintf(stderr,msg,ap);
      exit(1);
   }

   if (fprintf(dat,"%s\n",row) != strlen(row)+1)
   {
      fprintf(stderr,"\nvlogfile: The following message cannot be"
                     "written into the log file (fprintf fails):\n");
      vfprintf(stderr,msg,ap);
      exit(1);
   }

   fclose(dat);
}

void logfile(char *msg, ...)
{
   va_list ap;

   va_start(ap,msg);
   vlogfile(msg,ap);
   va_end(ap);
}

void vreport(char *msg,va_list ap)
{
   vfprintf(stderr,msg,ap);
   fprintf(stderr,"\n");
}

void report(char *msg, ...)
{
   va_list ap;

   va_start(ap,msg);
   vreport(msg,ap);
   va_end(ap);
}

void vinform(char *msg,va_list ap)
{
   va_list aq;

   va_copy(aq,ap);
   vreport(msg,aq);
   va_end(aq);

   vlogfile(msg,ap);
}

void inform(char *msg, ...)
{
   va_list ap;

   va_start(ap,msg);
   vinform(msg,ap);
   va_end(ap);
}

void inform_title(char *msg)
{
  char dash[256];

  memset(dash,(int)'-',255);
  dash[255] = 0;
  if (strlen(msg) < 255) dash[strlen(msg)] = 0;
  inform("");
  inform(msg);
  inform(dash);
}

void inform_table(char *msg)
{
  char dash[256];

  memset(dash,(int)'-',255);
  dash[255] = 0;
  if (strlen(msg) < 255) dash[strlen(msg)] = 0;
  inform(dash);
  inform(msg);
  inform(dash);
}

void get_error(char *msg, ...)
{
   va_list ap;

   fprintf(stderr,"\n");

   va_start(ap,msg);
   vinform(msg,ap);
   va_end(ap);

   exit(1);
}

int vget_system(int error_status,char *com,va_list ap)
{
   sys_com_type arg;
   int result;

   if (vbuilds(arg,sizeof(sys_com_type),com,ap) != 0)
      get_error("vget_system: System command too long:\n   %s\n",arg);

   logfile(arg);

   result = system(arg);

   if (result < 0 && error_status != 0)
      get_error("vget_system: system() function call failure!");

   if ((result>>8) >= error_status && error_status != 0)
      get_error("vget_system: Could not execute:\n   %s",arg);

   return(result);
}

void get_system(char *com, ...)
{
   va_list ap;

   va_start(ap,com);
   vget_system(1,com,ap);
   va_end(ap);
}

void get_system_err(int error_status,char *com, ...)
{
   va_list ap;

   va_start(ap,com);
   vget_system(error_status,com,ap);
   va_end(ap);
}

int go_system(char *com, ...)
{
   va_list ap;
   int result;

   va_start(ap,com);
   result = vget_system(0,com,ap);
   va_end(ap);
   return(result);
}

char *hrsp_dir(void)
{
   char *ptr;
   setenv("HRSP","./",1);
   ptr = getenv("HRSP");
   if (ptr == NULL)
      get_error("hrsp_dir: Environment variable HRSP not found!");

   return(ptr);
}
char *keep_tmpfiles(void)
{
   char *ptr;

   ptr = getenv("KEEP_TMPFILES");
   if (ptr == NULL) get_error
      ("keep_tmpfiles: Environment variable KEEP_TMPFILES not found!");

   return(ptr);
}

char *today(void)
{
   time_t time_since;
   struct tm time_today;

   time(&time_since);
   time_today = *localtime(&time_since);
   sprintf(date_today,"%04d-%02d-%02d",time_today.tm_year+1900,
                           time_today.tm_mon+1,time_today.tm_mday);

   return(date_today);
}

void split_path
  (path_name_type pathname,dir_name_type dirname,file_name_type filename)
{
  path_name_type fullpath;
  char *firstptr,*nextptr,*lastptr,*nameptr;

  strcpy(fullpath,pathname);

  if (strchr(fullpath,0x20) != NULL) get_error
    ("split_path: Blank character found in the path name '%s'!",fullpath);

  while ((firstptr = strstr(fullpath,"//")) != NULL)
  {
    nextptr = firstptr + 1;
    while (*firstptr != 0) *(firstptr++) = *(nextptr++);
  }

  lastptr = strrchr(fullpath,'/');

  if (lastptr == NULL)
  {
    nameptr = fullpath;
    strcpy(dirname,".");
  }
  else
  {
    nameptr = lastptr + 1;
    *lastptr = 0;
    if (strlen(fullpath) > MAX_DIR_NAME)
      get_error("split_path: Directory name too long: '%s'!",fullpath);
    strcpy(dirname,fullpath);
  }

  if (strlen(nameptr) > MAX_FILE_NAME)
    get_error("split_path: File name too long: '%s'!",nameptr);
  strcpy(filename,nameptr);
}

void parse_file_name(char *full,char *fnam,char *fext)
{
  char *fptr;

  strcpy(fnam,full);
  *fext = 0;
  fptr = strrchr(full,'.');
  if (fptr != NULL)
  {
    strcpy(fext,fptr);
    *fptr = 0;
  }
}

void vset_file_name(file_type *file,path_name_type path,va_list ap)
{
   char *slash;

   if (vbuilds(file->path,sizeof(path_name_type),path,ap) != 0)
      get_error("vset_file_name: File path name too long!\n   %s ...\n",
                 file->path);

   slash = strrchr(file->path,'/');

   if (slash == NULL)
   {
      strcpy(file->dir,".");
      if (strlen(file->path) > MAX_FILE_NAME)
         get_error("vset_file_name: File name too long!\n   %s\n",
                 file->path);
      strcpy(file->name,file->path);
   }
   else
   {
      if (strlen(slash) == 1)
         get_error("vset_file_name: Missing file name!\n   %s\n",file->path);
      *slash = 0;
      if (strlen(file->path) > MAX_DIR_NAME)
         get_error("vset_file_name: File directory name too long!\n   %s\n",
                 file->path);
      strcpy(file->dir,file->path);
      if (strlen(slash+1) > MAX_FILE_NAME)
         get_error("vset_file_name: File name too long!\n   %s\n",
                 file->path);
      strcpy(file->name,slash+1);
      *slash = '/';
   }

   if (builds(file->path,sizeof(path_name_type),"%s/%s",
                                               file->dir,file->name) != 0)
      get_error("vset_file_name: File path name too long!\n   %s/%s\n",
                 file->dir,file->name);

   file->dat = NULL;
}

void set_file_name(file_type *file,path_name_type path, ...)
{
   va_list ap;

   va_start(ap,path);
   vset_file_name(file,path,ap);
   va_end(ap);
}

void vremove_file(path_name_type path,va_list ap)
{
   file_type f;

   vset_file_name(&f,path,ap);
   get_system("/bin/rm -f %s",f.path);
}

void remove_file(path_name_type path, ...)
{
   va_list ap;

   va_start(ap,path);
   vremove_file(path,ap);
   va_end(ap);
}

int vfile_exists(path_name_type path,va_list ap)
{
   file_type f;

   vset_file_name(&f,path,ap);

   f.dat = fopen(f.path,"r");
   if (f.dat == NULL) return(0);
   fclose(f.dat);
   return(1);
}

int file_exists(path_name_type path, ...)
{
   va_list ap;
   int flag;

   va_start(ap,path);
   flag = vfile_exists(path,ap);
   va_end(ap);

   return(flag);
}

void get_fopen(file_type *file,char *mode)
{
   if (file->dat != NULL)
      get_error("get_file_open: File '%s' appears to be already open!",
                 file->path);
   file->dat = fopen(file->path,mode);
   if (file->dat == NULL)
      get_error("get_file_open: Cannot open '%s' for '%s'!",file->path,mode);
}

void check_file_open(file_type *file)
{
   if (file->dat == NULL)
      get_error("check_file_open: File '%s' is not open yet!",file->path);
}

void  get_fseek(file_type *file,int offset,int whence)
{
   check_file_open(file);
   if (fseek(file->dat,offset,whence))
      get_error("get_fseek: Cannot 'fseek' file '%s'!",file->path);
}

void  get_fread(byte *buf,int size,int count,file_type *file)
{
   check_file_open(file);
   if (fread(buf,size,count,file->dat) != count)
      get_error("get_fread: Cannot read %d bytes from '%s'!",
                 count*size,file->path);
}

void  get_fwrite(byte *buf,int size,int count,file_type *file)
{
   check_file_open(file);
   if (fwrite(buf,size,count,file->dat) != count)
      get_error("get_fwrite: Cannot write %d bytes into '%s'!",
                 count*size,file->path);
}

void get_fclose(file_type *file)
{
   check_file_open(file);
   if (fclose(file->dat))
      get_error("get_fclose: Cannot close '%s'!",file->path);
   file->dat = NULL;
}

void vopen_file(file_type *file,char *mode,path_name_type path,va_list ap)
{
   vset_file_name(file,path,ap);
   get_fopen(file,mode);
}

void open_file(file_type *file,char *mode,path_name_type path, ...)
{
   va_list ap;

   va_start(ap,path);
   vopen_file(file,mode,path,ap);
   va_end(ap);
}

int number_of_lines(path_name_type path, ...)
{
   va_list ap;
   file_type f;
   int n,r,c;

   va_start(ap,path);
   vopen_file(&f,"r",path,ap);
   va_end(ap);

   n = r = 0;
   while ((c=fgetc(f.dat)) != EOF)
   {
      if (c == 0x0A)
         n++,r=0;
      else
         r++;
   }
   get_fclose(&f);
   if (r != 0) n++;

   return(n);
}

int nint(double x)
{
   return((int)floor(x+0.5));
}

int upper_case_ok(char s)
{
   if (s>='A' && s<='Z')
      return(1);
   else
      return(0);
}

int lower_case_ok(char s)
{
   if (s>='a' && s<='z')
      return(1);
   else
      return(0);
}

void get_upper_case(char *s)
{
   char *p;

   for (p=s;*p!=0;p++) if (lower_case_ok(*p)) *p = *p & '\xDF';
}

double roundx(double x,int k)
{
   int n;
   double a,b;

   if (k < 1) get_error("roundx: Too few significant digits!");

   if (x == 0.0) return(x);
   a = fabs(x);
   n = (int)floor(log10(x));

   b = floor(a*pow(10.0,k-1-n)+0.5)*pow(10.0,n+1-k);

   if (x < 0.0)
      return(-b);
   else
      return(b);
}

int imin(int a,int b)
{
  if (b < a)
    return(b);
  else
    return(a);
}

int imax(int a,int b)
{
  if (b > a)
    return(b);
  else
    return(a);
}

int end_of_line_ok(char s)
{
   return(s == '\r' || s == '\n');
}

int white_spc_ok(char s)
{
   return(s == ' ' || s == '\t');
}

int empty_spc_ok(char s)
{
   return(white_spc_ok(s) || end_of_line_ok(s));
}

int quot_chr_ok(char s)
{
   return(s == '\x22' || s == '\x27');
}

int alpha_ok(char s)
{
   if (upper_case_ok(s) || lower_case_ok(s))
      return(1);
   else
      return(0);
}

int digit_ok(char s)
{
   if (s>='0' && s<='9')
      return(1);
   else
      return(0);
}

int alphanum_ok(char c)
{
   if (alpha_ok(c) || digit_ok(c))
      return(1);
   else
      return(0);
}

int ascii_ok(char s)
{
   if (s>=' ' && s<='~')
      return(1);
   else
      return(0);
}

int extended_ascii_ok(unsigned char s)
{
   if (s>=' ')
      return(1);
   else
      return(0);
}

int text_char_ok(char s)
{
   if (ascii_ok(s) || s=='\t')
      return(1);
   else
      return(0);
}

int ascii_string_ok(char *s,int n)
{
   int i;

   for (i=0;i<n;i++) if (!ascii_ok(s[i])) return(0);
   return(1);
}

int extended_ascii_string_ok(unsigned char *s,int n)
{
   int i;

   for (i=0;i<n;i++) if (!extended_ascii_ok(s[i])) return(0);
   return(1);
}

int text_string_ok(char *s,int n)
{
   int i;

   for (i=0;i<n;i++) if (!text_char_ok(s[i])) return(0);
   return(1);
}

int blank_string(char *s)
{
   while (*s == ' ') s++;
   if (*s == 0)
      return(1);
   else
      return(0);
}

int same_text(char *s,char *w)
{
   if (strncmp(s,w,strlen(w)) == 0)
      return(1);
   else
      return(0);
}

int same_case_text(char *s,char *w)
{
   if (strncasecmp(s,w,strlen(w)) == 0)
      return(1);
   else
      return(0);
}

int power_of_two(int n)
{
   if (n < 2) return(0);
   while ((n & 1) == 0) n = n>>1;
   if (n == 1)
      return(1);
   else
      return(0);
}

char sign_char(double x)
{
   if (x == 0.0) return(' ');
   if (x > 0.0) return('+');
   return('-');
}

double signed_double(char sgn,double x)
{
  if (sgn == '-')
    return(-fabs(x));
  else
    return(fabs(x));
}

double use_sign(double s,double x)
{
  if (s < 0.0)
    return(-fabs(x));
  else
    return(fabs(x));
}

void get_lower_case(char *s)
{
   while (*s != 0)
   {
      if (alpha_ok(*s)) *s |= '\x20';
      s++;
   }
}

void trim_right(char *s)
{
   char *ptr;

   ptr = s + strlen(s);
   while (*ptr <= ' ' && ptr >= s) *(ptr--) = 0;
}

int file_name_match_ok(file_name_type file_name,file_name_type template)
{
   file_name_type str;
   int k;

   if (*file_name == 0)
      get_error("file_name_match_ok: No characters found in the file name!");

   if (*template == 0)
      get_error("file_name_match_ok: No characters found in the template!");

   while (*template != 0)
   {
      switch(*template)
      {
         case '?': if (*file_name == 0) return(0);
                   template++;
                   file_name++;
                   if (*template == '*')
                      get_error("file_name_match_ok: '*' not expected "
                                "after '?' in the file name template!");
                   break;
         case '*': template++;
                   if (*template == '*' || *template == '?')
                      get_error("file_name_match_ok: '%c' not expected "
                                "after '*' in the file name template!",
                                 *template);
                   if (*template == 0) return(1);
                   for (k=0;template[k]!=0 && template[k]!='*'
                                           && template[k]!='?';k++);
                   memmove(str,template,k);
                   str[k] = 0;
                   file_name = strstr(file_name,str);
                   if (file_name == NULL) return(0);
                   file_name += k;
                   template += k;
                   break;
         default:  while (*template!=0 && *template!='*' && *template!='?')
                      if (*(template++) != *(file_name++)) return(0);
      }
   }

   if (*file_name == 0)
      return(1);
   else
      return(0);
}

char *next_tmp_file_name(file_name_type name)
{
   file_type f;
   char row[40];
   int kmax,k,next;

   get_system("touch HrspTmp000.fit");
   get_system("ls -1 HrspTmp*.fit > ls.out");
   open_file(&f,"r","ls.out");
   kmax = -1;
   while (fgets(row,40,f.dat) != NULL)
   {
     if (memcmp(row,"HrspTmp",7) != 0) continue;
     k = atoi(row+7);
     if (k > kmax) kmax = k;
   }
   get_fclose(&f);
   remove_file("HrspTmp000.fit");
   remove_file("ls.out");
   next = kmax + 1;
   if (next > 999) get_error("Too many temporary files!");
   sprintf(name,"HrspTmp%03d",next);
   return(name);
}

void remove_tmp_file(path_name_type path, ...)
{
   va_list ap;

   if (strcasecmp(keep_tmpfiles(),"YES") == 0) return;

   va_start(ap,path);
   vremove_file(path,ap);
   va_end(ap);
}

int count_items(char *s)
{
   int k;
   char *x;

   for (x=s,k=1;*x!=0;x++) if (*x == ',') k++;

   return(k);
}

void collect_int(char *s,int n,int *a)
{
   int k;
   char *p,*endptr;

   for (k=1,p=s;k<=n;k++,p++)
   {
      a[k] = (int)strtol(p,&endptr,10);
      if (endptr == p) get_error("An integer value expected!");
      p = endptr;
      if (k<n && *p!=',') get_error("A comma expected between the values!");
   }
}

void list_integers(int *a,int n,char *s,int max)
{
   int i;
   char buf[16];

   for (i=1,*s=0;i<=n;i++)
   {
      if (i == 1)
         sprintf(buf,"%d",a[i]);
      else
         sprintf(buf,",%d",a[i]);
      if (strlen(s)+strlen(buf) > max)
         get_error("list_integers: Output string too long!");
      strcat(s,buf);
   }
}

void copy_integer_array(int *a,int n,int *b)
{
   memmove(b,a,n*sizeof(int));
}

void copy_double_array(double *a,int n,double *b)
{
   memmove(b,a,n*sizeof(double));
}

void copy_string_array(char *a,int n,int m,char *b)
{
   memmove(b,a,n*m);
}

int parse_uint(char *intstr,char **endptr,int def)
{
  int32_t x;

  while (*intstr == ' ') intstr++;
  *endptr = NULL;
  errno = 0;
  x = strtol(intstr,endptr,10);
  if (*endptr == NULL) *endptr = intstr;
  if (errno != 0 || *endptr == intstr) return(def);
  return((int)x);
}

int str_to_int_def(char *a,int def)
{
   int n;
   char *eptr;

   errno = 0;
   n = strtol(a,&eptr,10);
   if (errno == 0 && eptr > a && *eptr == 0)
      return(n);
   else
      return(def);
}

double str_to_dbl_def(char *a,double def)
{
   double x;
   char *eptr;

   errno = 0;
   x = strtod(a,&eptr);
   if (errno == 0 && eptr > a && *eptr == 0)
      return(x);
   else
      return(def);
}

char str_to_log(char *a)
{
   if (*a == 'T' || *a == 't')
      return('T');
   else
      return('F');
}

char *unquote(char *s,char *u,int nmax)
{
   char *src,*dest;
   int nleft;
   char q;

   src = s;
   dest = u;
   nleft = nmax;
   q = *src;
   if (q == '\x22' || q == '\x27')
   {
      src++;
      while(*src != q && *src != 0 && nleft > 0) *(dest++) = *(src++),nleft--;
   }
   else
      strncpy(u,s,nmax+1);
   u[nmax] = 0;
   return(u);
}

int collect_single_integer(char *s,char **e,int n,char *id)
{
   int k;

   if (!digit_ok(*s)) get_error("%s expected!",id);
   k =  (int)strtol(s,e,10);
   if (k<1 || k>n) get_error("%s out of limits!",id);

   return(k);
}

int collect_integer_list(char *s,int *sel,int n,char *id)
{
   int first,last,step,count,seq,flag;

   memset(sel,0,(n+1)*sizeof(int));

   while (digit_ok(*s) || *s == '/')
   {
      if (*s == '/')
         flag = 0, s++;
      else
         flag = 1;

      first = collect_single_integer(s,&s,n,id);
      if (*s == '-')
         last = collect_single_integer(s+1,&s,n,id);
      else
         last = first;
      if (first > last) get_error("%s out of sequence!",id);

      if (last > first && *s == ':')
         step = collect_single_integer(s+1,&s,n,"Step");
      else
         step = 1;

      count = (last - first) / step;
      if (first + count * step != last)
         get_error("%s list: Bad range/step specification!",id);

      for (seq=first;seq<=last;seq+=step) sel[seq] = flag;

      if (flag == 0)
      {
         if (*s != '/') get_error("%s list: Character '/' expected!",id);
         s++;
      }

      if (*s == ',') s++;
   }
   if (*s != 0) get_error("%s list: Unexpected character %02Xh!",id,*s);

   for (seq=1,count=0;seq<=n;seq++) count += sel[seq];

   return(count);
}

int numeric_string_ok(char *str,int n)
{
  while (digit_ok(*str) && n > 0) str++,n--;
  return(n == 0);
}

int date_string_ok(char *datestr,int *year,int *month,int *day)
{
  if (strlen(datestr) < 10) return(0);
  if (!numeric_string_ok(datestr,4)) return(0);
  if (datestr[4] != '-') return(0);
  if (!numeric_string_ok(datestr+5,2)) return(0);
  if (datestr[7] != '-') return(0);
  if (!numeric_string_ok(datestr+8,2)) return(0);
  if (sscanf(datestr,"%d-%d-%d",year,month,day) != 3) return(0);
  return(1);
}

int time_string_ok(char *timestr,int dec,int *hour,int *min,double *sec)
{
  if (strlen(timestr) < 9+dec) return(0);
  if (!numeric_string_ok(timestr,2)) return(0);
  if (timestr[2] != ':') return(0);
  if (!numeric_string_ok(timestr+3,2)) return(0);
  if (timestr[5] != ':') return(0);
  if (!numeric_string_ok(timestr+6,2)) return(0);
  if (dec >= 0) if (timestr[8] != '.') return(0);
  if (dec > 0) if (!numeric_string_ok(timestr+9,dec)) return(0);
  if (sscanf(timestr,"%d:%d:%lf",hour,min,sec) != 3) return(0);
  return(1);
}
