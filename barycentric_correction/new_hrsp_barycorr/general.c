#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
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
 
   exit(1);
}

void *get_space(int n)
{
   void *a;

   if (memindx >= MAX_MEMPTR) get_abort("get_space: Too many memory blocks!");
   a = malloc(n);
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
   if (rem > 0) memcpy(&memptr[k],&memptr[k+1],rem*sizeof(void *));
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
   vsnprintf(temp,bufsize+1,msg,ap);
   memcpy(buf,temp,bufsize);
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
   logrow_type row;

   dat = fopen("hercules.log","a");
   if (dat == NULL)
   {
      fprintf(stderr,"\nvlogfile: Cannot open 'hercules.log'!\n");
      exit(1);
   }

   if (vbuilds(row,sizeof(logrow_type),msg,ap) != 0)
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
   vreport(msg,ap);
   vlogfile(msg,ap);
}

void inform(char *msg, ...)
{
   va_list ap;

   va_start(ap,msg);
   vinform(msg,ap);
   va_end(ap);
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

void vget_system(char *com,va_list ap)
{
   syscom_type arg;
 
   if (vbuilds(arg,sizeof(syscom_type),com,ap) != 0)
      get_error("vget_system: System command too long:\n   %s\n",arg);

   logfile(arg);

   if (system(arg) != 0)
      get_error("vget_system: Could not execute:\n   %s",arg);
}

void get_system(char *com, ...)
{
   va_list ap;

   va_start(ap,com);
   vget_system(com,ap);
   va_end(ap);
}

char *herc_dir(void)
{
   char *ptr;

   ptr = getenv("HERCULES");
   if (ptr == NULL)
      get_error("herc_dir: Environment variable HERCULES not found!");

   return(ptr);
}

char *keep_tmpfiles(void)
{
   char *ptr;

   ptr = getenv("KEEP_TMPFILES");
   if (ptr == NULL)
      get_error("herc_dir: Environment variable KEEP_TMPFILES not found!");

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

void vset_file_name(file_type *file,pathname_type path,va_list ap)
{
   char *slash;

   if (vbuilds(file->path,sizeof(pathname_type),path,ap) != 0)
      get_error("vset_file_name: File path name too long!\n   %s ...\n",
                 file->path);

   slash = strrchr(file->path,'/');

   if (slash == NULL)
   {
      strcpy(file->dir,".");
      if (strlen(file->path) > MAX_FILENAME)
         get_error("vset_file_name: File name too long!\n   %s\n",
                 file->path);
      strcpy(file->name,file->path);
   }
   else
   {
      if (strlen(slash) == 1)
         get_error("vset_file_name: Missing file name!\n   %s\n",file->path);
      *slash = 0;
      if (strlen(file->path) > MAX_DIRNAME)
         get_error("vset_file_name: File directory name too long!\n   %s\n",
                 file->path);
      strcpy(file->dir,file->path);
      if (strlen(slash+1) > MAX_FILENAME)
         get_error("vset_file_name: File name too long!\n   %s\n",
                 file->path);
      strcpy(file->name,slash+1);
      *slash = '/';
   }

   if (builds(file->path,sizeof(pathname_type),"%s/%s",
                                               file->dir,file->name) != 0)
      get_error("vset_file_name: File path name too long!\n   %s/%s\n",
                 file->dir,file->name);   

   file->dat = NULL;
}

void set_file_name(file_type *file,pathname_type path, ...)
{
   va_list ap;

   va_start(ap,path);
   vset_file_name(file,path,ap);
   va_end(ap);
}

void vremove_file(pathname_type path,va_list ap)
{
   file_type f;

   vset_file_name(&f,path,ap);
   get_system("/bin/rm -f %s",f.path);
}

void remove_file(pathname_type path, ...)
{
   va_list ap;

   va_start(ap,path);
   vremove_file(path,ap);
   va_end(ap);
}

int vfile_exists(pathname_type path,va_list ap)
{
   file_type f;
 
   vset_file_name(&f,path,ap);

   f.dat = fopen(f.path,"r");
   if (f.dat == NULL) return(0);
   fclose(f.dat);
   return(1);
}

int file_exists(pathname_type path, ...)
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

void  get_fread(char *buf,int size,int count,file_type *file)
{
   check_file_open(file);
   if (fread(buf,size,count,file->dat) != count)
      get_error("get_fread: Cannot read %d bytes from '%s'!",
                 count*size,file->path);
}

void  get_fwrite(char *buf,int size,int count,file_type *file)
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

void vopen_file(file_type *file,char *mode,pathname_type path,va_list ap)
{
   vset_file_name(file,path,ap);
   get_fopen(file,mode);
}

void open_file(file_type *file,char *mode,pathname_type path, ...)
{
   va_list ap;

   va_start(ap,path);
   vopen_file(file,mode,path,ap);
   va_end(ap);
}

int number_of_lines(pathname_type path, ...)
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
 
   if (k < 1) get_error("round: Too few significant digits!");
 
   if (x == 0.0) return(x);
   a = fabs(x);
   n = (int)floor(log10(x));
 
   b = floor(a*pow(10.0,k-1-n)+0.5)*pow(10.0,n+1-k);
 
   if (x < 0.0)
      return(-b);
   else
      return(b);
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

   for (i=0;i<n;i++) if (!ascii_ok(s[i])) return(1);
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

char sign_character(double x)
{
   if (x == 0.0) return(' ');
   if (x > 0.0) return('+');
   return('-');
}

void get_lower_case(char *s)
{
   while (*s != 0)
   {
      if (alpha_ok(*s)) *s |= '\x20';
      s++;
   }
}

void remove_trailing_blanks(char *s)
{
   char *ptr;

   ptr = s + strlen(s);
   while (*ptr <= ' ' && ptr >= s) *(ptr--) = 0;
}

int filename_match_ok(filename_type filename,filename_type template)
{
   filename_type str;
   int k;

   if (*filename == 0)
      get_error("filename_match_ok: No characters found in the file name!");

   if (*template == 0)
      get_error("filename_match_ok: No characters found in the template!");

   while (*template != 0)
   {
      switch(*template)
      {
         case '?': if (*filename == 0) return(0);
                   template++;
                   filename++;
                   if (*template == '*')
                      get_error("filename_match_ok: '*' not expected "
                                "after '?' in the filename template!");
                   break;
         case '*': template++;
                   if (*template == '*' || *template == '?')
                      get_error("filename_match_ok: '%c' not expected "
                                "after '*' in the filename template!",
                                 *template);
                   if (*template == 0) return(1);
                   for (k=0;template[k]!=0 && template[k]!='*' 
                                           && template[k]!='?';k++);
                   memcpy(str,template,k);
                   str[k] = 0;
                   filename = strstr(filename,str);
                   if (filename == NULL) return(0);
                   filename += k;
                   template += k;
                   break;
         default:  while (*template!=0 && *template!='*' && *template!='?')
                      if (*(template++) != *(filename++)) return(0);
      }
   }

   if (*filename == 0)
      return(1);
   else
      return(0);
}

char *next_tmp_filename(filename_type name)
{
   file_type f;
   char row[40];
   int kmax,k,next;

   get_system("touch HercTmp000.fit");
   get_system("ls -1 HercTmp*.fit > ls.out");
   open_file(&f,"r","ls.out");
   kmax = -1;
   while (fgets(row,40,f.dat) != NULL)
   {
     if (memcmp(row,"HercTmp",7) != 0) continue;
     k = atoi(row+7);
     if (k > kmax) kmax = k;
   }
   get_fclose(&f);
   remove_file("HercTmp000.fit");
   remove_file("ls.out");
   next = kmax + 1;
   if (next > 999) get_error("Too many temporary files!");
   sprintf(name,"HercTmp%03d",next);
   return(name);
}

void remove_tmp_file(pathname_type path, ...)
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
   memcpy(b,a,n*sizeof(int));
}

void copy_double_array(double *a,int n,double *b)
{
   memcpy(b,a,n*sizeof(double));
}

void copy_string_array(char *a,int n,int m,char *b)
{
   memcpy(b,a,n*m);
}
