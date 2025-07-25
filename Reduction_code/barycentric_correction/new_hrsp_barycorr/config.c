#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#include "general.h"
#include "angle.h"
#include "config.h"

static char *charptr;
static int white_brackets = 1;

static char ccd_flips[][MAX_VALSIZE+1] =
  {"NONE","X","Y","XY",""};

static char wave_list[][MAX_VALSIZE+1] =
  {"AIR","VACUUM",""};

static char back_methods[][MAX_VALSIZE+1] =
  {"CUBIC","SPLINE","POLYFIT",""};

static char flat_methods[][MAX_VALSIZE+1] =
  {"NONE","EXTRACTED","FULL",""};

static char rv_preps[][MAX_VALSIZE+1] =
  {"NONE","FLIP","MEAN",""};

static char rv_methods[][MAX_VALSIZE+1] =
  {"GAUSS","PARABOLA","SPLINE","PEAK","CCFMAX",""};

int cfg_white_space_ok(char c)
{
   if (c==' ' || c=='\t') return(1);
   if (white_brackets && (c=='[' || c==']')) return(1);
   return(0);
}

int cfg_quot_ok(char c)
{
   if (c=='\x27' || c=='\x22')
      return(1);
   else
      return(0);
}

int cfg_end_of_line_ok(char c)
{
   if (c=='\n' || c=='\r' || c=='%' || c==0)
      return(1);
   else
      return(0);
}

int keyname_char_ok(char c)
{
   if (alphanum_ok(c) || c=='_')
      return(1);
   else
      return(0);
}

int keyval_char_ok(char c)
{
   if (ascii_ok(c) && !cfg_end_of_line_ok(c) && !cfg_white_space_ok(c))
      return(1);
   else
      return(0);
}

char get_char(void)
{
   while (cfg_white_space_ok(*charptr)) charptr++;
   return(*charptr);
}

char next_char(void)
{
   charptr++;
   return(get_char());
}

void scan_config_line(char *s,keytype key,int *n,valtype val)
{
   cfglinetype cfg;
   char q;
   int k;

   if (strlen(s) > MAXCFGLINE)
      get_error("scan_config_line: Line too long!");
   strcpy(cfg,s);
   trim_right(cfg);
   charptr = cfg;
   if (cfg_end_of_line_ok(get_char()))
      get_error("scan_config_line: Line is empty!");

   if (!alpha_ok(get_char()))
      get_error("scan_config_line: A letter expected at the beginning "
                                  "of a keyword!\n%s",cfg);
   for (k=0;keyname_char_ok(charptr[k]);k++);
   if (k > MAX_KEYSIZE)
      get_error("scan_config_line: Keyword name too long!\n%s",cfg);
   memmove(key,charptr,k);
   key[k] = 0;
   charptr += k;
   if (get_char() == '=') charptr++;

   *n = 0;
   while (!cfg_end_of_line_ok(get_char()))
   {
      q = get_char();
      if (cfg_quot_ok(q))
      {
         for (k=1;charptr[k]!=q && charptr[k]!=0;k++);
         if (charptr[k] != q)
            get_error("scan_config_line: Quotation not closed!\n%s",cfg);
         k++;
      }
      else
         for (k=0;keyval_char_ok(charptr[k]);k++);
      if (k > MAX_VALSIZE)
         get_error("scan_config_line: Value string too long "
                   "for '%s'!\n%s",key,cfg);
      (*n)++;
      if (*n > MAX_NVALS)
         get_error("scan_config_line: Too many values found "
                   "for '%s'!\n%s",key,cfg);
      memmove(val[*n],charptr,k);
      val[*n][k] = 0;
      charptr+=k;
   }
}

void get_next_keyword(file_type *f,keytype key,int *n,valtype val)
{
   cfglinetype cfg;

   charptr = cfg;
   do
   {
      if (fgets(cfg,sizeof(cfglinetype),f->dat) == NULL) strcpy(cfg,"END");
      if (strlen(cfg) >= MAXCFGLINE)
         get_error("Reading '%s':\nData line too long!\n%s",f->path,cfg);
      trim_right(cfg);
   }
   while (cfg_end_of_line_ok(get_char()));

   scan_config_line(charptr,key,n,val);

   if (strcasecmp(key,"END") == 0 && *n != 0)
      get_error("Reading '%s':\nNo values expected after 'END'!\n%s",
                f->path,cfg);
}

void expect_next_keyword(file_type *f,char *s,int m,keytype key,valtype val)
{
   int n;

   if (m>MAX_NVALS)
      get_error("expect_next_keyword: Too many "
                "expected values (%d)!",m);
   get_next_keyword(f,key,&n,val);
   if (strcasecmp(key,s) != 0)
      get_error("Reading '%s':\nKeyword '%s' expected "
                "('%s' found)!",f->path,s,key);
   if (n != m)
      get_error("Reading '%s':\n %d value(s) for '%s' expected "
                "(%d found)!",f->path,m,key,n);
}

void get_string_value(char *s,int maxlen,valtype val,int k)
{
   char *ptr,q;

   if (strlen(val[k]) > maxlen)
      get_error("get_string_value: String value too long (%s)!",val[k]);
   ptr = val[k];
   q = *ptr;
   if (cfg_quot_ok(q))
   {
      ptr++;
      while (*ptr != q) *(s++) = *(ptr++);
      *s = 0;
   }
   else
      strcpy(s,ptr);
}

int get_integer_value(valtype val,int k)
{
   int a;
   char *eptr;

   a = strtol(val[k],&eptr,10);
   if (eptr == val[k])
      get_error("get_integer_value: No valid characters found!");
   if (*eptr != 0)
      get_error("get_integer_value: An invalid character detected "
                "after the integer value ('%s')!",val[k]);

   return(a);
}

void get_integer_value_array(int *a,valtype val,int k,int n)
{
   int i;

   for (i=1;i<=n;i++) a[i] = get_integer_value(val,k+i-1);
}

double get_double_value(valtype val,int k)
{
   double a;
   char *eptr;

   a = strtod(val[k],&eptr);
   if (eptr == val[k])
      get_error("get_double_value: No valid characters found!");
   if (*eptr != 0)
      get_error("get_double_value: An invalid character detected "
                "after the double value ('%s')!",val[k]);

   return(a);
}

void get_double_value_array(double *a,valtype val,int k,int n)
{
   int i;

   for (i=1;i<=n;i++) a[i] = get_double_value(val,k+i-1);
}

void expect_next_string_value(file_type *f,char *s,char *x,int maxlen)
{
   keytype key;
   valtype val;

   expect_next_keyword(f,s,1,key,val);
   get_string_value(x,maxlen,val,1);
}

int expect_next_integer_value(file_type *f,char *s)
{
   keytype key;
   valtype val;

   expect_next_keyword(f,s,1,key,val);
   return(get_integer_value(val,1));
}

double expect_next_double_value(file_type *f,char *s)
{
   keytype key;
   valtype val;

   expect_next_keyword(f,s,1,key,val);
   return(get_double_value(val,1));
}

void search_text_arg(char *arg,char list[][MAX_VALSIZE+1],int *seq,char *name)
{
  int i;

  *seq = 0;
  *name = 0;
  for (i=0;*list[i]!=0;i++)
    if (strcasecmp(list[i],arg) == 0) *seq = i+1;
  if (*seq > 0) strcpy(name,list[(*seq)-1]);
}

void expect_next_text_param
  (file_type *f,char *s,txtpar *p,char list[][MAX_VALSIZE+1])
{
   char txt[MAX_VALSIZE+1];

   expect_next_string_value(f,s,txt,MAX_VALSIZE);
   search_text_arg(txt,list,&p->id,p->txt);
   if (p->id < 1)
     get_error("Bad text parameter '%s = %s' in '%s'!",s,txt,f->path);
}

int expect_yes_or_no(file_type *file,char *key)
{
   char txt[MAX_VALSIZE+1];
   int result;

   expect_next_string_value(file,key,txt,MAX_VALSIZE);
   result = -1;
   if (strcasecmp(txt,"NO") == 0)
     result = 0;
   else
     if (strcasecmp(txt,"YES") == 0) result = 1;
   if (result < 0)
     get_error("YES or NO expected for '%s' in '%s'!",key,file->path);
   return(result);
}

void open_config_file(file_type *file,char *spectrograph,char *name,...)
{
   va_list ap;
   file_type f;

   if (strlen(spectrograph) > MAX_SPECNAME)
      get_error("open_config_file: Spectrograph name too long!");

   va_start(ap,name);
   vset_file_name(&f,name,ap);
   va_end(ap);

   if (file_exists("%s.cfg",f.name))
      open_file(file,"r","%s.cfg",f.name);
   else
      open_file(file,"r","%s/spec/%s/cfg/%s.cfg",
        hrsp_dir(),spectrograph,f.name);
}

void open_keydef_cfg(file_type *file,char *spectrograph)
{
   open_config_file(file,spectrograph,"keydef");
}

void get_regre_value(regre_type *reg,valtype val,int k,int n)
{
   if (n < 1 || n > MAX_REGRE) get_error("Bad number of axes!");
   reg->naxis = n;
   get_integer_value_array(reg->deg,val,k,n);
   list_integers(reg->deg,n,reg->degarg,MAX_DEGARG);
   reg->kappa = get_double_value(val,k+n);
}

void load_spec_info(char *specname,spec_type *spec)
{
   keytype key,ekey;
   valtype val;
   int k;
   file_type file;

   open_config_file(&file,specname,specname);

   strcpy(spec->specname,specname);

   expect_next_string_value(&file,"NAME",spec->fullname,MAX_SPECNAME);
   expect_next_keyword(&file,"BOXSIZE",2,key,val);
   spec->usize = get_double_value(val,1);
   spec->vsize = get_double_value(val,2);

   spec->uscale = 1.0/spec->usize;
   spec->vscale = 1.0/spec->vsize;

   expect_next_keyword(&file,"ORDERS",2,key,val);
   spec->firstord = get_integer_value(val,1);
   spec->lastord = get_integer_value(val,2);

   spec->nord = spec->lastord - spec->firstord + 1;
   spec->mscale = 1.0/(double)spec->nord;

   spec->bincount = expect_next_integer_value(&file,"BINCOUNT");
   spec->binsize =
      (double)expect_next_integer_value(&file,"BINSIZE") * 1.0e-3;

   expect_next_keyword(&file,"MLAM",2,key,val);
   spec->mlamstart = get_double_value(val,1);
   spec->mlamend = get_double_value(val,2);

   spec->mlamrange = spec->mlamend - spec->mlamstart;
   spec->mlamscale = 1.0/spec->mlamrange;

   expect_next_keyword(&file,"REFLINE",5,key,val);
   spec->refair = get_double_value(val,1);
   spec->refvac = get_double_value(val,2);
   spec->reford = get_integer_value(val,3);
   spec->uref = get_double_value(val,4);
   spec->vref = get_double_value(val,5);

   spec->nfib = expect_next_integer_value(&file,"NFIBRES");
   if (spec->nfib < 1 || spec->nfib > MAX_FIBRES)
      get_error("Reading '%s':\nBad number of fibres (%d)!",
                 file.path,spec->nfib);
   for (k=1;k<=spec->nfib;k++)
   {
      sprintf(ekey,"FIBRE%d",k);
      expect_next_keyword(&file,ekey,2,key,val);
      spec->fib[k].dx = get_double_value(val,1) / 1000.0;
      spec->fib[k].dy = get_double_value(val,2)/ 1000.0;
   }

   expect_next_keyword(&file,"ECHDEF",9,key,val);
   get_regre_value(&spec->echdef_vcen,val,1,2);
   get_regre_value(&spec->echdef_ord,val,4,2);
   get_regre_value(&spec->echdef_ucen,val,7,2);

   expect_next_keyword(&file,"END",0,key,val);

   get_fclose(&file);
}

void load_ccd_info(char *ccdname,char *spectrograph,ccd_type *ccd)
{
   keytype key,ekey;
   valtype val;
   int k;
   file_type file;

   if (strlen(ccdname) > MAX_CCDNAME) get_error("CCD name too long!");

   open_config_file(&file,spectrograph,ccdname);

   strcpy(ccd->ccdname,ccdname);
   expect_next_string_value(&file,"NAME",ccd->fullname,MAX_CCDNAME);
   expect_next_keyword(&file,"CHIPSIZE",2,key,val);
   get_integer_value_array(ccd->chipsize,val,1,2);
   expect_next_keyword(&file,"PIXSIZE",2,key,val);
   get_double_value_array(ccd->pixsize,val,1,2);
   for (k=1;k<=2;k++) ccd->pixsize[k] /= 1000.0;
   ccd->maxpix = expect_next_double_value(&file,"MAXPIX");
   ccd->ngain = expect_next_integer_value(&file,"NGAIN");
   for (k=1;k<=ccd->ngain;k++)
   {
      sprintf(ekey,"GAIN%d",k);
      expect_next_keyword(&file,ekey,3,key,val);
      ccd->gain[k].gainset = get_integer_value(val,1);
      ccd->gain[k].invgain = get_double_value(val,2);
      ccd->gain[k].bnoise = get_double_value(val,3);
   }
   ccd->rotate = expect_next_integer_value(&file,"ROTATE");

   expect_next_text_param(&file,"FLIP",&ccd->flip,ccd_flips);

   expect_next_keyword(&file,"END",0,key,val);

   get_fclose(&file);
}

void get_observatory(obs_type *obs,char *spectrograph)
{
   file_type file;
   keytype key;
   valtype val;
   char lonz,latz;
   int lonh,lonm,latd,latm;
   double lons,lats;

   open_config_file(&file,spectrograph,"observatory");

   expect_next_string_value(&file,"NAME",obs->name,MAX_OBSNAME);

   expect_next_keyword(&file,"LONGITUDE",3,key,val);
   if (*val[1] == '-')
     lonz = '-';
   else
     lonz = '+';
   lonh = abs(get_integer_value(val,1));
   lonm = get_integer_value(val,2);
   lons = get_double_value(val,3);
   set_angle_lambda(&obs->lon,lonz,lonh,lonm,lons);

   expect_next_keyword(&file,"LATITUDE",3,key,val);
   if (*val[1] == '-')
     latz = '-';
   else
     latz = '+';
   latd = abs(get_integer_value(val,1));
   latm = get_integer_value(val,2);
   lats = get_double_value(val,3);
   set_angle_delta(&obs->lat,latz,latd,latm,lats);

   obs->height = expect_next_integer_value(&file,"HEIGHT");

   get_fclose(&file);
}

int micro(double mm)
{
  return(nint(mm*1000.0));
}

void get_plot_info(plt_type *plt,char *spectrograph)
{
   file_type file;
   keytype key;
   valtype val;
   double minsep[3];
   int i,wtot,htot;
   int totcolgap,totrowgap;

   open_config_file(&file,spectrograph,"plot");

   plt->left_margin = micro(expect_next_double_value(&file,"LEFT_MARGIN"));
   if (plt->left_margin < 0)
     get_error("get_plot_info: Negative left margin!");

   plt->right_margin = micro(expect_next_double_value(&file,"RIGHT_MARGIN"));
   if (plt->right_margin < 0)
     get_error("get_plot_info: Negative right margin!");

   plt->top_margin = micro(expect_next_double_value(&file,"TOP_MARGIN"));
   if (plt->top_margin < 0)
     get_error("get_plot_info: Negative top margin!");

   plt->bottom_margin = micro(expect_next_double_value(&file,"BOTTOM_MARGIN"));
   if (plt->bottom_margin < 0)
     get_error("get_plot_info: Negative bottom margin!");

   plt->numrows = expect_next_integer_value(&file,"NUMROWS");
   if (plt->numrows < 1 || plt->numrows > 8)
     get_error("get_plot_info: Number of rows out of range!");

   plt->numcols = expect_next_integer_value(&file,"NUMCOLS");
   if (plt->numcols < 1 || plt->numcols > 3)
     get_error("get_plot_info: Number of columns out of range!");

   plt->rowgap = micro(expect_next_double_value(&file,"ROWGAP"));
   if (plt->rowgap < 0)
     get_error("get_plot_info: Negative row gap!");

   plt->colgap = micro(expect_next_double_value(&file,"COLGAP"));
   if (plt->colgap < 0)
     get_error("get_plot_info: Negative column gap!");

   plt->linewidth = micro(expect_next_double_value(&file,"LINEWIDTH"));
   if (plt->linewidth < 0)
     get_error("get_plot_info: Negative line width!");

   plt->fontsize = micro(expect_next_double_value(&file,"FONTSIZE"));
   if (plt->fontsize < 0)
     get_error("get_plot_info: Negative font size!");

   plt->ticksize = micro(expect_next_double_value(&file,"TICKSIZE"));
   if (plt->ticksize < 0)
     get_error("get_plot_info: Negative tick size!");

   expect_next_keyword(&file,"MINSEP",2,key,val);
   get_double_value_array(minsep,val,1,2);
   for (i=1;i<=2;i++)
   {
     plt->minsep[i] = micro(minsep[i]);
     if (plt->minsep[i] < 0)
       get_error("get_plot_info: Negative minimum separation!");
   }

   plt->global_limits = expect_yes_or_no(&file,"GLOBAL_LIMITS");
   plt->window_limits = expect_yes_or_no(&file,"WINDOW_LIMITS");

   get_fclose(&file);

   plt->plots_per_page = plt->numrows * plt->numcols;

   wtot = 210000 - plt->left_margin - plt->right_margin;
   if (wtot < 10000)
     get_error("get_plot_info: Left or right margin too wide!");

   htot = 297000 - plt->top_margin - plt->bottom_margin;
   if (htot < 10000)
     get_error("get_plot_info: Top or bottom margin too wide!");

   totcolgap = (plt->numcols - 1) * plt->colgap;
   totrowgap = (plt->numrows - 1) * plt->rowgap;

   plt->boxw = (wtot - totcolgap) / plt->numcols;
   if (plt->boxw < 10000)
     get_error("get_plot_info: Too much horizontal white space!");

   plt->boxh = (htot - totrowgap) / plt->numrows;
   if (plt->boxh < 10000)
     get_error("get_plot_info: Too much vertical white space!");

   plt->fontcent = (plt->fontsize * 7) / 20;
}

void get_parameters(param_type *par,char *spectrograph)
{
   file_type file;
   keytype key;
   valtype val;

   open_config_file(&file,spectrograph,"param");

   expect_next_keyword(&file,"TRANSFORM",5,key,val);
   get_integer_value_array(par->trans_rad,val,1,2);
   get_regre_value(&par->trans,val,3,2);

   expect_next_keyword(&file,"ECHDEF",10,key,val);
   par->echdef_ncols = get_integer_value(val,1);
   get_regre_value(&par->echdef_ycen,val,2,2);
   get_regre_value(&par->echdef_fwhm,val,5,2);
   get_regre_value(&par->echdef_ord,val,8,2);

   expect_next_text_param(&file,"BACK_METHOD",&par->back_method,back_methods);

   expect_next_keyword(&file,"ECHBACK",4,key,val);
   par->echback_ncols = get_integer_value(val,1);
   get_regre_value(&par->echback,val,2,2);

   expect_next_keyword(&file,"MEDFILT",2,key,val);
   par->medfilt_rx = get_integer_value(val,1);
   par->medfilt_ry = get_integer_value(val,2);

   expect_next_keyword(&file,"CRFILT",6,key,val);
   par->crfilt_nboxes = get_integer_value(val,1);
   par->crfilt_rslit = get_integer_value(val,2);
   par->crfilt_nsigma = get_double_value(val,3);
   get_regre_value(&par->crfilt,val,4,1);
   par->crfilt_npasses = get_integer_value(val,6);

   expect_next_keyword(&file,"FLAT_FILT",3,key,val);
   par->flat_rfilt = get_integer_value(val,1);
   get_regre_value(&par->flat,val,2,1);

   par->flat_smooth = expect_yes_or_no(&file,"FLAT_SMOOTH");

   expect_next_keyword(&file,"FLAT_THRES",2,key,val);
   par->flat_thres_full = get_double_value(val,1);
   par->flat_thres_ext = get_double_value(val,2);

   expect_next_text_param(&file,"FLAT_METHOD",&par->flat_method,flat_methods);

   par->extract_slit = expect_next_double_value(&file,"EXTRACT_SLIT");

   expect_next_keyword(&file,"END",0,key,val);

   get_fclose(&file);
}

void get_dispsol_info(dispsol_type *dsol,char *spectrograph)
{
   file_type file;
   keytype key;
   valtype val;

   open_config_file(&file,spectrograph,"dispsol");

   expect_next_keyword(&file,"ORD",2,key,val);
   dsol->firstord = get_integer_value(val,1);
   dsol->lastord = get_integer_value(val,2);

   expect_next_keyword(&file,"MLAM",2,key,val);
   dsol->mlamstart = get_double_value(val,1);
   dsol->mlamend = get_double_value(val,2);

   expect_next_text_param(&file,"WAVE",&dsol->wave,wave_list);

   expect_next_keyword(&file,"REG",3,key,val);
   get_regre_value(&dsol->reg,val,1,2);

   expect_next_keyword(&file,"END",0,key,val);

   get_fclose(&file);
}

void get_weights(file_type *f,int *sel,double *w)
{
   char weights[80];
   file_type dfile;
   char drow[80];
   int ord;
   double wval;

   expect_next_string_value(f,"WEIGHTS",weights,78);

   if (!alpha_ok(*weights)) get_error
     ("Bad WEIGHTS file name in '%s' (a letter expected)!",f->path);

   if (strcasecmp(weights,"NONE") == 0)
   {
      for (ord=1;ord<=MAXORD;ord++) w[ord] = 1.0;
   }
   else
   {
      if (!file_exists(weights)) get_error
       ("WEIGHTS file '%s' not found in '%s' !",weights,f->path);
      for (ord=1;ord<=MAXORD;ord++) w[ord] = -1.0;
      set_file_name(&dfile,weights);
      get_fopen(&dfile,"r");
      while(fgets(drow,64,dfile.dat) != NULL)
      {
         if (strlen(drow) > 60)
            get_error("Data row too long in '%s'!",dfile.path);
         if (sscanf(drow," %d %lf",&ord,&wval) != 2)
            get_error("Two numerical values per data row "
                      "expected in '%s'!",dfile.path);
         if (ord < 1 || ord > MAXORD)
            get_error("Order number out of range in '%s'!",dfile.path);
         if (wval <= 0.0)
            get_error("Invalid weight in '%s'!",dfile.path);
         w[ord] = wval;
      }
      get_fclose(&dfile);
   }

   for (ord=1;ord<=MAXORD;ord++)
      if (sel[ord]==1 && w[ord]<0.0)
            get_error("Missing weights for some orders in '%s'!",dfile.path);
}

void get_radius(file_type *f,int *sel,txtpar *method,int *r)
{
   char radius[80],*ptr;
   file_type dfile;
   char drow[80];
   int ord,rval,rcons;

   expect_next_string_value(f,"RADIUS",radius,78);

   if (!alphanum_ok(*radius)) get_error
     ("Bad RADIUS in '%s' (a digit or a letter expected)!",f->path);

   if (digit_ok(*radius))
   {
     for (ptr=radius;digit_ok(*ptr);ptr++);
     if (*ptr != 0)
       get_error("Bad RADIUS in '%s' (a digit expected)!",f->path);
   }

   memset(r,0,(MAXORD+1)*sizeof(int));

   if (strcasecmp(radius,"NONE") == 0)
   {
     if (method->id == METHOD_GAUSS ||
         method->id == METHOD_PARABOLA)
       get_error("Radius data missing for method '%s'!",method->txt);
   }
   else if (digit_ok(*radius))
   {
      rcons = atoi(radius);
      if (rcons < 1 || rcons > 32)
        get_error("RADIUS out of limits in '%s' (must be 1-32)!",f->path);
      for (ord=1;ord<=MAXORD;ord++) r[ord] = rcons;
   }
   else
   {
      if (!file_exists(radius)) get_error
       ("RADIUS file '%s' not found in '%s' !",radius,f->path);
      for (ord=1;ord<=MAXORD;ord++) r[ord] = -1;
      set_file_name(&dfile,radius);
      get_fopen(&dfile,"r");
      while(fgets(drow,64,dfile.dat) != NULL)
      {
         if (strlen(drow) > 60)
            get_error("Data row too long in '%s'!",dfile.path);
         if (sscanf(drow," %d %d",&ord,&rval) != 2)
            get_error("Two numerical values per data row "
                      "expected in '%s'!",dfile.path);
         if (ord < 1 || ord > MAXORD)
            get_error("Order number out of range in '%s'!",dfile.path);
         if (rval < 2) get_error("Radius too small in '%s'!",dfile.path);
         if (rval > MAXGR) get_error("Radius too big in '%s'!",dfile.path);
         r[ord] = rval;
      }
      get_fclose(&dfile);
      for (ord=1;ord<=MAXORD;ord++)
         if (sel[ord]==1 && r[ord]<0)
            get_error("Missing radius for some orders in '%s'!",dfile.path);
   }
}

void get_radvel_info(radvel_type *rvcfg,char *spectrograph)
{
   file_type file;
   char orders[80];
   int count;
   keytype key;
   valtype val;

   open_config_file(&file,spectrograph,"radvel");

   expect_next_string_value(&file,"ORDERS",orders,78);
   count = collect_integer_list(orders,rvcfg->sel,MAXORD,"Order number");
   if (count < 1) get_error("No orders selected in '%s'!",file.path);

   expect_next_keyword(&file,"MLAM",2,key,val);
   rvcfg->mlamstart = get_double_value(val,1);
   rvcfg->mlamend = get_double_value(val,2);
   if (rvcfg->mlamend < rvcfg->mlamstart)
     get_error("Bad MLAM values in '%s'!",file.path);

   expect_next_text_param(&file,"ORDPREP",&rvcfg->ordprep,rv_preps);

   rvcfg->cosbell = expect_next_integer_value(&file,"COSBELL");

   expect_next_text_param(&file,"METHOD",&rvcfg->method,rv_methods);

   rvcfg->checkrad = expect_next_integer_value(&file,"CHECKRAD");
   if (rvcfg->checkrad < 0 || rvcfg->checkrad > 8)
     get_error("Bad CHECKRAD value in '%s'!",file.path);

   get_weights(&file,rvcfg->sel,rvcfg->w);
   get_radius(&file,rvcfg->sel,&rvcfg->method,rvcfg->r);

   expect_next_keyword(&file,"END",0,key,val);

   get_fclose(&file);
}

void get_shift_info(shift_type *sft,char *spectrograph)
{
   file_type file;
   char orders[80];
   int count;
   keytype key;
   valtype val;

   open_config_file(&file,spectrograph,"shift");

   expect_next_string_value(&file,"ORDERS",orders,78);
   count = collect_integer_list(orders,sft->sel,MAXORD,"Order number");
   if (count < 1) get_error("No orders selected in '%s'!",file.path);

   expect_next_keyword(&file,"XBOUNDS",2,key,val);
   sft->xleft = get_double_value(val,1);
   sft->xright = get_double_value(val,2);
   if (sft->xleft < 1 || sft->xright < sft->xleft)
     get_error("Bad XBOUNDS in '%s'!",file.path);

   expect_next_text_param(&file,"ORDPREP",&sft->ordprep,rv_preps);

   sft->cosbell = expect_next_integer_value(&file,"COSBELL");

   expect_next_text_param(&file,"METHOD",&sft->method,rv_methods);

   sft->checkrad = expect_next_integer_value(&file,"CHECKRAD");
   if (sft->checkrad < 0 || sft->checkrad > 8)
     get_error("Bad CHECKRAD value in '%s'!",file.path);

   get_weights(&file,sft->sel,sft->w);
   get_radius(&file,sft->sel,&sft->method,sft->r);

   expect_next_keyword(&file,"END",0,key,val);

   get_fclose(&file);
}

void read_lsout(char * str)
{
   file_type lsout;
   char row[256];

   open_file(&lsout,"r","ls.out");
   if (fgets(row,256,lsout.dat) == NULL)
      get_error("A line of text expected in 'ls.out'!");
   if (strlen(row) == 255) get_error("Text line too long in 'ls.out'!");
   if (sscanf(row,"%255s",str) != 1)
      get_error("A character string expected in 'ls.out'!");
   if (fgets(row,256,lsout.dat) != NULL)
      get_error("A single line of text expected in 'ls.out'!");
   get_fclose(&lsout);
}

int matched_ok(char *f,char *t)
{
   char tmp[256],str[256];

   if (strlen(f) > MAX_FILE_NAME)
      get_error("matched_ok: File name '%s' too long!",f);
   sprintf(tmp,"TEMPORARY-%s",f);
   get_system("\\rm -f TEMPORARY-*.[Ff][Ii][Tt]");
   get_system("touch %s",tmp);
   get_system("\\rm -f ls.out");
   go_system("\\ls -1 TEMPORARY-%s >ls.out 2>&1",t);
   read_lsout(str);
   get_system("\\rm -f ls.out %s",tmp);
   return(strcmp(str,tmp) == 0);
}

void read_matching_cfg
  (char *cfgname,char *fname,char **keys,char **vals,int *n)
{
   file_type cfg;
   keytype akey;
   valtype aval;
   int nval,match;

   *n = 0;
   white_brackets = 0;
   open_file(&cfg,"r",cfgname);
   expect_next_keyword(&cfg,"file",1,akey,aval);
   while(strcmp(akey,"file") == 0)
   {
      match = matched_ok(fname,aval[1]);
      while(1)
      {
         get_next_keyword(&cfg,akey,&nval,aval);
         if (strcmp(akey,"file") == 0 || strcmp(akey,"END") == 0) break;
         if (match)
         {
            if (*n >= MAXCFGROWS) get_error
               ("Configuratioon file '%s' has too many rows!",cfg.name);
            strcpy(keys[*n],akey);
            strcpy(vals[*n],aval[1]);
            (*n)++;
         }
      }
   }
   get_fclose(&cfg);
   white_brackets = 1;
}

char **allocate_cfglines(void)
{
   char **a;
   int m,n,i;

   m = MAXCFGROWS;
   n = MAXCFGLINE + 1;

   a = (char **)get_space(m * sizeof(char *));
   a[0] = (char *)get_space(m * n);
   for (i=1;i<MAXCFGROWS;i++) a[i] = a[i-1] + n;

   return(a);
}

void free_cfglines(char **list)
{
   get_free((void *)list[0]);
   get_free((void *)list);
}
