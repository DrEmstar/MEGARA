#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#include "general.h"
#include "config.h"

static char *charptr;

int white_space_ok(char c)
{
   if (c==' ' || c=='\t' || c=='[' || c==']')
      return(1);
   else
      return(0);
}

int quot_ok(char c)
{
   if (c=='\x27' || c=='\x22')
      return(1);
   else
      return(0);
}

int end_of_line_ok(char c)
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
   if (ascii_ok(c) && !end_of_line_ok(c) && !white_space_ok(c))
      return(1);
   else
      return(0);
}

char get_char(void)
{
   while (white_space_ok(*charptr)) charptr++;
   return(*charptr);
}

char next_char(void)
{
   charptr++;
   return(get_char());
}

void scan_config_line(char *s,keytype key,int *n,valtype val)
{
   char cfg[MAXCFGLINE],q;
   int k;

   if (strlen(s) > MAXCFGLINE-1)
      get_error("scan_config_line: Line too long!");
   strcpy(cfg,s);
   remove_trailing_blanks(cfg);
   charptr = cfg;
   if (end_of_line_ok(get_char()))
      get_error("scan_config_line: Line is empty!");

   if (!alpha_ok(get_char()))
      get_error("scan_config_line: A letter expected at the beginning "
                                  "of a keyword!\n%s",cfg);
   for (k=0;keyname_char_ok(charptr[k]);k++);
   if (k > MAX_KEYSIZE)
      get_error("scan_config_line: Keyword name too long!\n%s",cfg);
   memcpy(key,charptr,k);
   key[k] = 0;
   charptr += k;
   if (get_char() == '=') charptr++;

   *n = 0;
   while (!end_of_line_ok(get_char()))
   {
      q = get_char();
      if (quot_ok(q))
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
      memcpy(val[*n],charptr,k);
      val[*n][k] = 0;
      charptr+=k;
   }

   if (*n > MAX_NVALS)
      get_error("scan_config_line: Too many values found "
                "for '%s'!\n%s",key,cfg);
}

void get_next_keyword(file_type *f,keytype key,int *n,valtype val)
{
   char cfg[MAXCFGLINE];

   charptr = cfg;
   do
   {
      if (fgets(cfg,MAXCFGLINE,f->dat) == NULL) strcpy(cfg,"END");
      if (strlen(cfg) >= MAXCFGLINE-1)
         get_error("Reading '%s':\nData line too long!\n%s",f->path,cfg);
      remove_trailing_blanks(cfg);
   }
   while (end_of_line_ok(get_char()));

   scan_config_line(charptr,key,n,val);

   if (strcasecmp(key,"END") == 0 && *n != 0)
      get_error("Reading '%s':\nNo values expected after 'END'!\n%s",
                f->path,cfg);
}

void expect_next_keyword(file_type *f,char *s,int m,keytype key,valtype val)
{
   int n;

   if (m>MAX_NVALS)
      get_error("expect_next_keyword: Bad "
                "expected number of values (%d)!",m);
   get_next_keyword(f,key,&n,val);
   if (strcasecmp(key,s) != 0)
      get_error("Reading '%s':\nKeyword '%s' expected "
                "('%s' found)!",f->path,s,key);
   if (n != m)
      get_error("Reading '%s':\n %d value(s) for '%s' expected "
                "(%d found)!",f->path,m,key,n);
}

int expect_next_keyword_list(file_type *f,key_info *s,keytype key,valtype val)
{
   int n,k;

   get_next_keyword(f,key,&n,val);
   for (k=1;*s[k].name!=0 && strcasecmp(s[k].name,key);k++);
   if (*s[k].name == 0)
      get_error("Reading '%s':\nKeyword '%s' not expected!",f->path,key);
   if (n != s[k].nval)
      get_error("Reading '%s':\n %d value(s) for '%s' expected "
                "(%d found)!",f->path,s[k].nval,key,n);
   return(k);
}

void get_string_value(char *s,int maxlen,valtype val,int k)
{
   char *ptr,q;

   if (strlen(val[k]) > maxlen)
      get_error("get_string_value: String value too long (%s)!",val[k]);
   ptr = val[k];
   q = *ptr;
   if (quot_ok(q))
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
      get_error("get_integer_value: A bad character detected "
                "after the integer value ('%s')!",val[k]);

   return(a);
}

void get_integer_value_array(int *a,valtype val,int k,int n)
{
   int i;

   for (i=1;i<=n;i++) a[i] = get_integer_value(val,k+i-1);
}

void get_box_type_value(box_type *box,valtype val)
{
   box->ax = get_integer_value(val,1);
   box->ay = get_integer_value(val,2);
   box->bx = get_integer_value(val,3);
   box->by = get_integer_value(val,4);
}

double get_double_value(valtype val,int k)
{
   double a;
   char *eptr;

   a = strtod(val[k],&eptr);
   if (eptr == val[k])
      get_error("get_double_value: No valid characters found!");
   if (*eptr != 0)
      get_error("get_double_value: A bad character detected "
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

void open_config_file(file_type *file,char *dir,char *name,...)
{
   va_list ap;
   file_type f;

   va_start(ap,name);
   vset_file_name(&f,name,ap);
   va_end(ap);

   if (file_exists("%s.cfg",f.name))
      open_file(file,"r","%s.cfg",f.name);
   else
      open_file(file,"r","%s/%s/%s.cfg",herc_dir(),dir,f.name);
}

void open_keydef_cfg(file_type *file)
{
   open_config_file(file,"key","keydef");
}

void get_regre_value(regre_type *reg,valtype val,int k,int n)
{
   if (n < 1 || n > MAX_REGRE) get_error("Bad number of axes!");
   reg->naxis = n;
   get_integer_value_array(reg->deg,val,k,n);
   list_integers(reg->deg,n,reg->degarg,MAX_DEGARG);
   reg->kappa = get_double_value(val,k+n);
}

void get_herc_info(herc_type *herc)
{
   keytype key,ekey;
   valtype val;
   int k;
   file_type file;

   open_config_file(&file,"def","hercules");

   expect_next_keyword(&file,"BOXSIZE",2,key,val);
   herc->usize = get_double_value(val,1);
   herc->vsize = get_double_value(val,2);

   herc->uscale = 1.0/herc->usize;
   herc->vscale = 1.0/herc->vsize;

   expect_next_keyword(&file,"ORDERS",2,key,val);
   herc->firstord = get_integer_value(val,1);
   herc->lastord = get_integer_value(val,2);

   herc->nord = herc->lastord - herc->firstord + 1;
   herc->mscale = 1.0/(double)herc->nord;

   expect_next_keyword(&file,"MLAM",2,key,val);
   herc->mlamstart = get_double_value(val,1);
   herc->mlamend = get_double_value(val,2);

   herc->mlamrange = herc->mlamend - herc->mlamstart;
   herc->mlamscale = 1.0/herc->mlamrange;

   expect_next_keyword(&file,"REFLINE",5,key,val);
   herc->refair = get_double_value(val,1);
   herc->refvac = get_double_value(val,2);
   herc->reford = get_integer_value(val,3);
   herc->uref = get_double_value(val,4);
   herc->vref = get_double_value(val,5);

   herc->nreg = expect_next_integer_value(&file,"NREGIONS");
   if (herc->nreg < 1 || herc->nreg > MAX_REGIONS)
      get_error("Reading '%s':\nBad number of regions (%d)!",
                 file.path,herc->nreg);
   for (k=1;k<=herc->nreg;k++)
   {
      sprintf(ekey,"REGCEN%d",k);
      expect_next_keyword(&file,ekey,4,key,val);
      herc->reg[k].ax = get_double_value(val,1);
      herc->reg[k].ay = get_double_value(val,2);
      herc->reg[k].bx = get_double_value(val,3);
      herc->reg[k].by = get_double_value(val,4);
   }

   herc->nfib = expect_next_integer_value(&file,"NFIBRES");
   if (herc->nfib < 1 || herc->nfib > MAX_FIBRES)
      get_error("Reading '%s':\nBad number of fibres (%d)!",
                 file.path,herc->nfib);
   for (k=1;k<=herc->nfib;k++)
   {
      sprintf(ekey,"FIBRE%d",k);
      expect_next_keyword(&file,ekey,2,key,val);
      herc->fib[k].dx = get_double_value(val,1) / 1000.0;
      herc->fib[k].dy = get_double_value(val,2)/ 1000.0;
   }

   expect_next_keyword(&file,"ECHDEF",6,key,val);
   get_regre_value(&herc->echdef_vcen,val,1,2);
   get_regre_value(&herc->echdef_ord,val,4,2);

   expect_next_keyword(&file,"END",0,key,val);

   get_fclose(&file);
}

void load_ccdinfo(char *ccdname,ccd_type *ccd)
{
   keytype key,ekey;
   valtype val;
   int k;
   file_type file;

   if (strlen(ccdname) > MAX_CCDNAME) get_error("CCD name too long!");

   open_config_file(&file,"ccd",ccdname);

   expect_next_string_value(&file,"NAME",ccd->name,MAX_CCDNAME);
   expect_next_keyword(&file,"CHIPSIZE",2,key,val);
   get_integer_value_array(ccd->chipsize,val,1,2);
   expect_next_keyword(&file,"PIXSIZE",2,key,val);
   get_double_value_array(ccd->pixsize,val,1,2);
   for (k=1;k<=2;k++) ccd->pixsize[k] /= 1000.0;
   ccd->maxpix = expect_next_double_value(&file,"MAXPIX");
   ccd->nbadlin = expect_next_integer_value(&file,"NBADLIN");
   for (k=1;k<=ccd->nbadlin;k++)
   {
      sprintf(ekey,"BADLIN%d",k);
      expect_next_keyword(&file,ekey,4,key,val);
      get_box_type_value(&ccd->badlin[k],val);
   }
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
   expect_next_string_value(&file,"FLIP",ccd->flip,2);

   expect_next_keyword(&file,"END",0,key,val);

   get_fclose(&file);
}

void get_parameters(param_type *par,int reg,int fib)
{
   file_type file;
   keytype key;
   valtype val;

   open_config_file(&file,"par","par-%d-%d",reg,fib);

   par->reg = expect_next_integer_value(&file,"REGION");
   if (par->reg != reg)
      get_error("Reading '%s':\nBad CCD position region number (%d)!",
                 file.path,par->reg);

   par->fib = expect_next_integer_value(&file,"FIBRE");
   if (par->fib != fib)
      get_error("Reading '%s':\nBad fibre number (%d)!",file.path,par->fib);

   expect_next_keyword(&file,"TRANSFORM",5,key,val);
   get_integer_value_array(par->trans_rad,val,1,2);
   get_regre_value(&par->trans,val,3,2);

   expect_next_keyword(&file,"ECHDEF",10,key,val);
   par->echdef_ncols = get_integer_value(val,1);
   get_regre_value(&par->echdef_ycen,val,2,2);
   get_regre_value(&par->echdef_fwhm,val,5,2);
   get_regre_value(&par->echdef_ord,val,8,2);

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

   par->extract_slit = expect_next_double_value(&file,"EXTRACT_SLIT");

   expect_next_keyword(&file,"MLAM",2,key,val);
   par->mlamstart = get_double_value(val,1);
   par->mlamend = get_double_value(val,2);

   par->mlamrange = par->mlamend - par->mlamstart;

   par->nbins = expect_next_integer_value(&file,"NBINS");

   expect_next_keyword(&file,"END",0,key,val);

   get_fclose(&file);
}

void get_dispsol_info(dispsol_type *dsol,int id)
{
   file_type file;
   keytype key;
   valtype val;

   open_config_file(&file,"dispsol","dispsol-%d",id);

   expect_next_keyword(&file,"ORD",2,key,val);
   dsol->firstord = get_integer_value(val,1);
   dsol->lastord = get_integer_value(val,2);

   expect_next_keyword(&file,"MLAM",2,key,val);
   dsol->mlamstart = get_double_value(val,1);
   dsol->mlamend = get_double_value(val,2);

   expect_next_string_value(&file,"WAVE",dsol->wave,6);

   expect_next_keyword(&file,"REG",3,key,val);
   get_regre_value(&dsol->reg,val,1,2);

   expect_next_keyword(&file,"END",0,key,val);

   get_fclose(&file);
}

void get_radvel_info(radvel_type *rvel,int id)
{
   file_type file;
   keytype key,ekey;
   valtype val;
   int k;
   file_type dfile;
   char drow[80];
   int ord,rval;
   double wval;

   open_config_file(&file,"radvel","radvel-%d",id);

   rvel->nblock = expect_next_integer_value(&file,"BLOCK");
   for (k=1;k<=rvel->nblock;k++)
   {
      sprintf(ekey,"BLOCK%d",k);
      expect_next_keyword(&file,ekey,2,key,val);
      rvel->block[k].first = get_integer_value(val,1);
      rvel->block[k].last = get_double_value(val,2);      
   }

   rvel->nomit = expect_next_integer_value(&file,"OMIT");
   for (k=1;k<=rvel->nomit;k++)
   {
      sprintf(ekey,"OMIT%d",k);
      expect_next_keyword(&file,ekey,2,key,val);
      rvel->omit[k].first = get_integer_value(val,1);
      rvel->omit[k].last = get_double_value(val,2);      
   }

   for (k=0;k<MAXORD;k++) rvel->rvok[k] = 0;

   for (k=1;k<=rvel->nblock;k++)
      for (ord=rvel->block[k].first;ord<=rvel->block[k].last;ord++)
         rvel->rvok[ord] = 1;

   for (k=1;k<=rvel->nomit;k++)
      for (ord=rvel->omit[k].first;ord<=rvel->omit[k].last;ord++)
         rvel->rvok[ord] = 0;

   expect_next_keyword(&file,"MLAM",2,key,val);
   rvel->mlamstart = get_double_value(val,1);
   rvel->mlamend = get_double_value(val,2);

   rvel->cosbell = expect_next_integer_value(&file,"COSBELL");

   expect_next_string_value(&file,"METHOD",rvel->method,12);

   if (strcasecmp(rvel->method,"GAUSS") != 0 && 
       strcasecmp(rvel->method,"PARABOLA") != 0 &&
       strcasecmp(rvel->method,"SPLINE") != 0 &&
       strcasecmp(rvel->method,"PEAK") != 0)
          get_error("Bad METHOD in RV configuration file!");

   expect_next_string_value(&file,"WEIGHTS",rvel->weights,78);

   if (strcasecmp(rvel->weights,"NONE") == 0)
   {
      for (k=0;k<MAXORD;k++) rvel->w[k] = 1.0;
   }
   else
   {
      for (k=0;k<MAXORD;k++) rvel->w[k] = -1.0;
      set_file_name(&dfile,rvel->weights);
      get_fopen(&dfile,"r");
      while(fgets(drow,64,dfile.dat) != NULL)
      {
         if (strlen(drow) > 60)
            get_error("Data row too long in `%s'!",dfile.path);
         if (sscanf(drow," %d %lf",&ord,&wval) != 2)
            get_error("Two numerical values per data row "
                      "expected in `%s'!",dfile.path);
         rvel->w[ord] = wval;
      }
      get_fclose(&dfile);
   }

   for (k=0;k<MAXORD;k++)
      if (rvel->rvok[k]==1 && rvel->w[k]<0.0)
            get_error("Missing weights for some orders in `%s'!",dfile.path);

   expect_next_string_value(&file,"RADIUS",rvel->radius,78);

   if (strcasecmp(rvel->radius,"NONE") == 0)
   {
      if (strcasecmp(rvel->method,"SPLINE") != 0 && 
          strcasecmp(rvel->method,"PEAK") != 0)
         get_error("Radius data missing for method `%s'!",rvel->method);
   }
   else
   {
      for (k=0;k<MAXORD;k++) rvel->r[k] = -1;
      set_file_name(&dfile,rvel->radius);
      get_fopen(&dfile,"r");
      while(fgets(drow,64,dfile.dat) != NULL)
      {
         if (strlen(drow) > 60)
            get_error("Data row too long in `%s'!",dfile.path);
         if (sscanf(drow," %d %d",&ord,&rval) != 2)
            get_error("Two numerical values per data row "
                      "expected in `%s'!",dfile.path);
         if (rval < 2) get_error("Radius too small in `%s'!",dfile.path);
         if (rval > MAXGR) get_error("Radius too big in `%s'!",dfile.path);
         rvel->r[ord] = rval;
      }
      get_fclose(&dfile);
   }

   for (k=0;k<MAXORD;k++)
      if (rvel->rvok[k]==1 && rvel->r[k]<0)
            get_error("Missing radius for some orders in `%s'!",dfile.path);

   expect_next_keyword(&file,"END",0,key,val);

   get_fclose(&file);
}
