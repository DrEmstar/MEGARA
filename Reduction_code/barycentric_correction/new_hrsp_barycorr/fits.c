#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#include "general.h"
#include "numeric.h"
#include "fits.h"

static char primary[] = "SIMPLE  =                    T ";
static char bintable[] = "XTENSION= 'BINTABLE' ";

/* --------------------------- FITS headers ---------------------------- */

void clear_fits_header(fits_head *head)
{
   memset(head->dat,0x20,HEAD_BYTES);
   strncpy(head->dat,"END",3);
   head->nrec = head->ncard = 1;
}

void clear_primary_header(fits_file *fits)
{
   clear_fits_header(&fits->head);
}

void clear_extension_header(fits_file *fits)
{
   clear_fits_header(&fits->extend);
}

/* --------------------------- FITS keywords --------------------------- */

int fits_key_char_ok(char c)
{
   if (upper_case_ok(c) || digit_ok(c) || c=='_' || c=='-')
      return(-1);
   else
      return(0);
}

int fits_keyword_name_ok(char *key)
{
  int i,n;

  n = strlen(key);
  if (n < 1 || n > MAX_FITS_KEYLEN) return(0);
  if (!upper_case_ok(*key)) return(0);
  for (i=0;i<n;i++)
     if (!fits_key_char_ok(key[i])) return(0);
  return(-1);
}

int extract_fits_keyname(char *s,fitskey_type key)
{
   int k,len;

   len = strlen(s);
   memset(key,0,sizeof(fitskey_type));
   for (k=0;k<MAX_FITS_KEYLEN && k<len;k++)
   {
      if (!fits_key_char_ok(s[k])) break;
      key[k] = s[k];
   }
   if (fits_key_char_ok(s[k]))
   {
      s[MAX_FITS_KEYLEN+8] = 0;
      return(-1);
   }

   return(0);
}

char *normalize_fits_keyname(fitskey_type key,fitskey_type normkey)
{
   int k;

   memmove(normkey,key,sizeof(fitskey_type));
   for (k=strlen(normkey);k<MAX_FITS_KEYLEN;k++) normkey[k] = ' ';
   normkey[MAX_FITS_KEYLEN] = 0;
   return(normkey);
}

char *collect_keyname(fits_file *fits,fits_head *head,int row,
                      fitskey_type key)
{
   char *card;

   card = head->dat + (row-1)*CARD_BYTES;
   if(extract_fits_keyname(card,key) != 0)
      get_error("FITS keyword name too long in '%s': '%s'",
                 fits->file.path,card);
   return(card);
}

int find_fits_card(fits_head *head,char *text)
{
   char *card;
   int row;

   for (row=1,card=head->dat;
        !same_case_text(card,text) && row<=head->ncard;
        row++,card+=CARD_BYTES);

   if (row > head->ncard) row=0;

   return(row);
}

void extract_fits_card(fits_head *head,int seq,char *card)
{
  memmove(card,head->dat+(seq-1)*CARD_BYTES,CARD_BYTES);
  card[CARD_BYTES] = 0;
}

int find_fits_comment(fits_head *head,char *com)
{
   char text[CARD_BYTES+1];

   if (strlen(com) > 71) get_error("find_fits_comment: Comment too long!");
   sprintf(text,"COMMENT  %s",com);
   return(find_fits_card(head,text));
}

int find_fits_keyword(fits_head *head,fitskey_type key)
{
   fitskey_type normkey;
   char text[MAX_FITS_KEYLEN+2];

   normalize_fits_keyname(key,normkey);
   sprintf(text,"%s=",normkey);
   return(find_fits_card(head,text));
}

int get_fits_keyword(fits_file *fits,fits_head *head,fitskey_type key)
{
   int row;

   row = find_fits_keyword(head,key);
   if (row < 1)
      get_error("Keyword '%s' not found in '%s'!",key,fits->file.path);

   return(row);
}

void read_keyword_value(fits_file *fits,fits_head *head,int row,
                        fitsval_type val)
{
   char *card,*limit,*first,*last;
   fitskey_type key;
   int len;

   card = head->dat + (row-1)*CARD_BYTES;
   limit = card + CARD_BYTES;
   collect_keyname(fits,head,row,key);
   memset(val,0,MAX_FITS_VALLEN+1);

   for (first=card+MAX_FITS_KEYLEN+2;first<limit;first++)
      if (*first != ' ') break;
   if (first == limit || *first == '/') return;

   if (*first == '\x27')
   {
      for (last=first+1;last<limit;last++) if (*last == '\x27') break;
      if (last == limit)
         get_error("Closing quotation mark not found for keyword '%s' "
                   "in '%s'!",key,fits->file.path);
   }
   else
   {
      for (last=first+1;last<limit;last++) if (*last == ' ') break;
      if (last < limit) last--;
   }
   len = (int)(last - first) + 1;
   memmove(val,first,len);
   val[len] = 0;
}

void read_keyword_textual(fits_file *fits,fits_head *head,int row,
                          char *dest,int maxlen)
{
   fitskey_type key;
   fitsval_type val;
   int len;

   collect_keyname(fits,head,row,key);
   read_keyword_value(fits,head,row,val);
   if (*val != '\x27')
      get_error("A character string value expected for '%s' in '%s'!",
                 key,fits->file.path);
   len = strlen(val)-2;
   if (len > maxlen)
      get_error("A character string too long for '%s' in '%s'!",
                 key,fits->file.path);
   memmove(dest,val+1,len);
   dest[len] = 0;
}

double read_keyword_double(fits_file *fits,fits_head *head,int row)
{
   fitskey_type key;
   fitsval_type val;
   char *eptr;
   double a;

   collect_keyname(fits,head,row,key);
   read_keyword_value(fits,head,row,val);
   a = strtod(val,&eptr);
   if (eptr == val || *eptr != 0)
      get_error("A numerical value expected for '%s' in '%s'!",
                key,fits->file.path);
   return(a);
}

int read_keyword_float(fits_file *fits,fits_head *head,int row)
{
   return((float)read_keyword_double(fits,head,row));
}

int read_keyword_integer(fits_file *fits,fits_head *head,int row)
{
   fitskey_type key;
   fitsval_type val;
   char *eptr;
   int a;

   collect_keyname(fits,head,row,key);
   read_keyword_value(fits,head,row,val);
   a = strtol(val,&eptr,10);
   if (eptr == val || *eptr != 0)
      get_error("An integer value expected for '%s' in '%s'!",
                 key,fits->file.path);
   return(a);
}

char read_keyword_logical(fits_file *fits,fits_head *head,int row)
{
   fitskey_type key;
   fitsval_type val;

   collect_keyname(fits,head,row,key);
   read_keyword_value(fits,head,row,val);
   if (strlen(val)!=1 || (*val!='T' && *val!='F'))
      get_error("A logical value expected for '%s' in '%s'!",
                 collect_keyname(fits,head,row,key),fits->file.path);
   return(*val);
}

void get_keyword_textual(fits_file *fits,fits_head *head,char *key,
                         char *str,int maxlen)
{
   int row;

   row = get_fits_keyword(fits,head,key);
   read_keyword_textual(fits,head,row,str,maxlen);
}

char get_keyword_logical(fits_file *fits,fits_head *head,char *key)
{
   int row;
   char a;

   row = get_fits_keyword(fits,head,key);
   a = read_keyword_logical(fits,head,row);

   return(a);
}

int get_keyword_integer(fits_file *fits,fits_head *head,char *key)
{
   int row;
   int a;

   row = get_fits_keyword(fits,head,key);
   a = read_keyword_integer(fits,head,row);

   return(a);
}

double get_keyword_double(fits_file *fits,fits_head *head,char *key)
{
   int row;
   double a;

   row = get_fits_keyword(fits,head,key);
   a = read_keyword_double(fits,head,row);

   return(a);
}

float get_keyword_float(fits_file *fits,fits_head *head,char *key)
{
   int row;
   float a;

   row = get_fits_keyword(fits,head,key);
   a = read_keyword_float(fits,head,row);

   return(a);
}

void get_keyword_filename(fits_file *fits,char *fname)
{
   get_keyword_textual(fits,&fits->head,"FILENAME",fname,MAX_FITS_VALLEN);
}

void check_keyword_logical(fits_file *fits,fits_head *head,char *key,char a)
{
   if (get_keyword_logical(fits,head,key) != a)
      get_error("'%s = %c' expected in '%s'!",key,a,fits->file.path);
}

void check_keyword_integer(fits_file *fits,fits_head *head,char *key,int a)
{
   if (get_keyword_integer(fits,head,key) != a)
      get_error("'%s = %d' expected in '%s'!",key,a,fits->file.path);
}

void ensure_not_original_file(fits_file *fits)
{
   fitsval_type filename;

   get_keyword_filename(fits,filename);

   if (strcasecmp(fits->file.name,filename) == 0)
      get_error("Fits file '%s' seems to contain an original image!\n"
                "Make a copy using a different file name before editing.",
                fits->file.name);
}

void remove_fits_cards(fits_file *fits,fits_head *head,int first,int last)
{
   int seq,count;
   char *src,*dest;
   div_t q;

   if (last < first) get_error("remove_fits_cards: Bad card range!");
   if (first<2 || last>=head->ncard)
      get_error("remove_fits_cards: Card number out of range for '%s'!",
                 fits->file.path);

   seq = last + 1;
   count = head->ncard - seq + 1;

   src = head->dat + (seq-1)*CARD_BYTES;
   dest = head->dat + (first-1)*CARD_BYTES;
   memmove(dest,src,count*CARD_BYTES);

   head->ncard -= last-first+1;
   q = div(head->ncard,RECORD_CARDS);
   head->nrec = q.quot;
   if (q.rem > 0) head->nrec++;
}

void overwrite_existing_card(fits_head *head,int seq,char *text)
{
   char *card;

   card = head->dat + (seq-1)*CARD_BYTES;
   memset(card,0x20,CARD_BYTES);
   strncpy(card,text,strlen(text));
}

void insert_new_card(fits_file *fits,fits_head *head,char *text)
{
   if (head->ncard == HEAD_CARDS)
      get_error("FITS header full in '%s'!",fits->file.path);
   overwrite_existing_card(head,head->ncard++,text);
   overwrite_existing_card(head,head->ncard,"END");
   if (head->ncard > head->nrec*RECORD_CARDS) head->nrec++;
}

void write_fits_card(fits_file *fits,fits_head *head,char *text)
{
   int seq;
   fitskey_type key;

   if (extract_fits_keyname(text,key) != 0)
      get_error("write_fits_card: FITS keyword name too long: '%s'",text);
   seq = find_fits_keyword(head,key);
   if (seq > 0)
      overwrite_existing_card(head,seq,text);
   else
      insert_new_card(fits,head,text);
}

void insert_blank_card(fits_file *fits,fits_head *head)
{
   insert_new_card(fits,head," ");
}

void insert_comment_card(fits_file *fits,fits_head *head,char *com)
{
   char text[CARD_BYTES+1];

   sprintf(text,"COMMENT  %s",com);
   insert_new_card(fits,head,text);
}

void add_comment(char *text,char *com)
{
   if (com == NULL) return;
   if (*com == 0) return;
   if (strlen(com) > 47)
      get_error("add_comment: FITS comment too long!\n   %s",com);
   if (strlen(text) + strlen(com) > 77)
      get_error("add_comment: No space for FITS comment!\n"
                "   %s / %s",text,com);
   while (strlen(text) < 30) strcat(text," ");
   strcat(text," / ");
   strcat(text,com);
}

void write_keyword_logical(fits_file *fits,fits_head *head,char *key,
                           char a,char *com)
{
   char text[CARD_BYTES+1];

   sprintf(text,"%-8s= %19s%c",key," ",a);
   add_comment(text,com);
   write_fits_card(fits,head,text);
}

void write_keyword_integer(fits_file *fits,fits_head *head,char *key,
                           int a,char *com)
{
   char text[CARD_BYTES+1];

   sprintf(text,"%-8s= %20d",key,a);
   add_comment(text,com);
   write_fits_card(fits,head,text);
}

void write_keyword_double(fits_file *fits,fits_head *head,char *key,
                          double a,char *com)
{
   char text[CARD_BYTES+1];

   sprintf(text,"%-8s= %20.13e",key,a);
   add_comment(text,com);
   write_fits_card(fits,head,text);
}

void write_keyword_float(fits_file *fits,fits_head *head,char *key,
                         float a,char *com)
{
   char text[CARD_BYTES+1];

   sprintf(text,"%-8s= %20.7e",key,a);
   add_comment(text,com);
   write_fits_card(fits,head,text);
}

void write_keyword_textual(fits_file *fits,fits_head *head,char *key,
                           char *a,char *com)
{
   char text[CARD_BYTES+1];

   sprintf(text,"%-8s= '%s'",key,a);
   add_comment(text,com);
   write_fits_card(fits,head,text);
}

void copy_fits_keyword
  (fits_head *src,char *key,fits_file *fits,fits_head *head)
{
  int row;
  char card[81];

  row = find_fits_keyword(src,key);
  if (row < 1) return;
  extract_fits_card(src,row,card);
  write_fits_card(fits,head,card);
}

void replace_value_field(fits_head *head,int seq,char *newval)
{
  char *card;

  card = head->dat + (seq - 1) * CARD_BYTES;
  memset(card+9,0x20,22);
  memmove(card+30-strlen(newval),newval,strlen(newval));
}

void comment_out_fits_card(fits_head *head,int seq)
{
  char *card;

  card = head->dat + (seq - 1) * CARD_BYTES;
  memset(card,0x20,10);
  memmove(card,"COMMENT",7);
}

void ensure_comment_card(fits_file *fits,fits_head *head,char *text)
{
  int seq;

  seq = find_fits_comment(head,text);
  if (seq < 1) insert_comment_card(fits,head,text);
}

/* ------------------------- Regression Blocks ------------------------- */

int find_regression_block(fits_head *head,char *key)
{
   char text[CARD_BYTES+1];

   sprintf(text,"HISTORY   REGRESSION %s START",key);

   return(find_fits_card(head,text));
}

int get_regression_block(fits_file *fits,fits_head *head,char *key)
{
   int row;

   row = find_regression_block(head,key);
   if (row < 1)
      get_error("Regression block '%s' not found in '%s'!",
                 key,fits->file.path);
   return(row);
}

void remove_regression_block(fits_file *fits,fits_head *head,char *key)
{
   int first,last;
   char *card,text[CARD_BYTES+1];

   first = find_regression_block(head,key);
   if (first < 1) return;
   sprintf(text,"HISTORY   REGRESSION %s END",key);
   for (last=first+1,card = head->dat + first*CARD_BYTES;
        last<=head->ncard && !same_text(card,text);
        last++,card+=CARD_BYTES);
   if (!same_text(card,text))
      get_error("Incomplete regression block '%s' in '%s' (END not found)!",
                 key,fits->file.path);
   remove_fits_cards(fits,head,first,last);
}

void insert_integer_value(fits_file *fits,fits_head *head,char *var,int val)
{
   char text[CARD_BYTES+1];

   sprintf(text,"HISTORY   %s %d",var,val);
   insert_new_card(fits,head,text);
}

void insert_double_value(fits_file *fits,fits_head *head,char *var,double val)
{
   char text[CARD_BYTES+1];

   sprintf(text,"HISTORY   %s %22.15e",var,val);
   insert_new_card(fits,head,text);
}

void insert_integer_list(fits_file *fits,fits_head *head,char *var,
                         int *x,int n)
{
   int i;
   char text[CARD_BYTES+1],val[16];

   sprintf(text,"HISTORY   %s",var);
   for (i=0;i<n;i++)
   {
      sprintf(val," %d",x[i]);
      if (strlen(text)+strlen(val) > CARD_BYTES) get_error("FITS card full!");
      strcat(text,val);
   }
   insert_new_card(fits,head,text);
}

int find_fits_array(fits_head *head,char *var)
{
   char text[CARD_BYTES+1];

   sprintf(text,"HISTORY   ARRAY %s",var);

   return(find_fits_card(head,text));
}

void remove_fits_array(fits_file *fits,fits_head *head,char *var)
{
   int n,count,first,last;
   div_t q;
   char *card,type;

   first = find_fits_array(head,var);
   if (first < 1) return;
   card = head->dat + (first-1)*CARD_BYTES;
   sscanf(card+16+strlen(var)," %d %c",&n,&type);
   switch(type)
   {
      case 'D': count = 3;
                break;
      default:  get_error("Bad array type (%c) for '%s' in '%s'!",
                           type,var,fits->file.path);
                count = 0;
   }
   q = div(n,count);
   last = first + q.quot;
   if (q.rem > 0) last++;

   remove_fits_cards(fits,head,first,last);
}

void insert_double_array(fits_file *fits,fits_head *head,char *var,
                         double *x,int n)
{
   char text[CARD_BYTES+1],val[32];
   int i,k;

   sprintf(text,"HISTORY   ARRAY %s %d D",var,n);
   for (i=0,k=3;i<n;i++,k++)
   {
      if (k == 3)
      {
         insert_new_card(fits,head,text);
         strcpy(text,"HISTORY   ");
         k = 0;
      }
      sprintf(val,"%23.15e",x[i]);
      strcat(text,val);
   }
   insert_new_card(fits,head,text);
}

void write_double_array(fits_file *fits,fits_head *head,char *var,
                        double *x,int n)
{
   remove_fits_array(fits,head,var);
   insert_double_array(fits,head,var,x,n);
}

void write_regression_block(fits_file *fits,fits_head *head,char *key,
                            regre_block *reg)
{
   char text[CARD_BYTES+1];

   remove_regression_block(fits,head,key);
   sprintf(text,"HISTORY   REGRESSION %s START",key);
   insert_new_card(fits,head,text);

   insert_integer_value(fits,head,"NDIM",reg->ndim);
   insert_integer_value(fits,head,"YCOL",reg->ycol);
   insert_integer_list(fits,head,"XCOL",&reg->xcol[1],reg->ndim);
   insert_integer_list(fits,head,"DEG",&reg->deg[1],reg->ndim);
   insert_integer_value(fits,head,"SCOL",reg->scol);
   insert_integer_value(fits,head,"FCOL",reg->fcol);
   insert_integer_value(fits,head,"RCOL",reg->rcol);
   insert_integer_value(fits,head,"NSEL",reg->nsel);
   insert_double_value(fits,head,"RMS",reg->rms);
   insert_double_value(fits,head,"KAPPA",reg->kappa);

   insert_double_array(fits,head,"COEFF",&reg->a[1],reg->ma);

   sprintf(text,"HISTORY   REGRESSION %s END",key);
   insert_new_card(fits,head,text);
}

char *get_regression_var(fits_file *fits,fits_head *head,char *key,
                         int row,char *var)
{
   char *card,*limit,*ptr;

   card = head->dat + (row-1)*CARD_BYTES;

   if (!same_text(card,"HISTORY   "))
      get_error("Bad regression block '%s' in '%s': "
                "keyword 'HISTORY' expected!",key,fits->file.path);
   if (!same_text(card+10,var))
      get_error("Bad regression block '%s' in '%s': "
                "keyword '%s' expected!",key,fits->file.path,var);

   limit = card + CARD_BYTES;
   for (ptr=card+10+strlen(var);*ptr==' ' && ptr<limit;ptr++);
   if (ptr >= limit)
      get_error("Bad regression block '%s' in '%s': "
                "Value not found for keyword '%s'!",key,fits->file.path,var);
   return(ptr);
}

int extract_integer_value(fits_file *fits,fits_head *head,char *key,
                          int row,char *var)
{
   char *ptr,*end;
   int a;

   ptr = get_regression_var(fits,head,key,row,var);
   a = strtol(ptr,&end,10);
   if (end == ptr)
      get_error("Bad regression block '%s' in '%s': "
                "Integer value expected for keyword '%s'!",
                 key,fits->file.path,var);
   return(a);
}

void extract_integer_list(fits_file *fits,fits_head *head,char *key,
                          int row,char *var,int *a,int n)
{
   char *ptr,*end;
   int i;

   ptr = get_regression_var(fits,head,key,row,var);

   for (i=0;i<n;i++)
   {
      a[i] = strtol(ptr,&end,10);
      if (end == ptr)
         get_error("Bad regression block '%s' in '%s': "
                   "Integer value expected for keyword '%s'!",
                    key,fits->file.path,var);
      ptr = end;
   }
}

double extract_double_value(fits_file *fits,fits_head *head,char *key,
                          int row,char *var)
{
   char *ptr,*end;
   double a;

   ptr = get_regression_var(fits,head,key,row,var);
   a = strtod(ptr,&end);
   if (end == ptr)
      get_error("Bad regression block '%s' in '%s': "
                "Numerical value expected for keyword '%s'!",
                 key,fits->file.path,var);
   return(a);
}

void extract_double_array(fits_file *fits,fits_head *head,char *key,
                          int row,char *var,double *a,int n)
{
   int dim;
   char type;
   int i;
   char *ptr,*eptr;

   ptr = get_regression_var(fits,head,key,row,"ARRAY");
   if (!same_text(ptr,var))
      get_error("Bad regression block '%s' in '%s': "
                "Double array '%s' expected!",
                 key,fits->file.path,var);
   sscanf(ptr+strlen(var),"%d %c",&dim,&type);
   if (type != 'D')
      get_error("Bad regression block '%s' in '%s': "
                "Double array type expected for keyword '%s'!",
                 key,fits->file.path,var);
   if (dim != n)
      get_error("Bad regression block '%s' in '%s': "
                "%d elements expected for keyword '%s'!",
                 key,fits->file.path,n,var);
   for (i=0;i<n;i++)
   {
      if ((i/3)*3 == i)
      {
         row++;
         ptr = head->dat + (row-1)*CARD_BYTES + 10;
      }
      a[i] = strtod(ptr,&eptr);
      if (eptr == ptr)
      get_error("Bad regression block '%s' in '%s': "
                "A double value expected for ellement element %d "
                "of array '%s'!",
                 key,fits->file.path,i+1,var);
      ptr = eptr;
   }
}

int get_array_dimension(fits_file *fits,fits_head *head,char *var)
{
   int row,dim;
   char *card,type;

   row = find_fits_array(head,var);
   if (row < 1)
      get_error("Array '%s' not found in '%s'!",var,fits->file.path);
   card = head->dat + (row-1)*CARD_BYTES;

   sscanf(card+16+strlen(var),"%d %c",&dim,&type);
   return(dim);
}

void get_double_array(fits_file *fits,fits_head *head,char *var,
                      double *a,int n)
{
   int row,dim,i;
   char *card,type,*ptr,*eptr;

   row = find_fits_array(head,var);
   if (row < 1)
      get_error("Array '%s' not found in '%s'!",var,fits->file.path);
   card = head->dat + (row-1)*CARD_BYTES;

   sscanf(card+16+strlen(var),"%d %c",&dim,&type);
   if (type != 'D')
      get_error("Double type expected for array '%s' in '%s'!",
                 var,fits->file.path);
   if (dim != n)
      get_error("%d elements expected for array '%s' in '%s' (%d found)!",
                 n,var,fits->file.path,dim);

   for (i=0;i<n;i++)
   {
      if ((i/3)*3 == i)
      {
         row++;
         ptr = head->dat + (row-1)*CARD_BYTES + 10;
      }
      a[i] = strtod(ptr,&eptr);
      if (eptr == ptr)
      get_error("A double value expected for element %d "
                "of array '%s' in '%s'!",
                 i+1,var,fits->file.path);
      ptr = eptr;
   }
}

void read_regression_block(fits_file *fits,fits_head *head,char *key,
                          regre_block *reg)
{
   int row;

   row = find_regression_block(head,key);
   if (row < 1)
      get_error("Regression block '%s' not found in '%s'!",
                 key,fits->file.path);
   row++;
   reg->ndim = extract_integer_value(fits,head,key,row++,"NDIM");
   if (reg->ndim > MAX_REGRESSION_NDIM)
      get_error("Bad regression block '%s' in '%s': NDIM too large!",
                 key,fits->file.path);
   reg->ycol = extract_integer_value(fits,head,key,row++,"YCOL");
   extract_integer_list(fits,head,key,row++,"XCOL",&reg->xcol[1],reg->ndim);
   extract_integer_list(fits,head,key,row++,"DEG",&reg->deg[1],reg->ndim);
   list_integers(reg->deg,reg->ndim,reg->arg,MAX_REGRESSION_ARG);
   reg->scol = extract_integer_value(fits,head,key,row++,"SCOL");
   reg->fcol = extract_integer_value(fits,head,key,row++,"FCOL");
   reg->rcol = extract_integer_value(fits,head,key,row++,"RCOL");
   reg->nsel = extract_integer_value(fits,head,key,row++,"NSEL");
   reg->rms = extract_double_value(fits,head,key,row++,"RMS");
   reg->kappa = extract_double_value(fits,head,key,row++,"KAPPA");
   reg->ma = number_of_coefficients(reg->deg,reg->ndim);
   extract_double_array(fits,head,key,row,"COEFF",&reg->a[1],reg->ma);
}

/* ---------------------------- Binary data ---------------------------- */

void set_totbinrec(fits_file *fits,int size)
{
   div_t q;

   q = div(size,RECORD_BYTES);
   fits->totbinrec = q.quot;
   if (q.rem != 0) fits->totbinrec++;
}

void set_img_totbinrec(fits_file *fits)
{
   set_totbinrec(fits,fits->totpix*fits->pixsize);
}

void set_tbl_totbinrec(fits_file *fits)
{
   set_totbinrec(fits,fits->nrows*fits->rowlen);
}

void allocate_whole_binary(fits_file *fits)
{
   fits->bin = (byte *)get_space(fits->totbinrec*RECORD_BYTES);
   memset(fits->bin,0,fits->totbinrec*RECORD_BYTES);
   fits->pix = (float *)fits->bin;
}

void free_binary_data(fits_file *fits)
{
   get_free((void *)fits->bin);
   fits->bin = NULL;
   fits->pix = NULL;
}

double read_pixel_double(byte *bptr,int bitpix)
{
   double val;

   switch(bitpix)
   {
      case   8: val = (double)*bptr;
                break;
      case  16: val = (double)*(short *)bptr;
                break;
      case  32: val = (double)*(int32_t *)bptr;
                break;
      case -32: val = (double)*(float *)bptr;
                break;
      case -64: val = *(double *)bptr;
                break;
      default:  get_error("read_pixel_double: Invalid BITPIX = %d!",bitpix);
   }

   return(val);
}

void write_pixel_double(byte *bptr,int bitpix,double val)
{
   switch(bitpix)
   {
      case   8: *bptr = (char)val;
                break;
      case  16: *(short *)bptr = (short)val;
                break;
      case  32: *(int32_t *)bptr = (int32_t)val;
                break;
      case -32: *(float *)bptr = (float)val;
                break;
      case -64: *(double *)bptr = val;
                break;
   }
}

float collect_pixel_value(fits_file *fits,int x,int y)
{
   if (x>=1 && x<=fits->npix[1] && y>=1 && y<=fits->npix[2])
      return(fits->pix[(y-1)*fits->npix[1]+x-1]);
   else
      return(0.0);
}

void collect_pixel_row(fits_file *fits,int row,double *val)
{
   float *pix;
   int col;

   if (row<1 || row>fits->npix[2])
      get_error("collect_pixel_row: Bad row number %d!",row);

   pix = fits->pix + (row - 1)*fits->npix[1];
   for (col=1;col<=fits->npix[1];col++,val++,pix++) *val = (double)(*pix);
}

void deposit_pixel_value(fits_file *fits,int x,int y,float val)
{
   if (x>=1 && x<=fits->npix[1] && y>=1 && y<=fits->npix[2])
      fits->pix[(y-1)*fits->npix[1]+x-1] = val;
}

void deposit_pixel_row(fits_file *fits,int row,double *val)
{
   float *pix;
   int col;

   if (row<1 || row>fits->npix[2])
      get_error("deposit_pixel_row: Bad row number %d!",row);

   pix = fits->pix + (row - 1)*fits->npix[1];
   for (col=1;col<=fits->npix[1];col++,val++,pix++) *pix = (float)(*val);
}

void invert_image_pixels(fits_file *fits)
{
   unsigned char *ptr,*src,*dest,t;
   int pix;

   for (pix=1,ptr=fits->bin;pix<=fits->totpix;pix++,ptr+=fits->pixsize)
   {
      for (src=ptr,dest=ptr+fits->pixsize-1;src<dest;src++,dest--)
         {t=*src; *src=*dest; *dest=t;}
   }
}

void copy_pixel_rows
   (fits_file *xfits,int xrow,fits_file *yfits,int yrow,int nrows)
{
   unsigned char *xptr,*yptr;
   int nbytes;

   if (nrows < 1) return;

   check_image_2d(xfits);
   check_image_2d(yfits);

   if (xfits->npix[1] != yfits->npix[1])
      get_error("copy_pixel_rows: NAXIS1 values must be the same!");
   if (xfits->bitpix != yfits->bitpix)
      get_error("copy_pixel_rows: BITPIX values must be the same!");
   if (xfits->bscale != yfits->bscale)
      get_error("copy_pixel_rows: BSCALE values must be the same!");
   if (xfits->bzero != yfits->bzero)
      get_error("copy_pixel_rows: BZERO values must be the same!");

   if (xrow < 1 || xrow > xfits->npix[2])
      get_error("Input row number out of range!");
   if (yrow < 1 || yrow > yfits->npix[2])
      get_error("Output row number out of range!");

   if (xrow + nrows - 1 > xfits->npix[2])
      get_error("No more input rows!");
   if (yrow + nrows - 1 > yfits->npix[2])
      get_error("No more output rows!");

   xptr = xfits->bin + (xrow - 1) * xfits->rowsize;
   yptr = yfits->bin + (yrow - 1) * yfits->rowsize;
   nbytes = nrows * xfits->rowsize;
   memmove(yptr,xptr,nbytes);
}

void copy_pixel_cols
   (fits_file *xfits,int xcol,fits_file *yfits,int ycol,int ncols)
{
   unsigned char *xptr,*yptr;
   int nbytes;
   int row;

   if (ncols < 1) return;

   check_image_2d(xfits);
   check_image_2d(yfits);

   if (xfits->npix[2] != yfits->npix[2])
      get_error("copy_pixel_cols: NAXIS2 values must be the same!");
   if (xfits->bitpix != yfits->bitpix)
      get_error("copy_pixel_cols: BITPIX values must be the same!");
   if (xfits->bscale != yfits->bscale)
      get_error("copy_pixel_cols: BSCALE values must be the same!");
   if (xfits->bzero != yfits->bzero)
      get_error("copy_pixel_cols: BZERO values must be the same!");

   if (xcol < 1 || xcol > xfits->npix[1])
      get_error("Input column number out of range!");
   if (ycol < 1 || ycol > yfits->npix[1])
      get_error("Output column number out of range!");

   if (xcol + ncols - 1 > xfits->npix[1])
      get_error("No more input columns!");
   if (ycol + ncols - 1 > yfits->npix[1])
      get_error("No more output columns!");

   xptr = xfits->bin + (xcol - 1) * xfits->pixsize;
   yptr = yfits->bin + (ycol - 1) * yfits->pixsize;
   nbytes = ncols * xfits->pixsize;
   for (row=1;row<=xfits->npix[2];row++)
   {
      memmove(yptr,xptr,nbytes);
      xptr += xfits->rowsize;
      yptr += yfits->rowsize;
   }
}

void get_raw_fits_table_data(fits_file *fits)
{
   unsigned char *bintbl;
   int col,row;
   unsigned char *coloffset,*colstart,*src,*dest;

   bintbl = fits->bin;
   allocate_whole_binary(fits);

   for (col=1,coloffset=fits->bin,src=bintbl;col<=fits->ncols;
                                    coloffset+=fits->colsize[col++])
   {
      for (row=1,colstart=coloffset;row<=fits->nrows;
                                    row++,colstart+=fits->rowlen)
      {
         if (fits->coltype[col] == 'C')
         {
            memmove(colstart,src,fits->colsize[col]);
            src += fits->colsize[col];
         }
         else
         {
            for (dest=colstart+fits->colsize[col]-1;dest>=colstart;
                 dest--,src++) *dest = *src;
         }
      }
   }

   get_free(bintbl);
}

void get_normal_table_data(fits_file *fits)
{
   unsigned char *rawtbl,*src,*dest,*coloffset,*colstart;
   int row,col;

   rawtbl = fits->bin;
   allocate_whole_binary(fits);

   for (col=1,coloffset=rawtbl,dest=fits->bin;col<=fits->ncols;
                                    coloffset+=fits->colsize[col++])
   {
      for (row=1,colstart=coloffset;row<=fits->nrows;
                                    row++,colstart+=fits->rowlen)
      {
         if (fits->coltype[col] == 'C')
         {
            memmove(dest,colstart,fits->colsize[col]);
            dest += fits->colsize[col];
         }
         else
         {
            for (src=colstart+fits->colsize[col]-1;src>=colstart;src--,dest++)
               *dest=*src;
         }
      }
   }

   get_free(rawtbl);
}

/* ---------------------------- FITS files ----------------------------- */

void vset_fits_name(fits_file *fits,char *name,va_list ap)
{
   file_type f;
   int fit_found;

   vset_file_name(&f,name,ap);
   fit_found = 0;
   if (strlen(f.path) > 4)
      fit_found = (strcasecmp(f.path+strlen(f.path)-4,".FIT") == 0);
   if (fit_found)
      set_file_name(&fits->file,f.path);
   else
      set_file_name(&fits->file,"%s.fit",f.path);
}

void set_fits_name(fits_file *fits,char *name, ...)
{
   va_list ap;

   va_start(ap,name);
   vset_fits_name(fits,name,ap);
   va_end(ap);
}

void open_fits_file(fits_file *fits)
{
   get_fopen(&fits->file,"r");
}

void create_fits_file(fits_file *fits)
{
   get_fopen(&fits->file,"w");
}

void seek_fits_file(fits_file *fits,int offset,int whence)
{
   get_fseek(&fits->file,offset,whence);
}

void read_fits_file(fits_file *fits,int size,int count,byte *buf)
{
  get_fread(buf,size,count,&fits->file);
}

void write_fits_file(fits_file *fits,int size,int count,byte *buf)
{
  get_fwrite(buf,size,count,&fits->file);
}

void close_fits_file(fits_file *fits)
{
   get_fclose(&fits->file);
}

/* --------- Accessing a FITS file that has already been open ---------- */

void load_fits_header(fits_file *fits,fits_head *head,char *signature)
{
   char *record,*card;
   char *headstop,*recstop;
   char *endcard;

   read_fits_file(fits,1,CARD_BYTES,(byte *)head->dat);
   if (!same_text(head->dat,signature))
      get_error("'%s' expected in '%s'!",signature,fits->file.path);
   seek_fits_file(fits,-CARD_BYTES,SEEK_CUR);

   memset(head->dat,0x20,HEAD_BYTES);
   endcard = NULL;
   headstop = head->dat + HEAD_BYTES;
   head->nrec = head->ncard = 0;
   for (record=head->dat;record<headstop && endcard==NULL;
        record+=RECORD_BYTES)
   {
      read_fits_file(fits,1,RECORD_BYTES,(byte *)record);
      if (!extended_ascii_string_ok((unsigned char *)record,RECORD_BYTES))
        get_error("An unexpected character found in '%s'!",fits->file.path);
      head->nrec++;
      recstop = record + RECORD_BYTES;
      for (card=record;card<recstop && endcard==NULL;card+=CARD_BYTES)
      {
         head->ncard++;
         if (same_text(card,"END       ")) endcard=card;
      }
   }
}

void load_primary_header(fits_file *fits)
{
   load_fits_header(fits,&fits->head,primary);
}

void load_extension_header(fits_file *fits)
{
   load_fits_header(fits,&fits->extend,bintable);
}

void save_primary_header(fits_file *fits)
{
   write_fits_file
      (fits,RECORD_BYTES,fits->head.nrec,(byte *)fits->head.dat);
}

void save_extension_header(fits_file *fits)
{
   write_fits_file
      (fits,RECORD_BYTES,fits->extend.nrec,(byte *)fits->extend.dat);
}

void load_whole_binary(fits_file *fits)
{
   read_fits_file(fits,RECORD_BYTES,fits->totbinrec,fits->bin);
}

void load_image_pixels(fits_file *fits)
{
   load_whole_binary(fits);
   invert_image_pixels(fits);
}

void load_table_data(fits_file *fits)
{
   load_whole_binary(fits);
   get_normal_table_data(fits);
}

void save_whole_binary(fits_file *fits)
{
   write_fits_file(fits,RECORD_BYTES,fits->totbinrec,fits->bin);
}

void save_image_pixels(fits_file *fits)
{
   invert_image_pixels(fits);
   save_whole_binary(fits);
}

void save_table_data(fits_file *fits)
{
   get_raw_fits_table_data(fits);
   save_whole_binary(fits);
}

void load_image_header(fits_file *fits)
{
   int axis,seq;
   fitskey_type key;

   load_primary_header(fits);

   fits->bitpix = get_keyword_integer(fits,&fits->head,"BITPIX");
   if (fits->bitpix!=8 && fits->bitpix!=16 && fits->bitpix!=32 &&
       fits->bitpix!=-32 && fits->bitpix!=-64)
      get_error("Invalid keyword BITPIX = %d in '%s'!",
                 fits->bitpix,fits->file.path);
   fits->pixsize = abs(fits->bitpix)/8;

   fits->naxis = get_keyword_integer(fits,&fits->head,"NAXIS");
   if (fits->naxis < 1 || fits->naxis > MAXNAXIS)
      get_error("Invalid keyword NAXIS = %d in '%s'!",
                 fits->naxis,fits->file.path);

   for (axis=1;axis<=MAXNAXIS;axis++) fits->npix[axis] = 1;

   for (axis=1,fits->totpix=1;axis<=fits->naxis;axis++)
   {
      sprintf(key,"NAXIS%d",axis);
      fits->npix[axis] = get_keyword_integer(fits,&fits->head,key);
      if (fits->npix[axis] < 1)
         get_error("Invalid keyword NAXIS%d = %d in '%s'!",
                    axis,fits->npix[axis],fits->file.path);
      fits->totpix *= fits->npix[axis];

      sprintf(key,"CRPIX%d",axis);
      seq = find_fits_keyword(&fits->head,key);
      if (seq > 0)
         fits->crpix[axis] = read_keyword_double(fits,&fits->head,seq);
      else
         fits->crpix[axis] = 1.0;

      sprintf(key,"CRVAL%d",axis);
      seq = find_fits_keyword(&fits->head,key);
      if (seq > 0)
         fits->crval[axis] = read_keyword_double(fits,&fits->head,seq);
      else
         fits->crval[axis] = 1.0;

      sprintf(key,"CDELT%d",axis);
      seq = find_fits_keyword(&fits->head,key);
      if (seq > 0)
         fits->cdelt[axis] = read_keyword_double(fits,&fits->head,seq);
      else
         fits->cdelt[axis] = 1.0;
   }

   fits->rowsize = fits->npix[1] * fits->pixsize;

   seq = find_fits_keyword(&fits->head,"BSCALE");
   if (seq > 0)
      fits->bscale = read_keyword_double(fits,&fits->head,seq);
   else
      fits->bscale = 1.0;

   seq = find_fits_keyword(&fits->head,"BZERO");
   if (seq > 0)
      fits->bzero = read_keyword_double(fits,&fits->head,seq);
   else
      fits->bzero = 0.0;

   set_img_totbinrec(fits);
   fits->bin = NULL;
   fits->pix = NULL;

   fits->binstart = fits->head.nrec * RECORD_BYTES;
}

void default_display_format(fits_file *fits,int col)
{
   switch(fits->coltype[col])
   {
      case 'I': strcpy(fits->coldisp[col],"I8");
                break;
      case 'R': strcpy(fits->coldisp[col],"E14.7");
                break;
      case 'D': strcpy(fits->coldisp[col],"E22.15");
                break;
      case 'C': sprintf(fits->coldisp[col],"A%d",fits->colsize[col]);
                break;
   }
}

void load_table_header(fits_file *fits)
{
   int col,total,seq;
   fitskey_type key;
   char fmt[9];

   load_primary_header(fits);

   check_keyword_integer(fits,&fits->head,"NAXIS",0);
   check_keyword_logical(fits,&fits->head,"EXTEND",'T');

   load_extension_header(fits);

   check_keyword_integer(fits,&fits->extend,"BITPIX",8);
   check_keyword_integer(fits,&fits->extend,"NAXIS",2);
   check_keyword_integer(fits,&fits->extend,"PCOUNT",0);
   check_keyword_integer(fits,&fits->extend,"GCOUNT",1);

   fits->rowlen = get_keyword_integer(fits,&fits->extend,"NAXIS1");
   fits->nrows = get_keyword_integer(fits,&fits->extend,"NAXIS2");
   fits->ncols = get_keyword_integer(fits,&fits->extend,"TFIELDS");

   if (fits->ncols > MAXFITSCOLS) get_error("Too many table columns!");

   set_tbl_totbinrec(fits);

   for (col=1,total=0;col<=fits->ncols;col++)
   {
      sprintf(key,"TFORM%d",col);
      get_keyword_textual(fits,&fits->extend,key,fmt,8);
      read_column_format(fmt,&fits->coltype[col],&fits->colsize[col]);
      total += fits->colsize[col];
      sprintf(key,"TDISP%d",col);
      seq = find_fits_keyword(&fits->extend,key);
      if (seq > 0)
         read_keyword_textual(fits,&fits->extend,seq,
                              fits->coldisp[col],MAXCOLDISP);
      else
         default_display_format(fits,col);
      sprintf(key,"TTYPE%d",col);
      get_keyword_textual(fits,&fits->extend,key,
                          fits->colname[col],MAXCOLNAME);
      trim_right(fits->colname[col]);
   }

   if (total != fits->rowlen)
      get_error("Bad total byte count per table row in '%s'!",
                 fits->file.path);

   fits->colstart[1]=1;
   fits->coloffset[1]=0;
   for (col=2;col<=fits->ncols;col++)
   {
      fits->colstart[col] = fits->colstart[col-1] + fits->colsize[col-1];
      fits->coloffset[col] = fits->coloffset[col-1] +
                                        fits->colsize[col-1]*fits->nrows;
   }

   fits->bin = NULL;
   fits->pix = NULL;

   fits->binstart = (fits->head.nrec+fits->extend.nrec) * RECORD_BYTES;
}

void save_image_header(fits_file *fits)
{
   save_primary_header(fits);
}

void save_table_header(fits_file *fits)
{
   save_primary_header(fits);
   save_extension_header(fits);
}

void load_image(fits_file *fits)
{
   load_image_header(fits);
   allocate_whole_binary(fits);
   load_image_pixels(fits);
}

void load_table(fits_file *fits)
{
   load_table_header(fits);
   allocate_whole_binary(fits);
   load_table_data(fits);
}

void load_raw_image(fits_file *fits)
{
   load_image_header(fits);
   allocate_whole_binary(fits);
   load_whole_binary(fits);
}

void load_raw_table(fits_file *fits)
{
   load_table_header(fits);
   allocate_whole_binary(fits);
   load_whole_binary(fits);
}

void write_image_date(fits_file *fits)
{
   write_keyword_textual(fits,&fits->head,"DATE",today(),
                         "Date of this FITS file creation");
}

void save_image(fits_file *fits)
{
   write_image_date(fits);
   save_image_header(fits);
   save_image_pixels(fits);
   free_binary_data(fits);
}

void save_raw_image(fits_file *fits)
{
   write_image_date(fits);
   save_image_header(fits);
   save_whole_binary(fits);
   free_binary_data(fits);
}

void save_table(fits_file *fits)
{
   write_image_date(fits);
   save_table_header(fits);
   save_table_data(fits);
   free_binary_data(fits);
}

void save_raw_table(fits_file *fits)
{
   write_image_date(fits);
   save_table_header(fits);
   save_whole_binary(fits);
   free_binary_data(fits);
}

/* ----------- Accessing a FITS file that has not been open ------------ */

int vfits_table_ok(char *name,va_list ap)
{
   fits_file fits;
   char card[CARD_BYTES];

   vset_fits_name(&fits,name,ap);
   open_fits_file(&fits);
   seek_fits_file(&fits,RECORD_BYTES,SEEK_SET);
   read_fits_file(&fits,1,CARD_BYTES,(byte *)card);
   close_fits_file(&fits);

   return(same_text(card,"XTENSION= 'BINTABLE'"));
}

int fits_table_ok(char *name, ...)
{
   int flag;
   va_list ap;

   va_start(ap,name);
   flag = vfits_table_ok(name,ap);
   va_end(ap);

   return(flag);
}

void vread_image_header(fits_file *fits,char *name,va_list ap)
{
   vset_fits_name(fits,name,ap);
   open_fits_file(fits);
   load_image_header(fits);
   close_fits_file(fits);
}

void read_image_header(fits_file *fits,char *name, ...)
{
   va_list ap;

   va_start(ap,name);
   vread_image_header(fits,name,ap);
   va_end(ap);
}

void vread_table_header(fits_file *fits,char *name,va_list ap)
{
   vset_fits_name(fits,name,ap);
   open_fits_file(fits);
   load_table_header(fits);
   close_fits_file(fits);
}

void read_table_header(fits_file *fits,char *name, ...)
{
   va_list ap;

   va_start(ap,name);
   vread_table_header(fits,name,ap);
   va_end(ap);
}

fits_head *vread_fits_header(fits_file *fits,char *name,va_list ap)
{
   int tab_ok;
   fits_head *head;
   va_list aq;

   va_copy(aq,ap);
   tab_ok = vfits_table_ok(name,aq);
   va_end(aq);
   vset_fits_name(fits,name,ap);
   open_fits_file(fits);
   if (tab_ok)
   {
      load_table_header(fits);
      head = &fits->extend;
   }
   else
   {
      load_image_header(fits);
      head = &fits->head;
   }
   close_fits_file(fits);

   return(head);
}

fits_head *read_fits_header(fits_file *fits,char *name, ...)
{
   va_list ap;
   fits_head *head;

   va_start(ap,name);
   head = vread_fits_header(fits,name,ap);
   va_end(ap);

   return(head);
}

void read_binary_data(fits_file *fits)
{
   open_fits_file(fits);
   seek_fits_file(fits,fits->binstart,SEEK_SET);
   load_whole_binary(fits);
   close_fits_file(fits);
}

void read_image_pixels(fits_file *fits)
{
   read_binary_data(fits);
   invert_image_pixels(fits);
}

void vread_image(fits_file *fits,char *name,va_list ap)
{
   vset_fits_name(fits,name,ap);
   open_fits_file(fits);
   load_image(fits);
   close_fits_file(fits);
}

void vread_table(fits_file *fits,char *name,va_list ap)
{
   vset_fits_name(fits,name,ap);
   open_fits_file(fits);
   load_table(fits);
   close_fits_file(fits);
}

void read_image(fits_file *fits,char *name,...)
{
   va_list ap;

   va_start(ap,name);
   vread_image(fits,name,ap);
   va_end(ap);
}

void read_table(fits_file *fits,char *name,...)
{
   va_list ap;

   va_start(ap,name);
   vread_table(fits,name,ap);
   va_end(ap);
}

void vread_raw_image(fits_file *fits,char *name,va_list ap)
{
   vset_fits_name(fits,name,ap);
   open_fits_file(fits);
   load_raw_image(fits);
   close_fits_file(fits);
}

void read_raw_image(fits_file *fits,char *name,...)
{
   va_list ap;

   va_start(ap,name);
   vread_raw_image(fits,name,ap);
   va_end(ap);
}

void vread_raw_table(fits_file *fits,char *name,va_list ap)
{
   vset_fits_name(fits,name,ap);
   open_fits_file(fits);
   load_raw_table(fits);
   close_fits_file(fits);
}

void read_raw_table(fits_file *fits,char *name,...)
{
   va_list ap;

   va_start(ap,name);
   vread_raw_table(fits,name,ap);
   va_end(ap);
}

void vwrite_image(fits_file *fits,char *name,va_list ap)
{
   vset_fits_name(fits,name,ap);
   create_fits_file(fits);
   save_image(fits);
   close_fits_file(fits);
}

void vwrite_table(fits_file *fits,char *name,va_list ap)
{
   vset_fits_name(fits,name,ap);
   create_fits_file(fits);
   save_table(fits);
   close_fits_file(fits);
}

void write_image(fits_file *fits,char *name, ...)
{
   va_list ap;

   va_start(ap,name);
   vwrite_image(fits,name,ap);
   va_end(ap);
}

void write_table(fits_file *fits,char *name, ...)
{
   va_list ap;

   va_start(ap,name);
   vwrite_table(fits,name,ap);
   va_end(ap);
}

void vwrite_raw_image(fits_file *fits,char *name,va_list ap)
{
   vset_fits_name(fits,name,ap);
   create_fits_file(fits);
   save_raw_image(fits);
   close_fits_file(fits);
}

void write_raw_image(fits_file *fits,char *name, ...)
{
   va_list ap;

   va_start(ap,name);
   vwrite_raw_image(fits,name,ap);
   va_end(ap);
}

void vwrite_raw_table(fits_file *fits,char *name,va_list ap)
{
   vset_fits_name(fits,name,ap);
   create_fits_file(fits);
   save_raw_table(fits);
   close_fits_file(fits);
}

void write_raw_table(fits_file *fits,char *name, ...)
{
   va_list ap;

   va_start(ap,name);
   vwrite_raw_table(fits,name,ap);
   va_end(ap);
}

void check_image_1d(fits_file *fits)
{
   if (fits->naxis != 1)
      get_error("One-dimensional image expected (%s)!",fits->file.path);
}

void check_image_2d(fits_file *fits)
{
   if (fits->naxis != 2)
      get_error("Two-dimensional image expected (%s)!",fits->file.path);
}

void check_same_image_size(fits_file *x,fits_file *y)
{
   int axis;

   if (x->naxis != y->naxis)
      get_error("Same number of axes expected (%s %s)!",
                 x->file.path,y->file.path);

   for (axis=1;axis<=x->naxis;axis++)
      if (x->npix[axis] != y->npix[axis])
         get_error("Same number of pixels along axis %d expected(%s %s)!",
                 axis,x->file.path,y->file.path);
}

void check_same_image_step(fits_file *x,fits_file *y)
{
   int axis;

   if (x->naxis != y->naxis)
      get_error("Same number of axes expected (%s %s)!",
                 x->file.path,y->file.path);

   for (axis=1;axis<=x->naxis;axis++)
      if (x->cdelt[axis] != y->cdelt[axis])
         get_error("Same step along axis %d expected(%s %s)!",
                 axis,x->file.path,y->file.path);
}

void check_fft_image(fits_file *x)
{
   int axis;

   for (axis=1;axis<=x->naxis;axis++)
   {
      if (!power_of_two(x->npix[axis]))
         get_error("check_fft_image: Invalid length (N=%d) "
                 "for axis No. %d in '%s'!",
                 x->npix[axis],axis,x->file.path);
   }
}

void check_real_image(fits_file *fits)
{
   if (fits->bitpix != -32)
      get_error("A 32-bit real image expected (%s)!",fits->file.path);
}

int pixel_inside(fits_file *fits,int x,int y)
{
   if (x < 1 || x > fits->npix[1]) return(0);
   if (y < 1 || y > fits->npix[2]) return(0);
   return(1);
}

double world_coordinate(fits_file *x,int axis,double pix)
{
   return(x->crval[axis] + (pix-x->crpix[axis])*x->cdelt[axis]);
}

double pixel_coordinate(fits_file *x,int axis,double val)
{
   return(x->crpix[axis] + (val-x->crval[axis])/x->cdelt[axis]);
}

double get_axis_start(fits_file *x,int axis)
{
   return(world_coordinate(x,axis,1.0));
}

void add_constant(char *xname,double c,char *yname)
{
   fits_file fits;
   int k;

   read_image(&fits,xname);
   check_real_image(&fits);
   for (k=0;k<fits.totpix;k++) fits.pix[k] += (float)c;
   write_image(&fits,yname);
}

void write_image_descriptor_integer(char *img,char *key,int a,char *com)
{
   fits_file fits;

   read_raw_image(&fits,img);
   write_keyword_integer(&fits,&fits.head,key,a,com);
   write_raw_image(&fits,img);
}

void write_image_descriptor_double(char *img,char *key,double a,char *com)
{
   fits_file fits;

   read_raw_image(&fits,img);
   write_keyword_double(&fits,&fits.head,key,a,com);
   write_raw_image(&fits,img);
}

void write_image_descriptor_float(char *img,char *key,float a,char *com)
{
   fits_file fits;

   read_raw_image(&fits,img);
   write_keyword_float(&fits,&fits.head,key,a,com);
   write_raw_image(&fits,img);
}

void write_image_descriptor_logical(char *img,char *key,char a,char *com)
{
   fits_file fits;

   read_raw_image(&fits,img);
   write_keyword_logical(&fits,&fits.head,key,a,com);
   write_raw_image(&fits,img);
}

void write_image_descriptor_textual(char *img,char *key,char *a,char *com)
{
   fits_file fits;

   read_raw_image(&fits,img);
   write_keyword_textual(&fits,&fits.head,key,a,com);
   write_raw_image(&fits,img);
}

void update_fits_minmax(fits_file *fits)
{
   int k;
   float min,max;

   check_real_image(fits);
   for (k=0,min=max=*fits->pix;k<fits->totpix;k++)
   {
      if (fits->pix[k] < min)
         min = fits->pix[k];
      else
         if (fits->pix[k] > max) max = fits->pix[k];
   }
   write_keyword_float(fits,&fits->head,"DATAMIN",min,"Data minimum value");
   write_keyword_float(fits,&fits->head,"DATAMAX",max,"Data maximum value");
}

void update_image_minmax(char *name)
{
   fits_file fits;

   read_image(&fits,name);
   update_fits_minmax(&fits);
   write_image(&fits,name);
}

void locate_fits_max(fits_file *fits,float *maxpix,int *maxcol,int *maxrow)
{
   int row,col;
   float *pixptr;

   *maxpix = *fits->pix;
   *maxrow = *maxcol = 1;
   for (row=1,pixptr=fits->pix;row<=fits->npix[2];row++)
   {
      for (col=1;col<=fits->npix[1];col++,pixptr++)
      {
         if (*pixptr > *maxpix)
         {
            *maxpix = *pixptr;
            *maxrow = row;
            *maxcol = col;
         }
      }
   }
}

void normalize_pixel_sum(fits_file *fits)
{
   int k;
   float tot;

   check_real_image(fits);
   for (k=0,tot=0.0;k<fits->totpix;k++) tot += fits->pix[k];
   if (tot == 0.0)
      get_error("Cannot normalize a zero image (%s)!",fits->file.path);
   for (k=0;k<fits->totpix;k++) fits->pix[k] /= tot;
   update_fits_minmax(fits);
}

void normalize_image(char *xname,char *yname)
{
   fits_file fits;

   read_image(&fits,xname);
   normalize_pixel_sum(&fits);
   write_image(&fits,yname);
}

void locate_pixel(fits_file *fits,int k,int *pix)
{
   int block,rest;
   div_t q;
   int axis;

   block = fits->totpix;
   rest = k-1;

   for (axis=fits->naxis;axis>=1;axis--)
   {
      block /= fits->npix[axis];
      q = div(rest,block);
      pix[axis] = q.quot+1;
      rest = q.rem;
   }
}

void xshift_image_fits(fits_file *fits,double q)
{
   int w,h;
   int row,col,iq,acol,bcol;
   float a,b;
   double uq,vq;

   if (q == 0) return;

   w = fits->npix[1];
   h = fits->npix[2];

   iq = (int)floor(fabs(q));
   uq = fabs(q) - iq;
   vq = 1.0 - uq;

   if (q > 0)
   {
      for (row=1;row<=h;row++)
      {
         for (col=w,bcol=w-iq;col>=1;col--,bcol--)
         {
            b = collect_pixel_value(fits,bcol,row);
            a = collect_pixel_value(fits,bcol-1,row);
            deposit_pixel_value(fits,col,row,a*uq+b*vq);
         }
      }
   }
   else
   {
      for (row=1;row<=h;row++)
      {
         for (col=1,acol=iq+1;col<=w;col++,acol++)
         {
            a = collect_pixel_value(fits,acol,row);
            b = collect_pixel_value(fits,acol+1,row);
            deposit_pixel_value(fits,col,row,a*vq+b*uq);
         }
      }
   }
}

void yshift_image_fits(fits_file *fits,double q)
{
   int h,w;
   int col,row,iq,arow,brow;
   float a,b;
   double uq,vq;

   if (q == 0) return;

   h = fits->npix[2];
   w = fits->npix[1];

   iq = (int)floor(fabs(q));
   uq = fabs(q) - iq;
   vq = 1.0 - uq;

   if (q > 0)
   {
      for (col=1;col<=w;col++)
      {
         for (row=h,brow=h-iq;row>=1;row--,brow--)
         {
            b = collect_pixel_value(fits,col,brow);
            a = collect_pixel_value(fits,col,brow-1);
            deposit_pixel_value(fits,col,row,a*uq+b*vq);
         }
      }
   }
   else
   {
      for (col=1;col<=w;col++)
      {
         for (row=1,arow=iq+1;row<=h;row++,arow++)
         {
            a = collect_pixel_value(fits,col,arow);
            b = collect_pixel_value(fits,col,arow+1);
            deposit_pixel_value(fits,col,row,a*vq+b*uq);
         }
      }
   }
}

void stat_image_fits(fits_file *fits,imgstat *s)
{
   int k,kmax,kmin;
   float *ptr,max,min;

   max=min=*(fits->pix);
   kmax=kmin=1;
   s->tot = 0.0;

   for (k=1,ptr=fits->pix;k<=fits->totpix;k++,ptr++)
   {
      s->tot += (double)(*ptr);
      if (*ptr < min)
         min = *ptr, kmin = k;
      else
         if (*ptr > max) max = *ptr, kmax = k;
   }

   s->n = fits->totpix;
   s->min = (double)min;
   s->max = (double)max;

   locate_pixel(fits,kmin,s->pixmin);
   locate_pixel(fits,kmax,s->pixmax);
}

void stat_image_file(char *name,imgstat *s)
{
   fits_file fits;

   read_image(&fits,name);
   check_real_image(&fits);

   stat_image_fits(&fits,s);

   free_binary_data(&fits);
}

void multiply_image_fits(fits_file *xfits,fits_file *yfits,fits_file *zfits)
{
   float *xptr,*yptr,*zptr;
   int pix;

   check_same_image_size(xfits,yfits);

   xptr = xfits->pix;
   yptr = yfits->pix;
   zptr = zfits->pix;

   for (pix=1;pix<=xfits->totpix;pix++,xptr++,yptr++,zptr++)
      *zptr = (*xptr) * (*yptr);
}

void add_image_fits(fits_file *xfits,fits_file *yfits,fits_file *zfits)
{
   float *xptr,*yptr,*zptr;
   int pix;

   check_same_image_size(xfits,yfits);

   xptr = xfits->pix;
   yptr = yfits->pix;
   zptr = zfits->pix;

   for (pix=1;pix<=xfits->totpix;pix++,xptr++,yptr++,zptr++)
      *zptr = *xptr + *yptr;
}

void subtract_image_fits(fits_file *xfits,fits_file *yfits,fits_file *zfits)
{
   float *xptr,*yptr,*zptr;
   int pix;

   check_same_image_size(xfits,yfits);

   xptr = xfits->pix;
   yptr = yfits->pix;
   zptr = zfits->pix;

   for (pix=1;pix<=xfits->totpix;pix++,xptr++,yptr++,zptr++)
      *zptr = *xptr - *yptr;
}

void add_image_file(char *xname,char *yname,char *zname,int pmod)
{
   fits_file xfits,yfits,zfits;

   read_image(&xfits,xname);
   read_image(&yfits,yname);

   zfits = xfits;
   allocate_whole_binary(&zfits);

   add_image_fits(&xfits,&yfits,&zfits);

   write_image(&zfits,zname);
   if (pmod) inform("Image '%s' created.",zfits.file.path);

   free_binary_data(&xfits);
   free_binary_data(&yfits);
}

void subtract_image_file(char *xname,char *yname,char *zname,int pmod)
{
   fits_file xfits,yfits,zfits;

   read_image(&xfits,xname);
   read_image(&yfits,yname);

   zfits = xfits;
   allocate_whole_binary(&zfits);

   subtract_image_fits(&xfits,&yfits,&zfits);

   write_image(&zfits,zname);
   if (pmod) inform("Image '%s' created.",zfits.file.path);

   free_binary_data(&xfits);
   free_binary_data(&yfits);
}

void image_div_scalar_fits(fits_file *xfits,fits_file *yfits,double val)
{
   float *xptr,*yptr;
   int pix;

   xptr = xfits->pix;
   yptr = yfits->pix;

   for (pix=1;pix<=xfits->totpix;pix++,xptr++,yptr++)
      *yptr = *xptr / val;
}

void image_div_scalar_file(char *xname,char *yname,double val,int pmod)
{
   fits_file xfits,yfits;

   read_image(&xfits,xname);

   yfits = xfits;
   allocate_whole_binary(&yfits);

   image_div_scalar_fits(&xfits,&yfits,val);

   write_image(&yfits,yname);
   if (pmod) inform("Image '%s' created.",yfits.file.path);

   free_binary_data(&xfits);
}

void copy_keyword_fits(fits_file *afits,fits_head *ahead,
  fits_file *bfits,fits_head *bhead,char *key)
{
  int seq;
  char card[CARD_BYTES+1];

  seq = get_fits_keyword(afits,ahead,key);
  extract_fits_card(ahead,seq,card);
  write_fits_card(bfits,bhead,card);
}

void copy_keyword_file(char *src,char *dest,char *key,int pmod)
{
  fits_file afits,bfits;
  fits_head *head;

  head = read_fits_header(&afits,src);

  if (fits_table_ok(dest))
  {
    read_raw_table(&bfits,dest);
    copy_keyword_fits(&afits,head,&bfits,&bfits.extend,key);
    write_raw_table(&bfits,bfits.file.path);
  }
  else
  {
    read_raw_image(&bfits,dest);
    copy_keyword_fits(&afits,head,&bfits,&bfits.head,key);
    write_raw_image(&bfits,bfits.file.path);
  }

  if (pmod) inform("\nKeyword '%s' copied from '%s' to '%s'.",key,src,dest);
}

void create_image(char *name,fits_head *head,int bitpix,int naxis,...)
{
   fits_file fits;
   va_list ap;
   int axis;
   fitskey_type key;
   char *card,text[CARD_BYTES+1];

   set_fits_name(&fits,name);

   fits.bitpix = bitpix;
   fits.naxis = naxis;
   fits.pixsize = abs(bitpix)/8;

   clear_primary_header(&fits);

   write_keyword_logical(&fits,&fits.head,"SIMPLE",'T',"FITS format");
   write_keyword_integer(&fits,&fits.head,"BITPIX",bitpix,"Data type");
   write_keyword_integer(&fits,&fits.head,"NAXIS",naxis,"Number of axes");

   va_start(ap,naxis);
   for (axis=1,fits.totpix=1;axis<=naxis;axis++)
   {
      fits.npix[axis] = va_arg(ap,int);
      fits.totpix *= fits.npix[axis];
   }
   for (axis=1;axis<=naxis;axis++)
   {
      fits.crpix[axis] = va_arg(ap,double);
      fits.crval[axis] = va_arg(ap,double);
      fits.cdelt[axis] = va_arg(ap,double);
   }
   va_end(ap);

   for (axis=1;axis<=naxis;axis++)
   {
      sprintf(key,"NAXIS%d",axis);
      write_keyword_integer(&fits,&fits.head,key,fits.npix[axis],
                                                      "Number of pixels");
   }

   for (axis=1;axis<=naxis;axis++)
   {
      sprintf(key,"CRPIX%d",axis);
      write_keyword_double(&fits,&fits.head,key,fits.crpix[axis],
                                                      "Reference pixel");
      sprintf(key,"CRVAL%d",axis);
      write_keyword_double(&fits,&fits.head,key,fits.crval[axis],
                                        "Coordinate at reference pixel");
      sprintf(key,"CDELT%d",axis);
      write_keyword_double(&fits,&fits.head,key,fits.cdelt[axis],
                                        "Coordinate increment per pixel");
      sprintf(key,"CTYPE%d",axis);
      write_keyword_textual(&fits,&fits.head,key,"PIXEL",
                                        "Units of coordinate");
   }

   write_keyword_integer(&fits,&fits.head,"DATAMIN",0,
                                     "Data minimum not defined");
   write_keyword_integer(&fits,&fits.head,"DATAMAX",0,
                                     "Data maximum not defined");

   if (head != NULL)
   {
      for (card=head->dat;!same_text(card,"END       ");card+=CARD_BYTES)
      {
         if (same_text(card,"SIMPLE  =")) continue;
         if (same_text(card,"BITPIX  =")) continue;
         if (same_text(card,"NAXIS   =")) continue;
         for (axis=1;axis<=MAXNAXIS;axis++)
         {
            sprintf(text,"NAXIS%-3d=",axis);
            if (same_text(card,text)) break;
            sprintf(text,"CRPIX%-3d=",axis);
            if (same_text(card,text)) break;
            sprintf(text,"CRVAL%-3d=",axis);
            if (same_text(card,text)) break;
            sprintf(text,"CDELT%-3d=",axis);
            if (same_text(card,text)) break;
            sprintf(text,"CTYPE%-3d=",axis);
            if (same_text(card,text)) break;
         }
         if (same_text(card,text)) continue;
         if (same_text(card,"BUNIT   =")) continue;
         if (same_text(card,"DATE    =")) continue;
         if (same_text(card,"DATAMIN =")) continue;
         if (same_text(card,"DATAMAX =")) continue;
         memmove(text,card,CARD_BYTES);
         text[CARD_BYTES] = 0;
         insert_new_card(&fits,&fits.head,text);
      }
   }

   set_img_totbinrec(&fits);
   fits.bin = NULL;
   fits.pix = NULL;

   fits.binstart = fits.head.nrec * RECORD_BYTES;

   allocate_whole_binary(&fits);
   write_raw_image(&fits,name);
}

/* ------------------------ FITS tables -------------------------- */

void read_column_format(char *form,char *type,int *size)
{
   char *ptr;
   int count;

   if (!digit_ok(*form))
      get_error("read_column_format: Bad format '%s'!",form);
   count = strtol(form,&ptr,10);

   switch(*ptr)
   {
      case 'J': *type = 'I';
                *size = 4;
                break;
      case 'E': *type = 'R';
                *size = 4;
                break;
      case 'D': *type = 'D';
                *size = 8;
                break;
      case 'A': *type = 'C';
                *size = count;
                break;
      default:  *type = '?';
   }

   if (*type == '?')
      get_error("read_column_format: Bad format '%s'!",form);

   if (!blank_string(ptr+1))
      get_error("read_column_format: Bad format '%s'!",form);

   if (*type!='C' && count!=1)
      get_error("read_column_format: Bad format '%s'!",form);

   if (*type=='C' && count<1)
      get_error("read_column_format: Bad format '%s'!",form);
}

void vcreate_fits_table(char **coldef,int nrows,char *name,va_list ap)
{
   fits_file fits;
   int k,col,seq,*val;
   fitskey_type key;
   va_list aq;

   va_copy(aq,ap);
   vset_fits_name(&fits,name,aq);
   va_end(aq);

   clear_primary_header(&fits);
   clear_extension_header(&fits);

   fits.nrows = nrows;
   for (k=0,fits.ncols=0;*coldef[k]!='#';k+=4,fits.ncols++);

   for (col=1,k=0,fits.rowlen=0;col<=fits.ncols;col++,k+=4)
   {
      read_column_format(coldef[k+2],&fits.coltype[col],&fits.colsize[col]);
      fits.rowlen += fits.colsize[col];
   }

   write_keyword_logical(&fits,&fits.head,"SIMPLE",'T',
                                        "FITS format");
   write_keyword_integer(&fits,&fits.head,"BITPIX",8,
                                        "Required for a binary table");
   write_keyword_integer(&fits,&fits.head,"NAXIS",0,
                                        "No primary data");
   write_keyword_logical(&fits,&fits.head,"EXTEND",'T',
                                        "There are extensions");

   write_keyword_textual(&fits,&fits.extend,"XTENSION","BINTABLE",
                                        "Binary table");
   write_keyword_integer(&fits,&fits.extend,"BITPIX",8,
                                        "Character format");
   write_keyword_integer(&fits,&fits.extend,"NAXIS",2,
                                        "Two axes");
   write_keyword_integer(&fits,&fits.extend,"NAXIS1",fits.rowlen,
                                        "Number of bytes per row");
   write_keyword_integer(&fits,&fits.extend,"NAXIS2",fits.nrows,
                                        "Number of rows");
   write_keyword_integer(&fits,&fits.extend,"PCOUNT",0,
                                        "Required for a binary table");
   write_keyword_integer(&fits,&fits.extend,"GCOUNT",1,
                                        "Required for a binary table");
   write_keyword_integer(&fits,&fits.extend,"TFIELDS",fits.ncols,
                                        "Number of columns");

   for (k=0,col=1;*coldef[k]!='#';k+=4,col++)
   {
      sprintf(key,"TTYPE%d",col);
      write_keyword_textual(&fits,&fits.extend,key,coldef[k],
                                           "Column type");
      sprintf(key,"TUNIT%d",col);
      if (*coldef[k+1] != 0)
         write_keyword_textual(&fits,&fits.extend,key,coldef[k+1],
                                           "Units");
      sprintf(key,"TFORM%d",col);
      write_keyword_textual(&fits,&fits.extend,key,coldef[k+2],
                                           "Data type");
      sprintf(key,"TDISP%d",col);
      write_keyword_textual(&fits,&fits.extend,key,coldef[k+3],
                                           "Display format");
   }

   set_tbl_totbinrec(&fits);

   allocate_whole_binary(&fits);

   val = (int *)fits.bin;
   for (seq=0;seq<fits.nrows;seq++) val[seq] = seq+1;

   vwrite_table(&fits,name,ap);
}

int find_table_column(fits_file *t,char *c)
{
   int col;

   for (col=1;col<=t->ncols && strcmp(t->colname[col],c);col++);
   if (strcmp(t->colname[col],c)) col = -1;

   return(col);
}

int get_table_column(fits_file *t,char *c)
{
   int col;

   col = find_table_column(t,c);
   if (col < 1)
      get_error("Table column '%s' not found in '%s'!",c,t->file.name);

   return(col);
}

void get_several_table_columns(fits_file *t,char *list,int *c,int n)
{
   int i,k;
   char cname[MAXCOLNAME+1],*src,*dest;

   for (i=1,src=list;i<=n;i++)
   {
      memset(cname,0,MAXCOLNAME+1);
      for (k=0,dest=cname;k<MAXCOLNAME && *src>',';k++,src++,dest++)
         *dest = *src;
      if (k >= MAXCOLNAME) get_error("Table column name too long!");
      c[i] = get_table_column(t,cname);
      if (i < n)
      {
         if (*src == ',')
            src++;
         else
            get_error("A comma expected after a table column name!");
      }
   }
   if (*src != 0) get_error("End of character string expected "
                            "after the last table column name!");
}

void read_table_column_double(fits_file *t,int col,double *x)
{
   int row;
   int *iptr;
   float *fptr;
   char *offset;

   offset = (char *)t->bin + t->coloffset[col];
   switch(t->coltype[col])
   {
      case 'I': iptr = (int *)offset;
                for (row=1,x++;row<=t->nrows;row++,iptr++,x++)
                   *x = (double)(*iptr);
                break;
      case 'R': fptr = (float *)offset;
                for (row=1,x++;row<=t->nrows;row++,fptr++,x++)
                   *x = (double)(*fptr);
                break;
      case 'D': memmove(&x[1],offset,t->nrows*sizeof(double));
                break;
      default:  get_error("Bad column type '%c'!",t->coltype[col]);
                x=0;
   }
}

void read_table_column_integer(fits_file *t,int col,int *x)
{
   int row;
   float *fptr;
   double *dptr;
   char *offset;

   offset = (char *)t->bin + t->coloffset[col];
   switch(t->coltype[col])
   {
      case 'D': dptr = (double *)offset;
                for (row=1,x++;row<=t->nrows;row++,dptr++,x++)
                   *x = (int)(*dptr);
                break;
      case 'R': fptr = (float *)offset;
                for (row=1,x++;row<=t->nrows;row++,fptr++,x++)
                   *x = (int)(*fptr);
                break;
      case 'I': memmove(&x[1],offset,t->nrows*sizeof(int));
                break;
      default:  get_error("Bad column type '%c'!",t->coltype[col]);
                x=0;
   }
}

void write_table_column_double(fits_file *t,int col,double *x)
{
   int row;
   int *iptr;
   float *fptr;
   char *offset;

   offset = (char *)t->bin + t->coloffset[col];
   switch(t->coltype[col])
   {
      case 'I': iptr = (int *)offset;
                for (row=1,x++;row<=t->nrows;row++,iptr++,x++)
                   *iptr = (int)(*x);
                break;
      case 'R': fptr = (float *)offset;
                for (row=1,x++;row<=t->nrows;row++,fptr++,x++)
                   *fptr = (float)(*x);
                break;
      case 'D': memmove(offset,&x[1],t->nrows*sizeof(double));
                break;
      default:  get_error("Bad column type '%c'!",t->coltype[col]);
                x=0;
   }
}

void write_table_column_integer(fits_file *t,int col,int *x)
{
   int row;
   float *fptr;
   double *dptr;
   char *offset;

   offset = (char *)t->bin + t->coloffset[col];
   switch(t->coltype[col])
   {
      case 'D': dptr = (double *)offset;
                for (row=1,x++;row<=t->nrows;row++,dptr++,x++)
                   *dptr = (double)(*x);
                break;
      case 'R': fptr = (float *)offset;
                for (row=1,x++;row<=t->nrows;row++,fptr++,x++)
                   *fptr = (float)(*x);
                break;
      case 'I': memmove(offset,&x[1],t->nrows*sizeof(int));
                break;
      default:  get_error("Bad column type '%c'!",t->coltype[col]);
                x=0;
   }
}
