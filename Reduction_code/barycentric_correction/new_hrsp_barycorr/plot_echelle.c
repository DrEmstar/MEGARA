/* ------------------------------------------------------------------------

   Program:  plot_echelle
   Purpose:  Generate a PostScript plot of echelle orders in a given image.

   ------------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#include "general.h"
#include "angle.h"
#include "datetime.h"
#include "vector.h"
#include "numeric.h"
#include "nutation.h"
#include "astro.h"
#include "config.h"
#include "fits.h"
#include "process.h"

#define MAX_ORD 1000

#define PIX_SCALE 1
#define WAV_SCALE 2
#define LOG_SCALE 3
#define CCF_SCALE 4

typedef struct
{
  int axis;
  int size;
  int first_tick,last_tick;
  int mark_freq,tick_freq;
  double first_mark,last_mark;
  char fmt[32];
  char title[80];
  int expo;
}
plot_axis;

fits_file fits;
plt_type plt;
double *startval,*stepval;
double xstart[MAX_ORD+1],xstep[MAX_ORD+1];
int sel[MAX_ORD+1];
file_type ps;
double ymin,ymax;

void pps(char *fmt,...)
{
  va_list ap;

  va_start(ap,fmt);
  vfprintf(ps.dat,fmt,ap);
  va_end(ap);
  fprintf(ps.dat,"\n");
}

void bps(void)
{
  pps("");
}

void init_ps(int npages)
{
  pps("%%!PS-Adobe-3.0");
  pps("%%%%BoundingBox: 0 0 595 842");
  pps("%%%%Pages: %d",npages);
  pps("%%EndComments");
  bps();
  pps("/rshow");
  pps("{");
  pps("  dup stringwidth pop neg 0 rmoveto show");
  pps("}");
  pps("def");
  bps();
  pps("/cshow");
  pps("{");
  pps("  dup stringwidth pop 2 div neg 0 rmoveto show");
  pps("}");
  pps("def");
  bps();
}

void endpage(void)
{
  pps("showpage");
}

void newpage(int page)
{
  if (page > 1) bps();
  pps("%%%%Page: %d %d",page,page);
}

void draw_axis(plot_axis *a)
{
  int n,i,t,s;
  double val;
  char mark[16],template[32];

  pps("gsave");
  if (a->axis == YAXIS) pps("90 rotate");
  pps("newpath 0 0 moveto %d 0 rlineto stroke",a->size);
  n = a->mark_freq * a->tick_freq;
  for (i=0;i<=n;i++)
  {
    t = ((n-i)*a->first_tick + i*a->last_tick) / n;
    if (t < 0 || t > a->size) continue;
    s = plt.ticksize;
    if (a->axis == XAXIS) s = -s;
    if ((i % a->tick_freq) == 0)
    {
      s = (s * 3) / 2;
      val = ((n-i)*a->first_mark + i*a->last_mark) / (double)n;
      sprintf(mark,a->fmt,val);
      if (a->axis == YAXIS)
      {
        pps("gsave");
        pps("%d %d translate",t-plt.fontcent,s+700);
        pps("-90 rotate");
        pps("0 0 moveto (%s) rshow",mark);
        pps("grestore");
      }
      else
      {
        pps("gsave");
        pps("%d %d translate",t,s-plt.fontsize);
        pps("0 0 moveto (%s) cshow",mark);
        pps("grestore");
      }
    }
    pps("newpath %d 0 moveto 0 %d rlineto stroke",t,s);
  }
  if (a->axis == XAXIS)
  {
    pps("%d %d moveto",a->size/2,-(plt.ticksize*3)/2-2*plt.fontsize);
    pps("(%s) cshow",a->title);
  }
  else
  {
    strcpy(template,a->title);
    if (a->expo > 0) strcat(template,"u(x10u)");
    if (a->expo < 0) strcat(template,"u(x10-u)");
    pps("%d %d moveto",a->size/2,(plt.ticksize*3)/2);
    pps("(%s) stringwidth pop 2 div neg 0 rmoveto",template);
    pps("(-0.0) stringwidth pop 0 exch 700 add rmoveto");
    pps("0 %d rmoveto",plt.fontsize / 2);
    pps("(%s) show",a->title);
    if (a->expo != 0)
    {
      pps("%d 0 rmoveto",plt.fontsize / 2);
      pps("(\\() show");
      pps("(x10) show");
      pps("0 %d rmoveto",(plt.fontsize*5)/10);
      pps("/Helvetica findfont %d scalefont setfont",(plt.fontsize*8)/10);
      pps("(%d) show",a->expo);
      pps("0 %d rmoveto",-(plt.fontsize*5)/10);
      pps("/Helvetica findfont %d scalefont setfont",plt.fontsize);
      pps("(\\)) show");
    }
  }
  pps("grestore");
}

void newrange(double *yrange,double *ymag)
{
  *yrange = ymax - ymin;
  if (fabs(ymin) > fabs(ymax))
    *ymag = fabs(ymin);
  else
    *ymag = fabs(ymax);
}

void set_ylim(float *fptr,int n,int m,int k)
{
  int i,j;
  double val,g,u,yrange,ymag;
  float *p;

  p = fptr;
  ymin = ymax = (double)(*p);
  for (i=0;i<k;i++)
  {
    for (j=0;j<n;j++)
    {
      val = (double)(*(p++));
      if (val < ymin) ymin = val;
      if (val > ymax) ymax = val;
    }
    p += m;
  }

  newrange(&yrange,&ymag);

  if (ymag == 0.0)
  {
    ymin = -1.0;
    ymax = 1.0;
    newrange(&yrange,&ymag);
  }

  if (yrange/ymag < 0.1)
  {
    g = floor(log10(ymag));
    u = pow(10.0,g);
    ymin -= u;
    ymax += u;
    newrange(&yrange,&ymag);
  }
}

void set_xaxis(plot_axis *a,int size,double start,double step,int count)
{
  int scale;
  double h,s,r,g,u,b,c;
  double pa,pb;
  int f;

  a->axis = XAXIS;
  a->size = size;

  scale = PIX_SCALE;
  if (step < 1.0)
  {
    scale = WAV_SCALE;
    if (start < 10.0)
    {
      scale = LOG_SCALE;
      if (start < 1.0) scale = CCF_SCALE;
    }
  }

  switch(scale)
  {
    case WAV_SCALE:
      if (step >= 1.0)
        get_error("Wavelength scale: Unexpected STEPVAL!");
      if (start < 2000.0 || start >= 10000.0)
        get_error("Wavelength scale: Unexpected STARTVAL!");
      strcpy(a->title,"Wavelength (Angstrom)");
      break;

    case LOG_SCALE:
      if (step >= 1.0)
        get_error("Logarithmic scale: Unexpected STEPVAL!");
      if (start >= 10.0)
        get_error("Logarithmic scale: Unexpected STARTVAL!");
      strcpy(a->title,"Logarithm of Wavelength");
      break;

    case CCF_SCALE:
      if (step >= 1.0)
        get_error("CCF scale: Unexpected STEPVAL!");
      if (start >= 1.0)
        get_error("CCF scale: Unexpected STARTVAL!");
      strcpy(a->title,"Logarithmic Shift");
      break;

    default:
      if (step != 1.0)
        get_error("Pixel scale: Unexpected STEPVAL!");
      if (floor(start) != start)
        get_error("Pixel scale: Unexpected STARTVAL!");
      strcpy(a->title,"Pixel number");
  }

  h = (double)(count - 1) * step;
  s = (double)a->size / h;
  r = (double)plt.minsep[1] / s;
  g = floor(log10(r));
  u = pow(10.0,g);
  b = r / u;
  c = 1.0, f = 5;
  if (c < b) c = 2.0, f = 4;
  if (c < b) c = 5.0, f = 5;
  if (c < b) c = 10.0, f = 5;
  c *= u;

  pa = floor(start/c);
  pb = ceil((start+h)/c);

  if (g < 0.0)
    sprintf(a->fmt,"%%1.%df",(int)fabs(g));
  else
    strcpy(a->fmt,"%1.0f");

  a->first_mark = pa * c;
  a->last_mark = pb * c;
  a->mark_freq = (int)(pb - pa);

  a->first_tick = nint((a->first_mark - start) * s);
  a->last_tick = nint((a->last_mark - start) * s);
  a->tick_freq = f;
}

void set_yaxis(plot_axis *a,int size)
{
  double s,r,g,u,b,c;
  double pa,pb;
  int f;

  a->axis = YAXIS;
  a->size = size;

  s = (double)a->size / (ymax-ymin);
  r = (double)plt.minsep[2] / s;
  g = floor(log10(r));
  u = pow(10.0,g);
  b = r / u;
  c = 1.0, f = 5;
  if (c < b) c = 2.0, f = 4;
  if (c < b) c = 5.0, f = 5;
  if (c < b) c = 10.0, f = 5;

  c *= u;

  pa = floor(ymin/c);
  pb = ceil(ymax/c);

  if (u != 1.0) u *= 10.0;

  a->first_mark = pa * c/u;
  a->last_mark = pb * c/u;
  a->mark_freq = (int)(pb - pa);

  a->first_tick = nint((a->first_mark*u - ymin) * s);
  a->last_tick = nint((a->last_mark*u - ymin) * s);
  a->tick_freq = f;
  strcpy(a->fmt,"%3.1f");

  strcpy(a->title,"Intensity");
  a->expo = (int)log10(u);
}

void plot_order
  (int ord,int row,int col,float *fptr,int n,double start,double step)
{
  int x,y;
  int i;
  char ordtxt[16];
  plot_axis x_axis,y_axis;

  x = plt.left_margin + (col - 1) * (plt.boxw + plt.colgap);
  y = plt.bottom_margin + (plt.numrows - row) * (plt.boxh + plt.rowgap);

  sprintf(ordtxt,"Order %d",ord);

  set_xaxis(&x_axis,plt.boxw,start,step,n);
  set_yaxis(&y_axis,plt.boxh);

  pps("gsave");
  pps("72 25400 div dup scale");
  pps("%d %d translate",x,y);
  pps("/Helvetica findfont %d scalefont setfont",plt.fontsize);
  pps("%d setlinewidth",plt.linewidth);

  draw_axis(&x_axis);
  draw_axis(&y_axis);

  pps("newpath %d 0 moveto 0 %d rlineto %d 0 rlineto stroke",
    plt.boxw,plt.boxh,-plt.boxw);
  pps("%d %d moveto (%s) rshow",plt.boxw,plt.boxh+700,ordtxt);
  pps("%d %d moveto (%s) show",0,plt.boxh+700,fits.file.name);
  pps("0 0 1 setrgbcolor");
  pps("newpath");
  for (i=0;i<n;i++)
  {
    x = (plt.boxw * i) / (n - 1);
    y = nint(((double)fptr[i]-ymin)/(ymax-ymin) * (double)plt.boxh);
    if (i == 0)
      pps("%d %d moveto",x,y);
    else
      pps("%d %d lineto",x,y);
  }
  pps("stroke");
  pps("grestore");
}

int main(int argc,char **argv)
{
  char spectrograph[MAX_SPECNAME+1];
  int first_bin,last_bin,numbins;
  int count,numord,firstord,lastord,ord,seq;
  double axis_start;
  int npages,len;
  path_name_type pspath;
  char *eptr;
  int page,row,col,empty;
  float *fptr,*rptr;
  int nsel;

  if (argc != 5) get_error("Usage:  plot_echelle <img> <start> <end> <ord>");

  read_image(&fits,argv[1]);
  check_image_2d(&fits);

  get_spectrograph(&fits,&fits.head,spectrograph);

  first_bin = atoi(argv[2]);
  last_bin = atoi(argv[3]);

  if (first_bin < 1 || first_bin > fits.npix[1])
    get_error("Start pixel out of range!");

  if (last_bin < 1 || last_bin > fits.npix[1])
    get_error("End pixel out of range!");

  numbins = last_bin - first_bin + 1;

  if (numbins < 8) get_error("Invalid pixel selection!");

  numord = fits.npix[2];
  firstord = nint(world_coordinate(&fits,2,(double)numord));
  lastord = nint(world_coordinate(&fits,2,1.0));
  if (firstord < 1 || firstord > MAX_ORD)
    get_error("First order number out of limits in '%s'!",fits.file.path);
  if (lastord < 1 || lastord > MAX_ORD)
    get_error("Last order number out of limits in '%s'!",fits.file.path);
  if (lastord-firstord != numord-1)
    get_error("Unexpected order numbers in '%s'!",fits.file.path);

  startval = dvector(numord);
  stepval = dvector(numord);

  memset(xstart,0,sizeof(xstart));
  memset(xstep,0,sizeof(xstep));

  get_double_array(&fits,&fits.head,"STARTVAL",&startval[1],numord);
  get_double_array(&fits,&fits.head,"STEPVAL",&stepval[1],numord);

  for (seq=1,ord=lastord;seq<=numord;seq++,ord--)
  {
    xstart[ord] = startval[seq];
    xstep[ord] = stepval[seq];
  }

  count = collect_integer_list(argv[4],sel,MAX_ORD,"Order number");
  if (count < 1) get_error("No orders selected!");

  for (ord=1,nsel=0;ord<=MAX_ORD;ord++)
  {
    if (!sel[ord]) continue;
    if (ord < firstord || ord > lastord)
      get_error("Order number out of range!");
    nsel++;
  }

  if (nsel != count) get_error("Order count mismatch!");

  get_plot_info(&plt,spectrograph);

  npages = (nsel - 1) / plt.plots_per_page + 1;

  strcpy(pspath,fits.file.path);
  len = strlen(pspath);
  if (len < 4) get_error("File name too short: '%s'!",fits.file.path);
  eptr = pspath + len - 4;
  if (strcmp(eptr,".fit") != 0)
    get_error("Extension '.fit' expected in '%s'!",fits.file.path);
  strcpy(eptr,".ps");
  set_file_name(&ps,pspath);

  if (plt.global_limits)
  {
    if (plt.window_limits)
      set_ylim(fits.pix+first_bin-1,numbins,fits.npix[1]-numbins,fits.npix[2]);
    else
      set_ylim(fits.pix,fits.npix[1],0,fits.npix[2]);
  }

  get_fopen(&ps,"w");

  init_ps(npages);

  rptr = fits.pix + fits.totpix - fits.npix[1];
  fptr = rptr + first_bin - 1;
  page = 1;
  row = col = 1;
  empty = -1;
  for (ord=firstord;ord<=lastord;ord++)
  {
    if (sel[ord])
    {
      if (empty) newpage(page);
      axis_start = xstart[ord]+(first_bin-1)*xstep[ord];
      if (!plt.global_limits)
      {
        if (plt.window_limits)
          set_ylim(fptr,numbins,0,1);
        else
          set_ylim(rptr,fits.npix[1],0,1);
      }
      plot_order(ord,row,col,fptr,numbins,axis_start,xstep[ord]);
      empty = 0;
      col++;
      if (col > plt.numcols)
      {
        col = 1;
        row++;
        if (row > plt.numrows)
        {
          row = 1;
          page++;
          endpage();
          empty=-1;
        }
      }
    }
    rptr -= fits.npix[1];
    fptr -= fits.npix[1];
  }
  if (!empty) endpage();

  get_fclose(&ps);

  free_binary_data(&fits);

  free_dvector(startval);
  free_dvector(stepval);

  check_memory();

  return(0);
}
