#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#include "general.h"
#include "vector.h"
#include "angle.h"
#include "numeric.h"
#include "astro.h"

#define HIP_TYPE      1
#define HD_TYPE       2
#define HR_TYPE       3
#define FLM_TYPE      4
#define BYR_TYPE      5

#define HRLIMIT   10000
#define HDLIMIT  360000
#define HIPLIMIT 130000

#define HIPTOT   118218
#define STARLEN  68

#define INTERPOLATION_RADIUS   5

#define MEAN_POSITION 0
#define TRUE_POSITION 1

#define ZETA   0
#define THETA  1
#define ZEE    2

#define arg_l      0
#define arg_Omega  4

#define RA_GAL    192.85948
#define DEC_GAL    27.12825
#define NODE_GAL   32.93192

#define arc_sec  4.848136811095359936e-6
#define arc_tsec 7.272205216643039904e-5

#define parsec_km      3.08567758e13

#define unit_distance 149597870.66
#define unit_velocity 1731.4568363425926

#define LIGHT_TIME 499.004782

char *constell[] = 
   { "And","Ant","Aps","Aqr","Aql","Ara","Ari","Aur","Boo","Cae","Cam",
     "Cnc","CVn","CMa","CMi","Cap","Car","Cas","Cen","Cep","Cet","Cha",
     "Cir","Col","Com","CrA","CrB","Crv","Crt","Cru","Cyg","Del","Dor",
     "Dra","Equ","Eri","For","Gem","Gru","Her","Hor","Hya","Hyi","Ind",
     "Lac","Leo","LMi","Lep","Lib","Lup","Lyn","Lyr","Men","Mic","Mon",
     "Mus","Nor","Oct","Oph","Ori","Pav","Peg","Per","Phe","Pic","Psc",
     "PsA","Pup","Pyx","Ret","Sge","Sgr","Sco","Scl","Sct","Ser","Sex",
     "Tau","Tel","Tri","TrA","Tuc","UMa","UMi","Vel","Vir","Vol","Vul",
     "#" };

char *greek[] =
   { "alp","bet","gam","del","eps","zet","eta","the",
     "iot","kap","lam","mu","nu","xi","omi","pi",
     "rho","sig","tau","ups","phi","chi","psi","ome",
     "#" };

double c_prec[3][3]= { {2306.2181,  0.30188,  0.017998},
                       {2004.3109, -0.42665, -0.041833},
                       {2306.2181,  1.09468,  0.018203} };

double ec[4] = { 84381.448, -46.8150, -0.00059, 0.001813 };

double frq[5][5] =
   { { 485866.733,  715922.633,  31.310,  0.064, 1325},
     {1287099.804, 1292581.224,  -0.577, -0.012,   99},
     { 335778.877,  295263.137, -13.257,  0.011, 1342},
     {1072261.307, 1105601.328,  -6.891,  0.019, 1236},
     { 450160.280, -482890.539,   7.455,  0.008,   -5} };

double st[4] = { 24110.54841, 8640184.812866, 0.093104, -6.2e-6 };

double equatorial_radius = 6378.137;
double inverse_flattening = 298.257223563;
double angular_velocity = 7.292115e-5;

double nut[106][9]=
   { { 0,  0,  0,  0,  1, -171996, -174.2, 92025,  8.9},
     { 0,  0,  0,  0,  2,    2062,    0.2,  -895,  0.5},
     {-2,  0,  2,  0,  1,      46,    0.0,   -24,  0.0},
     { 2,  0, -2,  0,  0,      11,    0.0,     0,  0.0},
     {-2,  0,  2,  0,  2,      -3,    0.0,     1,  0.0},
     { 1, -1,  0, -1,  0,      -3,    0.0,     0,  0.0},
     { 0, -2,  2, -2,  1,      -2,    0.0,     1,  0.0},
     { 2,  0, -2,  0,  1,       1,    0.0,     0,  0.0},
     { 0,  0,  2, -2,  2,  -13187,   -1.6,  5736, -3.1},
     { 0,  1,  0,  0,  0,    1426,   -3.4,    54, -0.1},
     { 0,  1,  2, -2,  2,    -517,    1.2,   224, -0.6},
     { 0, -1,  2, -2,  2,     217,   -0.5,   -95,  0.3},
     { 0,  0,  2, -2,  1,     129,    0.1,   -70,  0.0},
     { 2,  0,  0, -2,  0,      48,    0.0,     1,  0.0},
     { 0,  0,  2, -2,  0,     -22,    0.0,     0,  0.0},
     { 0,  2,  0,  0,  0,      17,   -0.1,     0,  0.0},
     { 0,  1,  0,  0,  1,     -15,    0.0,     9,  0.0},
     { 0,  2,  2, -2,  2,     -16,    0.1,     7,  0.0},
     { 0, -1,  0,  0,  1,     -12,    0.0,     6,  0.0},
     {-2,  0,  0,  2,  1,      -6,    0.0,     3,  0.0},
     { 0, -1,  2, -2,  1,      -5,    0.0,     3,  0.0},
     { 2,  0,  0, -2,  1,       4,    0.0,    -2,  0.0},
     { 0,  1,  2, -2,  1,       4,    0.0,    -2,  0.0},
     { 1,  0,  0, -1,  0,      -4,    0.0,     0,  0.0},
     { 2,  1,  0, -2,  0,       1,    0.0,     0,  0.0},
     { 0,  0, -2,  2,  1,       1,    0.0,     0,  0.0},
     { 0,  1, -2,  2,  0,      -1,    0.0,     0,  0.0},
     { 0,  1,  0,  0,  2,       1,    0.0,     0,  0.0},
     {-1,  0,  0,  1,  1,       1,    0.0,     0,  0.0},
     { 0,  1,  2, -2,  0,      -1,    0.0,     0,  0.0},
     { 0,  0,  2,  0,  2,   -2274,   -0.2,   977, -0.5},
     { 1,  0,  0,  0,  0,     712,    0.1,    -7,  0.0},
     { 0,  0,  2,  0,  1,    -386,   -0.4,   200,  0.0},
     { 1,  0,  2,  0,  2,    -301,    0.0,   129, -0.1},
     { 1,  0,  0, -2,  0,    -158,    0.0,    -1,  0.0},
     {-1,  0,  2,  0,  2,     123,    0.0,   -53,  0.0},
     { 0,  0,  0,  2,  0,      63,    0.0,    -2,  0.0},
     { 1,  0,  0,  0,  1,      63,    0.1,   -33,  0.0},
     {-1,  0,  0,  0,  1,     -58,   -0.1,    32,  0.0},
     {-1,  0,  2,  2,  2,     -59,    0.0,    26,  0.0},
     { 1,  0,  2,  0,  1,     -51,    0.0,    27,  0.0},
     { 0,  0,  2,  2,  2,     -38,    0.0,    16,  0.0},
     { 2,  0,  0,  0,  0,      29,    0.0,    -1,  0.0},
     { 1,  0,  2, -2,  2,      29,    0.0,   -12,  0.0},
     { 2,  0,  2,  0,  2,     -31,    0.0,    13,  0.0},
     { 0,  0,  2,  0,  0,      26,    0.0,    -1,  0.0},
     {-1,  0,  2,  0,  1,      21,    0.0,   -10,  0.0},
     {-1,  0,  0,  2,  1,      16,    0.0,    -8,  0.0},
     { 1,  0,  0, -2,  1,     -13,    0.0,     7,  0.0},
     {-1,  0,  2,  2,  1,     -10,    0.0,     5,  0.0},
     { 1,  1,  0, -2,  0,      -7,    0.0,     0,  0.0},
     { 0,  1,  2,  0,  2,       7,    0.0,    -3,  0.0},
     { 0, -1,  2,  0,  2,      -7,    0.0,     3,  0.0},
     { 1,  0,  2,  2,  2,      -8,    0.0,     3,  0.0},
     { 1,  0,  0,  2,  0,       6,    0.0,     0,  0.0},
     { 2,  0,  2, -2,  2,       6,    0.0,    -3,  0.0},
     { 0,  0,  0,  2,  1,      -6,    0.0,     3,  0.0},
     { 0,  0,  2,  2,  1,      -7,    0.0,     3,  0.0},
     { 1,  0,  2, -2,  1,       6,    0.0,    -3,  0.0},
     { 0,  0,  0, -2,  1,      -5,    0.0,     3,  0.0},
     { 1, -1,  0,  0,  0,       5,    0.0,     0,  0.0},
     { 2,  0,  2,  0,  1,      -5,    0.0,     3,  0.0},
     { 0,  1,  0, -2,  0,      -4,    0.0,     0,  0.0},
     { 1,  0, -2,  0,  0,       4,    0.0,     0,  0.0},
     { 0,  0,  0,  1,  0,      -4,    0.0,     0,  0.0},
     { 1,  1,  0,  0,  0,      -3,    0.0,     0,  0.0},
     { 1,  0,  2,  0,  0,       3,    0.0,     0,  0.0},
     { 1, -1,  2,  0,  2,      -3,    0.0,     1,  0.0},
     {-1, -1,  2,  2,  2,      -3,    0.0,     1,  0.0},
     {-2,  0,  0,  0,  1,      -2,    0.0,     1,  0.0},
     { 3,  0,  2,  0,  2,      -3,    0.0,     1,  0.0},
     { 0, -1,  2,  2,  2,      -3,    0.0,     1,  0.0},
     { 1,  1,  2,  0,  2,       2,    0.0,    -1,  0.0},
     {-1,  0,  2, -2,  1,      -2,    0.0,     1,  0.0},
     { 2,  0,  0,  0,  1,       2,    0.0,    -1,  0.0},
     { 1,  0,  0,  0,  2,      -2,    0.0,     1,  0.0},
     { 3,  0,  0,  0,  0,       2,    0.0,     0,  0.0},
     { 0,  0,  2,  1,  2,       2,    0.0,    -1,  0.0},
     {-1,  0,  0,  0,  2,       1,    0.0,    -1,  0.0},
     { 1,  0,  0, -4,  0,      -1,    0.0,     0,  0.0},
     {-2,  0,  2,  2,  2,       1,    0.0,    -1,  0.0},
     {-1,  0,  2,  4,  2,      -2,    0.0,     1,  0.0},
     { 2,  0,  0, -4,  0,      -1,    0.0,     0,  0.0},
     { 1,  1,  2, -2,  2,       1,    0.0,    -1,  0.0},
     { 1,  0,  2,  2,  1,      -1,    0.0,     1,  0.0},
     {-2,  0,  2,  4,  2,      -1,    0.0,     1,  0.0},
     {-1,  0,  4,  0,  2,       1,    0.0,     0,  0.0},
     { 1, -1,  0, -2,  0,       1,    0.0,     0,  0.0},
     { 2,  0,  2, -2,  1,       1,    0.0,    -1,  0.0},
     { 2,  0,  2,  2,  2,      -1,    0.0,     0,  0.0},
     { 1,  0,  0,  2,  1,      -1,    0.0,     0,  0.0},
     { 0,  0,  4, -2,  2,       1,    0.0,     0,  0.0},
     { 3,  0,  2, -2,  2,       1,    0.0,     0,  0.0},
     { 1,  0,  2, -2,  0,      -1,    0.0,     0,  0.0},
     { 0,  1,  2,  0,  1,       1,    0.0,     0,  0.0},
     {-1, -1,  0,  2,  1,       1,    0.0,     0,  0.0},
     { 0,  0, -2,  0,  1,      -1,    0.0,     0,  0.0},
     { 0,  0,  2, -1,  2,      -1,    0.0,     0,  0.0},
     { 0,  1,  0,  2,  0,      -1,    0.0,     0,  0.0},
     { 1,  0, -2, -2,  0,      -1,    0.0,     0,  0.0},
     { 0, -1,  2,  0,  1,      -1,    0.0,     0,  0.0},
     { 1,  1,  0, -2,  1,      -1,    0.0,     0,  0.0},
     { 1,  0, -2,  2,  0,      -1,    0.0,     0,  0.0},
     { 2,  0,  0,  2,  0,       1,    0.0,     0,  0.0},
     { 0,  0,  2,  4,  2,      -1,    0.0,     0,  0.0},
     { 0,  1,  0,  1,  0,       1,    0.0,     0,  0.0} };

double ephemeris_correction[]=
   { -9999 ,   0.00,
      1801, 13.40, 1802, 13.10, 1803, 12.90, 1804, 12.70, 1805, 12.60,
      1806, 12.50, 1807, 12.50, 1808, 12.50, 1809, 12.50, 1810, 12.50,
      1811, 12.50, 1812, 12.50, 1813, 12.50, 1814, 12.50, 1815, 12.50,
      1816, 12.50, 1817, 12.40, 1818, 12.30, 1819, 12.20, 1820, 12.00,
      1821, 11.70, 1822, 11.40, 1823, 11.10, 1824, 10.60, 1825, 10.20,
      1826,  9.60, 1827,  9.10, 1828,  8.60, 1829,  8.00, 1830,  7.50,
      1831,  7.00, 1832,  6.60, 1833,  6.30, 1834,  6.00, 1835,  5.80,
      1836,  5.70, 1837,  5.60, 1838,  5.60, 1839,  5.60, 1840,  5.70,
      1841,  5.80, 1842,  5.90, 1843,  6.10, 1844,  6.20, 1845,  6.30,
      1846,  6.50, 1847,  6.60, 1848,  6.80, 1849,  6.90, 1850,  7.10,
      1851,  7.20, 1852,  7.30, 1853,  7.40, 1854,  7.50, 1855,  7.60,
      1856,  7.70, 1857,  7.70, 1858,  7.80, 1859,  7.80, 1860,  7.88,
      1861,  7.82, 1862,  7.54, 1863,  6.97, 1864,  6.40, 1865,  6.02,
      1866,  5.41, 1867,  4.10, 1868,  2.92, 1869,  1.82, 1870,  1.61,
      1871,  0.10, 1872, -1.02, 1873, -1.28, 1874, -2.69, 1875, -3.24,
      1876, -3.64, 1877, -4.54, 1878, -4.71, 1879, -5.11, 1880, -5.40,
      1881, -5.42, 1882, -5.20, 1883, -5.46, 1884, -5.46, 1885, -5.79,
      1886, -5.63, 1887, -5.64, 1888, -5.80, 1889, -5.66, 1890, -5.87,
      1891, -6.01, 1892, -6.19, 1893, -6.64, 1894, -6.44, 1895, -6.47,
      1896, -6.09, 1897, -5.76, 1898, -4.66, 1899, -3.74, 1900, -2.72,
      1901, -1.54, 1902, -0.02, 1903,  1.24, 1904,  2.64, 1905,  3.86,
      1906,  5.37, 1907,  6.14, 1908,  7.75, 1909,  9.13, 1910, 10.46,
      1911, 11.53, 1912, 13.36, 1913, 14.65, 1914, 16.01, 1915, 17.20,
      1916, 18.24, 1917, 19.06, 1918, 20.25, 1919, 20.95, 1920, 21.16,
      1921, 22.25, 1922, 22.41, 1923, 23.03, 1924, 23.49, 1925, 23.62,
      1926, 23.86, 1927, 24.49, 1928, 24.34, 1929, 24.08, 1930, 24.02,
      1931, 24.00, 1932, 23.87, 1933, 23.95, 1934, 23.86, 1935, 23.93,
      1936, 23.73, 1937, 23.92, 1938, 23.96, 1939, 24.02, 1940, 24.33,
      1941, 24.83, 1942, 25.30, 1943, 25.70, 1944, 26.24, 1945, 26.77,
      1946, 27.28, 1947, 27.78, 1948, 28.25, 1949, 28.71, 1950, 29.15,
      1951, 29.57, 1952, 29.97, 1953, 30.36, 1954, 30.72, 1955, 31.07,
      1956, 31.35, 1957, 31.68, 1958, 32.18, 1959, 32.68, 1960, 33.15,
      1961, 33.59, 1962, 34.00, 1963, 34.47, 1964, 35.03, 1965, 35.73,
      1966, 36.54, 1967, 37.43, 1968, 38.29, 1969, 39.20, 1970, 40.18,
      1971, 41.17, 1972, 42.23, 1973, 43.37, 1974, 44.49, 1975, 45.48,
      1976, 46.46, 1977, 47.52, 1978, 48.53, 1979, 49.59, 1980, 50.54,
      1981, 51.38, 1982, 52.17, 1983, 52.96, 1984, 53.79, 1985, 54.34,
      1986, 54.87, 1987, 55.32, 1988, 55.82, 1989, 56.30, 1990, 56.86,
      1991, 57.57, 1992, 58.31, 1993, 59.12, 1994, 59.98, 1995, 60.78,
      1996, 61.63, 1997, 62.29, 1998, 62.97, 1999, 63.47, 2000, 63.83,
      2001, 64.00, 2002, 65.00, 2003, 66.00, 2004, 67.00, 2005, 68.00,
      2006, 69.00, 9999, 0.00 };

int monthlen[] = { 0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
double J2000JD = 2451545.0;
double sider = 1.002737909350795;
double vt_factor = 4.74047;
double julian_year = 365.25;
mat EQ2000_to_G,G_to_EQ2000;
double first_jd,last_jd;
int num_rec;

double *earth;

int leap_year(int year)
{
   div_t cent,quad;

   cent = div(year,100);

   if (year>1582 && cent.rem==0)
      quad = div(cent.quot,4);
   else
      quad = div(year,4);

   if (quad.rem == 0)
      return(1);
   else
      return(0);
}

void prep_feb(int year)
{
   if (leap_year(year))
      monthlen[2]=29;
   else
      monthlen[2]=28;
}

int year_length(int year)
{
   if (leap_year(year))
      return(366);
   else
      return(365);
}

long date_order(long d,long m,long y)
{
   return(d+(m+y*100L)*100L);
}

double julian_century(double jd)
{
   return((jd-J2000JD)/36525.0);
}

double get_dt(double year)
{
   int i;
   double *dt;

   dt = ephemeris_correction;
   for (i=0;year>dt[i];i+=2);
   
   return(dt[i-1]+(dt[i+1]-dt[i-1])/(dt[i]-dt[i-2])*(year-dt[i-2]));
}

long gregorian_date_to_jd(long d,long m,long y)
{
   long jd;

   jd = (1461L*(y+4800L+(m-14L)/12L))/4L+(367L*(m-2L-12L*((m-14L)/12L)))/12L-
        (3L*((y+4900L+(m-14L)/12L)/100L))/4L+d-32075L;
   return(jd);
}

long julian_date_to_jd(long d,long m,long y)
{
   long jd;

   jd = 367L*y-(7L*(y+5001L+(m-9L)/7L))/4L+(275L*m)/9L+d+1729777L;

   return(jd);
}

long date_to_jd(long d,long m,long y)
{
   if (date_order(d,m,y) >= date_order(15,10,1582))
      return(gregorian_date_to_jd(d,m,y));
   else
      return(julian_date_to_jd(d,m,y));
}

void jd_to_gregorian_date(long jd,long *d,long *m,long *y)
{
   long l,n,i,j;

   l=jd+68569L;
   n=(4L*l)/146097L;
   l=l-(146097L*n+3L)/4L;
   i=(4000L*(l+1L))/1461001L;
   l=l-(1461L*i)/4L+31L;
   j=(80L*l)/2447L;
   *d=l-(2447L*j)/80L;
   l=j/11L;
   *m=j+2L-12L*l;
   *y=100L*(n-49L)+i+l;
}

void jd_to_julian_date(long jd,long *d,long *m,long *y)
{
   long l,n,i,j,k;

   j=jd+1402L;
   k=(j-1L)/1461L;
   l=j-1461L*k;
   n=(l-1L)/365L-l/1461L;
   i=l-365L*n+30L;
   j=(80L*i)/2447L;
   *d=i-(2447L*j)/80L;
   i=j/11L;
   *m=j+2L-12L*i;
   *y=4L*k+n+i-4716L;
}

void jd_to_date(long jd,long *d,long *m,long *y)
{
   if (jd >= 2299161L)
      jd_to_gregorian_date(jd,d,m,y);
   else
      jd_to_julian_date(jd,d,m,y);
}

momentype make_moment(double day,double month,double year,
                      double hour,double min,double sec)
{
   momentype w;

   w.date.day=day, w.date.month=month, w.date.year=year;

   w.date.day += glue_angle(' ',hour,min,sec)/24;
   
   while (w.date.month < 1.0) w.date.month+=12.0,w.date.year--;
   while (w.date.month >= 13.0) w.date.month-=12.0,w.date.year++;
   year=floor(w.date.year);
   w.date.day+=(w.date.year-year)*year_length((int)year);
   w.date.year=year;
   prep_feb((int)w.date.year);
   month=floor(w.date.month);
   w.date.day+=(w.date.month-month)*monthlen[(int)month];
   w.date.month=month;
   
   while (w.date.day < 1.0)
   {
      w.date.month--;
      if (w.date.month < 1.0)
      {
         w.date.month=12.0,w.date.year--;
         prep_feb((int)w.date.year);
      }
      w.date.day += (double)monthlen[(int)w.date.month];
   }

   while (w.date.day >= (double)monthlen[(int)w.date.month]+1.0)
   {
      w.date.day -= (double)monthlen[(int)w.date.month];
      w.date.month++;
      if (w.date.month > 12.0)
      {
         w.date.month=1.0,w.date.year++;
         prep_feb(w.date.year);
      }
   }

   day=floor(w.date.day);
   w.time = form_angle((w.date.day-day)*TWOPI);
   w.date.day=day;

   w.date.jd = (double)date_to_jd((long)w.date.day,(long)w.date.month,
                       (long)w.date.year)-0.5;
   w.JD = w.date.jd+w.time.alpha.frac/24;
   w.T = julian_century(w.JD);
   w.yr = w.date.year + 
         (w.JD - date_to_jd(1,1,w.date.year))/year_length((int)w.date.year);

   return(w);
}

void get_precarg(double precarg[3],double t)
{
   int arg,deg;

   for (arg=ZETA; arg<=ZEE; arg++)
   {
      for (deg=2,precarg[arg]=0.0; deg>=0; deg--)
         precarg[arg] = (precarg[arg]+c_prec[arg][deg])*t;
      precarg[arg] *= arc_sec;
   }
}

void prerot(double precarg[3],mat prec)
{
   mat rzeta,rtheta,rzee,c;

   rotate(ZAXIS,-precarg[ZETA],rzeta);
   rotate(YAXIS,precarg[THETA],rtheta);
   rotate(ZAXIS,-precarg[ZEE],rzee);
   matrix_product(rtheta,rzeta,c);
   matrix_product(rzee,c,prec);
}

double obliquity(double t)
{
   int deg;
   double e;

   for (deg=3,e=0.0; deg>=0; deg--) e = e*t+ec[deg];
   e *= arc_sec;

   return(e);
}

double fundamental_argument(int arg,double t)
{
   double r,argval;
   int deg;

   r = range(frq[arg][4]*t,1)*1296000.0;
   for (deg=3,argval=0.0; deg>=0; deg--) argval = argval*t+frq[arg][deg];
   argval = range(argval+r,1296000.0);
   return(argval*arc_sec);
}

double nutarg(int row,double t)
{
   double r,argval,sumarg;
   int arg,deg;

   for (arg=arg_l,sumarg=0.0; arg<=arg_Omega; arg++)
   {
      r = range(frq[arg][4]*t,1)*1296000.0;
      for (deg=3,argval=0.0; deg>=0; deg--)
         argval = argval*t+frq[arg][deg];
      argval = range(argval,1296000.0)+r;
      sumarg += argval*nut[row][arg];
   }
   sumarg = range(sumarg,1296000.0)*arc_sec;
   return(sumarg);
}

void nutate(double *delpsi,double *deleps,double t)
{
   double argval;
   int row;

   for (row=0,*delpsi=*deleps=0.0; row<106; row++)
   {
      argval = nutarg(row,t);
      *delpsi += (nut[row][5] + nut[row][6]*t)*sin(argval);
      *deleps += (nut[row][7] + nut[row][8]*t)*cos(argval);
   }
   *delpsi *= arc_sec/10000;
   *deleps *= arc_sec/10000;
}

void nutrot(mat nut,double epsilon,double delpsi,double deleps)
{
   mat reps,rdpsi,rdeps,c;

   rotate(XAXIS,epsilon,reps);
   rotate(ZAXIS,-delpsi,rdpsi);
   rotate(XAXIS,-epsilon-deleps,rdeps);
   matrix_product(rdpsi,reps,c);
   matrix_product(rdeps,c,nut);
}

void compute_precnut(epochtype *w)
{
   double precarg[3];

   get_precarg(precarg,w->tdt.T);
   prerot(precarg,w->prec);
   w->epsilon = obliquity(w->tdt.T);
   nutate(&w->delpsi,&w->deleps,w->tdt.T);
   nutrot(w->nut,w->epsilon,w->delpsi,w->deleps);
   matrix_product(w->nut,w->prec,w->precnut);
   transpose_matrix(w->prec,w->invprec);
   inverse_matrix(w->nut,w->invnut);
   matrix_product(w->invprec,w->invnut,w->invprecnut);
}

void compute_gst(epochtype *w)
{
   int deg;
   double gmst,gast,t;

   t = julian_century(w->ut.date.jd);

   for (deg=3,gmst=0.0;deg>=0;deg--) gmst = gmst*t + st[deg];
   gmst = range(gmst,86400)*arc_tsec;
   gmst = range(gmst+w->ut.time.alpha.rad*sider,TWOPI);
   gast = range(gmst + w->delpsi*cos(w->epsilon+w->deleps),TWOPI);
   w->gmst = form_angle(gmst);
   w->gast = form_angle(gast);
}

void compute_lst(topotype *t,epochtype *w)
{
   double lmst,last;

   lmst = range(w->gmst.alpha.rad+t->longitude.alpha.rad,TWOPI);
   last = range(w->gast.alpha.rad+t->longitude.alpha.rad,TWOPI);
   t->lmst = form_angle(lmst);
   t->last = form_angle(last);
}

void galactic_matrices(void)
{
   mat rx,rz,t;

   rotate(ZAXIS,HALFPI+getrad(RA_GAL),rz);
   rotate(XAXIS,HALFPI-getrad(DEC_GAL),rx);
   matrix_product(rx,rz,t);
   rotate(ZAXIS,-getrad(NODE_GAL),rz);
   matrix_product(rz,t,EQ2000_to_G);

   inverse_matrix(EQ2000_to_G,G_to_EQ2000);
}

void compute_galactic(epochtype *w)
{
   matrix_product(EQ2000_to_G,w->invprecnut,w->true_to_gal);
   matrix_product(EQ2000_to_G,w->invprec,w->mean_to_gal);

   matrix_product(w->precnut,G_to_EQ2000,w->gal_to_true);
   matrix_product(w->prec,G_to_EQ2000,w->gal_to_mean);
}

epochtype set_ut(double day,double month,double year,
                 double hour,double min,double sec)
{
   epochtype w;

   w.ut = make_moment(day,month,year,hour,min,sec);
   w.tdt = make_moment(day,month,year,hour,min,sec+get_dt(w.ut.yr));

   compute_precnut(&w);
   compute_galactic(&w);
   compute_gst(&w);

   return(w);
}

epochtype set_tdt(double day,double month,double year,
                  double hour,double min,double sec)
{
   epochtype w;

   w.tdt = make_moment(day,month,year,hour,min,sec);
   w.ut = make_moment(day,month,year,hour,min,sec-get_dt(w.tdt.yr));

   compute_precnut(&w);
   compute_galactic(&w);
   compute_gst(&w);

   return(w);
}

epochtype set_jd(double jd)
{
   double frac,h,m,s;
   long jday,day,month,year;
   char z;

   jday = (long)floor(jd+0.5);
   frac = jd-(double)jday+0.5;
   jd_to_date(jday,&day,&month,&year);
   break_angle(frac*24,&z,&h,&m,&s);
   return(set_tdt((double)day,(double)month,(double)year,h,m,s));
}

epochtype julian_epoch(double year)
{
   double jd;

   jd = (year-2000.0)*julian_year + J2000JD;
   return(set_jd(jd));
}

double jd_to_julian_year(double jd)
{
   double year;

   year = 2000.0 + (jd - J2000JD) / julian_year;
   return(year);
}

topotype set_topocentre(char l_sgn,double l_hr,double l_min,double l_sec,
                        char f_sgn,double f_dg,double f_min,double f_sec,
                        double height,epochtype *w)
{
   topotype t;
   double p,g,q,u,c,s;
   vec omega;

   t.longitude = form_lambda(l_sgn,l_hr,l_min,l_sec);
   t.latitude = form_delta(f_sgn,f_dg,f_min,f_sec);
   t.height = height;

   compute_lst(&t,w);

   set_vector(omega,0,0,angular_velocity);

   p = cos(t.latitude.delta.rad);
   g = sin(t.latitude.delta.rad);
   u = 1 - 1/inverse_flattening;
   q = u*g;
   c = 1/sqrt(p*p+q*q);
   s = u*u*c;

   t.r[0] = (equatorial_radius*c+t.height/1000)*p*cos(t.last.alpha.rad);
   t.r[1] = (equatorial_radius*c+t.height/1000)*p*sin(t.last.alpha.rad);
   t.r[2] = (equatorial_radius*s+t.height/1000)*g;

   vector_product(omega,t.r,t.v);

   return(t);
}

void astro_initialize(void)
{
   galactic_matrices();
}

void equatorial_matrices(startype *s)
{
   if (s->position == MEAN_POSITION)
   {
      copy_matrix(s->equinox.mean_to_gal,s->eq_to_gal);
      copy_matrix(s->equinox.gal_to_mean,s->gal_to_eq);
   }
   else
   {
      copy_matrix(s->equinox.true_to_gal,s->eq_to_gal);
      copy_matrix(s->equinox.gal_to_true,s->gal_to_eq);
   }
}

void local_matrices(startype *s)
{
   mat p,r;

   rotate(XAXIS,-(HALFPI-s->dec.delta.rad),p);
   rotate(ZAXIS,-(HALFPI+s->ra.alpha.rad),r);
   matrix_product(r,p,s->loc_to_eq);

   rotate(XAXIS,HALFPI-s->dec.delta.rad,p);
   rotate(ZAXIS,HALFPI+s->ra.alpha.rad,r);
   matrix_product(p,r,s->eq_to_loc);

   matrix_product(s->eq_to_gal,s->loc_to_eq,s->loc_to_gal);
   matrix_product(s->eq_to_loc,s->gal_to_eq,s->gal_to_loc);
}

void fix_star(startype *s)
{
   double kp,dist,l,b;
   vec r,v;

   if (s->par <= 0) s->par=0.001;
   dist = 1/s->par;
   kp = vt_factor/s->par;
   s->absmag = s->vmag+5-5*log10(dist);

   equatorial_matrices(s);
   local_matrices(s);

   sphere_to_xyz(s->ra.alpha.rad,s->dec.delta.rad,dist,r);
   matrix_times_vector(s->eq_to_gal,r,s->r);
   xyz_to_sphere(s->r,&l,&b,&dist);
   s->l = form_angle(l);
   s->b = form_angle(b);
   set_vector(v,kp*s->mura,kp*s->mudec,s->rv);
   matrix_times_vector(s->loc_to_gal,v,s->v);
}

void unfix_star(startype *s)
{
   double kp,dist,ra,dec;
   vec r,v;

   dist = vector_length(s->r);
   s->par = 1/dist;
   kp = vt_factor/s->par;
   s->vmag = s->absmag-5+5*log10(dist);

   equatorial_matrices(s);

   matrix_times_vector(s->gal_to_eq,s->r,r);
   xyz_to_sphere(r,&ra,&dec,&dist);
   s->ra = form_angle(ra);
   s->dec = form_angle(dec);

   local_matrices(s);

   matrix_times_vector(s->gal_to_loc,s->v,v);
   s->mura = v[0]/kp;
   s->mudec = v[1]/kp;
   s->rv = v[2];
}

startype new_position(startype s,epochtype equinox,int pos,epochtype epoch)
{
   startype z;
   double dt,l,b,dist;
   vec dr;

   z = s;
   fix_star(&z);

   z.equinox = equinox;
   z.epoch = epoch;
   z.position = pos;

   dt = (z.epoch.tdt.JD-s.epoch.tdt.JD)*86400/parsec_km;
   vector_times_scalar(z.v,dt,dr);
   vector_sum(z.r,dr,z.r);
   xyz_to_sphere(z.r,&l,&b,&dist);
   z.l = form_angle(l);
   z.b = form_angle(b);

   unfix_star(&z);

   return(z);
}

void load_earth(void)
{
   file_type f;
   long len;
   ldiv_t q;

   open_file(&f,"r","./astro/earth.bin");

   get_fseek(&f,0,SEEK_END);
   len = ftell(f.dat);
   if (len < 1040) get_error("File 'earth.bin' too short!");

   q = ldiv(len,104L);
   if (q.rem != 0L)
      get_error("An integer number of records expected in 'earth.bin'!");

   get_fseek(&f,0,SEEK_SET);

   earth = (double *)get_space(len);

   get_fread((void *)earth,104,q.quot,&f);

   get_fclose(&f);

   num_rec = q.quot;
   first_jd = earth[0];
   last_jd = earth[(num_rec-1)*13];
}

void free_earth(void)
{
   get_free((void *)earth);
}

void interpolate_earth(double jd,vec rh,vec vh,vec rb,vec vb)
{
   int icen,i,k,p;
   double x[2*INTERPOLATION_RADIUS+1],y[2*INTERPOLATION_RADIUS+1];
   double *source;
   double val[12];

   icen = nint(jd-first_jd);

   if (icen < INTERPOLATION_RADIUS || icen+INTERPOLATION_RADIUS >= num_rec)
      get_error("JD out of limits!");

   for (i=-INTERPOLATION_RADIUS,k=0,
        source=earth+(icen-INTERPOLATION_RADIUS)*13;
        i<=INTERPOLATION_RADIUS;i++,k++,source+=13) x[k] = *source;

   for (p=0;p<12;p++)
   {
      for (i=-INTERPOLATION_RADIUS,k=0,
           source=earth+(icen-INTERPOLATION_RADIUS)*13+p+1;
           i<=INTERPOLATION_RADIUS;i++,k++,source+=13) y[k] = *source;
      val[p] =  interpolation(x,y,2*INTERPOLATION_RADIUS+1,jd);
   }

   copy_vector(val,rh);
   copy_vector(val+3,vh);
   copy_vector(val+6,rb);
   copy_vector(val+9,vb);
}

int hd_to_hip(int hd)
{
   int hip;
   file_type file;

   if (hd < 1 || hd >= HDLIMIT) return(0);

   open_file(&file,"r","./astro/hdhip.bin");

   get_fseek(&file,hd*sizeof(int),SEEK_SET);
   get_fread((void *)&hip,sizeof(int),1,&file);
   get_fclose(&file);

   return(hip);
}

int hr_to_hd(int hr)
{
   int hd;
   file_type file;

   if (hr < 1 || hr >= HRLIMIT) return(0);

   open_file(&file,"r","./astro/hrhd.bin");

   get_fseek(&file,hr*sizeof(int),SEEK_SET);
   get_fread((void *)&hd,sizeof(int),1,&file);
   get_fclose(&file);

   return(hd);
}

int find_bright_star(char *star,char *mask)
{
   int hr;
   file_type file;
   char *bsc,*ptr;
   int k,found;

   bsc = get_space(10*HRLIMIT);

   open_file(&file,"r","./astro/hrstar.bin");
   get_fread(bsc,10,HRLIMIT,&file);
   get_fclose(&file);

   for (hr=1,ptr=bsc+10;hr<HRLIMIT;hr++,ptr+=10)
   {
      for (k=0,found=1;k<10;k++)
      {
         if (mask[k] == '0') continue;
         if ((star[k] | ' ') != (ptr[k] | ' ')) found=0,k=9;
      }
      if (found) break;
   }

   get_free(bsc);

   if (found)
      return(hr);
   else
      return(0);
}

int get_hipparcos_number(char *starname)
{
   char *ptr,*eptr;
   int type,flm,cmp,cat,hip,hd,hr;
   char byr[4],con[4],star[11],mask[11];

   int number_of_digits(char *s)
   {
      int k;

      for (k=0;digit_ok(s[k]);k++);

      return(k);
   }

   int number_of_letters(char *s)
   {
      int k;

      for (k=0;alpha_ok(s[k]);k++);

      return(k);
   }

   char *skip_star_name_gap(char *p)
   {
      if (*p == 0) return(p);

      if (*p == ' ')
      {
         while (*p == ' ') p++;
         return(p);
      }

      if (*p=='-' || *p=='_')
         return(p+1);
      else
         return(p);
   }

   void collect_constellation(void)
   {
      int n;

      if (!alpha_ok(*ptr)) get_error("Bad star name (%s)!",starname);
      n = number_of_letters(ptr);
      if (n < 3) get_error("Bad star name (%s)!",starname);
      memcpy(con,ptr,3);
      con[3] = 0;
      ptr += 3;
      for (n=0;*constell[n]!='#' && strcasecmp(con,constell[n]);n++);
      if (strcasecmp(con,constell[n]))
         get_error("Bad star name (%s)!",starname);
   }

   void collect_bayer(void)
   {
      int n;

      n = number_of_letters(starname);
      if (n != 3 && n != 2) get_error("Bad star name (%s)!",starname);
      memcpy(byr,starname,n);
      byr[n] = 0;
      ptr = starname+n;
      for (n=0;*greek[n]!='#' && strcasecmp(byr,greek[n]);n++);
      if (strcasecmp(byr,greek[n])) get_error("Bad star name (%s)!",starname);
   }

   int catalog_ok(void)
   {
      if (strncasecmp(starname,"HD",2) == 0)
      {
         type = HD_TYPE;
         ptr = starname+2;
         return(1);
      }

      if (strncasecmp(starname,"HR",2) == 0)
      {
         type = HR_TYPE;
         ptr = starname+2;
         return(1);
      }

       if (strncasecmp(starname,"HIP",3) == 0)
      {
         type = HIP_TYPE;
         ptr = starname+3;
         return(1);
      }

     return(0);
   }

   int end_of_text(char *p)
   {
      while (*p == ' ') p++;
      if (*p == 0)
         return(1);
      else
         return(0);
   }

   type = 0;

   if (!alphanum_ok(*starname)) get_error("Bad star name (%s)!",starname);

   if (digit_ok(*starname))
   {
      type = FLM_TYPE;
      if (*starname == '0') get_error("Bad star name (%s)!",starname);
      if (number_of_digits(starname) > 4)
         get_error("Bad star name (%s)!",starname);
      flm = strtol(starname,&eptr,10);
      ptr = skip_star_name_gap(eptr);
      collect_constellation();
      if (!end_of_text(ptr)) get_error("Bad star name (%s)!",starname);
   }
   else
   {
      if (catalog_ok())
      {
         ptr = skip_star_name_gap(ptr);
         if (*ptr == '0')
            get_error("Leading zero not allowed in catalogue number (%s)!",
                       starname);
         if (!digit_ok(*ptr)) get_error("Bad star name (%s)!",starname);
         if (*ptr == '0') get_error("Bad star name (%s)!",starname);
         if (number_of_digits(ptr) > 6)
            get_error("Bad star name (%s)!",starname);
         cat = strtol(ptr,&eptr,10);
         if (!end_of_text(eptr)) get_error("Bad star name (%s)!",starname);
      }
      else
      {
         type = BYR_TYPE;
         collect_bayer();
         if (digit_ok(*ptr))
         {
            cmp = (int)(*ptr - '0');
            if (cmp < 1) get_error("Bad star name (%s)!",starname);
            ptr++;
         }
         else
            cmp = 0;
         ptr = skip_star_name_gap(ptr);
         collect_constellation();
         if (!end_of_text(ptr)) get_error("Bad star name (%s)!",starname);
      }
   }

   switch(type)
   {
      case HIP_TYPE: hip = cat;
                     break;
      case HD_TYPE:  hd = cat;
                     hip = hd_to_hip(hd);
                     break;
      case HR_TYPE:  hr = cat;
                     hd = hr_to_hd(hr);
                     hip = hd_to_hip(hd);
                     break;
      case FLM_TYPE: sprintf(star,"%3d    %s",flm,con);
                     sprintf(mask,"1110000111");
                     hr = find_bright_star(star,mask);
                     hd = hr_to_hd(hr);
                     hip = hd_to_hip(hd);
                     break;
      case BYR_TYPE: if (cmp == 0)
                        sprintf(star,"   %-3s %s",byr,con);
                     else
                        sprintf(star,"   %-3s%d%s",byr,cmp,con);
                     sprintf(mask,"0001111111");
                     hr = find_bright_star(star,mask);
                     hd = hr_to_hd(hr);
                     hip = hd_to_hip(hd);
                     break;
      default: get_error("Bad star name (%s)!",starname);
   }

   return(hip);
}

double rounded(float x)
{
   return(floor((double)x*100+0.5)/100);
}

startype find_star(int hip)

{
   char starbuf[STARLEN];
   startype star;
   file_type f;
   int seq;
   int *iptr;
   short *sptr;
   float *fptr;
   double *dptr;

   if (hip < 1 || hip >= HIPLIMIT) get_error("HIP number out of limits!");

   open_file(&f,"r","./astro/index.bin");
   get_fseek(&f,hip*sizeof(int),SEEK_SET);
   get_fread((void *)&seq,sizeof(int),1,&f);
   get_fclose(&f);

   if (seq < 0 || seq >= HIPTOT) get_error("Sequence number out of limits!");
   if (seq == 0) get_error("Star not found!");

   memset(star.id,0,11);
   memset(star.sp,0,13);

   iptr = (int *)(&starbuf[0]);
   sptr = (short *)(&starbuf[8]);
   dptr = (double *)(&starbuf[20]);
   fptr = (float *)(&starbuf[36]);

   open_file(&f,"r","./astro/stars.bin");
   get_fseek(&f,(seq-1)*STARLEN,SEEK_SET);
   get_fread(starbuf,1,STARLEN,&f);
   get_fclose(&f);

   if (hip != iptr[0]) get_error("HIP number mismatch!");

   star.hd = iptr[1];
   star.hr = (int)sptr[0];
   memcpy(star.id,starbuf+10,10);
   star.ra = form_angle(getrad(dptr[0]));
   star.dec = form_angle(getrad(dptr[1]));
   
   star.mura = rounded(fptr[0])/1000;
   star.mudec = rounded(fptr[1])/1000;
   star.par = rounded(fptr[2])/1000;
   star.rv = rounded(fptr[3]);
   star.vmag = rounded(fptr[4]);
   memcpy(star.sp,starbuf+56,12);

   star.equinox = julian_epoch(2000.0);
   star.epoch = julian_epoch(1991.25);

   return(star);
}

int get_object_info(char *starname,int *hipnum,int *usrnum,startype *star)
{
   file_type catfile;
   char row[256];
   int ndat,i;
   char label[256],dectxt[256],z;
   double eq,ep,ah,am,as,dd,dm,ds,mura,mudec,par,rv;

   *hipnum = *usrnum = 0;

   if (file_exists("starcat.dat"))
   {
      open_file(&catfile,"r","starcat.dat");
      for (i=1;fgets(row,250,catfile.dat) != NULL;i++)
      {
         if (strlen(row) > 200)
            get_error("Data line too long in `%s'!",catfile.path);
         ndat = sscanf(row,"%s %lf %lf %lf %lf %lf %s %lf %lf "
                           "%lf %lf %lf %lf",
                       label,&eq,&ep,&ah,&am,&as,dectxt,&dm,&ds,
                             &mura,&mudec,&par,&rv);
         if (ndat != 13)
            get_error("Missing data in `%s'!",catfile.path);
         if (strcasecmp(starname,label) != 0) continue;
         *usrnum = i;
         z = *dectxt;
         if (z!='-' && z!='+')
            get_error("Bad declination sign in `%s'!",catfile.path);
         dd = fabs(atof(dectxt));
         star->equinox = julian_epoch(eq);
         star->epoch = julian_epoch(ep);
         star->ra = form_alpha(ah,am,as);
         star->dec = form_delta(z,dd,dm,ds);
         star->mura = mura/1000;
         star->mudec = mudec/1000;
         star->par = par/1000;
         star->rv = rv;
         break;
      }
      get_fclose(&catfile);
   }

   if (*usrnum < 1)
   {
      *hipnum = get_hipparcos_number(starname);
      if (*hipnum > 0) *star = find_star(*hipnum);
   }
   
   return(1);
}

void get_barycentric_correction(epochtype epoch,startype star1,topotype topo,
      double *drv,double *djd)
{
   startype star2;
   vec rs,rhel_2000,vhel_2000,rbar_2000,vbar_2000;
   vec rbar,vbar;

   star2 = new_position(star1,epoch,TRUE_POSITION,epoch);

   sphere_to_xyz(star2.ra.alpha.rad,star2.dec.delta.rad,1.0,rs);

   interpolate_earth(epoch.tdt.JD,rhel_2000,vhel_2000,rbar_2000,vbar_2000);

   matrix_times_vector(epoch.precnut,rbar_2000,rbar);
   matrix_times_vector(epoch.precnut,vbar_2000,vbar);   

   vector_times_scalar(vbar,unit_velocity,vbar);

   *drv = scalar_product(vbar,rs) + scalar_product(topo.v,rs);

   *djd = scalar_product(rbar,rs)*LIGHT_TIME/86400.0;
}

void get_sun_correction(epochtype epoch,topotype topo,double *drv,double *djd)
{
   vec rhel_2000,vhel_2000,rbar_2000,vbar_2000;
   vec rhel,vhel,ahel,ehel,rho;

   interpolate_earth(epoch.tdt.JD,rhel_2000,vhel_2000,rbar_2000,vbar_2000);

   matrix_times_vector(epoch.precnut,rhel_2000,rhel);
   matrix_times_vector(epoch.precnut,vhel_2000,vhel);

   vector_times_scalar(vhel,unit_velocity,vhel);

   vector_times_scalar(topo.r,1.0/unit_distance,rho);
   vector_sum(rhel,rho,ahel);

   unit_vector(ahel,ehel);

   *drv = -(scalar_product(vhel,ehel)+scalar_product(topo.v,ehel));

   *djd = -scalar_product(ahel,ehel)*LIGHT_TIME/86400.0;
}
