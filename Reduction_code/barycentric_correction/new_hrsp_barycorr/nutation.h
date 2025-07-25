#ifndef NUTATION_DEFINED

  #define NUTATION_DEFINED

  #define LS_COUNT 678
  #define PL_COUNT 687

  #define D2PI 6.283185307179586476925287
  #define DAS2R 4.848136811095359935899141e-6
  #define DC2R (DAS2R / 1e7)
  #define TURNAS  1296000.0

  #define PL_MERCURY  5
  #define PL_VENUS    6
  #define PL_EARTH    7
  #define PL_ACCUM   13

#endif

double get_fund_arg(double t,int j);
double get_planet_arg(double t,int j);
void nut00a(double t,double *dpsi,double *deps);
void nut06a(double t,double *dpsi,double *deps);
