#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#include "general.h"
#include "angle.h"
#include "datetime.h"
#include "vector.h"
#include "nutation.h"
#include "astro.h"

#define BRIGHT_COUNT     9110
#define CCDM_COUNT     105838
#define BARBIER_COUNT   36145
#define HIP_MAIN_COUNT 118218
#define HIP_COM_COUNT   24588

#define BRIGHT_RECORD   32
#define CCDM_RECORD     24
#define BARBIER_RECORD  32
#define HIP_RECORD      80

#define MAX_STAR_NAME 80

#define YEAR_SEC     31557600.0
#define AU_KM       149597870.691
#define PC_AU       R2ASEC

#define LIGHT_TIME  499.004783806135643345637
#define PC_KM       30856775813057.2895329155
#define UNIT_VEL    1731.45683670138888888889
#define TRANS_VEL   4.74047046324815575328922
#define DEFLECT_CON 9.87062871570064612827663e-9

#define EARTH_EQU_RAD   6378.137
#define EARTH_INV_FLAT   298.257223563

#define TAI_OFFSET 32.184 /* Time difference: TDT - TAI */

static constellationtype constell[] =
 {{ "Andromeda",           "Andromedae",          "And" },
  { "Antlia",              "Antliae",             "Ant" },
  { "Apus",                "Apodis",              "Aps" },
  { "Aquarius",            "Aquarii",             "Aqr" },
  { "Aquila",              "Aquilae",             "Aql" },
  { "Ara",                 "Arae",                "Ara" },
  { "Aries",               "Arietis",             "Ari" },
  { "Auriga",              "Aurigae",             "Aur" },
  { "Bootes",              "Bootis",              "Boo" },
  { "Caelum",              "Caeli",               "Cae" },
  { "Camelopardalis",      "Camelopardalis",      "Cam" },
  { "Cancer",              "Cancri",              "Cnc" },
  { "Canes Venatici",      "Canum Venaticorum",   "CVn" },
  { "Canis Major",         "Canis Majoris",       "CMa" },
  { "Canis Minor",         "Canis Minoris",       "CMi" },
  { "Capricornus",         "Capricorni",          "Cap" },
  { "Carina",              "Carinae",             "Car" },
  { "Cassiopeia",          "Cassiopeiae",         "Cas" },
  { "Centaurus",           "Centauri",            "Cen" },
  { "Cepheus",             "Cephei",              "Cep" },
  { "Cetus",               "Ceti",                "Cet" },
  { "Chamaeleon",          "Chamaeleontis",       "Cha" },
  { "Circinus",            "Circini",             "Cir" },
  { "Columba",             "Columbae",            "Col" },
  { "Coma Berenices",      "Comae Berenices",     "Com" },
  { "Corona Australis",    "Coronae Australis",   "CrA" },
  { "Corona Borealis",     "Coronae Borealis",    "CrB" },
  { "Corvus",              "Corvi",               "Crv" },
  { "Crater",              "Crateris",            "Crt" },
  { "Crux",                "Crucis",              "Cru" },
  { "Cygnus",              "Cygni",               "Cyg" },
  { "Delphinus",           "Delphini",            "Del" },
  { "Dorado",              "Doradus",             "Dor" },
  { "Draco",               "Draconis",            "Dra" },
  { "Equuleus",            "Equulei",             "Equ" },
  { "Eridanus",            "Eridani",             "Eri" },
  { "Fornax",              "Fornacis",            "For" },
  { "Gemini",              "Geminorum",           "Gem" },
  { "Grus",                "Gruis",               "Gru" },
  { "Hercules",            "Herculis",            "Her" },
  { "Horologium",          "Horologii",           "Hor" },
  { "Hydra",               "Hydrae",              "Hya" },
  { "Hydrus",              "Hydri",               "Hyi" },
  { "Indus",               "Indi",                "Ind" },
  { "Lacerta",             "Lacertae",            "Lac" },
  { "Leo",                 "Leonis",              "Leo" },
  { "Leo Minor",           "Leonis Minoris",      "LMi" },
  { "Lepus",               "Leporis",             "Lep" },
  { "Libra",               "Librae",              "Lib" },
  { "Lupus",               "Lupi",                "Lup" },
  { "Lynx",                "Lyncis",              "Lyn" },
  { "Lyra",                "Lyrae",               "Lyr" },
  { "Mensa",               "Mensae",              "Men" },
  { "Microscopium",        "Microscopii",         "Mic" },
  { "Monoceros",           "Monocerotis",         "Mon" },
  { "Musca",               "Muscae",              "Mus" },
  { "Norma",               "Normae",              "Nor" },
  { "Octans",              "Octantis",            "Oct" },
  { "Ophiuchus",           "Ophiuchi",            "Oph" },
  { "Orion",               "Orionis",             "Ori" },
  { "Pavo",                "Pavonis",             "Pav" },
  { "Pegasus",             "Pegasi",              "Peg" },
  { "Perseus",             "Persei",              "Per" },
  { "Phoenix",             "Phoenicis",           "Phe" },
  { "Pictor",              "Pictoris",            "Pic" },
  { "Pisces",              "Piscium",             "Psc" },
  { "Piscis Austrinus",    "Piscis Austrini",     "PsA" },
  { "Puppis",              "Puppis",              "Pup" },
  { "Pyxis",               "Pyxidis",             "Pyx" },
  { "Reticulum",           "Reticuli",            "Ret" },
  { "Sagitta",             "Sagittae",            "Sge" },
  { "Sagittarius",         "Sagittarii",          "Sgr" },
  { "Scorpius",            "Scorpii",             "Sco" },
  { "Sculptor",            "Sculptoris",          "Scl" },
  { "Scutum",              "Scuti",               "Sct" },
  { "Serpens",             "Serpentis",           "Ser" },
  { "Sextans",             "Sextantis",           "Sex" },
  { "Taurus",              "Tauri",               "Tau" },
  { "Telescopium",         "Telescopii",          "Tel" },
  { "Triangulum",          "Trianguli",           "Tri" },
  { "Triangulum Australe", "Trianguli Australis", "TrA" },
  { "Tucana",              "Tucanae",             "Tuc" },
  { "Ursa Major",          "Ursae Majoris",       "UMa" },
  { "Ursa Minor",          "Ursae Minoris",       "UMi" },
  { "Vela",                "Velorum",             "Vel" },
  { "Virgo",               "Virginis",            "Vir" },
  { "Volans",              "Volantis",            "Vol" },
  { "Vulpecula",           "Vulpeculae",          "Vul" },
  { "",                    "",                    ""    }};

static greektype greek[] =
 {{ "Alpha",   "Alp" },
  { "Beta",    "Bet" },
  { "Gamma",   "Gam" },
  { "Delta",   "Del" },
  { "Epsilon", "Eps" },
  { "Zeta",    "Zet" },
  { "Eta",     "Eta" },
  { "Theta",   "The" },
  { "Iota",    "Iot" },
  { "Kappa",   "Kap" },
  { "Lambda",  "Lam" },
  { "Mu",      "Mu"  },
  { "Nu",      "Nu"  },
  { "Xi",      "Xi"  },
  { "Omicron", "Omi" },
  { "Pi",      "Pi"  },
  { "Rho",     "Rho" },
  { "Sigma",   "Sig" },
  { "Tau",     "Tau" },
  { "Upsilon", "Ups" },
  { "Phi",     "Phi" },
  { "Chi",     "Chi" },
  { "Psi",     "Psi" },
  { "Omega",   "Ome" },
  { "",        ""    }};

  static char *bad_star[] = 
    { "",
      "parse_greek_letter: A letter expected",
      "parse_greek_letter: Bad Greek letter",
      "parse_constellation: A letter expected",
      "parse_constellation: Bad constellation name",
      "parse_flamsteed_number: Bad Flamsteed number",
      "parse_bayer_seq: Bad Bayer sequence number",
      "parse_catalog_number: Bad catalog number",
      "parse_multiple_component: Bad double/multiple component",
      "parse_multiple_component: End of string expected",
      "parse_base_name: A letter or a digit expected" };

static double fukushima_williams[4][6] =
 {{   -0.052928,   10.556378,  0.4932044, -0.00031238, -2.788e-6,   2.60e-8},
  {84381.412819,  -46.811016,  0.0511268,  0.00053289, -4.40e-7,   -1.76e-8},
  {   -0.041775, 5038.481484,  1.5584175, -0.00018522, -2.6452e-5, -1.48e-8},
  {84381.406,     -46.836769, -0.0001831,  0.00200340, -5.76e-7,   -4.34e-8}};

static int cio_ia[5] = {  0, 33, 36, 61, 65 };
static int cio_ib[5] = { 32, 35, 60, 64, 65 };

static int cio_ks[66][8] =
 {{  0,  0,  0,  0,  1,  0,  0,  0 },
  {  0,  0,  0,  0,  2,  0,  0,  0 },
  {  0,  0,  2, -2,  3,  0,  0,  0 },
  {  0,  0,  2, -2,  1,  0,  0,  0 },
  {  0,  0,  2, -2,  2,  0,  0,  0 },
  {  0,  0,  2,  0,  3,  0,  0,  0 },
  {  0,  0,  2,  0,  1,  0,  0,  0 },
  {  0,  0,  0,  0,  3,  0,  0,  0 },
  {  0,  1,  0,  0,  1,  0,  0,  0 },
  {  0,  1,  0,  0, -1,  0,  0,  0 },
  {  1,  0,  0,  0, -1,  0,  0,  0 },
  {  1,  0,  0,  0,  1,  0,  0,  0 },
  {  0,  1,  2, -2,  3,  0,  0,  0 },
  {  0,  1,  2, -2,  1,  0,  0,  0 },
  {  0,  0,  4, -4,  4,  0,  0,  0 },
  {  0,  0,  1, -1,  1, -8, 12,  0 },
  {  0,  0,  2,  0,  0,  0,  0,  0 },
  {  0,  0,  2,  0,  2,  0,  0,  0 },
  {  1,  0,  2,  0,  3,  0,  0,  0 },
  {  1,  0,  2,  0,  1,  0,  0,  0 },
  {  0,  0,  2, -2,  0,  0,  0,  0 },
  {  0,  1, -2,  2, -3,  0,  0,  0 },
  {  0,  1, -2,  2, -1,  0,  0,  0 },
  {  0,  0,  0,  0,  0,  8,-13, -1 },
  {  0,  0,  0,  2,  0,  0,  0,  0 },
  {  2,  0, -2,  0, -1,  0,  0,  0 },
  {  0,  1,  2, -2,  2,  0,  0,  0 },
  {  1,  0,  0, -2,  1,  0,  0,  0 },
  {  1,  0,  0, -2, -1,  0,  0,  0 },
  {  0,  0,  4, -2,  4,  0,  0,  0 },
  {  0,  0,  2, -2,  4,  0,  0,  0 },
  {  1,  0, -2,  0, -3,  0,  0,  0 },
  {  1,  0, -2,  0, -1,  0,  0,  0 },
  {  0,  0,  0,  0,  2,  0,  0,  0 },
  {  0,  0,  0,  0,  1,  0,  0,  0 },
  {  0,  0,  2, -2,  3,  0,  0,  0 },
  {  0,  0,  0,  0,  1,  0,  0,  0 },
  {  0,  0,  2, -2,  2,  0,  0,  0 },
  {  0,  0,  2,  0,  2,  0,  0,  0 },
  {  0,  0,  0,  0,  2,  0,  0,  0 },
  {  0,  1,  0,  0,  0,  0,  0,  0 },
  {  1,  0,  0,  0,  0,  0,  0,  0 },
  {  0,  1,  2, -2,  2,  0,  0,  0 },
  {  0,  0,  2,  0,  1,  0,  0,  0 },
  {  1,  0,  2,  0,  2,  0,  0,  0 },
  {  0,  1, -2,  2, -2,  0,  0,  0 },
  {  1,  0,  0, -2,  0,  0,  0,  0 },
  {  0,  0,  2, -2,  1,  0,  0,  0 },
  {  1,  0, -2,  0, -2,  0,  0,  0 },
  {  0,  0,  0,  2,  0,  0,  0,  0 },
  {  1,  0,  0,  0,  1,  0,  0,  0 },
  {  1,  0, -2, -2, -2,  0,  0,  0 },
  {  1,  0,  0,  0, -1,  0,  0,  0 },
  {  1,  0,  2,  0,  1,  0,  0,  0 },
  {  2,  0,  0, -2,  0,  0,  0,  0 },
  {  2,  0, -2,  0, -1,  0,  0,  0 },
  {  0,  0,  2,  2,  2,  0,  0,  0 },
  {  2,  0,  2,  0,  2,  0,  0,  0 },
  {  2,  0,  0,  0,  0,  0,  0,  0 },
  {  1,  0,  2, -2,  2,  0,  0,  0 },
  {  0,  0,  2,  0,  0,  0,  0,  0 },
  {  0,  0,  0,  0,  1,  0,  0,  0 },
  {  0,  0,  2, -2,  2,  0,  0,  0 },
  {  0,  0,  2,  0,  2,  0,  0,  0 },
  {  0,  0,  0,  0,  2,  0,  0,  0 },
  {  0,  0,  0,  0,  1,  0,  0,  0 }};

static double cio_ss[66][2] =
  {{  -2640.73e-6,      0.39e-6 },
  {    -63.53e-6,      0.02e-6 },
  {    -11.75e-6,     -0.01e-6 },
  {    -11.21e-6,     -0.01e-6 },
  {      4.57e-6,      0.00e-6 },
  {     -2.02e-6,      0.00e-6 },
  {     -1.98e-6,      0.00e-6 },
  {      1.72e-6,      0.00e-6 },
  {      1.41e-6,      0.01e-6 },
  {      1.26e-6,      0.01e-6 },
  {      0.63e-6,      0.00e-6 },
  {      0.63e-6,      0.00e-6 },
  {     -0.46e-6,      0.00e-6 },
  {     -0.45e-6,      0.00e-6 },
  {     -0.36e-6,      0.00e-6 },
  {      0.24e-6,      0.12e-6 },
  {     -0.32e-6,      0.00e-6 },
  {     -0.28e-6,      0.00e-6 },
  {     -0.27e-6,      0.00e-6 },
  {     -0.26e-6,      0.00e-6 },
  {      0.21e-6,      0.00e-6 },
  {     -0.19e-6,      0.00e-6 },
  {     -0.18e-6,      0.00e-6 },
  {      0.10e-6,     -0.05e-6 },
  {     -0.15e-6,      0.00e-6 },
  {      0.14e-6,      0.00e-6 },
  {      0.14e-6,      0.00e-6 },
  {     -0.14e-6,      0.00e-6 },
  {     -0.14e-6,      0.00e-6 },
  {     -0.13e-6,      0.00e-6 },
  {      0.11e-6,      0.00e-6 },
  {     -0.11e-6,      0.00e-6 },
  {     -0.11e-6,      0.00e-6 },
  {     -0.07e-6,      3.57e-6 },
  {      1.73e-6,     -0.03e-6 },
  {      0.00e-6,      0.48e-6 },
  {    743.52e-6,     -0.17e-6 },
  {     56.91e-6,      0.06e-6 },
  {      9.84e-6,     -0.01e-6 },
  {     -8.85e-6,      0.01e-6 },
  {     -6.38e-6,     -0.05e-6 },
  {     -3.07e-6,      0.00e-6 },
  {      2.23e-6,      0.00e-6 },
  {      1.67e-6,      0.00e-6 },
  {      1.30e-6,      0.00e-6 },
  {      0.93e-6,      0.00e-6 },
  {      0.68e-6,      0.00e-6 },
  {     -0.55e-6,      0.00e-6 },
  {      0.53e-6,      0.00e-6 },
  {     -0.27e-6,      0.00e-6 },
  {     -0.27e-6,      0.00e-6 },
  {     -0.26e-6,      0.00e-6 },
  {     -0.25e-6,      0.00e-6 },
  {      0.22e-6,      0.00e-6 },
  {     -0.21e-6,      0.00e-6 },
  {      0.20e-6,      0.00e-6 },
  {      0.17e-6,      0.00e-6 },
  {      0.13e-6,      0.00e-6 },
  {     -0.13e-6,      0.00e-6 },
  {     -0.12e-6,      0.00e-6 },
  {     -0.11e-6,      0.00e-6 },
  {      0.30e-6,    -23.42e-6 },
  {     -0.03e-6,     -1.46e-6 },
  {     -0.01e-6,     -0.25e-6 },
  {      0.00e-6,      0.23e-6 },
  {     -0.26e-6,     -0.01e-6 }};

static double cio_sp[6] =
  {94.00e-6, 3808.65e-6, -122.68e-6, -72574.11e-6, 27.98e-6, 15.62e-6};

static double era_coeff[2] = { 0.7790572732640,  0.00273781191135448};

static double gmstp_coeff[6] =
  {0.014506, 4612.156534, 1.3915817, -4.4e-7, -2.9956e-5, -3.68e-8};

static byte *jpl = NULL;
static jpl_header jplh;

static leap_type *leap_list = NULL;
static int leap_count = 0;

static deltat_type *deltat_list = NULL;
static int deltat_count = 0;

static brighttype *bright = NULL;
static ccdmtype *ccdm = NULL;
static barbiertype *barbier = NULL;
static hiptype *hip_main = NULL;
static hiptype *hip_com = NULL;
static char *comp_buff = NULL;
static char **comp_list = NULL;
static binary_type *binary = NULL;
static int binary_count = 0;

/* ----------------------------- Coordinates  ----------------------------- */

double get_fw_arg(int argseq,double t)
{
  int i;
  double arg;

  for (i=5,arg=0.0;i>=0;i--)
    arg = arg * t + fukushima_williams[argseq][i];

  return(asec2rad(arg));
}

void get_fw_all(double *fw,double t)
{
  int seq;

  for (seq=0;seq<4;seq++) fw[seq] = get_fw_arg(seq,t);
}

void fw_to_matrix(double *fw,mat r)
{
  mat rgam,rphi,rpsi,reps;

  rotate(ZAXIS,fw[FW_GAM],rgam);
  rotate(XAXIS,fw[FW_PHI],rphi);
  rotate(ZAXIS,-fw[FW_PSI],rpsi);
  rotate(XAXIS,-fw[FW_EPS],reps);

  matrix_quad_product(reps,rpsi,rphi,rgam,r);
}

double cio_sum(int k,double *fa)
{
  int i,j;
  double a,s;

  for (i=cio_ib[k],s=0.0;i>=cio_ia[k];i--)
  {
    for (j=0,a=0.0;j<8;j++) a += (double)cio_ks[i][j] * fa[j];
    s += cio_ss[i][0] * sin(a) + cio_ss[i][1] * cos(a);
  }

  return(s);
}

double cio_locator(double t,double x,double y)
{
  double fa[8],sk[6],s;
  int k;

  for (k=0;k<5;k++) fa[k] = get_fund_arg(t,k);

  fa[5] = get_planet_arg(t,PL_VENUS);
  fa[6] = get_planet_arg(t,PL_EARTH);
  fa[7] = get_planet_arg(t,PL_ACCUM);

  for (k=0;k<6;k++) sk[k] = cio_sp[k];
  for (k=0;k<5;k++) sk[k] += cio_sum(k,fa);

  for (k=5,s=0.0;k>=0;k--) s = s * t + sk[k];

  return(asec2rad(s) - x * y / 2.0);
}

double earth_rot_angle(datetime_type *ut1)
{
  double frac;

  frac = era_coeff[0]
        + fmod(era_coeff[1] * ut1->djd,1.0)
        + ut1->jdf + 0.5; /* The result is never negative */

  return(fmod(frac,1.0) * TWOPI);
}

double gmstp(datetime_type *tdt)
{
  int k;
  double s;

  for (k=5,s=0.0;k>=0;k--) s = s * tdt->jcen + gmstp_coeff[k];

  return(asec2rad(s));
}

double get_gsd(datetime_type *ut1,datetime_type *tdt)
{
  double gsd;

  gsd = (ut1->djd + 2458257.0)
      + (era_coeff[0] + era_coeff[1] * ut1->djd + gmstp(tdt)/TWOPI);

  return(gsd);
}

void deflection_get_apparent(vec p,vec r_hel,vec app)
{
  double gam,pe,g;
  vec e,a,b,d,s;

  gam = 2.0 * DEFLECT_CON * scalar_product(r_hel,r_hel);
  unit_vector(r_hel,e);
  pe = scalar_product(p,e);
  g = gam / (1.0 + pe);
  vector_times_scalar(p,pe,a);
  vector_difference(e,a,b);
  vector_times_scalar(b,g,d);
  vector_sum(p,d,s);
  unit_vector(s,app);
}

void deflection_get_true(vec app,vec r_hel,vec p)
{
  vec a,del;

  copy_vector(app,p);
  do
  {
    deflection_get_apparent(p,r_hel,a);
    vector_difference(app,a,del);
    vector_sum(p,del,p);
    unit_vector(p,p);
  }
  while(vector_length(del) > 1.0e-12);
}

void aberration_get_apparent(vec p,vec v_obs,vec app)
{
  vec v,a,b,s;
  double gam,pv,g;

  vector_times_scalar(v_obs,UNIT_VEL/LIGHT_SPEED,v);
  gam = sqrt(1.0 - scalar_product(v,v));
  pv = scalar_product(p,v);
  g = 1.0 + pv / (1.0 + gam);
  vector_times_scalar(p,gam,a);
  vector_times_scalar(v,g,b);
  vector_sum(a,b,s);
  unit_vector(s,app);
}

void aberration_get_true(vec app,vec v_obs,vec p)
{
  vec a,del;

  copy_vector(app,p);
  do
  {
    aberration_get_apparent(p,v_obs,a);
    vector_difference(app,a,del);
    vector_sum(p,del,p);
    unit_vector(p,p);
  }
  while(vector_length(del) > 1.0e-12);
}

/* ----------------------------- Leap Seconds ----------------------------- */

void load_leap(void)
{
  file_type f;
  int32_t flen;
  int line_count,seq;
  char row[80];
  int day,month,year,adj,leap_old,leap_new;

  if (leap_list != NULL) get_error("load_leap: Leap seconds already loaded!");

  open_file(&f,"r","%s/astro/leap.dat",hrsp_dir());
  get_fseek(&f,0,SEEK_END);
  flen = ftell(f.dat);
  if (flen < 336L) get_error("load_leap: Leap seconds file too short!");
  line_count = (int)(flen / 12L);
  leap_list = (leap_type *)get_space(line_count * sizeof(leap_type));
  get_fseek(&f,0,SEEK_SET);
  seq = 0;
  leap_old = 10; /* Leap seconds accumulated up to 1 Jan 1972 */
  while (fgets(row,80,f.dat) != NULL)
  {
    if (seq >= line_count)
      get_error("load_leap: Unexpected leap seconds data line!");
    if (strlen(row) > 64)
      get_error("load_leap: Leap seconds data line too long!");
    if (strlen(row) < 7)
      get_error("load_leap: Leap seconds data line too short!");
    if (sscanf(row," %d %d %d %d",&year,&month,&day,&adj) != 4)
      get_error("load_leap: Bad leap seconds data line!");
    if (abs(adj) != 1)
      get_error("load_leap: Invalid leap second adjustment!");
    leap_new = leap_old + adj;
    set_date_dmy(&leap_list[seq].utc,day,month,year);
    if (seq > 0)
    {
      if (leap_list[seq].utc.jd <= leap_list[seq-1].utc.jd)
        get_error("load_leap: Date out of order!");
    }
    leap_list[seq].leap_old = leap_old;
    leap_list[seq].adj = adj;
    leap_list[seq].leap_new = leap_new;
    if (leap_new >= 0)
    {
      next_day(&leap_list[seq].utc,&leap_list[seq].tai.date);
      set_time_decsec(&leap_list[seq].tai.time,(double)leap_new);
    }
    else
    {
      leap_list[seq].tai.date = leap_list[seq].utc;
      set_time_decsec(&leap_list[seq].tai.time,(double)(leap_new+86400));
    }
    leap_old = leap_list[seq].leap_new;
    seq++;
  }
  leap_count = seq;
  get_fclose(&f);
}

void free_leap(void)
{
  get_free((void *)leap_list);
  leap_list = NULL;
  leap_count = 0;
}

int get_leap_count(void)
{
  return(leap_count);
}

void get_leap_item(leap_type *a,int seq)
{
  if (leap_count < 1)
    get_error("get_leap_item: Leap seconds not loaded!");
  if (seq < 0 || seq >= leap_count)
    get_error("get_leap_item: Leap item sequence number out of range!");
  *a = leap_list[seq];
}

int get_leap_total_utc(int32_t utc_jd)
{
  int a,b,c;

  if (leap_count < 1)
    get_error("get_leap_total_utc: Leap seconds not loaded!");

  a = 0;
  b = leap_count - 1;

  while (b - a > 1)
  {
    c = (a + b) / 2;
    if (utc_jd < leap_list[c].utc.jd)
      b = c;
    else
      a = c;
  }

  if (utc_jd > leap_list[b].utc.jd)
    return(leap_list[b].leap_new);

  if (utc_jd > leap_list[a].utc.jd)
    return(leap_list[a].leap_new);
  else
    return(leap_list[a].leap_old);
}

int get_leap_adj_utc(int32_t utc_jd)
{
  int i,adj;

  if (leap_count < 1)
    get_error("get_leap_adj_utc: Leap seconds not loaded!");

  for (i=adj=0;i<leap_count;i++)
  {
    if (leap_list[i].utc.jd == utc_jd) adj = leap_list[i].adj;
  }
  return(adj);
}

int get_leap_total_tai(datetime_type *tai)
{
  int a,b,c;

  if (leap_count < 1)
    get_error("get_leap_total_tai: Leap seconds not loaded!");

  a = 0;
  b = leap_count - 1;

  while (b - a > 1)
  {
    c = (a + b) / 2;
    if (compare_datetime(tai,&leap_list[c].tai) < 0)
      b = c;
    else
      a = c;
  }

  if (compare_datetime(tai,&leap_list[b].tai) >= 0)
    return(leap_list[b].leap_new);

  if (compare_datetime(tai,&leap_list[a].tai) >= 0)
    return(leap_list[a].leap_new);
  else
    return(leap_list[a].leap_old);
}

int leap_seq(datetime_type *tai)
{
  int i;
  int32_t jd;
  double decsec;
  leap_type *leap;

  if (leap_count < 1)
    get_error("leap_tai: Leap seconds not loaded!");

  jd = tai->date.jd;
  decsec = tai->time.decsec + 1.0;
  if (decsec >= FULL_TSEC) decsec -= FULL_TSEC, jd++;

  for (i=leap_count-1;i>=0;i--)
  {
    leap = &leap_list[i];
    if (leap->adj < 0) continue;
    if (leap->tai.date.jd != jd) continue;
    if (leap->tai.time.decsec != decsec) continue;
    return(i);
  }

  return(-1);
}

/* ------------------------------- Delta T  ------------------------------- */

void load_deltat(void)
{
  file_type f;
  int32_t flen;
  int line_count,seq;
  char row[80];
  int day,month,year;
  double deltat;

  if (deltat_list != NULL) get_error("load_deltat: Delta T already loaded!");

  open_file(&f,"r","%s/astro/deltat.dat",hrsp_dir());
  get_fseek(&f,0,SEEK_END);
  flen = ftell(f.dat);
  if (flen < 12L) get_error("load_deltat: Delta T file too short!");
  line_count = (int)(flen / 12L);
  if (line_count < 8) get_error("load_deltat: Too few Delta T data lines!");
  deltat_list = (deltat_type *)get_space(line_count * sizeof(deltat_type));
  get_fseek(&f,0,SEEK_SET);
  seq = 0;
  while (fgets(row,80,f.dat) != NULL)
  {
    if (seq >= line_count)
      get_error("load_deltat: Unexpected Delta T data line!");
    if (strlen(row) > 64)
      get_error("load_deltat: Delta T data line too long!");
    if (strlen(row) < 12)
      get_error("load_deltat: Delta T data line too short!");
    if (sscanf(row," %d %d %d %lf",&year,&month,&day,&deltat) != 4)
      get_error("load_deltat: Bad Delta T data line!");
    if (deltat < 0.0) get_error("load_deltat: Invalid Delta T value!");
    set_date_dmy(&deltat_list[seq].utc,day,month,year);
    if (seq > 0)
    {
      if (deltat_list[seq].utc.jd <= deltat_list[seq-1].utc.jd)
        get_error("load_deltat: Date out of order!");
    }
    deltat_list[seq].deltat = deltat;
    seq++;
  }
  deltat_count = seq;
  get_fclose(&f);
}

void free_deltat(void)
{
  get_free((void *)deltat_list);
  deltat_list = NULL;
  deltat_count = 0;
}

int get_deltat_count(void)
{
  return(deltat_count);
}

void get_deltat_item(deltat_type *a,int seq)
{
  if (deltat_count < 1)
    get_error("get_deltat_item: Delta T values not loaded!");
  if (seq < 0 || seq >= deltat_count)
    get_error("get_deltat_item: Delta T item sequence number out of range!");
  *a = deltat_list[seq];
}

double get_deltat_utc(datetime_type *utc)
{
  int a,b,c;
  double dx,dy,x;

  if (deltat_count < 1)
    get_error("get_deltat_utc: Delta T values not loaded!");

  a = 0;
  b = deltat_count - 1;

  while (b - a > 1)
  {
    c = (a + b) / 2;
    if (utc->date.jd < deltat_list[c].utc.jd)
      b = c;
    else
      a = c;
  }

  if (utc->date.jd < deltat_list[a].utc.jd)
    return(deltat_list[a].deltat);
  if (utc->date.jd >= deltat_list[b].utc.jd)
    return(deltat_list[b].deltat);
  dx = deltat_list[b].utc.jd0 - deltat_list[a].utc.jd0;
  dy = deltat_list[b].deltat - deltat_list[a].deltat;
  x = (utc->jd0 - deltat_list[a].utc.jd0) + utc->jdf;
  return(deltat_list[a].deltat + (dy / dx) * x);
}

/* --------------------------------- JPL  --------------------------------- */

void load_jpl(void)
{
  file_type f;
  int32_t flen;
  int i;
  byte *src;
  int32_t kmax,nd;
  int imax;

  if (jpl != NULL) get_error("load_jpl: JPL ephemeris already loaded!");

  open_file(&f,"r","%s/astro/de405",hrsp_dir());
  get_fseek(&f,0,SEEK_END);
  flen = ftell(f.dat);
  if (flen < 2856) get_error("load_jpl: JPL ephemeris file too short!");
  jpl = (byte *)get_space(flen);
  get_fseek(&f,0,SEEK_SET);
  get_fread((void *)jpl,1,flen,&f);
  get_fclose(&f);

  memset(&jplh,0,sizeof(jpl_header));

  memmove(jplh.ephem_id,src=jpl,EPHEM_ID_LENGTH);
  memmove(jplh.ephem_start_epoch,src+=EPHEM_ID_LENGTH,EPHEM_ID_LENGTH);
  memmove(jplh.ephem_final_epoch,src+=EPHEM_ID_LENGTH,EPHEM_ID_LENGTH);

  jplh.first_epoch = *((double *)(jpl+2652));
  jplh.last_epoch = *((double *)(jpl+2660));
  jplh.days_per_record = *((double *)(jpl+2668));

  jplh.ccount = *((int32_t *)(jpl+2676));

  for (i=0,src=jpl+252;i<jplh.ccount;i++,src+=CNAMLEN)
    memmove(jplh.cnam[i],src,CNAMLEN);

  jplh.au = *((double *)(jpl+2680));
  jplh.invau = 1.0 / jplh.au;
  jplh.emrat = *((double *)(jpl+2688));
  jplh.gamma = 1.0 + jplh.emrat;

  for (i=MERCURY,src=jpl+2696;i<=NUTATION;i++,src+=12)
  {
    jplh.pt[i] = *((int32_t *)(src));
    jplh.ncf[i] = *((int32_t *)(src+4));
    jplh.na[i] = *((int32_t *)(src+8));
  }

  jplh.numde = *((int32_t *)(jpl+2840));

  jplh.pt[LIBRATION] = *((int32_t *)(jpl+2844));
  jplh.ncf[LIBRATION] = *((int32_t *)(jpl+2848));
  jplh.na[LIBRATION] = *((int32_t *)(jpl+2852));

  for (i=MERCURY,kmax=0,imax=0;i<=LIBRATION;i++)
    if (jplh.pt[i] > kmax) kmax=jplh.pt[i],imax=i;

  if (imax == NUTATION)
    nd = 2L;
  else
    nd = 3L;

  jplh.ksize =
    2L * (kmax + nd * jplh.ncf[imax] * jplh.na[imax] - 1L);

  jplh.recsize = jplh.ksize * 4L;

  memmove(jplh.cval,jpl+jplh.recsize,8*jplh.ccount);
}

void free_jpl(void)
{
  get_free((void *)jpl);
  jpl = NULL;
}

void jpl_interpolate
  (double *buf,double *t,int ncf,int ncm,int na,int ifl,double *pv)
{
  static double pc[18],vc[18];
  static int np=2,nv=3,first_call=1;
  static double twot=0.0;
  double dna,dt0,temp,tc,vfac;
  int i,j,k,ncl,skip;

  if (first_call)
  {
    pc[0] = 1.0;
    pc[1] = 0.0;
    vc[1] = 1.0;
    first_call = 0;
  }

  memset(pv,0,6*sizeof(double));

  dna = (double)na;
  dt0 = (double)(int)t[0];
  temp = dna * t[0];
  ncl = (int)(temp - dt0);

  tc = 2.0 * (fmod(temp,1.0) + dt0) - 1.0;

  if (tc != pc[1])
  {
    np = 2;
    nv = 3;
    pc[1] = tc;
    twot = tc + tc;
  }

  if (np < ncf)
  {
    for (i=np;i<ncf;i++) pc[i] = twot * pc[i-1] - pc[i-2];
    np = ncf;
  }

  for (i=0,skip=ncl*ncf*ncm;i<ncm;i++,skip+=ncf)
  {
    pv[i] = 0.0;
    for (j=ncf-1;j>=0;j--) pv[i] += pc[j] * buf[skip+j];
  }

  if (ifl != POSVEL) return;

  vfac = (dna + dna) / t[1];
  vc[2] = twot + twot;

  if (nv < ncf)
  {
    for (i=nv;i<ncf;i++) vc[i] = twot * vc[i-1] + 2.0 * pc[i-1] - vc[i-2];
    nv = ncf;
  }

  for (i=0,k=ncm,skip=ncl*ncf*ncm;i<ncm;i++,k++,skip+=ncf)
  {
    pv[k] = 0.0;
    for (j=ncf-1;j>=0;j--) pv[k] += vc[j] * buf[skip+j];
    pv[k] *= vfac;
  }
}

void jpl_target(datetime_type *tdb,int targ,int fl,double *pv)
{
  double u[2],*cbuf;
  int32_t nr;
  double jd0;
  int i;

  if (jpl == NULL) get_error("jpl_target: JPL ephemeris data not loaded!");

  if (tdb->jd < jplh.first_epoch || tdb->jd > jplh.last_epoch)
    get_error("jpl_target: Interpolation time out of limits!");

  nr = (int)((tdb->jd0 - jplh.first_epoch) / jplh.days_per_record) + 2L;
  if (tdb->jd0 == jplh.last_epoch) nr--;

  jd0 = (double)(nr - 2L) * jplh.days_per_record + jplh.first_epoch;
  u[0] = ((tdb->jd0 - jd0) + tdb->jdf) / jplh.days_per_record;
  u[1] = jplh.days_per_record;

  cbuf = (double *)(jpl + nr * jplh.recsize);

  switch(targ)
  {
    case NUTATION:
      if (jplh.ncf[targ] <= 0)
        get_error
          ("jpl_target: JPL ephemeris file does not include nutation!");
      jpl_interpolate
        (cbuf+jplh.pt[targ]-1,u,jplh.ncf[targ],2,jplh.na[targ],fl,pv);
      break;
    case LIBRATION:
      if (jplh.ncf[targ] <= 0)
        get_error
          ("jpl_target: JPL ephemeris file does not include libration!");
      jpl_interpolate
        (cbuf+jplh.pt[targ]-1,u,jplh.ncf[targ],3,jplh.na[targ],fl,pv);
      break;
    default:
      jpl_interpolate
        (cbuf+jplh.pt[targ]-1,u,jplh.ncf[targ],3,jplh.na[targ],fl,pv);
      for (i=0;i<6;i++) pv[i] *= jplh.invau;
  }
}

void jpl_posvel(datetime_type *tdb,int targ,vec r,vec v)
{
  double pv[6];

  jpl_target(tdb,targ,POSVEL,pv);
  memmove(&r[0],&pv[0],3*sizeof(double));
  memmove(&v[0],&pv[3],3*sizeof(double));
}

void jpl_earth(datetime_type *tdb,vec rh,vec vh,vec rb,vec vb)
{
  vec rs,vs,re,ve,rm,vm,rc,vc;

  jpl_posvel(tdb,SUN,rs,vs);
  jpl_posvel(tdb,EMBARY,re,ve);
  jpl_posvel(tdb,MOON,rm,vm);

  vector_div_scalar(rm,jplh.gamma,rc);
  vector_div_scalar(vm,jplh.gamma,vc);

  vector_difference(re,rc,rb);
  vector_difference(ve,vc,vb);

  vector_difference(rb,rs,rh);
  vector_difference(vb,vs,vh);
}

/* -------------------------------- Astro  -------------------------------- */

void load_astro(void)
{
  load_jpl();
  load_leap();
  load_deltat();
}

void free_astro(void)
{
  free_jpl();
  free_leap();
  free_deltat();
}

/* -------------------------------- Epoch  -------------------------------- */

void clear_epoch(epoch_type *w)
{
  memset(w,0,sizeof(epoch_type));
}

void make_epoch_data(epoch_type *w)
{
  double fw[4];
  double x,y,z,rr,e,d,s,a,p,q,stp;
  mat re,rd,res,gst;
  vec rs,mi,mj;

  get_fw_all(fw,w->tdt.jcen);
  fw_to_matrix(fw,w->gcrs_to_mean_eq);
  transpose_matrix(w->gcrs_to_mean_eq,w->mean_eq_to_gcrs);
  w->epsilon = fw[FW_EPS];
  nut06a(w->tdt.jcen,&w->delpsi,&w->deleps);
  fw[FW_PSI] += w->delpsi;
  fw[FW_EPS] += w->deleps;
  fw_to_matrix(fw,w->gcrs_to_true_eq);
  transpose_matrix(w->gcrs_to_true_eq,w->true_eq_to_gcrs);

  x = w->gcrs_to_true_eq[2][0];
  y = w->gcrs_to_true_eq[2][1];
  z = w->gcrs_to_true_eq[2][2];

  rr = x*x + y*y;
  if (rr > 0.0)
    e = atan2(y,x);
  else
    e = 0.0;
  d = atan(sqrt(rr/(1.0-rr)));
  s = cio_locator(w->tdt.jcen,x,y);
  w->xpol = x;
  w->ypol = y;
  w->orig = s;
  rotate(ZAXIS,-(e+s),res);
  rotate(YAXIS,d,rd);
  rotate(ZAXIS,e,re);
  matrix_triple_product(res,rd,re,w->gcrs_to_cio);

  a = 1.0 / (1.0 + z);
  set_vector(rs,1.0 - a * x * x,-a * x * y,-x);
  matrix_row(w->gcrs_to_true_eq,0,mi);
  matrix_row(w->gcrs_to_true_eq,1,mj);
  p = scalar_product(mi,rs);
  q = scalar_product(mj,rs);
  w->equorig = s;
  if (p!=0.0 || q!=0.0) w->equorig -= atan2(q,p);

  stp = gmstp(&w->tdt);
  w->equeq = -(w->equorig + stp);
  w->erot = earth_rot_angle(&w->ut1);

  set_alpha_rad(&w->gmst,w->erot+stp);
  set_alpha_rad(&w->gast,w->erot-w->equorig);

  rotate(ZAXIS,w->gast.rad,gst);
  matrix_product(gst,w->gcrs_to_true_eq,w->gcrs_to_earth);
  transpose_matrix(w->gcrs_to_earth,w->earth_to_gcrs);

  jpl_earth(&w->tdb,w->reh,w->veh,w->reb,w->veb);
}

double get_tdb_correction(datetime_type *tdt)
{
  double g,u,sing,sinu;

  g = 357.53 + 0.98560028 * tdt->djd;
  u = 246.11 + 0.90251792 * tdt->djd;

  sing = sin(deg2rad(g));
  sinu = sin(deg2rad(u));

  return(0.001657 * sing + 0.000022 * sinu);
}

void tai_to_tdt(datetime_type *tai,datetime_type *tdt)
{
  shift_datetime(tai,tdt,TAI_OFFSET);
}

void tdt_to_tai(datetime_type *tdt,datetime_type *tai)
{
  shift_datetime(tdt,tai,-TAI_OFFSET);
}

void tdt_to_tdb(datetime_type *tdt,datetime_type *tdb)
{
  double delta;

  delta = get_tdb_correction(tdt);
  shift_datetime(tdt,tdb,delta);
}

void tdb_to_tdt(datetime_type *tdb,datetime_type *tdt)
{
  double delta;

  delta = get_tdb_correction(tdb);
  shift_datetime(tdb,tdt,-delta);
}

void utc_to_tai(datetime_type *utc,datetime_type *tai)
{
  int leap;

  leap = get_leap_total_utc(utc->date.jd);
  shift_datetime(utc,tai,(double)leap);
}

void tai_to_utc(datetime_type *tai,datetime_type *utc)
{
  double fsec,isec;
  datetime_type tmp;
  int seq,leap,adj;

  fsec = modf(tai->time.decsec,&isec);

  tmp.date = tai->date;
  set_time_decsec(&tmp.time,isec);
  seq = leap_seq(&tmp);
  if (seq >= 0)
  {
    utc->date = leap_list[seq].utc;
    set_time_decsec_adj(&utc->time,FULL_TSEC + fsec,1);
  }
  else
  {
    leap = get_leap_total_tai(tai);
    shift_datetime(tai,utc,-(double)leap);
    adj = get_leap_adj_utc(utc->date.jd);
    set_time_decsec_adj(&utc->time,utc->time.decsec,adj);
  }
  update_time_args(utc);
}

void set_datetime_utc(datetime_type *utc,
  int day,int month,int year,int hour,int min,double sec)
{
  date_type date;
  int adj;

  set_date_dmy(&date,day,month,year);
  adj = get_leap_adj_utc(date.jd);
  set_datetime_adj(utc,day,month,year,hour,min,sec,adj);
}

void use_epoch_utc(epoch_type *w)
{
  utc_to_tai(&w->utc,&w->tai);
  tai_to_tdt(&w->tai,&w->tdt);
  tdt_to_tdb(&w->tdt,&w->tdb);
}

void use_epoch_tai(epoch_type *w)
{
  tai_to_tdt(&w->tai,&w->tdt);
  tdt_to_tdb(&w->tdt,&w->tdb);
  tai_to_utc(&w->tai,&w->utc);
}

void use_epoch_tdt(epoch_type *w)
{
  tdt_to_tai(&w->tdt,&w->tai);
  tdt_to_tdb(&w->tdt,&w->tdb);
  tai_to_utc(&w->tai,&w->utc);
}

void use_epoch_tdb(epoch_type *w)
{
  tdb_to_tdt(&w->tdb,&w->tdt);
  tdt_to_tai(&w->tdt,&w->tai);
  tai_to_utc(&w->tai,&w->utc);
}

void derive_ut1(epoch_type *w)
{
  w->leap = get_leap_total_utc(w->utc.date.jd);
  w->deltat = get_deltat_utc(&w->utc);
  shift_datetime(&w->tdt,&w->ut1,-w->deltat);
}

void finalize_epoch(epoch_type *w)
{
  derive_ut1(w);
  make_epoch_data(w);
}

void set_epoch_utc(epoch_type *w,
  int day,int month,int year,int hour,int min,double sec)
{
  set_datetime_utc(&w->utc,day,month,year,hour,min,sec);
  use_epoch_utc(w);
  finalize_epoch(w);
}

void set_epoch_tdt(epoch_type *w,
  int day,int month,int year,int hour,int min,double sec)
{
  set_datetime(&w->tdt,day,month,year,hour,min,sec);
  use_epoch_tdt(w);
  finalize_epoch(w);
}

void set_epoch_tdb(epoch_type *w,
  int day,int month,int year,int hour,int min,double sec)
{
  set_datetime(&w->tdb,day,month,year,hour,min,sec);
  use_epoch_tdb(w);
  finalize_epoch(w);
}

void new_epoch_ut1(epoch_type *w)
{
  w->deltat = get_deltat_utc(&w->ut1);
  shift_datetime(&w->ut1,&w->tdt,w->deltat);
  use_epoch_tdt(w);
  w->leap = get_leap_total_utc(w->utc.date.jd);
  make_epoch_data(w);
}

void set_epoch_ut1(epoch_type *w,
  int day,int month,int year,int hour,int min,double sec)
{
  set_datetime(&w->ut1,day,month,year,hour,min,sec);
  new_epoch_ut1(w);
}

void set_julian_year_tdt(epoch_type *w,double year)
{
  set_julian_year(&w->tdt,year);
  use_epoch_tdt(w);
  finalize_epoch(w);
}

void add_tdt_seconds(epoch_type *a,epoch_type *b,double secs)
{
  shift_datetime(&a->tdt,&b->tdt,secs);
  use_epoch_tdt(b);
  finalize_epoch(b);
}

void add_tdb_seconds(epoch_type *a,epoch_type *b,double secs)
{
  shift_datetime(&a->tdb,&b->tdb,secs);
  use_epoch_tdb(b);
  finalize_epoch(b);
}

void add_ut1_seconds(epoch_type *a,epoch_type *b,double secs)
{
  shift_datetime(&a->ut1,&b->ut1,secs);
  new_epoch_ut1(b);
}

/* -------------------------------- Orbits -------------------------------- */

void set_orbit(orbit_t *orb,double a,double i,double ome,double Ome,
  double e,interval_type *P,datetime_type *T)
{
  mat rlong,rincl,rnode,rswap;

  orb->a = a;
  orb->i = deg2rad(i);
  orb->ome = deg2rad(ome);
  orb->Ome = deg2rad(Ome);
  orb->e = e;
  orb->P = P->days;
  orb->T = T->jd;

  orb->b = a * sqrt(1 - e * e);

  rotate(ZAXIS,-orb->ome,rlong);
  rotate(XAXIS,-orb->i,rincl);
  rotate(ZAXIS,-orb->Ome-HALFPI,rnode);
  rotate(YAXIS,PI,rswap);

  matrix_quad_product(rswap,rnode,rincl,rlong,orb->r);
}

void get_apparent_position(orbit_t *orb,datetime_type *t,vec r)
{
  double dt,M,E,E1;
  vec h;

  dt = t->jd - orb->T;
  M = full_range(TWOPI / orb->P * dt,TWOPI);
  E = M;
  do
  {
    E1 = E;
    E = M + orb->e * sin(E);
  }
  while(fabs(E - E1) >= 1.0e-8);
  set_vector(h,orb->a * (cos(E) - orb->e),orb->b * sin(E),0);
  matrix_times_vector(orb->r,h,r);
}

/* ------------------------------ Topocentre ------------------------------ */

void clear_topo(topo_type *t)
{
  memset(t,0,sizeof(topo_type));
}

void set_topo(topo_type *topo,angle_type *lon,angle_type *lat,
  double height,epoch_type *epoch)
{
  double cosphi,sinphi,coslam,sinlam;
  double a,f,u,q,h;
  double C,S;
  double x,y,z,angvel;
  vec omega,re,ve;

  topo->lon = *lon;
  topo->lat = *lat;
  topo->height = height;
  topo->epoch = *epoch;

  set_alpha_frac(&topo->lmst,topo->epoch.gmst.frac + topo->lon.lambda.frac);
  set_alpha_frac(&topo->last,topo->epoch.gast.frac + topo->lon.lambda.frac);

  cosphi = cos(topo->lat.delta.rad);
  sinphi = sin(topo->lat.delta.rad);
  coslam = cos(topo->lon.lambda.rad);
  sinlam = sin(topo->lon.lambda.rad);

  a = EARTH_EQU_RAD;
  f = 1.0 / EARTH_INV_FLAT;

  u = (1.0 - f);
  q = u * u;

  C = 1.0 / sqrt(cosphi * cosphi + q * sinphi * sinphi);
  S = q * C;

  h = topo->height * 1.0e-3;

  x = (a * C + h) * cosphi * coslam;
  y = (a * C + h) * cosphi * sinlam;
  z = (a * S + h) * sinphi;

  angvel = (1.0 + era_coeff[1]) * TWOPI / FULL_TSEC;
  set_vector(omega,0,0,angvel);

  set_vector(topo->r,x,y,z);
  vector_product(omega,topo->r,topo->v);

  matrix_times_vector(epoch->earth_to_gcrs,topo->r,topo->re);
  matrix_times_vector(epoch->earth_to_gcrs,topo->v,topo->ve);
  vector_div_scalar(topo->re,AU_KM,re);
  vector_div_scalar(topo->ve,UNIT_VEL,ve);
  vector_sum(epoch->reb,re,topo->rb);
  vector_sum(epoch->veb,ve,topo->vb);
  vector_sum(epoch->reh,re,topo->rh);
  vector_sum(epoch->veh,ve,topo->vh);
}

/* -------------------------------- Stars  -------------------------------- */

void clear_star(star_type *s)
{
  memset(s,0,sizeof(star_type));
}

void local_matrices(star_type *s)
{
  mat rx,rz;

  rotate(ZAXIS,HALFPI+s->ra.alpha.rad,rz);
  rotate(XAXIS,HALFPI-s->dec.delta.rad,rx);
  matrix_product(rx,rz,s->xyz_to_loc);

  rotate(XAXIS,s->dec.delta.rad-HALFPI,rx);
  rotate(ZAXIS,-s->ra.alpha.rad-HALFPI,rz);
  matrix_product(rz,rx,s->loc_to_xyz);
}

void space_vectors(star_type *s)
{
  double rho,vx,vy,vz;
  vec r,v;

  local_matrices(s);

  if (s->par > 0.0)
  {
    rho = 1.0 / s->par;
    vx = TRANS_VEL * s->mura / s->par;
    vy = TRANS_VEL * s->mudec / s->par;
    vz = s->rv;
    set_vector(r,0,0,rho);
    set_vector(v,vx,vy,vz);
    s->model = FULL_SPACE_MOTION;
  }
  else
  {
    set_ort(ZAXIS,r);
    set_vector(v,s->mura,s->mudec,0);
    s->model = PROPER_MOTION_ONLY;
  }

  matrix_times_vector(s->loc_to_xyz,r,s->r);
  matrix_times_vector(s->loc_to_xyz,v,s->v);
}

void sky_angles(star_type *s)
{
  double ra,dec,rho;
  vec v;

  xyz_to_sphere(s->r,&ra,&dec,&rho);
  set_angle_rad(&s->ra,ra);
  set_angle_rad(&s->dec,dec);

  local_matrices(s);

  switch(s->model)
  {
    case FULL_SPACE_MOTION:
      s->par = 1.0 / rho;
      matrix_times_vector(s->xyz_to_loc,s->v,v);
      s->mura = v[0] * s->par / TRANS_VEL;
      s->mudec = v[1] * s->par / TRANS_VEL;
      s->rv = v[2];
      break;
    case PROPER_MOTION_ONLY:
      s->par = 0.0;
      matrix_times_vector(s->xyz_to_loc,s->v,v);
      s->mura = v[0];
      s->mudec = v[1];
      break;
    default: get_error("sky_angles: Invalid space motion model!");
  }
}

void cent_bary(star_type *s)
{
  vec robs;

  switch(s->cent)
  {
    case BARY_CENT:
      return;
    case GEO_CENT:
      vector_div_scalar(s->epoch.reb,PC_AU,robs);
      vector_sum(robs,s->r,s->r);
      break;
    case TOPO_CENT:
      vector_div_scalar(s->topo.rb,PC_AU,robs);
      vector_sum(robs,s->r,s->r);
      break;
    default: get_error("cent_bary: Invalid observer location!");
  }
}

void cent_obs(star_type *s)
{
  vec robs;

  switch(s->cent)
  {
    case BARY_CENT:
      return;
    case GEO_CENT:
      vector_div_scalar(s->epoch.reb,PC_AU,robs);
      vector_difference(s->r,robs,s->r);
      break;
    case TOPO_CENT:
      vector_div_scalar(s->topo.rb,PC_AU,robs);
      vector_difference(s->r,robs,s->r);
      break;
    default: get_error("cent_obs: Invalid observer location!");
  }
}

void apply_deflection(star_type *s)
{
  vec p,app;

  unit_vector(s->r,p);
  switch(s->cent)
  {
    case GEO_CENT:
      deflection_get_apparent(p,s->epoch.reh,app);
      break;
    case TOPO_CENT:
      deflection_get_apparent(p,s->topo.rh,app);
      break;
    default: get_error("apply_deflection: Invalid observer location!");
  }
  vector_times_scalar(app,vector_length(s->r),s->r);
}

void remove_deflection(star_type *s)
{
  vec p,app;

  unit_vector(s->r,app);
  switch(s->cent)
  {
    case GEO_CENT:
      deflection_get_true(app,s->epoch.reh,p);
      break;
    case TOPO_CENT:
      deflection_get_true(app,s->topo.rh,p);
      break;
    default: get_error("remove_deflection: Invalid observer location!");
  }
  vector_times_scalar(p,vector_length(s->r),s->r);
}

void apply_aberration(star_type *s)
{
  vec p,app;

  unit_vector(s->r,p);
  switch(s->cent)
  {
    case GEO_CENT:
      aberration_get_apparent(p,s->epoch.veb,app);
      break;
    case TOPO_CENT:
      aberration_get_apparent(p,s->topo.vb,app);
      break;
    default: get_error("apply_aberration: Invalid observer location!");
  }
  vector_times_scalar(app,vector_length(s->r),s->r);
}

void remove_aberration(star_type *s)
{
  vec p,app;

  unit_vector(s->r,p);
  switch(s->cent)
  {
    case GEO_CENT:
      aberration_get_true(p,s->epoch.veb,app);
      break;
    case TOPO_CENT:
      aberration_get_true(p,s->topo.vb,app);
      break;
    default: get_error("remove_aberration: Invalid observer location!");
  }
  vector_times_scalar(app,vector_length(s->r),s->r);
}

void rotate_to_gcrs(star_type *s)
{
  mat rot;
  vec r,v;

  switch(s->ref)
  {
    case ICRS:
      return;
    case TRUE_EQ:
      copy_matrix(s->equinox.true_eq_to_gcrs,rot);
      break;
    case MEAN_EQ:
      copy_matrix(s->equinox.mean_eq_to_gcrs,rot);
      break;
    default: get_error("rotate_to_gcrs: Invalid reference frame!");
  }

  matrix_times_vector(rot,s->r,r);
  matrix_times_vector(rot,s->v,v);

  copy_vector(r,s->r);
  copy_vector(v,s->v);
}

void rotate_to_obs(star_type *s)
{
  mat rot;
  vec r,v;

  switch(s->ref)
  {
    case ICRS:
      return;
    case TRUE_EQ:
      copy_matrix(s->equinox.gcrs_to_true_eq,rot);
      break;
    case MEAN_EQ:
      copy_matrix(s->equinox.gcrs_to_mean_eq,rot);
      break;
    default: get_error("rotate_to_obs: Invalid reference frame!");
  }

  matrix_times_vector(rot,s->r,r);
  matrix_times_vector(rot,s->v,v);

  copy_vector(r,s->r);
  copy_vector(v,s->v);
}

void space_motion(star_type *s,double dt)
{
  double mu;
  vec h,r,v;
  mat rot;

  copy_vector(s->r,r);
  copy_vector(s->v,v);

  switch(s->model)
  {
    case FULL_SPACE_MOTION:
      vector_times_scalar(v,dt*YEAR_SEC/PC_KM,h);
      vector_sum(r,h,s->r);
      break;
    case PROPER_MOTION_ONLY:
      mu = vector_length(v);
      if (mu == 0.0) break;
      vector_product(r,v,h);
      arbitrary_rotation(h,asec2rad(mu*dt),rot);
      matrix_times_vector(rot,r,s->r);
      matrix_times_vector(rot,v,s->v);
      break;
    default: get_error("space_motion: Invalid space motion model!");
  }
}

void locate_components(star_type *s,vec ab,vec ac,vec cb)
{
  vec u;

  if (s->mass_ratio > 0.0)
  {
    get_apparent_position(&s->orb,&s->epoch.tdt,u);
    matrix_times_vector(s->loc_to_xyz,u,ab);
    vector_times_scalar(ab,ASEC2R*vector_length(s->r),ab);
    vector_times_scalar(ab,s->mass_ratio,ac);
    vector_difference(ab,ac,cb);
  }
  else
  {
    zero_vector(ab);
    zero_vector(ac);
    zero_vector(cb);
  }
}

void select_primary_star(star_type *s)
{
  vec ab,ac,cb,u;

  locate_components(s,ab,ac,cb);

  switch(s->component)
  {
    case PRIMARY_STAR:
      copy_vector(s->r,u);
      break;
    case SECONDARY_STAR:
      vector_difference(s->r,ab,u);
      break;
    case CENTER_OF_MASS:
      vector_difference(s->r,ac,u);
      break;
    default: get_error("select_primary_star: Invalid component!");
  }
  unit_vector(u,u);
  vector_times_scalar(u,vector_length(s->r),s->r);
  s->component = PRIMARY_STAR;
}

void select_secondary_star(star_type *s)
{
  vec ab,ac,cb,u;

  locate_components(s,ab,ac,cb);

  switch(s->component)
  {
    case PRIMARY_STAR:
      vector_sum(s->r,ab,u);
      break;
    case SECONDARY_STAR:
      copy_vector(s->r,u);
      break;
    case CENTER_OF_MASS:
      vector_sum(s->r,cb,u);
      break;
    default: get_error("select_secondary_star: Invalid component!");
  }
  unit_vector(u,u);
  vector_times_scalar(u,vector_length(s->r),s->r);
  s->component = SECONDARY_STAR;
}

void select_center_of_mass(star_type *s)
{
  vec ab,ac,cb,u;

  locate_components(s,ab,ac,cb);

  switch(s->component)
  {
    case PRIMARY_STAR:
      vector_sum(s->r,ac,u);
      break;
    case SECONDARY_STAR:
      vector_difference(s->r,cb,u);
      break;
    case CENTER_OF_MASS:
      copy_vector(s->r,u);
      break;
    default: get_error("select_center_of_mass: Invalid component!");
  }
  unit_vector(u,u);
  vector_times_scalar(u,vector_length(s->r),s->r);
  s->component = CENTER_OF_MASS;
}

void select_component(star_type *s,int comp)
{
  switch(comp)
  {
    case PRIMARY_STAR:
      select_primary_star(s);
      break;
    case SECONDARY_STAR:
      select_secondary_star(s);
      break;
    case CENTER_OF_MASS:
      select_center_of_mass(s);
      break;
    default: get_error("select_component: Invalid component!");
  }
}

void move_star(star_type *a,star_type *b,int ref,int cent,int pos,
  epoch_type *eq,epoch_type *ep,topo_type *topo)
{
  *b = *a;
  space_vectors(b);
  rotate_to_gcrs(b);
  if (a->pos == APPARENT)
  {
    remove_aberration(b);
    remove_deflection(b);
  }
  cent_bary(b);
  select_center_of_mass(b);

  b->ref = ref;
  b->cent = cent;
  b->pos = pos;
  if (eq != NULL)
    b->equinox = *eq;
  else
    clear_epoch(&b->equinox);
  b->epoch = *ep;
  if (topo != NULL)
    b->topo = *topo;
  else
    clear_topo(&b->topo);

  space_motion(b,b->epoch.tdb.jyr - a->epoch.tdb.jyr);

  select_component(b,a->component);
  cent_obs(b);
  if (b->pos == APPARENT)
  {
    apply_deflection(b);
    apply_aberration(b);
  }
  rotate_to_obs(b);
  sky_angles(b);
}

void full_corr(star_type *astar,star_type *bstar,epoch_type *epoch,
  topo_type *topo,double *vorb,double *vrot,double *dv,double *dt)
{
  vec p;
  double dr;

  move_star(astar,bstar,ICRS,TOPO_CENT,GEOMETRIC,epoch,epoch,topo);

  space_vectors(bstar);
  unit_vector(bstar->r,p);

  *vorb = scalar_product(epoch->veb,p) * UNIT_VEL;
  *vrot = scalar_product(topo->ve,p);
  dr = scalar_product(topo->rb,p);

  *dv = *vorb + *vrot;
  *dt = (dr * LIGHT_TIME) / 86400.0;
}

void star_corr(star_type *astar,star_type *bstar,epoch_type *epoch,
  topo_type *topo,double *dv,double *dt)
{
  double vorb,vrot;

  full_corr(astar,bstar,epoch,topo,&vorb,&vrot,dv,dt);
}

void sun_corr(epoch_type *epoch,topo_type *topo,double *dv,double *dt)
{
  vec r,p;
  double vorb,vrot,dr;

  opposite_vector(topo->rh,r);
  unit_vector(r,p);

  vorb = scalar_product(epoch->veh,p);
  vrot = scalar_product(topo->ve,p);
  dr = scalar_product(topo->rh,p);

  *dv = vorb * UNIT_VEL + vrot;
  *dt = (dr * LIGHT_TIME) / 86400.0;
}

/* ----------------------------- Star labels ------------------------------ */

void clear_label(starlabel_type *lab)
{
  memset(lab,0,sizeof(starlabel_type));
}

void clear_lab_bayer(starlabel_type *lab)
{
  memset(lab->bayer,0,4);
}

void clear_lab_cons(starlabel_type *lab)
{
  memset(lab->cons,0,4);
}

void clear_lab_comp(starlabel_type *lab)
{
  memset(lab->comp,0,4);
}

void set_hip_label(starlabel_type *lab,int hip,char *comp)
{
  clear_label(lab);
  lab->kind = HIP_STAR;
  lab->catnum = hip;
  strcpy(lab->comp,comp);
}

void set_hr_label(starlabel_type *lab,int hr,char *comp)
{
  clear_label(lab);
  lab->kind = HR_STAR;
  lab->catnum = hr;
  strcpy(lab->comp,comp);
}

void set_usr_label(starlabel_type *lab,int usr)
{
  clear_label(lab);
  lab->kind = USR_STAR;
  lab->catnum = usr;
}

/* ------------------------------ Star names ------------------------------ */

void compact_star_name(char *instr,char *outstr,int maxlen)
{
  int outlen;

  memset(outstr,0,maxlen+1);
  outlen = 0;
  while (*instr != 0 && outlen < maxlen)
  {
    while (*instr != 0 && !alphanum_ok(*instr)) instr++;
    if (*instr != 0 && outlen > 0 && outlen < maxlen)
      *(outstr++) = ' ', outlen++;
    while (alphanum_ok(*instr) && outlen < maxlen)
      *(outstr++) = *(instr++), outlen++;
  }
}

int get_starname_kind(char *starname)
{
  if (digit_ok(*starname)) return(FLM_STAR);
  if (strncasecmp(starname,"HD",2) == 0) return(HD_STAR);
  if (strncasecmp(starname,"HR",2) == 0) return(HR_STAR);
  if (strncasecmp(starname,"HIP",3) == 0) return(HIP_STAR);
  if (strncasecmp(starname,"USR",3) == 0) return(USR_STAR);
  return(BYR_STAR);
}

int parse_greek_letter(char *greekname,char **endptr,starlabel_type *lab)
{
  int seq,n;

  clear_lab_bayer(lab);
  while(*greekname == ' ') greekname++;
  if (!alpha_ok(*greekname)) return(GREEK_LETTER_EXPECTED);
  *endptr = greekname;
  for (seq=0;*greek[seq].full != 0;seq++)
  {
    n = strlen(greek[seq].full);
    if (strncasecmp(greekname,greek[seq].full,n) == 0)
    {
      strcpy(lab->bayer,greek[seq].abbr);
      *endptr = greekname + n;
      break;
    }
    n = strlen(greek[seq].abbr);
    if (strncasecmp(greekname,greek[seq].abbr,n) == 0)
    {
      strcpy(lab->bayer,greek[seq].abbr);
      *endptr = greekname + n;
      break;
    }
  }
  if (*endptr == greekname) return(BAD_GREEK_LETTER);
  return(0);
}

int parse_constellation(char *consname,char **endptr,starlabel_type *lab)
{
  int seq,n;

  clear_lab_cons(lab);
  while(*consname == ' ') consname++;
  if (!alpha_ok(*consname)) return(CONSTELL_LETTER_EXPECTED);
  *endptr = consname;
  for (seq=0;*constell[seq].gen != 0;seq++)
  {
    n = strlen(constell[seq].gen);
    if (strncasecmp(consname,constell[seq].gen,n) == 0)
    {
      strcpy(lab->cons,constell[seq].abbr);
      *endptr = consname + n;
      break;
    }
    if (strncasecmp(consname,constell[seq].abbr,3) == 0)
    {
      strcpy(lab->cons,constell[seq].abbr);
      *endptr = consname + 3;
      break;
    }
  }
  if (*endptr == consname) return(BAD_CONSTELL_NAME);
  return(0);
}

int parse_flamsteed_number(char *flmno,char **endptr,starlabel_type *lab)
{
  lab->flam = parse_uint(flmno,endptr,0);
  if (lab->flam < 1) return(BAD_FLAMSTEED_NUMBER);
  return(0);
}

int parse_bayer_seq(char *bayerseq,char **endptr,starlabel_type *lab)
{
  lab->seq = parse_uint(bayerseq,endptr,0);
  if (lab->seq < 0 || lab->seq > 9) return(BAD_BAYER_SEQ);
  return(0);
}

int parse_catalog_number(char *catno,char **endptr,starlabel_type *lab)
{
  lab->catnum = parse_uint(catno,endptr,0);
  if (lab->catnum < 1) return(BAD_CATALOG_NUMBER);
  return(0);
}

int parse_bayer_star(char *starname,char **endptr,starlabel_type *lab)
{
  int parse_code;

  parse_code = parse_greek_letter(starname,endptr,lab);
  if (parse_code == 0) 
    parse_code = parse_bayer_seq(*endptr,endptr,lab);
  if (parse_code == 0)
    parse_code = parse_constellation(*endptr,endptr,lab);
  return(parse_code);
}

int parse_flamsteed_star(char *starname,char **endptr,starlabel_type *lab)
{
  int parse_code;

  parse_code = parse_flamsteed_number(starname,endptr,lab);
  if (parse_code == 0)
    parse_code = parse_constellation(*endptr,endptr,lab);
  return(parse_code);
}

int parse_catalog_star(char *starname,char **endptr,starlabel_type *lab)
{
  if (lab->kind == HD_STAR || lab->kind == HR_STAR)
    return(parse_catalog_number(starname+2,endptr,lab));
  else
    return(parse_catalog_number(starname+3,endptr,lab));
}

int parse_multiple_component(char *comptr,starlabel_type *lab)
{
  clear_lab_comp(lab);
  while (*comptr == ' ') comptr++;
  if (*comptr == 0) return(0);
  if (!alpha_ok(*comptr)) return(BAD_MULTI_COMP);
  *lab->comp = *(comptr++) & '\xDF';
  while (*comptr == ' ') comptr++;
  if (*comptr != 0) return(MULTI_COMP_END_STR);
  return(0);
}

int parse_base_name(char *starname,char **endptr,starlabel_type *lab)
{
  char newname[MAX_STAR_NAME+1];

  compact_star_name(starname,newname,MAX_STAR_NAME);
  if (!alphanum_ok(*newname)) return(STAR_NAME_ALPHANUM);
  clear_label(lab);
  lab->kind = get_starname_kind(newname);
  switch (lab->kind)
  {
    case BYR_STAR:
      return(parse_bayer_star(newname,endptr,lab));
    case FLM_STAR:
      return(parse_flamsteed_star(newname,endptr,lab));
    default:
      return(parse_catalog_star(newname,endptr,lab));
  }
}

int parse_star_name(char *starname,starlabel_type *lab)
{
  int parse_code;
  char *chrptr;

  parse_code = parse_base_name(starname,&chrptr,lab);
  if (parse_code == 0)
    parse_code = parse_multiple_component(chrptr,lab);
  return(parse_code);
}

void bad_star_name(char *starname,int parse_code)
{
  get_error("%s in '%s'!",bad_star[parse_code],starname);
}

void parse_star_name_err(char *starname,starlabel_type *lab)
{
  int parse_code;

  parse_code = parse_star_name(starname,lab);
  if (parse_code != 0) bad_star_name(starname,parse_code);
}

/* ---------------------------- Star catalogs  ---------------------------- */

int32_t load_data_long(byte **dptr)
{
  int32_t x;

  x = *((int32_t *)(*dptr));
  (*dptr) += 4;

  return(x);
}

int load_data_int(byte **dptr)
{
  return((int)load_data_long(dptr));
}

double load_data_double(byte **dptr)
{
  double x;

  x = *((double *)(*dptr));
  (*dptr) += 8;

  return(x);
}

void load_data_string(char *a,int n,byte **dptr)
{
  memmove(a,*dptr,n);
  (*dptr) += n;
}

void load_bright_data(brighttype *buff,byte *data,int n)
{
  int i;
  brighttype *bptr;
  byte *dptr;

  dptr = data;
  for (i=0,bptr=buff;i<n;i++,bptr++)
  {
    bptr->hr = load_data_int(&dptr);
    bptr->flam = load_data_int(&dptr);
    load_data_string(bptr->bayer,4,&dptr);
    bptr->seq = load_data_int(&dptr);
    load_data_string(bptr->cons,4,&dptr);
    bptr->hd = load_data_int(&dptr);
    load_data_string(bptr->comp,4,&dptr);
    bptr->valid = load_data_int(&dptr);
  }
}

void load_ccdm_data(ccdmtype *buff,byte *data,int n)
{
  int i;
  ccdmtype *bptr;
  byte *dptr;

  dptr = data;
  for (i=0,bptr=buff;i<n;i++,bptr++)
  {
    load_data_string(bptr->ccdm,12,&dptr);
    load_data_string(bptr->comp,4,&dptr);
    bptr->hd = load_data_int(&dptr);
    bptr->dblhd = load_data_int(&dptr);
  }
}

void load_barbier_data(barbiertype *buff,byte *data,int n)
{
  int i;
  barbiertype *bptr;
  byte *dptr;

  dptr = data;
  for (i=0,bptr=buff;i<n;i++,bptr++)
  {
    bptr->seq = load_data_int(&dptr);
    load_data_string(bptr->ccdm,12,&dptr);
    load_data_string(bptr->comp,4,&dptr);
    bptr->hip = load_data_int(&dptr);
    bptr->rv = load_data_double(&dptr);
  }
}

void load_hip_data(hiptype *buff,byte *data,int n)
{
  int i;
  hiptype *bptr;
  byte *dptr;

  dptr = data;
  for (i=0,bptr=buff;i<n;i++,bptr++)
  {
    bptr->hip = load_data_int(&dptr);
    bptr->ra = load_data_double(&dptr);
    bptr->dec = load_data_double(&dptr);
    bptr->par = load_data_double(&dptr);
    bptr->mura = load_data_double(&dptr);
    bptr->mudec = load_data_double(&dptr);
    load_data_string(bptr->aref,4,&dptr);
    load_data_string(bptr->ccdm,12,&dptr);
    bptr->ncomp = load_data_int(&dptr);
    load_data_string(bptr->annex,4,&dptr);
    bptr->hd = load_data_int(&dptr);
    load_data_string(bptr->comp,4,&dptr);
    bptr->valid = load_data_int(&dptr);
  }
}

void load_bright(void)
{
  file_type f;
  int32_t flen,blen;
  byte *data;

  if (bright != NULL)
    get_error("load_bright: Bright stars already loaded!");

  blen = BRIGHT_COUNT * sizeof(brighttype);
  flen = BRIGHT_COUNT * BRIGHT_RECORD;

  open_file(&f,"r","%s/astro/bright-stars",hrsp_dir());
  get_fseek(&f,0,SEEK_END);
  if (ftell(f.dat) != flen)
    get_error("load_bright: Unexpected file size for 'bright-stars'!");
  get_fseek(&f,0,SEEK_SET);
  data = (byte *)get_space(flen);
  get_fread(data,1,flen,&f);
  get_fclose(&f);
  bright = (brighttype *)get_space(blen);
  load_bright_data(bright,data,BRIGHT_COUNT);
  get_free(data);
}

void free_bright(void)
{
  get_free((void *)bright);
  bright = NULL;
}

int bright_loaded(void)
{
  return(bright != NULL);
}

void load_ccdm(void)
{
  file_type f;
  int32_t flen,blen;
  byte *data;

  if (ccdm != NULL)
    get_error("load_ccdm: CCDM stars already loaded!");

  blen = CCDM_COUNT * sizeof(ccdmtype);
  flen = CCDM_COUNT * CCDM_RECORD;

  open_file(&f,"r","%s/astro/ccdm-stars",hrsp_dir());
  get_fseek(&f,0,SEEK_END);
  if (ftell(f.dat) != flen)
    get_error("load_ccdm: Unexpected file size for 'ccdm-stars'!");
  get_fseek(&f,0,SEEK_SET);
  data = (byte *)get_space(flen);
  get_fread(data,1,flen,&f);
  get_fclose(&f);
  ccdm = (ccdmtype *)get_space(blen);
  load_ccdm_data(ccdm,data,CCDM_COUNT);
  get_free(data);
}

void free_ccdm(void)
{
  get_free((void *)ccdm);
  ccdm = NULL;
}

int ccdm_loaded(void)
{
  return(ccdm != NULL);
}

void load_barbier(void)
{
  file_type f;
  int32_t flen,blen;
  byte *data;

  if (barbier != NULL)
    get_error("load_barbier: Barbier-Brossat stars already loaded!");

  blen = BARBIER_COUNT * sizeof(barbiertype);
  flen = BARBIER_COUNT * BARBIER_RECORD;

  open_file(&f,"r","%s/astro/barbier-brossat",hrsp_dir());
  get_fseek(&f,0,SEEK_END);
  if (ftell(f.dat) != flen)
    get_error("load_barbier: Unexpected file size for 'barbier-brossat'!");
  get_fseek(&f,0,SEEK_SET);
  data = (byte *)get_space(flen);
  get_fread(data,1,flen,&f);
  get_fclose(&f);
  barbier = (barbiertype *)get_space(blen);
  load_barbier_data(barbier,data,BARBIER_COUNT);
  get_free(data);
}

void free_barbier(void)
{
  get_free((void *)barbier);
  barbier = NULL;
}

int barbier_loaded(void)
{
  return(barbier != NULL);
}

void load_hip_main(void)
{
  file_type f;
  int32_t flen,blen;
  byte *data;

  if (hip_main != NULL)
    get_error("load_hip_main: Hipparcos stars already loaded!");

  blen = HIP_MAIN_COUNT * sizeof(hiptype);
  flen = HIP_MAIN_COUNT * HIP_RECORD;

  open_file(&f,"r","%s/astro/hip-main",hrsp_dir());
  get_fseek(&f,0,SEEK_END);
  if (ftell(f.dat) != flen)
    get_error("load_hip_main: Unexpected file size for 'hip-main'!");
  get_fseek(&f,0,SEEK_SET);
  data = (byte *)get_space(flen);
  get_fread(data,1,flen,&f);
  get_fclose(&f);
  hip_main = (hiptype *)get_space(blen);
  load_hip_data(hip_main,data,HIP_MAIN_COUNT);
  get_free(data);
}

void free_hip_main(void)
{
  get_free((void *)hip_main);
  hip_main = NULL;
}

void load_hip_com(void)
{
  file_type f;
  int32_t flen,blen;
  byte *data;

  if (hip_com != NULL)
    get_error("load_hip_com: Hipparcos components already loaded!");

  blen = HIP_COM_COUNT * sizeof(hiptype);
  flen = HIP_COM_COUNT * HIP_RECORD;

  open_file(&f,"r","%s/astro/hip-com",hrsp_dir());
  get_fseek(&f,0,SEEK_END);
  if (ftell(f.dat) != flen)
    get_error("load_hip_com: Unexpected file size for 'hip-com'!");
  get_fseek(&f,0,SEEK_SET);
  data = (byte *)get_space(flen);
  get_fread(data,1,flen,&f);
  get_fclose(&f);
  hip_com = (hiptype *)get_space(blen);
  load_hip_data(hip_com,data,HIP_COM_COUNT);
  get_free(data);
}

void free_hip_com(void)
{
  get_free((void *)hip_com);
  hip_com = NULL;
}

void make_comp_list(void)
{
  int i,mhip,hip;
  hiptype *hptr;
  char *comp;

  if (comp_list != NULL)
    get_error("make_comp_list: Hipparcos component list already created!");
  if (comp_buff != NULL)
    get_error("make_comp_list: Hipparcos component buffer not empty!");
  for (i=mhip=0,hptr=hip_main;i<HIP_MAIN_COUNT;i++,hptr++)
    if (hptr->hip > mhip) mhip = hptr->hip;
  comp_list = (char **)get_space((mhip+1)*sizeof(char **));
  comp_buff = (char *)get_space((mhip+1)*8);
  memset(comp_list,0,(mhip+1)*sizeof(char **));
  memset(comp_buff,0,(mhip+1)*8);

  for (hip=0,comp=comp_buff;hip<=mhip;hip++,comp+=8) comp_list[hip] = comp;

  for (i=0,hptr=hip_com;i<HIP_COM_COUNT;i++,hptr++)
    strcat(comp_list[hptr->hip],hptr->comp);
}

void free_comp_list(void)
{
  get_free((void *)comp_list);
  comp_list = NULL;
  get_free((void *)comp_buff);
  comp_buff = NULL;
}

void load_hip(void)
{
  load_hip_main();
  load_hip_com();
  make_comp_list();
}

void free_hip(void)
{
  free_hip_main();
  free_hip_com();
  free_comp_list();
}

int hip_loaded(void)
{
  if (hip_main == NULL) return(0);
  if (hip_com == NULL) return(0);
  if (comp_buff == NULL) return(0);
  if (comp_list == NULL) return(0);

  return(-1);
}

void load_binary(void)
{
  file_type f;
  char row[256],aname[256],bname[256];
  int seq;
  binary_type *x;

  if (binary != NULL) get_error("load_binary: Binary stars already loaded!");

  binary_count = 0;
  open_file(&f,"r","%s/astro/binary.dat",hrsp_dir());
  while (fgets(row,sizeof(row),f.dat) != NULL) binary_count++;
  binary = (binary_type *)get_space(binary_count*sizeof(binary_type));
  get_fseek(&f,0,SEEK_SET);
  for(seq=0,x=binary;seq<binary_count;seq++,x++)
  {
    memset(x,0,sizeof(binary_type));
    if (fgets(row,sizeof(row),f.dat) == NULL)
      get_error("load_binary: Cannot read data line!");
    if (strlen(row) > 250) get_error("load_binary: Data line too long!");
    if (strlen(row) < 16) get_error("load_binary: Data line too short!");
    if (sscanf(row," %s %s %lf "
      "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
      aname,bname,&x->epoch,
      &x->ra,&x->dec,&x->par,&x->mura,&x->mudec,&x->rv,
      &x->mass_ratio,
      &x->a,&x->i,&x->ome,&x->Ome,&x->e,&x->P,&x->T) != 17)
      get_error("load_binary: Missing data!");
    parse_star_name_err(aname,&x->prim);
    parse_star_name_err(bname,&x->seco);
  }
  get_fclose(&f);
}

void free_binary(void)
{
  get_free((void *)binary);
  binary = NULL;
  binary_count = 0;
}

int binary_loaded(void)
{
  return(binary != NULL);
}

void load_cats(void)
{
  load_hip();
  load_bright();
  load_ccdm();
  load_barbier();
  load_binary();
}

void free_cats(void)
{
  free_hip();
  free_bright();
  free_ccdm();
  free_barbier();
  free_binary();
}

double barbier_rv(int hip,char *comp)
{
  int i;
  barbiertype *bptr;

  if (hip < 1) return(0.0);

  for (i=0,bptr=barbier;i<BARBIER_COUNT;i++,bptr++)
  {
    if (bptr->hip != hip) continue;
    if (*comp != 0) if (strstr(bptr->comp,comp) == NULL) continue;
    return(bptr->rv);
  }

  return(0.0);
}

/* ------------------------- Star identification  ------------------------- */

char *print_star_label(starlabel_type *lab,char *buf)
{
  char seqtxt[8];

  if (lab->seq > 0)
    sprintf(seqtxt,"%d",lab->seq);
  else
    memset(seqtxt,0,8);

  switch(lab->kind)
  {
    case BYR_STAR:
      strcpy(buf,lab->bayer);
      strcat(buf,seqtxt);
      strcat(buf," ");
      strcat(buf,lab->cons);
      break;
    case FLM_STAR:
      sprintf(buf,"%d %s",lab->flam,lab->cons);
      break;
    case HD_STAR:
      sprintf(buf,"HD %d",lab->catnum);
      break;
    case HR_STAR:
      sprintf(buf,"HR %d",lab->catnum);
      break;
    case HIP_STAR:
      sprintf(buf,"HIP %d",lab->catnum);
      break;
    case USR_STAR:
      sprintf(buf,"USR %d",lab->catnum);
      break;
  }
  if (*lab->comp != 0) strcat(buf," ");
  strcat(buf,lab->comp);
  return(buf);
}

char *print_comp_choice(char *comp,char *buf)
{
  int n,i;
  char item[8];

  n = strlen(comp);
  *buf = 0;
  for (i=0;i<n;i++)
  {
    if (i > 0)
    {
      if (i < n-1)
        sprintf(item,", ");
      else
        sprintf(item," or ");
      strcat(buf,item);
    }
    sprintf(item,"%c",comp[i]);
    strcat(buf,item);
  }
  return(buf);
}

void check_bsc_entry(brighttype *fptr,int n,starlabel_type *lab)
{
  char btxt[80];

  print_star_label(lab,btxt);

  if (n == 0) get_error("check_bsc_entry: Unknown star '%s'!",btxt);

  if (n > 1) get_error
    ("check_bsc_entry: Ambiguous star '%s' (%d entries found)!",btxt,n);

  if (lab->kind == BYR_STAR && lab->seq == 0) lab->seq = fptr->seq;
}

void check_hip_entry(hiptype *fptr,int n,starlabel_type *lab)
{
  char stxt[80],ctxt[80];

  print_star_label(lab,stxt);

  if (n < 1 || fptr == NULL) get_error
    ("check_hip_entry: Unknown star '%s'!",stxt);

  if (n > 1) get_error
    ("check_hip_entry: Ambiguous star '%s' (%d entries found)!",stxt,n);

  print_comp_choice(comp_list[fptr->hip],ctxt);

  if (*lab->comp == 0)
  {
    if (strcmp(fptr->annex,"C") == 0 && fptr->ncomp > 1)
      get_error("check_hip_entry: Ambiguous star '%s' "
        "(%s component required)!",stxt,ctxt);
  }
  else
  {
    if (strchr(comp_list[fptr->hip],(int)*lab->comp) == NULL)
      get_error("check_hip_entry: Invalid component '%s'!",stxt);
  }
}

void bayer_to_hd(starlabel_type *lab)
{
  int i,n;
  brighttype *bptr,*fptr;
  char stxt[80];

  print_star_label(lab,stxt);

  if (lab->kind != BYR_STAR)
    get_error("bayer_to_hd: Bayer designation expected (%s)!",stxt);

  if (!bright_loaded())
    get_error("bayer_to_hd: Bright stars not loaded!");

  for (i=n=0,bptr=bright,fptr=NULL;i<BRIGHT_COUNT;i++,bptr++)
  {
    if (!bptr->valid) continue;
    if (strcmp(lab->bayer,bptr->bayer) != 0) continue;
    if (strcmp(lab->cons,bptr->cons) != 0) continue;
    if (lab->seq != 0 && lab->seq != bptr->seq) continue;
    if (*lab->comp != 0 && strchr(bptr->comp,(int)*lab->comp) == NULL)
      continue;
    fptr = bptr;
    n++;
  }

  check_bsc_entry(fptr,n,lab);

  lab->kind = HD_STAR;
  lab->catnum = fptr->hd;
}

void flamsteed_to_hd(starlabel_type *lab)
{
  int i,n;
  brighttype *bptr,*fptr;
  char stxt[80];

  print_star_label(lab,stxt);

  if (lab->kind != FLM_STAR)
    get_error("flamsteed_to_hd: Flamsteed designation expected (%s)!",stxt);

  if (!bright_loaded())
    get_error("flamsteed_to_hd: Bright stars not loaded!");

  for (i=n=0,bptr=bright,fptr=NULL;i<BRIGHT_COUNT;i++,bptr++)
  {
    if (!bptr->valid) continue;
    if (lab->flam != bptr->flam) continue;
    if (strcmp(lab->cons,bptr->cons) != 0) continue;
    if (*lab->comp != 0 && strchr(bptr->comp,(int)(*lab->comp)) == NULL)
      continue;
    fptr = bptr;
    n++;
  }

  check_bsc_entry(fptr,n,lab);

  lab->kind = HD_STAR;
  lab->catnum = fptr->hd;
}

void hr_to_hd(starlabel_type *lab)
{
  brighttype *fptr;
  char stxt[80];

  print_star_label(lab,stxt);

  if (lab->kind != HR_STAR)
    get_error("hr_to_hd: HR designation expected (%s)!",stxt);

  if (lab->catnum < 1 || lab->catnum > BRIGHT_COUNT)
    get_error("hr_to_hd: HR number out of range (%s)!",stxt);

  if (!bright_loaded())
    get_error("hr_to_hd: Bright stars not loaded!");

  fptr = &bright[lab->catnum - 1];

  check_bsc_entry(fptr,1,lab);

  lab->kind = HD_STAR;
  lab->catnum = fptr->hd;
}

void hip_to_hip(starlabel_type *lab)
{
  int i,nh;
  hiptype *hptr,*fptr;
  char stxt[80];

  if (lab->hip_ok) return;

  print_star_label(lab,stxt);

  if (lab->kind != HIP_STAR)
    get_error("hip_to_hip: HIP designation expected (%s)!",stxt);

  if (strlen(lab->comp) > 1)
    get_error("hip_to_hip: Too many components (%s)!",stxt);

  if (!hip_loaded())
    get_error("hip_to_hip: Hipparcos catalogue not loaded!");

  for (i=nh=0,hptr=hip_main,fptr=NULL;i<HIP_MAIN_COUNT;i++,hptr++)
  {
    if (!hptr->valid) continue;
    if (hptr->hip != lab->catnum) continue;
    if (*lab->comp != 0 &&
      strchr(comp_list[hptr->hip],(int)*lab->comp) == NULL) continue;
    fptr = hptr, nh++;
  }

  check_hip_entry(fptr,nh,lab);

  if (strlen(comp_list[fptr->hip]) < 2) *lab->comp = 0;

  lab->hip_ok = -1;
}

void hd_to_hip(starlabel_type *lab)
{
  int i,nc,nhc,nh;
  ccdmtype *cptr,*gptr;
  hiptype *hptr,*fptr;
  char stxt[80];

  print_star_label(lab,stxt);

  if (lab->kind != HD_STAR)
    get_error("hd_to_hip: HD designation expected (%s)!",stxt);

  if (strlen(lab->comp) > 1)
    get_error("hd_to_hip: Too many components (%s)!",stxt);

  if (!ccdm_loaded())
    get_error("hd_to_hip: CCDM catalogue not loaded!");

  if (!hip_loaded())
    get_error("hd_to_hip: Hipparcos catalogue not loaded!");

  nc = nhc = nh = 0;
  cptr = gptr = NULL;
  hptr = fptr = NULL;

  for (i=nc=0,cptr=ccdm,gptr=NULL;i<CCDM_COUNT;i++,cptr++)
  {
    if (cptr->hd < 1) continue;
    if (cptr->hd != lab->catnum)
    {
      if (!cptr->dblhd) continue;
      if (cptr->hd != (lab->catnum+1)) continue;
    }
    if (*lab->comp != 0 && *lab->comp != cptr->comp[1])
      get_error("hd_to_hip: Unexpected component (%s)!",stxt);
    gptr = cptr, nc++;
  }

  if (nc > 1) get_error
    ("hd_to_hip: Ambigous star in CCDM catalogue (%s)!",stxt);

  if (nc == 1)
  {
    for (i=nhc=0,hptr=hip_com,fptr=NULL;i<HIP_COM_COUNT;i++,hptr++)
    {
      if (!hptr->valid) continue;
      if (strcmp(hptr->ccdm,gptr->ccdm) != 0) continue;
      if (hptr->comp[0] != gptr->comp[1]) continue;
      fptr = hptr, nhc++;
    }

    if (nhc > 1) get_error
      ("hd_to_hip: Ambigous component %s%s in Hipparcos catalogue (%s)!",
      gptr->ccdm,gptr->comp,stxt);
    if (nhc == 1)
    {
      check_hip_entry(fptr,nhc,lab);
      lab->kind = HIP_STAR;
      lab->catnum = fptr->hip;
      *lab->comp = *fptr->comp;
      hip_to_hip(lab);
      return;
    }
  }

  for (i=nh=0,hptr=hip_main,fptr=NULL;i<HIP_MAIN_COUNT;i++,hptr++)
  {
    if (!hptr->valid) continue;
    if (hptr->hd < 1) continue;
    if (hptr->hd != lab->catnum) continue;
    if (*lab->comp != 0 &&
      strchr(comp_list[hptr->hip],(int)*lab->comp) == NULL) continue;
    fptr = hptr, nh++;
  }
  check_hip_entry(fptr,nh,lab);
  lab->kind = HIP_STAR;
  lab->catnum = fptr->hip;
  hip_to_hip(lab);
}

void any_to_hip(starlabel_type *lab)
{
  switch(lab->kind)
  {
    case FLM_STAR:
      flamsteed_to_hd(lab);
      hd_to_hip(lab);
      break;
    case BYR_STAR:
      bayer_to_hd(lab);
      hd_to_hip(lab);
      break;
    case HR_STAR:
      hr_to_hd(lab);
      hd_to_hip(lab);
      break;
    case HD_STAR:
      hd_to_hip(lab);
      break;
    case HIP_STAR:
      hip_to_hip(lab);
      break;
    default: get_error("any_to_hip: Invalid star label type!");
  }
}

int same_star(starlabel_type *a,starlabel_type *b)
{
  starlabel_type p,q;

  if (a->kind == USR_STAR || b->kind == USR_STAR)
  {
    if (a->kind != USR_STAR || b->kind != USR_STAR) return(0);
    return(a->catnum == b->catnum);
  }

  p = *a;
  q = *b;

  any_to_hip(&p);
  any_to_hip(&q);

  return(p.catnum == q.catnum && *p.comp == *q.comp);
}

void load_hip_star(starlabel_type *lab,star_type *star)
{
  int i,n,k;
  hiptype *hptr,*fptr,*cptr;
  binary_type *bptr;
  datetime_type pass;
  interval_type per;

  if (lab->kind != HIP_STAR)
    get_error("load_hip_star: HIP designation expected!");

  if (!binary_loaded())
    get_error("load_hip_star: Binary stars not loaded!");

  clear_star(star);

  if (binary_count > 0)
  {
    for (i=0,bptr=binary;i<binary_count;i++,bptr++)
    {
      k = 0;
      if (same_star(lab,&bptr->prim))
        k = 1;
      else
        if (same_star(lab,&bptr->seco)) k = 2;
      if (k > 0)
      {
        set_interval_years(&per,bptr->P);
        set_besselian_year(&pass,bptr->T);
        star->ref = ICRS;
        star->cent = BARY_CENT;
        star->pos = GEOMETRIC;
        set_julian_year_tdt(&star->equinox,2000.0);
        set_julian_year_tdt(&star->epoch,bptr->epoch);
        set_angle_deg(&star->ra,bptr->ra);
        set_angle_deg(&star->dec,bptr->dec);
        star->par = bptr->par * 1.0e-3;
        star->mura = bptr->mura * 1.0e-3;
        star->mudec = bptr->mudec * 1.0e-3;
        star->rv = bptr->rv;
        star->mass_ratio = bptr->mass_ratio;
        set_orbit
          (&star->orb,bptr->a,bptr->i,bptr->ome,bptr->Ome,bptr->e,&per,&pass);
        star->component = CENTER_OF_MASS;
        space_vectors(star);
        if (k == 1)
          select_primary_star(star);
        else
          select_secondary_star(star);
        sky_angles(star);
        return;
      }
    }
  }

  if (!hip_loaded())
    get_error("load_hip_star: Hipparcos catalogue not loaded!");

  for (i=n=0,hptr=hip_main,fptr=NULL;i<HIP_MAIN_COUNT;i++,hptr++)
  {
    if (!hptr->valid) continue;
    if (hptr->hip != lab->catnum) continue;
    if (*lab->comp != 0)
      if (strchr(comp_list[hptr->hip],(int)*lab->comp) == NULL) continue;
    fptr = hptr;
    n++;
  }

  check_hip_entry(fptr,n,lab);

  if (strcmp(fptr->annex,"C") == 0 && *lab->comp != 0)
  {
    cptr = NULL;
    for (i=n=0,hptr=hip_com;i<HIP_COM_COUNT;i++,hptr++)
    {
      if (strcmp(fptr->ccdm,hptr->ccdm) != 0) continue;
      if (fptr->hip != hptr->hip) continue;
      if (strcmp(lab->comp,hptr->comp) != 0) continue;
      cptr = hptr;
      n++;
    }
    check_hip_entry(cptr,n,lab);
    fptr = cptr;
  }

  star->ref = ICRS;
  star->cent = BARY_CENT;
  star->pos = GEOMETRIC;
  set_julian_year_tdt(&star->equinox,2000.0);
  set_julian_year_tdt(&star->epoch,1991.25);
  set_angle_deg(&star->ra,fptr->ra);
  set_angle_deg(&star->dec,fptr->dec);
  star->par = fptr->par * 1.0e-3;
  star->mura = fptr->mura * 1.0e-3;
  star->mudec = fptr->mudec * 1.0e-3;
  star->rv = barbier_rv(fptr->hip,fptr->comp);
  star->mass_ratio = 0.0;
  star->component = CENTER_OF_MASS;
}

void load_usr_star(starlabel_type *lab,star_type *star)
{
  file_type file;
  char row[1024];
  int nvals;
  int usr,seq,found;
  char eqtxt[24];
  int ref;
  double eq,ep;
  int rah,ram;
  double ras;
  char dect[16],decz;
  int decd,decm;
  double decs;
  double par,mura,mudec,rv;

  if (lab->kind != USR_STAR)
    get_error("load_usr_star: USR designation expected!");

  clear_star(star);

  set_file_name(&file,"stars.dat");
  get_fopen(&file,"r");
  usr = 1;
  found = 0;
  while (fgets(row,sizeof(row),file.dat) != NULL)
  {
    if (strlen(row) > 1000)
      get_error("Catalogue line too long in '%s'!",file.path);
    nvals = sscanf(row,"%d %s %lf %d %d %lf %s %d %lf %lf %lf %lf %lf",
      &seq,eqtxt,&ep,&rah,&ram,&ras,dect,&decm,&decs,&par,&mura,&mudec,&rv);
    if (nvals != 13)
      get_error("13 data values per line expected in '%s'!",file.path);
    if (seq != usr)
      get_error("Sequence number out of order in '%s'!",file.path);
    if (seq == lab->catnum) found = 1;
    if (found) break;
    usr++;
  }
  get_fclose(&file);

  if (!found) get_error("Star not found in '%s'!",file.path);

  if (strcasecmp(eqtxt,"ICRS") == 0)
  {
    ref = ICRS;
    eq = 2000.0;
  }
  else
  {
    ref = MEAN_EQ;
    eq = str_to_dbl_def(eqtxt,0.0);
    if (eq == 0.0) get_error("Invalid catalogue equinox in '%s'!",file.path);
  }

  if (*dect == '-')
    decz = '-';
  else
    decz = '+';
  decd = abs(atoi(dect));

  star->ref = ref;
  star->cent = BARY_CENT;
  star->pos = GEOMETRIC;
  set_julian_year_tdt(&star->equinox,eq);
  set_julian_year_tdt(&star->epoch,ep);
  set_angle_alpha(&star->ra,rah,ram,ras);
  set_angle_delta(&star->dec,decz,decd,decm,decs);
  star->par = par * 1.0e-3;
  star->mura = mura * 1.0e-3;
  star->mudec = mudec * 1.0e-3;
  star->rv = rv;
  star->mass_ratio = 0.0;
  star->component = CENTER_OF_MASS;
}

void load_star(starlabel_type *lab,star_type *star)
{
  switch(lab->kind)
  {
    case HIP_STAR: load_hip_star(lab,star);
                   break;
    case USR_STAR: load_usr_star(lab,star);
                   break;
    default: get_error("load_star: Invalid star label kind!");
  }
}
