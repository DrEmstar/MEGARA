#ifndef DATETIME_DEFINED

  #define DATETIME_DEFINED

  typedef struct
  {
    char full[12];
    char abbr[4];
  }
  monthtype;

  typedef struct
  {
    int day,month,year;
    int32_t jd;
    double jd0;
  }
  date_type;

  typedef struct
  {
    int hour,min,adj;
    double sec;
    double full_day,decsec,frac,rad;
  }
  time_type;

  typedef struct
  {
    date_type date;
    time_type time;
    double jd0,jdf,jd,djd,jcen,jyr;
  }
  datetime_type;

  typedef struct
  {
    double julian_cents;
    double julian_years;
    double days;
    double seconds;
  }
  interval_type;

#endif

char *full_month_name(int m);
char *short_month_name(int m);
int month_number(char *mname);
int32_t gregorian_date_to_jd(int d,int m,int y);
int32_t julian_date_to_jd(int d,int m,int y);
int32_t date_to_jd(int d,int m,int y);
void jd_to_gregorian_date(int32_t jd,int *d,int *m,int *y);
void jd_to_julian_date(int32_t jd,int *d,int *m,int *y);
void jd_to_date(int32_t jd,int *d,int *m,int *y);
void set_date_dmy(date_type *u,int d,int m,int y);
void set_date_jd(date_type *u,int32_t jd);
void add_days(date_type *a,date_type *b,int d);
void next_day(date_type *a,date_type *b);
void set_time_decsec_adj(time_type *t,double decsec,int adj);
void set_time_decsec(time_type *t,double decsec);
void set_time_hms_adj(time_type *t,int hour,int min,double sec,int adj);
void set_time_hms(time_type *t,int hour,int min,double sec);
void set_time_frac_adj(time_type *t,double frac,int adj);
void set_time_frac(time_type *t,double frac);
void update_time_args(datetime_type *u);
void set_datetime_adj(datetime_type *u,int day,int month,int year,
  int hour,int min,double dsec,int adj);
void set_datetime(datetime_type *u,int day,int month,int year,
  int hour,int min,double dsec);
void split_jd(double jda,double jdb,double *jd0,double *jdf);
void set_datetime_jd(datetime_type *u,double jda,double jdb);
void set_julian_year(datetime_type *u,double year);
void set_besselian_year(datetime_type *u,double year);
void shift_datetime(datetime_type *a,datetime_type *b,double sec);
int compare_datetime(datetime_type *a,datetime_type *b);
void set_interval_years(interval_type *t,double years);
void set_interval_days(interval_type *t,double days);
