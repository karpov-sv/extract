#ifndef TIME_STR_H
#define TIME_STR_H

#include <sys/types.h>

typedef struct time_str {
    u_int32_t year __attribute__ ((packed));
    u_int32_t month __attribute__ ((packed));
    u_int32_t day __attribute__ ((packed));
    u_int32_t hour __attribute__ ((packed));
    u_int32_t minute __attribute__ ((packed));
    u_int32_t second __attribute__ ((packed));
    u_int32_t microsecond __attribute__ ((packed));
} time_str;

/* Time interval in msec */
int64_t time_interval(time_str , time_str );
time_str time_current();
time_str time_zero();
int time_is_zero(time_str );
time_str time_make(int , int , int , int , int , int , int );

/* time_str to string */
char *time_str_get_date(time_str );
char *time_str_get_time(time_str );
char *time_str_get_date_time(time_str );
char *time_str_get_date_time_static(time_str );

char *time_str_get_evening_date(time_str );
time_str time_str_get_morning_time(time_str );

/* Julianic date with one day precision (always *.5) */
long double time_str_get_JD(time_str );
time_str time_str_from_JD(long double );

/* Standard unix time representation (one second precision) */
time_t time_unix(time_str );
time_str time_from_unix(time_t );
void time_increment(time_str *, double );
time_str time_incremented(time_str , double );

/* String to time_str */
void time_str_set_date(time_str *, char *);
void time_str_set_time(time_str *, char *);
time_str time_str_from_date_time(char *);
time_str time_str_from_filename(char *);

/* Static timestamp - no need to free() */
const char *timestamp();

/* Sidereal time - Greenwich and local */
double time_str_get_sidereal_time(time_str );
double time_str_get_local_sidereal_time(time_str , double );

/* UUID, sort of */
u_int64_t time_str_get_uuid(time_str );
time_str time_str_from_uuid(u_int64_t);
time_str time_str_from_uuid_string(char *);

#endif /* TIME_STR_H */
