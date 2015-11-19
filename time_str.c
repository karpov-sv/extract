#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>

#include "setup.h"
#include "time_str.h"
#include "utils.h"

static int days_in_months[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

static __thread char *timestamp_string = NULL;
static __thread char *date_time_string = NULL;

/* Time interval in millisecods */
int64_t time_interval(time_str t1, time_str t2)
{
    long long res = 0;

    res = (time_str_get_JD(t2) - time_str_get_JD(t1))*1000*60*60*24;
    /*
    res += t2.millisecond - t1.millisecond;
    res += (t2.second - t1.second)*1000;
    res += (t2.minute - t1.minute)*1000*60;
    res += (t2.hour - t1.hour)*1000*60*60;
    res += (time_str_get_JD(t2) - time_str_get_JD(t1))*1000*60*60*24;
    */

    return res;
}

/* Return current time with one second precision */
time_str time_current()
{
    struct tm tm;
    struct timeval tv;
    time_str time;

    gettimeofday(&tv, NULL);

    gmtime_r(&tv.tv_sec, &tm);

    time.year = tm.tm_year + 1900;
    time.month = tm.tm_mon + 1;
    time.day = tm.tm_mday;
    time.hour = tm.tm_hour;
    time.minute = tm.tm_min;
    time.second = tm.tm_sec;
    time.microsecond = tv.tv_usec;

    return time;
}

/* Return zero time */
time_str time_zero()
{
    time_str time = {1900, 1, 1, 0, 0, 0, 0};

    return time;
}

int time_is_zero(time_str time)
{
    if(!time_interval(time, time_zero()))
        return TRUE;
    else
        return FALSE;
}

time_str time_make(int year, int month, int day, int hour, int minute, int second, int microsecond)
{
    time_str time = {year, month, day, hour, minute, second, microsecond};

    return time;
}

/* Format DD.MM.YYYY */
char *time_str_get_date(time_str time)
{
    return make_string("%02d.%02d.%04d", time.day, time.month, time.year);
}

/* format HH:MM:SS.mmm */
char *time_str_get_time(time_str time)
{
    return make_string("%02d:%02d:%02d.%03d", time.hour, time.minute,
                       time.second, time.microsecond/1000);
}

/* Format DD.MM.YYYY HH:MM:SS.mmm */
char *time_str_get_date_time(time_str time)
{
    return make_string("%02d.%02d.%04d %02d:%02d:%02d.%03d",
                       time.day, time.month, time.year,
                       time.hour, time.minute, time.second, time.microsecond/1000);
}

char *time_str_get_date_time_static(time_str time)
{
    if(date_time_string)
        free(date_time_string);

    date_time_string = time_str_get_date_time(time);

    return date_time_string;
}

/* format DD_MM_YYYY, for last evening - for logs etc */
/* This date is the same from 12am till 12am next day in LOCAL time */
char *time_str_get_evening_date(time_str time_in)
{
    time_str time = time_in;

    tzset();
    /* Convert to local time */
    time_increment(&time, -timezone);
    /* Shift the time back accordingly to let the local time prior to 12am next
     * day be on the same date as the start of observations */
    time_increment(&time, -12*60*60);

    return make_string("%04d_%02d_%02d", time.year, time.month, time.day);
}

/* Returns time of the end of observations */
time_str time_str_get_morning_time(time_str time_in)
{
    time_str time = time_in;

    time.minute = time.second = time.microsecond = 0;

    if(time.hour > 12)
        time.day ++;

    time.hour = MORNING_TIME_HOUR;

    if(time.day > days_in_months[time.month - 1]){
        if(time.month < 11)
            time.month ++;
        else {
            time.month = 1;
            time.year ++;
        }

        if(time.month == 3 && time.day == 29 &&
           (time.year % 4 == 0 && (time.year % 100 != 0 ||  time.year % 400 == 0))){
            time.day = 29;
            time.month = 2;
        } else
            time.day = 1;
    }

    return time;
}

/*
  Julianic date, with ~2 millisecond precision
 */

/* Stupid math.h doesn't contain it... */
double trunc(double );

long double time_str_get_JD(time_str time)
{
    /* Version adapted from libAstronomy */
    double year = time.year;
    double month = time.month;
    double day = time.day;
    double hour = time.hour;
    double minute = time.minute;
    double second = time.second;
    double millisecond = time.microsecond/1000;
    double a;
    int b;

    a = 10000.0*year + 100.0*time.month + time.day;
    if (month <= 2) {
        month = month + 12;
        year = year - 1;
    }

    if (a <= 15821004.1)
        b = -2 + (int)((year + 4716)/4) - 1179;
    else
        b = (int)(year/400) - (int)(year/100) + (int)(year/4);
    a = 365.0*year - 679004.0;

    return (a + b + (int)(30.6001*(month + 1)) + day + (hour + minute/60.0 + second/3600.0 + millisecond/3600.0/1000.0)/24.0) + 2400000.5;

    /* This version has some kind of 1-month systematic error
    double a, b, c, e, f;
    double yr, mo, day;

    if (time.month < 3) {
        yr = (double)(time.year - 1);
        mo = (double)(time.month + 12);
    }
    else {
        yr = (double)time.year;
        mo = (double)(time.month);
    }
    day = (double)time.day;

    a = trunc(yr / 100);
    b = trunc(a / 4);
    c = 2 - a + b;
    e = trunc(365.25 * (yr + 4716));
    f = trunc(30.6001*mo);

    return c + day + e + f - 1524.5;
    */

    /*
    long int jd12h;
    double tjd;

    jd12h = (long) time.day - 32075L + 1461L * ((long) time.year + 4800L
                                           + ((long) time.month - 14L) / 12L) / 4L
        + 367L * ((long) time.month - 2L - ((long) time.month - 14L) / 12L * 12L)
        / 12L - 3L * (((long) time.year + 4900L + ((long)time. month - 14L) / 12L)
                      / 100L) / 4L;
    tjd = (double) jd12h - 0.5 + time.hour / 24.0;

    return (tjd);
    */
}

time_str time_str_from_JD(long double tjd)
{
    time_str time = time_zero();
    long int k, m, n;
    long double djd = tjd + 0.5;
    long int jd = (long int) djd;

    time.hour = floorl( fmodl (djd, 1.0) * 24.0);
    time.minute = floorl( fmodl (djd*24.0, 1.0) * 60);
    time.second = floorl( fmodl (djd*24.0*60.0, 1.0) * 60);
    time.microsecond = floorl( fmodl (djd*24.0*60.0*60, 1.0) * 1e6);

    k     = jd + 68569L;
    n     = 4L * k / 146097L;

    k     = k - (146097L * n + 3L) / 4L;
    m     = 4000L * (k + 1L) / 1461001L;
    k     = k - 1461L * m / 4L + 31L;

    time.month = (short int) (80L * k / 2447L);
    time.day   = (short int) (k - 2447L * (long int) time.month / 80L);
    k      = (long int) time.month / 11L;

    time.month = (short int) ((long int) time.month + 2L - 12L * k);
    time.year  = (short int) (100L * (n - 49L) + m + k);

    return time;
}

time_t time_unix(time_str time)
{
    time_t unix_time;
    struct tm t;

    t.tm_sec = time.second;
    t.tm_min = time.minute;
    t.tm_hour = time.hour;
    t.tm_mday = time.day;
    t.tm_mon = time.month - 1;
    t.tm_year = time.year - 1900;
    t.tm_isdst = -1;
#ifndef __CYGWIN__
    t.tm_gmtoff = 0;
#endif

    unix_time = mktime(&t);

    return unix_time;
}

time_str time_from_unix(time_t t)
{
    struct tm tm;
    time_str time;

    localtime_r(&t, &tm);

    time.year = tm.tm_year + 1900;
    time.month = tm.tm_mon + 1;
    time.day = tm.tm_mday;
    time.hour = tm.tm_hour;
    time.minute = tm.tm_min;
    time.second = tm.tm_sec;
    time.microsecond = 0;

    return time;
}

/* increments time by delta sec. */
/* FIXME: works for dates after the Epoch only! */
void time_increment(time_str *time, double delta)
{
    int idelta = floor(delta);
    double usec_new = time->microsecond + 1e6*(delta - idelta);
    time_str new;

    while(usec_new > 1e6){
        usec_new -= 1e6;
        idelta ++;
    }

    new = time_from_unix(time_unix(*time) + idelta);

    new.microsecond = usec_new;

    *time = new;
}

time_str time_incremented(time_str time, double delta)
{
    time_str copy = time;

    time_increment(&copy, delta);

    return copy;
}

/* Format DD.MM.YYYY or DD/MM/YYYY or DD-MM-YYYY */
void time_str_set_date(time_str *time, char *string)
{
    char *tmp = make_string("%s", string);
    char *pos = tmp;
    int day = 0;
    int month = 0;
    int year = 0;

    /* replace all dots with spaces */
    while(*pos){
        if(*pos == '.' || *pos == '/' || *pos == '-')
            *pos = ' ';
        pos ++;
    }

    sscanf(tmp, "%d %d %d", &day, &month, &year);

    if(day > 31){
        int tmp = day;

        day = year;
        year = tmp;
    }

    time->year = year;
    time->month = month;
    time->day = day;

    free(tmp);
}

/* format HH:MM:SS.mmm */
void time_str_set_time(time_str *time, char *string)
{
    int hour = 0;
    int min = 0;
    double sec = 0;

    if(string)
        sscanf(string, "%2d:%2d:%lf", &hour, &min, &sec);

    time->hour = hour;
    time->minute = min;
    time->second = floor(sec);
    time->microsecond = 1000000*(sec - time->second);
}

/* Format DD.MM.YYYY HH:MM:SS.mmm */
time_str time_str_from_date_time(char *string)
{
    time_str time = time_zero();

    if(string){
        char *tmp = make_string("%s", string);
        char *pos = tmp;
        int day = 0;
        int month = 0;
        int year = 0;
        int hour = 0;
        int min = 0;
        int sec = 0;
        int msec = 0;

        /* replace all dots with spaces */
        while(*pos){
            if(*pos == '.' || *pos == ':' || *pos == '-')
                *pos = ' ';
            pos ++;
        }

        /* FIXME: here we incorrectly parse milliseconds if not all three digits are provided
           - e.g. 10.12 instead of 10.120 */
        sscanf(tmp, "%d %d %d %d %d %d %d",
               &day, &month, &year, &hour, &min, &sec, &msec);

        if(day > 31){
            int tmp = day;

            day = year;
            year = tmp;
        }

        time.year = year;
        time.month = month;
        time.day = day;
        time.hour = hour;
        time.minute = min;
        time.second = sec;
        time.microsecond = msec*1000;

        free(tmp);
    }

    return time;
}

const char *timestamp()
{
    if(timestamp_string)
        free(timestamp_string);

    timestamp_string = time_str_get_date_time(time_current());

    return timestamp_string;
}

double time_str_get_sidereal_time(time_str time)
{
    double mjd = time_str_get_JD(time) - 2400000.5;
    double mjd0, t, ut, g;

    mjd0 = (double)(long)mjd;
    ut = (mjd - mjd0)*24.0;
    t = (mjd0 - 51544.5)/36525.0;
    g = (6.697374558 + 1.0027379093*ut +
         (8640184.812866 + (0.093104 - 6.2E-6*t)*t)*t/3600.0);

    return 24.0*mod(g/24.0, 1);
}

double time_str_get_local_sidereal_time(time_str time, double lambda)
{
    return 24.0*mod((time_str_get_sidereal_time(time) - lambda/15.0)/24.0, 1);
}

time_str time_str_from_filename(char *filename)
{
    time_str time = time_zero();

    sscanf(filename, "%4d%2d%2d_%2d%2d%2d",
           &time.year, &time.month, &time.day, &time.hour, &time.minute, &time.second);

    return time;
}

/* u_int64_t of YYYYMMDDHHMMSSMMM form - sort of unique ID for frames */
u_int64_t time_str_get_uuid(time_str time)
{
    u_int64_t result = 0;

    result = result*1000 + time.year;
    result = result*100 + time.month;
    result = result*100 + time.day;

    result = result*100 + time.hour;
    result = result*100 + time.minute;
    result = result*100 + time.second;

    result = result*1000 + time.microsecond/1000;

    return result;
}

time_str time_str_from_uuid(u_int64_t uuid)
{
    time_str time = time_zero();

    time.microsecond = 1000*(uuid % 1000);
    uuid /= 1000;

    time.second = uuid % 100;
    uuid /= 100;
    time.minute = uuid % 100;
    uuid /= 100;
    time.hour = uuid % 100;
    uuid /= 100;

    time.day = uuid % 100;
    uuid /= 100;
    time.month = uuid % 100;
    uuid /= 100;
    time.year = uuid % 10000;
    uuid /= 1000;

    return time;
}

time_str time_str_from_uuid_string(char *string)
{
    u_int64_t uuid = 0;

    sscanf(string, "%017lld", &uuid);

    return time_str_from_uuid(uuid);
}
