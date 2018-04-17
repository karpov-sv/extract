#ifndef EXTRACT_H
#define EXTRACT_H

#include "lists.h"
#include "image.h"

#define PEAK_INITIAL 0x1
#define PEAK_MEASURED 0x2
#define PEAK_FAILED 0x10

/* Flags for peaks and for the mask frame */
#define FLAG_NORMAL 0x0
#define FLAG_SATURATED 0x1
#define FLAG_TRUNCATED 0x2
#define FLAG_DEBLENDED 0x4
#define FLAG_ELONGATED 0x8
#define FLAG_UNMEASURED 0x10
#define FLAG_BAD 0x100

/* Single peak */
typedef struct peak_str {
    LIST_HEAD(struct peak_str);

    int id;

    double x;
    double dx;

    double y;
    double dy;

    double bg;
    double dbg;

    /* Measured model flux */
    double flux;
    double dflux;

    double minflux;
    double maxflux;

    /* Initial shape parameters */
    double a;
    double b;
    double theta;

    /* Some PSF statistics */
    double params[10];

    double fwhm;
    double ellipticity;
    double beta;
    double egg;
    double beta2;

    double excess;

    double chisq;

    double ra;
    double dec;

    int flags;
    int state;
} peak_str;

void find_peaks(image_str *, image_str *, image_str *, image_str *, double , struct list_head *);
void load_peaks(char *, image_str *, image_str *, image_str *, int , struct list_head *);

void dump_peaks_to_file(struct list_head *, char *, int );

#endif /* EXTRACT_H */
