#ifndef EXTRACT_H
#define EXTRACT_H

#include "lists.h"
#include "image.h"

#define PEAK_INITIAL 0x1
#define PEAK_MEASURED 0x2
#define PEAK_FAILED 0x10

/* Flags for peaks and for the mask frame */
#define FLAG_NORMAL 0x0
#define FLAG_SATURATED 0x2
#define FLAG_TRUNCATED 0x4
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

    double A; /* Peak amplitude */

    double bg;
    double dbg;

    /* Measured model flux */
    double flux;
    double dflux;

    /* Initial shape parameters */
    double a;
    double b;
    double theta;

    double excess;

    double chisq;

    int flags;
    int state;
} peak_str;

void find_peaks(image_str *, image_str *, image_str *, image_str *, double , struct list_head *);

void dump_peaks_to_file(struct list_head *, char *, int );

#endif /* EXTRACT_H */
