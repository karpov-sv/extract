#ifndef EXTRACT_H
#define EXTRACT_H

#include "lists.h"
#include "image.h"

#define PEAK_INITIAL 0
#define PEAK_MEASURED 1
#define PEAK_FAILED 100

/* Single peak */
typedef struct peak_str {
    LIST_HEAD(struct peak_str);

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

    int state;
} peak_str;

void find_peaks(image_str *, image_str *, struct list_head *);

void dump_peaks_to_file(struct list_head *, char *);
void dump_peaks_to_file_measured(struct list_head *, char *);
void dump_peaks_to_file_failed(struct list_head *, char *);
void dump_peaks_to_file_full(struct list_head *, char *);

#endif /* EXTRACT_H */
