#ifndef COORDS_H
#define COORDS_H

#include "lists.h"
#include "time_str.h"

#define MAX_SIP_ORDER 3

typedef struct {
    double CRPIX1;
    double CD11;
    double CD12;
    double CRPIX2;
    double CD21;
    double CD22;

    /* SIP polynomial */
    int A_ORDER;
    double A[(MAX_SIP_ORDER + 1)*(MAX_SIP_ORDER + 1)];
    int B_ORDER;
    double B[(MAX_SIP_ORDER + 1)*(MAX_SIP_ORDER + 1)];
    int AP_ORDER;
    double AP[(MAX_SIP_ORDER + 1)*(MAX_SIP_ORDER + 1)];
    int BP_ORDER;
    double BP[(MAX_SIP_ORDER + 1)*(MAX_SIP_ORDER + 1)];

    double ra0; /* CRVAL1 */
    double dec0; /* CRVAL2 */

    int filter; /* Filter used for photometric calibration, 0 - Clear, 1 - B, 2 - V, 3 - R */
    /* double mag_scale; /\* Cat = mag_scale*Instr + mag0 *\/ */
    /* double mag_scale_err; */
    /* double mag0; /\* Cat = mag_scale*Instr + mag0 *\/ */
    /* double mag0_err; */
    /* double mag_covar; /\* Covariance of errors *\/ */
    double mag_sigma; /* Sttdev of Cat - Instr */

    double mag_C[6];
} coords_str;

/* x, y -> ra, dec */
void coords_get_ra_dec(coords_str *, double , double , double *, double *);
double coords_get_ra(coords_str *, double , double );
double coords_get_dec(coords_str *, double , double );
inline double coords_sky_distance(double , double , double , double );
double coords_distance(coords_str *, double , double , double , double );
/* Spherical rotation */
void coords_get_ra_dec_shifted(double, double, double, double, double *, double *);
void coords_get_ra_dec_shift(double, double, double, double, double *, double *);

/* xi, eta - standard coordinates */
void coords_get_x_y_from_xi_eta(coords_str *, double , double , double *, double *);
void coords_get_ra_dec_from_xi_eta(coords_str *, double , double , double *, double *);
void coords_get_xi_eta_from_ra_dec(coords_str *, double , double , double *, double *);
void coords_get_xi_eta(coords_str *, double , double , double *, double *);
inline double coords_get_xi(coords_str *, double , double );
inline double coords_get_eta(coords_str *, double , double );

/* ra, dec -> x, y */
void coords_get_x_y(coords_str *, double , double , double *, double *);
double coords_get_x(coords_str *, double , double );
double coords_get_y(coords_str *, double , double );

/* instr mag -> catalogue mag */
double coords_get_mag(coords_str *, double , double , double , double *);

void coords_clear(coords_str *);
coords_str coords_empty();
int coords_is_empty(coords_str *);

double coords_x0(coords_str *);
double coords_y0(coords_str *);

double coords_get_pixscale(coords_str *);

#endif /* COORDS_H */
