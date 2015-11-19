#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "setup.h"

#include "utils.h"
#include "coords.h"
/* #include "objects.h" */

void coords_get_ra_dec_from_xi_eta(coords_str *coords, double xi, double eta,
                                   double *ra_ptr, double *dec_ptr)
{
    double sin_dec0 = sin(coords->dec0*M_PI/180);
    double cos_dec0 = cos(coords->dec0*M_PI/180);

    double delta_ra;
    double ra;
    double dec;

    /* For FITS WCS standard the xi/eta is in degrees, not radians */
    xi *= M_PI/180;
    eta *= M_PI/180;

    delta_ra = atan2(xi, cos_dec0 - sin_dec0*eta);
    ra = mod(coords->ra0 + delta_ra*180./M_PI, 360);

    dec = atan2((sin_dec0 + eta*cos_dec0)*cos(delta_ra),
                  cos_dec0 - eta*sin_dec0) * 180./M_PI;

    /* FIXME: Polar exceptions */
    if(dec < -90){
        dec = dec + 180;
    } else if(dec > 90){
        dec = dec - 180;
    }

    if(ra_ptr)
        *ra_ptr = ra;
    if(dec_ptr)
        *dec_ptr = dec;
}

/* x, y -> ra, dec
 * ra in hours, dec in degrees */
void coords_get_ra_dec(coords_str *coords, double x, double y,
                       double *ra_ptr, double *dec_ptr)
{
    double xi;
    double eta;

    coords_get_xi_eta(coords, x, y, &xi, &eta);

    coords_get_ra_dec_from_xi_eta(coords, xi, eta, ra_ptr, dec_ptr);
}

inline double coords_get_ra(coords_str *coords, double x, double y)
{
    double ra = 0;

    coords_get_ra_dec(coords, x, y, &ra, NULL);

    return ra;
}

inline double coords_get_dec(coords_str *coords, double x, double y)
{
    double dec = 0;

    coords_get_ra_dec(coords, x, y, NULL, &dec);

    return dec;
}

/* Angular distance in degrees between points on the sky */
/* implementation from q3c */
inline double coords_sky_distance(double ra1, double dec1, double ra2, double dec2)
{
    double x, y, z;

    x = sin((ra1 - ra2)/2*M_PI/180.0);
    x *= x;
    y = sin((dec1 - dec2)/2*M_PI/180.0);
    y *= y;

    /* Seem to be more precise :) */
    z = cos((dec1 + dec2)/2*M_PI/180.0);
    z*=z;

    return 2*asin(sqrt(x*(z - y) + y))*180.0/M_PI;
}

/* Angular distance in degrees between points on frame */
double coords_distance(coords_str *coords, double x1, double y1,
                       double x2, double y2)
{
    double ra1 = 0;
    double ra2 = 0;
    double dec1 = 0;
    double dec2 = 0;

    coords_get_ra_dec(coords, x1, y1, &ra1, &dec1);
    coords_get_ra_dec(coords, x2, y2, &ra2, &dec2);

    return coords_sky_distance(ra1, dec1, ra2, dec2);
}

void coords_get_x_y_from_xi_eta(coords_str *coords, double xi, double eta, double *x_ptr, double *y_ptr)
{
    double xx = (xi*coords->CD22 - eta*coords->CD12)/
        (coords->CD11*coords->CD22 - coords->CD21*coords->CD12);
    double yy = (eta*coords->CD11 - xi*coords->CD21)/
        (coords->CD11*coords->CD22 - coords->CD21*coords->CD12);

    if(coords->AP_ORDER){
        double dx = 0;
        double dy = 0;
        int i;
        int j;

        for(i = 0; i < MAX_SIP_ORDER; i++)
            for(j = 0; j < MAX_SIP_ORDER; j++){
                if(i + j <= coords->AP_ORDER)
                    dx += coords->AP[i + j*MAX_SIP_ORDER]*pow(xx, i)*pow(yy, j);
                if(i + j <= coords->BP_ORDER)
                    dy += coords->BP[i + j*MAX_SIP_ORDER]*pow(xx, i)*pow(yy, j);
            }

        xx += dx;
        yy += dy;
    }

    xx += coords->CRPIX1;
    yy += coords->CRPIX2;

    if(x_ptr)
        *x_ptr = xx;
    if(y_ptr)
        *y_ptr = yy;
}

void coords_get_xi_eta_from_ra_dec(coords_str *coords, double ra, double dec, double *xi_ptr, double *eta_ptr)
{
    double xi = 0;
    double eta = 0;
    double dec_center = coords->dec0*M_PI/180;
    double delta_ra;
    double xx;
    double yy;

    if((ra < 10) && (coords->ra0 > 350))
        delta_ra = (ra + 360) - coords->ra0;
    else if ((ra > 350) && (coords->ra0 < 10))
        delta_ra = (ra - 360) - coords->ra0;
    else
        delta_ra = ra - coords->ra0;

    delta_ra *= M_PI/180;

    xx = cos(dec*M_PI/180)*sin(delta_ra);
    yy = sin(dec_center)*sin(dec*M_PI/180) +
        cos(dec_center)*cos(dec*M_PI/180)*cos(delta_ra);
    xi = (xx/yy);

    xx = cos(dec_center)*sin(dec*M_PI/180) -
        sin(dec_center)*cos(dec*M_PI/180)*cos(delta_ra);
    eta = (xx/yy);

    /* For FITS WCS standard the xi/eta is in degrees, not radians */
    xi *= 180./M_PI;
    eta *= 180./M_PI;

    if(xi_ptr)
        *xi_ptr = xi;
    if(eta_ptr)
        *eta_ptr = eta;
}

void coords_get_xi_eta(coords_str *coords, double x, double y, double *xi_ptr, double *eta_ptr)
{
    double xi = coords_get_xi(coords, x, y);
    double eta = coords_get_eta(coords, x, y);

    if(xi_ptr)
        *xi_ptr = xi;
    if(eta_ptr)
        *eta_ptr = eta;
}

inline double coords_get_xi(coords_str *coords, double x, double y)
{
    double xx = x - coords->CRPIX1;
    double yy = y - coords->CRPIX2;

    if(coords->A_ORDER){
        double dx = 0;
        double dy = 0;
        int i;
        int j;

        for(i = 0; i < MAX_SIP_ORDER; i++)
            for(j = 0; j < MAX_SIP_ORDER; j++){
                if(i + j <= coords->A_ORDER)
                    dx += coords->A[i + j*MAX_SIP_ORDER]*pow(xx, i)*pow(yy, j);
                if(i + j <= coords->B_ORDER)
                    dy += coords->B[i + j*MAX_SIP_ORDER]*pow(xx, i)*pow(yy, j);
            }

        xx += dx;
        yy += dy;
    }

    return coords->CD11*xx + coords->CD12*yy;
}

inline double coords_get_eta(coords_str *coords, double x, double y)
{
    double xx = x - coords->CRPIX1;
    double yy = y - coords->CRPIX2;

    if(coords->A_ORDER){
        double dx = 0;
        double dy = 0;
        int i;
        int j;

        for(i = 0; i < MAX_SIP_ORDER; i++)
            for(j = 0; j < MAX_SIP_ORDER; j++){
                if(i + j <= coords->A_ORDER)
                    dx += coords->A[i + j*MAX_SIP_ORDER]*pow(xx, i)*pow(yy, j);
                if(i + j <= coords->B_ORDER)
                    dy += coords->B[i + j*MAX_SIP_ORDER]*pow(xx, i)*pow(yy, j);
            }

        xx += dx;
        yy += dy;
    }

    return coords->CD21*xx + coords->CD22*yy;
}

/* ra, dec -> x, y */
void coords_get_x_y(coords_str *coords, double ra, double dec, double *x_ptr, double *y_ptr)
{
    double xi = 0;
    double eta = 0;
    double x;
    double y;

    coords_get_xi_eta_from_ra_dec(coords, ra, dec, &xi, &eta);

    coords_get_x_y_from_xi_eta(coords, xi, eta, &x, &y);

    if(x_ptr)
        *x_ptr = x;
    if(y_ptr)
        *y_ptr = y;
}

inline double coords_get_x(coords_str *coords, double ra, double dec)
{
    double x = 0;

    coords_get_x_y(coords, ra, dec, &x, NULL);

    return x;
}

inline double coords_get_y(coords_str *coords, double ra, double dec)
{
    double y = 0;

    coords_get_x_y(coords, ra, dec, NULL, &y);

    return y;
}

void coords_clear(coords_str *coords)
{
    coords->CRPIX1 = 0;
    coords->CD11 = 0;
    coords->CD12 = 0;

    coords->CRPIX2 = 0;
    coords->CD21 = 0;
    coords->CD22 = 0;

    coords->ra0 = 0;
    coords->dec0 = 0;

    coords->A_ORDER = 0;
    coords->B_ORDER = 0;
    coords->AP_ORDER = 0;
    coords->BP_ORDER = 0;

    memset(coords->A, 0, sizeof(double)*(MAX_SIP_ORDER + 1)*(MAX_SIP_ORDER + 1));
    memset(coords->B, 0, sizeof(double)*(MAX_SIP_ORDER + 1)*(MAX_SIP_ORDER + 1));
    memset(coords->AP, 0, sizeof(double)*(MAX_SIP_ORDER + 1)*(MAX_SIP_ORDER + 1));
    memset(coords->BP, 0, sizeof(double)*(MAX_SIP_ORDER + 1)*(MAX_SIP_ORDER + 1));

    coords->filter = 0;
    /* coords->mag_scale = 1; */
    /* coords->mag_scale_err = 0; */
    /* coords->mag0 = 0; */
    /* coords->mag0_err = 0; */
    /* coords->mag_covar = 0; */
    /* coords->mag_sigma = 0; */

    memset(coords->mag_C, 0, sizeof(double)*6);
}

coords_str coords_empty()
{
    coords_str coords;

    coords_clear(&coords);

    return coords;
}

int coords_is_empty(coords_str *coords)
{
    if(!coords->CRPIX1 &&
       !coords->CD11 &&
       !coords->CD12 &&
       !coords->CRPIX2 &&
       !coords->CD21 &&
       !coords->CD22)
        return TRUE;
    else
        return FALSE;
}

double coords_x0(coords_str *coords)
{
    return coords_get_x(coords, coords->ra0, coords->dec0);
}

double coords_y0(coords_str *coords)
{
    return coords_get_y(coords, coords->ra0, coords->dec0);
}

double coords_get_pixscale(coords_str *coords)
{
    double delta = 0.5*(sqrt(coords->CD11*coords->CD11 + coords->CD21*coords->CD21) +
                        sqrt(coords->CD12*coords->CD12 + coords->CD22*coords->CD22));
    double scale = delta;

    return scale;
}

/* Rotate the coordinates */
void coords_get_ra_dec_shifted(double ra0, double dec0, double delta_ra, double delta_dec, double *ra_ptr, double *dec_ptr)
{
    coords_str coords;

    coords.ra0 = ra0;
    coords.dec0 = dec0;

    /* FIXME: is this correct? */
    coords_get_ra_dec_from_xi_eta(&coords, delta_ra, delta_dec, ra_ptr, dec_ptr);
}

/* Get the rotation between coordinates */
void coords_get_ra_dec_shift(double ra0, double dec0, double ra, double dec, double *delta_ra_ptr, double *delta_dec_ptr)
{
    coords_str coords;

    coords.ra0 = ra0;
    coords.dec0 = dec0;

    /* FIXME: is this correct? */
    coords_get_xi_eta_from_ra_dec(&coords, ra, dec, delta_ra_ptr, delta_dec_ptr);
}

/* instr mag -> catalogue mag */
/* Error is systematic error of this calibration, to be added to the flux error */
double coords_get_mag(coords_str *coords, double instr, double x, double y, double *err_ptr)
{
    if(err_ptr)
        *err_ptr = 0;//sqrt(coords->mag_scale_err*coords->mag_scale_err*instr*instr + coords->mag0_err*coords->mag0_err + 2*coords->mag_covar*instr);

    //return coords->mag_scale*instr + coords->mag0;

    return instr + coords->mag_C[0] + x*x*coords->mag_C[1] + x*coords->mag_C[2] + x*y*coords->mag_C[3] + y*coords->mag_C[4] + y*y*coords->mag_C[5];
}
