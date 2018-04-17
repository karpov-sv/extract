#ifndef PSF_H
#define PSF_H

#include "image.h"

int psf_var_degree;
int psf_vignet_size;
int psf_size;
int psf_subtract_bg;
double psf_gain;
double psf_satur_level;
double psf_phot_aper;

typedef struct {
    /* Scaling */
    double x0;
    double y0;
    double sx;
    double sy;

    double z0;
    double sz;

    double pix_step;
    double fwhm;

    int degree;
    int zdegree;
    int width;
    int height;
    int ncoeffs;

    double *C;
} psf_str;

psf_str *psf_create(char *);
void psf_delete(psf_str *);

inline double psf_sampled_value(psf_str *, int , int , double , double , double );
image_str *psf_sampled_image_3d(psf_str *, double , double , double );
image_str *psf_sampled_image(psf_str *, double , double );

void psf_image_fill(psf_str *psf, image_str *, image_str *, double , double );
image_str *psf_image(psf_str *, image_str *, double , double , int );

psf_str *psf_create_from_fits_and_save(char *, char *, char *);
psf_str *psf_create_from_fits(char *, char *);

image_str *image_smooth_psf(image_str *, psf_str *);
image_str *image_unsharp_psf(image_str *, psf_str *);

double image_ellipticity(image_str *);
double image_fwhm(image_str *);

int psf_is_positive(psf_str *, image_str *);

int psf_photometry(image_str *, char *, char *, int );

#endif /* PSF_H */
