#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <fitsio2.h>

#include "utils.h"

#include "psf.h"

psf_str *psf_create(char *filename)
{
    psf_str *psf = (psf_str *)malloc(sizeof(psf_str));
    fitsfile *fits = NULL;
    int status = 0;
    int tmp_i = 0;

    if(!filename || !file_exists_and_normal(filename)){
        dprintf("Can't read PSFex file: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    fits_open_file(&fits, filename, READONLY, &status);
    fits_movabs_hdu(fits, 2, &tmp_i, &status);

    fits_read_key(fits, TDOUBLE, "POLZERO1", &psf->x0, NULL, &status);
    fits_read_key(fits, TDOUBLE, "POLSCAL1", &psf->sx, NULL, &status);
    fits_read_key(fits, TDOUBLE, "POLZERO2", &psf->y0, NULL, &status);
    fits_read_key(fits, TDOUBLE, "POLSCAL2", &psf->sy, NULL, &status);

    fits_read_key(fits, TINT, "POLDEG1", &psf->degree, NULL, &status);
    fits_read_key(fits, TDOUBLE, "PSF_SAMP", &psf->pix_step, NULL, &status);
    fits_read_key(fits, TDOUBLE, "PSF_FWHM", &psf->fwhm, NULL, &status);

    fits_read_key(fits, TINT, "PSFAXIS1", &psf->width, NULL, &status);
    fits_read_key(fits, TINT, "PSFAXIS2", &psf->height, NULL, &status);
    fits_read_key(fits, TINT, "PSFAXIS3", &psf->ncoeffs, NULL, &status);

    psf->C = malloc(sizeof(double)*psf->ncoeffs*psf->width*psf->height);
    fits_read_col(fits, TDOUBLE, 1, 1, 1, psf->ncoeffs*psf->width*psf->height, NULL, psf->C, NULL, &status);

    fits_close_file(fits, &status);

    if(status != 0){
        dprintf("Can't parse PSFex file: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    /* Shift the origin by 1 pixel, as PSFEx origin is 1,1 */
    psf->x0 -= 1;
    psf->y0 -= 1;

    return psf;
}

void psf_delete(psf_str *psf)
{
    free(psf);
}

inline double psf_sampled_value(psf_str *psf, int x, int y, int x0, int y0)
{
    double dx = (x0 - psf->x0)/psf->sx;
    double dy = (y0 - psf->y0)/psf->sy;
    int i1;
    int i2;
    double value = 0;
    double pi2 = 1;
    int idx =  + y*psf->width + x;

    for(i2 = 0; i2 <= psf->degree; i2++){
        double pi1 = 1;

        for(i1 = 0; i1 <= psf->degree - i2; i1++){
            value += pi1*pi2*psf->C[idx];
            pi1 *= dx;
            idx += psf->width*psf->height;
        }

        pi2 *= dy;
    }

    return value;
}

image_str *psf_sampled_image(psf_str *psf, double x0, double y0)
{
    image_str *image = image_create_double(psf->width, psf->height);
    int x;
    int y;

    for(x = 0; x < psf->width; x++)
        for(y = 0; y < psf->height; y++)
            PIXEL_DOUBLE(image, x, y) = psf_sampled_value(psf, x, y, x0, y0);

    return image;
}

/* SExtractor stuff starts */

#define         INTERPW		8	/* Interpolation function range */
#define PI      	3.1415926535898

static const double INTERPFAC = 3.0;
static const double IINTERPFAC = .3333333333333333333333333333;

#define	INTERPF(x)	(x<1e-5 && x>-1e-5? 1.0 \
			:(x>INTERPFAC?0.0:(x<-INTERPFAC?0.0 \
			:sinf(PI*x)*sinf(PI/INTERPFAC*x)/(PI*PI/INTERPFAC*x*x))))
				/* Lanczos approximation */

static int _psfex_vignet_resample(double *pix1, int w1, int h1,
                                  double *pix2, int w2, int h2,
                                  double dx, double dy, double step2)
{
    double	*mask,*maskt, xc1,xc2,yc1,yc2, xs1,ys1, x1,y1, x,y, dxm,dym,
            val, norm,
            *pix12, *pixin,*pixin0, *pixout,*pixout0;
    int		i,j,k,n,t, *start,*startt, *nmask,*nmaskt,
            ixs2,iys2, ix2,iy2, dix2,diy2, nx2,ny2, iys1a, ny1, hmw,hmh,
            ix,iy, ix1,iy1;


    /* Initialize destination buffer to zero */
    memset(pix2, 0, w2*h2*sizeof(double));

    xc1 = (double)(w1/2);	/* Im1 center x-coord*/
    xc2 = (double)(w2/2);	/* Im2 center x-coord*/
    xs1 = xc1 + dx - xc2*step2;	/* Im1 start x-coord */

    if ((int)xs1 >= w1)
        return -1;
    ixs2 = 0;			/* Int part of Im2 start x-coord */
    if (xs1<0.0)
    {
        dix2 = (int)(1-xs1/step2);
        /*-- Simply leave here if the images do not overlap in x */
        if (dix2 >= w2)
            return -1;
        ixs2 += dix2;
        xs1 += dix2*step2;
    }
    nx2 = (int)((w1-1-xs1)/step2+1);/* nb of interpolated Im2 pixels along x */
    if (nx2>(ix2=w2-ixs2))
        nx2 = ix2;
    if (nx2<=0)
        return -1;
    yc1 = (double)(h1/2);	/* Im1 center y-coord */
    yc2 = (double)(h2/2);	/* Im2 center y-coord */
    ys1 = yc1 + dy - yc2*step2;	/* Im1 start y-coord */
    if ((int)ys1 >= h1)
        return -1;
    iys2 = 0;			/* Int part of Im2 start y-coord */
    if (ys1<0.0)
    {
        diy2 = (int)(1-ys1/step2);
        /*-- Simply leave here if the images do not overlap in y */
        if (diy2 >= h2)
            return -1;
        iys2 += diy2;
        ys1 += diy2*step2;
    }
    ny2 = (int)((h1-1-ys1)/step2+1);/* nb of interpolated Im2 pixels along y */
    if (ny2>(iy2=h2-iys2))
        ny2 = iy2;
    if (ny2<=0)
        return -1;

    /* Set the yrange for the x-resampling with some margin for interpolation */
    iys1a = (int)ys1;		/* Int part of Im1 start y-coord with margin */
    hmh = INTERPW/2 - 1;		/* Interpolant start */
    if (iys1a<0 || ((iys1a -= hmh)< 0))
        iys1a = 0;
    ny1 = (int)(ys1+ny2*step2)+INTERPW-hmh;	/* Interpolated Im1 y size */
    if (ny1>h1)					/* with margin */
        ny1 = h1;
    /* Express everything relative to the effective Im1 start (with margin) */
    ny1 -= iys1a;
    ys1 -= (double)iys1a;

    /* Allocate interpolant stuff for the x direction */
    if ((mask = (double *) malloc(sizeof(double) * nx2 * INTERPW)) == NULL) /* Interpolation masks */
        return -1;
    if ((nmask = (int *) malloc(sizeof(int) * nx2)) == NULL) /* Interpolation mask sizes */
        return -1;
    if ((start = (int *) malloc(sizeof(int) * nx2)) == NULL) /* Int part of Im1 conv starts */
        return -1;

    /* Compute the local interpolant and data starting points in x */
    hmw = INTERPW/2 - 1;
    x1 = xs1;
    maskt = mask;
    nmaskt = nmask;
    startt = start;
    for (j=nx2; j--; x1+=step2)
    {
        ix = (ix1=(int)x1) - hmw;
        dxm = ix1 - x1 - hmw;	/* starting point in the interpolation func */
        if (ix < 0)
        {
            n = INTERPW+ix;
            dxm -= (double)ix;
            ix = 0;
        }
        else
            n = INTERPW;
        if (n>(t=w1-ix))
            n=t;
        *(startt++) = ix;
        *(nmaskt++) = n;
        norm = 0.0;
        for (x=dxm, i=n; i--; x+=1.0)
            norm += (*(maskt++) = INTERPF(x));
        norm = norm>0.0? 1.0/norm : 1.0;
        maskt -= n;
        for (i=n; i--;)
            *(maskt++) *= norm;
    }

    if ((pix12 = (double *) calloc(nx2*ny1, sizeof(double))) == NULL) { /* Intermediary frame-buffer */
        return -1;
    }

    /* Make the interpolation in x (this includes transposition) */
    pixin0 = pix1+iys1a*w1;
    pixout0 = pix12;
    for (k=ny1; k--; pixin0+=w1, pixout0++)
    {
        maskt = mask;
        nmaskt = nmask;
        startt = start;
        pixout = pixout0;
        for (j=nx2; j--; pixout+=ny1)
        {
            pixin = pixin0+*(startt++);
            val = 0.0;
            for (i=*(nmaskt++); i--;)
                val += *(maskt++)**(pixin++);
            *pixout = val;
        }
    }

    /* Reallocate interpolant stuff for the y direction */
    if ((mask = (double *) realloc(mask, sizeof(double) * ny2 * INTERPW)) == NULL) { /* Interpolation masks */
        return -1;
    }
    if ((nmask = (int *) realloc(nmask, sizeof(int) * ny2)) == NULL) { /* Interpolation mask sizes */
        return -1;
    }
    if ((start = (int *) realloc(start, sizeof(int) * ny2)) == NULL) { /* Int part of Im1 conv starts */
        return -1;
    }

    /* Compute the local interpolant and data starting points in y */
    hmh = INTERPW/2 - 1;
    y1 = ys1;
    maskt = mask;
    nmaskt = nmask;
    startt = start;
    for (j=ny2; j--; y1+=step2)
    {
        iy = (iy1=(int)y1) - hmh;
        dym = iy1 - y1 - hmh;	/* starting point in the interpolation func */
        if (iy < 0)
        {
            n = INTERPW+iy;
            dym -= (double)iy;
            iy = 0;
        }
        else
            n = INTERPW;
        if (n>(t=ny1-iy))
            n=t;
        *(startt++) = iy;
        *(nmaskt++) = n;
        norm = 0.0;
        for (y=dym, i=n; i--; y+=1.0)
            norm += (*(maskt++) = INTERPF(y));
        norm = norm>0.0? 1.0/norm : 1.0;
        maskt -= n;
        for (i=n; i--;)
            *(maskt++) *= norm;
    }

    /* Make the interpolation in y  and transpose once again */
    pixin0 = pix12;
    pixout0 = pix2+ixs2+iys2*w2;
    for (k=nx2; k--; pixin0+=ny1, pixout0++)
    {
        maskt = mask;
        nmaskt = nmask;
        startt = start;
        pixout = pixout0;
        for (j=ny2; j--; pixout+=w2)
        {
            pixin = pixin0+*(startt++);
            val = 0.0;
            for (i=*(nmaskt++); i--;)
                val += *(maskt++)**(pixin++);
            *pixout = val;
        }
    }

    /* Free memory */
    free(pix12);
    free(mask);
    free(nmask);
    free(start);

    return 0;
}

/* SExtractor stuff ends */

void psf_image_fill(psf_str *psf, image_str *sampled, image_str *image, double dx, double dy)
{
    _psfex_vignet_resample(sampled->double_data, sampled->width, sampled->height,
                           image->double_data, image->width, image->height,
                           -dx/psf->pix_step, -dy/psf->pix_step, 1.0/psf->pix_step);
}

image_str *psf_image(psf_str *psf, image_str *sampled, double dx, double dy, int size)
{
    image_str *image = image_create_double(size, size);

    psf_image_fill(psf, sampled, image, dx, dy);

    return image;
}

psf_str *psf_create_from_fits_and_save(char *filename, char *outname)
{
    char *dirname = make_temp_dirname("/tmp/psfex_XXXXXX");

    char *empty_filename = make_string("%s/default.empty", dirname);
    char *filter_filename = make_string("%s/default.filter", dirname);
    char *param_filename = make_string("%s/default.params", dirname);
    char *cat_filename = make_string("%s/out.psf.cat", dirname);
    char *psf_name = make_string("%s/out.psf.psf", dirname);

    psf_str *psf = NULL;

    FILE *file = NULL;

    /* SExtractor part */
    file = fopen(empty_filename, "w");
    fclose(file);

    file = fopen(filter_filename, "w");
    fprintf(file, "CONV NORM\n"
            "# 5x5 convolution mask of a gaussian PSF with FWHM = 2.0 pixels.\n"
            "0.006319 0.040599 0.075183 0.040599 0.006319\n"
            "0.040599 0.260856 0.483068 0.260856 0.040599\n"
            "0.075183 0.483068 0.894573 0.483068 0.075183\n"
            "0.040599 0.260856 0.483068 0.260856 0.040599\n"
            "0.006319 0.040599 0.075183 0.040599 0.006319\n");
    fclose(file);

    file = fopen(param_filename, "w");
    fprintf(file,
            "VIGNET(40,40)\n"
            "X_IMAGE\nY_IMAGE\n"
            "FLUX_RADIUS\nFLUX_APER\nFLUXERR_APER\n"
            "ELONGATION\nFLAGS\nFLUX_MAX\nSNR_WIN\n");
    fclose(file);

    system_run_silently("sex %s -c %s"
                        " -CATALOG_TYPE FITS_LDAC -PARAMETERS_NAME %s"
                        " -FILTER Y -FILTER_NAME %s"
                        " -CLEAN Y -VERBOSE_TYPE QUIET"
                        " -DETECT_THRESH 3"
                        " -CATALOG_NAME %s",
                        filename, empty_filename, param_filename, filter_filename, cat_filename);

    /* PSFEx part */
    system_run_silently("psfex %s -c %s -OUTCAT_NAME %s"
                        " -PSFVAR_DEGREES 4 -PSF_RECENTER Y -PSF_SIZE 25,25 -PSFVAR_NSNAP 20"
                        " -NTHREADS 0 -CHECKIMAGE_TYPE NONE -VERBOSE_TYPE QUIET -WRITE_XML N",
                        cat_filename, empty_filename, psf_name);

    if(file_exists_and_normal(psf_name)){
        psf = psf_create(psf_name);

        if(psf && outname)
            /* psf_name contains valid PSF file - we shall copy it to outname */
            copy_file(psf_name, outname);
    }

    remove_dir(dirname);

    free(psf_name);
    free(cat_filename);
    free(param_filename);
    free(filter_filename);
    free(empty_filename);

    if(!psf)
        exit_with_error("Error creating PSF for %s", filename);

    return psf;
}

psf_str *psf_create_from_fits(char *filename)
{
    return psf_create_from_fits_and_save(filename, NULL);
}
