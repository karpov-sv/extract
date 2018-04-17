#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_sf_erf.h>

#include "utils.h"

#include "image.h"
#include "csa/csa.h"
#include "mpfit.h"

/* /\* Delta arrays - possible steps from given pixel *\/ */
/* static int dx[] = {0, 0,  0, 1, 1,  1, -1, -1, -1}; */
/* static int dy[] = {0, 1, -1, 1, 0, -1,  1,  0, -1}; */
/* static int dN = 9; */

#define IS_PIXEL_VALID(image, saturation, x, y) ((x) >= 0 && (x) < (image)->width && (y) >= 0 && (y) < (image)->height && PIXEL_DOUBLE(image, x, y) < saturation)

/* /\* Get the direction to the nearest local maximum as an index in delta */
/*    arrays *\/ */
/* static inline int get_max_n_idx(image_str *image, double saturation, int x, int y) */
/* { */
/*     int d; */
/*     int max_value = 0; */
/*     int max_id = 0; */

/*     for(d = 0; d < dN; d++){ */
/*         if(IS_PIXEL_VALID(image, saturation, x + dx[d], y + dy[d])){ */
/*             int value = PIXEL_DOUBLE(image, x + dx[d], y + dy[d]); */

/*             if(d == 0 || value > max_value){ */
/*                 max_value = value; */
/*                 max_id = d; */
/*             } */
/*         } */
/*     } */

/*     return max_id; */
/* } */

/* /\* Finds the locally reachable maximum position *\/ */
/* static inline int get_reachable_maximum_xy(image_str *image, double saturation, int x0, int y0, int *x_ptr, int *y_ptr) */
/* { */
/*     int x = x0; */
/*     int y = y0; */
/*     int id = 0; */
/*     int Nsteps = 0; */

/*     while((id = get_max_n_idx(image, saturation, x, y))){ */
/*         x += dx[id]; */
/*         y += dy[id]; */
/*         Nsteps ++; */
/*     } */

/*     if(x_ptr) */
/*         *x_ptr = x; */
/*     if(y_ptr) */
/*         *y_ptr = y; */

/*     return Nsteps; */
/* } */

/* static inline int is_minimum(image_str *image, double saturation, int x, int y) */
/* { */
/*     int d; */
/*     int value = PIXEL_DOUBLE(image, x, y); */
/*     int result = TRUE; */

/*     for(d = 1; d < dN; d++) */
/*         if(IS_PIXEL_VALID(image, saturation, x + dx[d], y + dy[d]) && */
/*            PIXEL_DOUBLE(image, x + dx[d], y + dy[d]) < value){ */
/*             result = FALSE; */
/*             break; */
/*         } */

/*     return result; */
/* } */

/* static inline double get_local_value(image_str *image, double saturation, int x, int y) */
/* { */
/*     int d; */
/*     double sum = 0; */
/*     int N = 0; */

/*     for(d = 0; d < dN; d++) */
/*         if(IS_PIXEL_VALID(image, saturation, x + dx[d], y + dy[d])){ */
/*             sum += PIXEL_DOUBLE(image, x + dx[d], y + dy[d]); */
/*             N += 1; */
/*         } */

/*     return sum/N; */
/* } */

static int cmpfn(const void *a, const void *b)
{
    double res = (*(double*)a - *(double*)b);

    return (res > 0 ? 1 : (res < 0 ? -1 : 0));
}

static double get_P(double I, double bg, double sigma, double a, double p)
{
    /* Skewed Gaussian, convolved with exponential distribution exp(-delta/a)/a */
    /* FIXME: this will not work in low-count limit! Should use Poissonian one */
    double arg1 = (2.0*a*bg - 2*a*I + sigma*sigma)/2/a/a;
    double arg2 = 0.5*(a*bg - a*I + sigma*sigma)*sqrt(2)/(sigma*a);
    double value = exp(arg1 + gsl_sf_log_erfc(arg2))/2/a;

    /* Normal Gaussian */
    double value2 = exp(-(I - bg)*(I - bg)/2/sigma/sigma)/sigma/sqrt(2*M_PI);

    /* if(!isfinite(value)){ */
    /*     dprintf("%g %g %g %g - %g\n", I, bg, sigma, a, value); */
    /*     exit(1); */
    /* } */

    /* Mixture model - Gaussian + skewed Gaussian */
    return value*(1 - p) + value2*p;
}

struct model_str {
    double *I;
    double *H;
    double *dH;
    int nbins;
};

static int fn_mpfit(int Npoints, int Nparams, double *params, double *residuals, double **derivatives, void *data)
{
    struct model_str *model = (struct model_str *)data;
    double S = params[0];
    double sigma = params[1];
    double a = params[2];
    double p = params[3];

    int d;

    for(d = 0; d < model->nbins; d++){
        residuals[d] = (model->H[d] - get_P(model->I[d], S, sigma, a, p))/model->dH[d];
    }

    return 0;
}

static double get_bijaoui_bg(double *data, int N, double *error_ptr)
{
    int nbins = 100;
    double *I = (double*)calloc(nbins, sizeof(double));
    double *H = (double*)calloc(nbins, sizeof(double));
    double *dH = (double*)calloc(nbins, sizeof(double));
    double mad;
    double med = get_median_mad(data, N, &mad);
    double min = med - 3.0*mad;
    double max = med + 5.0*mad;
    int N1 = 0;
    int d;

    double params[4];
    double parerrs[4];
    mp_par paropts[4];

    mp_result result;
    struct model_str model;
    struct mp_config_struct conf;

    memset(&conf, 0, sizeof(conf));
    /* conf.ftol = 1e-14; */
    /* conf.xtol = 1e-14; */
    /* conf.gtol = 1e-16; */
    /* conf.maxiter = 1000; */
    /* conf.douserscale = 1; */

    /* conf.nofinitecheck = 1; */

    /* Bin size should not be smaller than 1-2 ADU */
    if((max - min)/nbins < 2.0){
        nbins = floor((max - min)/2.0);
    }

    for(d = 0; d < nbins; d++)
        I[d] = min + (max - min)*d/nbins;

    for(d = 0; d < N; d++){
        int bin = round((data[d] - min)*nbins/(max - min));

        if(!isfinite(data[d]))
            continue;

        if(bin >= 0 && bin < nbins){
            H[bin] ++;
        }

        N1 ++;
    }

    if(N1 < 100)
        return NAN;

    for(d = 0; d < nbins; d++){
        double step = I[1] - I[0];

        dH[d] = sqrt(MAX(1, H[d]))/N1/step;
        H[d] = H[d]/N1/step;
    }

    model.I = I;
    model.H = H;
    model.dH = dH;
    model.nbins = nbins;

    memset(paropts, 0, sizeof(mp_par)*4);
    memset(&result, 0, sizeof(mp_result));
    result.xerror = parerrs;

    paropts[0].step = 0.1;
    paropts[1].relstep = 0.01;
    paropts[2].relstep = 0.01;

    paropts[0].limited[0] = 1;
    paropts[0].limits[0] = min;
    paropts[0].limited[1] = 1;
    paropts[0].limits[1] = max;

    paropts[1].fixed = 0;
    paropts[1].limited[0] = 1;
    paropts[1].limits[0] = 1;
    paropts[1].limited[1] = 1;
    paropts[1].limits[1] = 5.0*mad;

    paropts[2].limited[0] = 1;
    paropts[2].limits[0] = 1;
    paropts[2].limited[1] = 1;
    paropts[2].limits[1] = 1000;

    paropts[3].limited[0] = 1;
    paropts[3].limits[0] = 0;
    paropts[3].limited[1] = 1;
    paropts[3].limits[1] = 1;

    params[0] = med;
    params[1] = 1.4*mad;
    params[2] = 100;
    params[3] = 0.0;

    mpfit((mp_func)fn_mpfit, nbins, 4, params, paropts, &conf, &model, &result);

    if(error_ptr)
        *error_ptr = params[1];

    return params[0];
}

static double get_mode_bg(double *data, int N)
{
    int iter;
    double mad = 0;
    double med = 0;
    double *data1 = data;
    int N1 = N;

    qsort(data, N, sizeof(double), cmpfn);

    for(iter = 0; iter < 10; iter++){
        med = get_median_mad(data1, N1, &mad);

        while(data1[N1-1] > med + 2.0*mad)
            N1 --;

        if(iter > 7)
            while(data1[0] < med - 3.0*mad){
                data1 ++;
                N1 --;
            }
    }

    return 2.5*med - 1.5*get_mean(data1, N1, NULL);
}

image_str *image_background(image_str *image, image_str *errors, image_str *mask, int step)
{
    image_str *bgimage = image_create_double(image->width, image->height);
    double saturation;
    int i;
    int j;

    int nx = ceil(3.0*image->width/step);
    int ny = ceil(3.0*image->height/step);

    csa *vcsa = csa_create();
    csa *ecsa = csa_create();
    point *points = NULL;
    point *epoints = NULL;
    int Npoints = 0;

    if(image_keyword_find(image, "SATURATE"))
        saturation = image_keyword_get_double(image, "SATURATE");
    else
        saturation = image_max_value(image);

    for(i = 0; i <= nx; i++)
        for(j = 0; j <= ny; j++){
            int x0 = image->width*i/nx;
            int y0 = image->height*j/ny;
            int x;
            int y;

            int N = 0;
            double *data = NULL;

            for(x = MAX(0, x0 - 0.5*step); x < MIN(image->width, x0 + 0.5*step); x++)
                for(y = MAX(0, y0 - 0.5*step); y < MIN(image->height, y0 + 0.5*step); y++){
                    if(!IS_PIXEL_VALID(image, saturation, x, y))
                        continue;
                   if(mask && PIXEL(mask, x, y))
                        continue;

                    data = realloc(data, sizeof(double)*(N + 1));
                    data[N] = PIXEL_DOUBLE(image, x, y);
                    N += 1;
                }

            if(N > 3){
                double error = 0;
                //double value = get_mode_bg(data, N);
                double value = get_bijaoui_bg(data, N, &error);

                points = realloc(points, sizeof(point)*(Npoints + 1));
                epoints = realloc(epoints, sizeof(point)*(Npoints + 1));

                points[Npoints].x = x0;
                points[Npoints].y = y0;
                points[Npoints].z = value;

                epoints[Npoints].x = x0;
                epoints[Npoints].y = y0;
                epoints[Npoints].z = error;

                Npoints ++;
            }

            free(data);
        }

    csa_setk(vcsa, 80);
    csa_addpoints(vcsa, Npoints, points);
    csa_calculatespline(vcsa);

    csa_setk(ecsa, 80);
    csa_addpoints(ecsa, Npoints, epoints);
    if(errors){
        csa_calculatespline(ecsa);
    }

    free(points);
    free(epoints);

    for(j = 0; j < image->height; j++)
        for(i = 0; i < image->width; i++){
            point p;

            p.x = i;
            p.y = j;
            p.z = 0;

            csa_approximatepoint(vcsa, &p);
            PIXEL_DOUBLE(bgimage, i, j) = p.z;

            if(errors){
                csa_approximatepoint(ecsa, &p);
                PIXEL_DOUBLE(errors, i, j) = p.z;
            }
        }

    csa_destroy(vcsa);
    csa_destroy(ecsa);

    return bgimage;
}
