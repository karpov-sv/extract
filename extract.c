#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_cdf.h>

#include "utils.h"

#include "extract.h"
#include "psf.h"
#include "mpfit.h"
#include "kdtree.h"
#include "random.h"

static int debug = FALSE;
static int peak_to_stop = -1;

static image_str *image0 = NULL;

typedef struct {
    image_str *image;
    image_str *errors;
    image_str *mask;
    psf_str *psf;

    int x0;
    int y0;

    int size;

    double threshold; /* Minimal signal/noise for fitted flux */

    int Ngood;
    int *good;

    int Nglobals;
    int Npeakpars;

    int Npeaks;
    peak_str **peaks;

    double *params;
    double *parerrs;
    mp_par *paropts;

    double *dx;
    double *dy;
    double *A;

    image_str **psfs;
    image_str **stamps;

    mp_result result;
    double *residuals;

    /* Deblended peaks */
    struct list_head newpeaks;
} model_str;

model_str *model_create(image_str *image, image_str *errors, image_str *mask, psf_str *psf, int x0, int y0, int size)
{
    model_str *model = (model_str *)malloc(sizeof(model_str));
    int d;

    model->image = image;
    model->errors = errors;
    model->mask = mask;
    model->psf = psf;

    model->x0 = x0;
    model->y0 = y0;
    model->size = size;

    model->threshold = 2;

    model->Nglobals = 3; // 6;
    model->Npeakpars = 3;
    model->Npeaks = 0;

    model->peaks = NULL;
    model->params = (double *)calloc(model->Nglobals, sizeof(double));
    model->parerrs = (double *)calloc(model->Nglobals, sizeof(double));
    model->paropts = (mp_par *)calloc(model->Nglobals, sizeof(mp_par));

    model->dx = NULL;
    model->dy = NULL;
    model->A = NULL;
    model->psfs = NULL;
    model->stamps = NULL;

    model->residuals = calloc(model->image->width*model->image->height, sizeof(double));

    init_list(model->newpeaks);

    /* Background */
    for(d = 0; d < model->Nglobals; d++){
        model->params[d] = 0;
        model->paropts[d].fixed = 1;
        model->paropts[d].side = 0; /* Analytical derivative available */
        //model->paropts[d].deriv_debug = 1; /* Analytical derivative debug */
    }

    model->paropts[0].fixed = 0;
    /* model->params[0] = image_min_value(model->image); */

    /* Mask */
    model->good = calloc(model->image->width*model->image->height, sizeof(int));
    model->Ngood = 0;
    for(d = 0; d < mask->width*mask->height; d++)
        if(!(mask->data[d] & FLAG_SATURATED) && !(mask->data[d] & FLAG_BAD)){
            model->good[d] = model->Ngood;
            model->Ngood ++;
        } else
            model->good[d] = -1; /* Fake index which will not be used */

    return model;
}

void model_delete(model_str *model)
{
    if(model->params)
        free(model->params);
    if(model->parerrs)
        free(model->parerrs);
    if(model->paropts)
        free(model->paropts);

    if(model->residuals)
        free(model->residuals);
    if(model->good)
        free(model->good);

    if(model->dx)
        free(model->dx);
    if(model->dy)
        free(model->dy);
    if(model->A)
        free(model->A);

    if(model->psfs){
        int d;

        for(d = 0; d < model->Npeaks; d++){
            image_delete(model->psfs[d]);
            image_delete(model->stamps[d]);
        }

        free(model->psfs);
        free(model->stamps);
    }

    free_list(model->newpeaks);

    free(model);
}

void model_add_peak(model_str *model, peak_str *peak)
{
    /* double delta = 1.0; /\* Positional uncertainty *\/ */
    double delta = MAX(1.0, 0.5*model->psf->fwhm); /* Positional uncertainty */
    int idx = model->Nglobals + model->Npeakpars*model->Npeaks;

    if(peak->flags & FLAG_SATURATED)
        delta *= 3;

    /* dprintf("Adding peak %d at %g %g state %d\n", peak->id, peak->x - model->x0, peak->y - model->y0, peak->state); */

    model->params = realloc(model->params, sizeof(double)*(idx + model->Npeakpars));
    model->parerrs = realloc(model->parerrs, sizeof(double)*(idx + model->Npeakpars));
    model->paropts = realloc(model->paropts, sizeof(mp_par)*(idx + model->Npeakpars));
    memset(model->paropts + idx, 0, model->Npeakpars*sizeof(mp_par));

    model->peaks = realloc(model->peaks, sizeof(peak_str *)*(model->Npeaks + 1));
    model->dx = realloc(model->dx, sizeof(double)*(model->Npeaks + 1));
    model->dy = realloc(model->dy, sizeof(double)*(model->Npeaks + 1));
    model->A = realloc(model->A, sizeof(double)*(model->Npeaks + 1));
    model->psfs = realloc(model->psfs, sizeof(image_str *)*(model->Npeaks + 1));
    model->stamps = realloc(model->stamps, sizeof(image_str *)*(model->Npeaks + 1));

    model->peaks[model->Npeaks] = peak;

    model->psfs[model->Npeaks] = psf_sampled_image(model->psf, peak->x, peak->y);
    model->stamps[model->Npeaks] = image_create_double(model->size, model->size);
    psf_image_fill(model->psf, model->psfs[model->Npeaks], model->stamps[model->Npeaks], 0.0, 0.0);
    model->dx[model->Npeaks] = 0.0;
    model->dy[model->Npeaks] = 0.0;
    model->A[model->Npeaks] = 0.0;

    model->params[idx + 0] = peak->x - model->x0;
    model->params[idx + 1] = peak->y - model->y0;
    model->params[idx + 2] = peak->flux;

    if(model->params[idx + 2] < 0)
        model->params[idx + 2] = 0;

    model->paropts[idx + 0].step = 0.01;
    model->paropts[idx + 1].step = 0.01;
    model->paropts[idx + 2].relstep = 0.001;

    /* model->paropts[idx + 0].fixed = 1; */
    /* model->paropts[idx + 1].fixed = 1; */

    if(peak->state == PEAK_MEASURED){
        /* For already measured peaks we fix all parameters */
        model->paropts[idx + 0].fixed = 1;
        model->paropts[idx + 1].fixed = 1;
        model->paropts[idx + 2].fixed = 1;
    } else {
        if(model->params[idx + 0] < 0 || model->params[idx + 0] >= model->size ||
           model->params[idx + 1] < 0 || model->params[idx + 1] >= model->size){
            /* For peaks outside the subimage we fix the positions */
            model->paropts[idx + 0].fixed = 1;
            model->paropts[idx + 1].fixed = 1;
        } else {
            model->paropts[idx + 0].limited[0] = 1;
            model->paropts[idx + 0].limited[1] = 1;
            model->paropts[idx + 0].limits[0] = peak->x - model->x0 - delta;
            model->paropts[idx + 0].limits[1] = peak->x - model->x0 + delta;

            model->paropts[idx + 1].limited[0] = 1;
            model->paropts[idx + 1].limited[1] = 1;
            model->paropts[idx + 1].limits[0] = peak->y - model->y0 - delta;
            model->paropts[idx + 1].limits[1] = peak->y - model->y0 + delta;
        }

        model->paropts[idx + 2].limited[0] = 1;
        model->paropts[idx + 2].limited[1] = 0;
        model->paropts[idx + 2].limits[0] = 0;
        model->paropts[idx + 2].side = 0; /* Analytical derivative available */
        //model->paropts[idx + 2].deriv_debug = 1;
    }

    model->Npeaks ++;
}

void model_fill_image(model_str *model, double *params, double *imdata, double **derivatives)
{
    int Npoints = model->image->width*model->image->height;
    int d;

    for(d = 0; d < Npoints; d++){
        int x = d % model->image->width;
        int y = (d - x) / model->image->width;

        imdata[d] = params[0] + x*params[1] + y*params[2];// + x*y*params[3] + x*x*params[4] + y*y*params[5];
    }

    /* Imdata has full model size, but derivatives should contain only 'good' points! */
    if(derivatives){
        if(derivatives[0])
            for(d = 0; d < Npoints; d++)
                /* Derivatives should be normalized inplace, as they are used for residuals only */
                if(model->good[d] >= 0)
                    derivatives[0][model->good[d]] = -1.0/model->errors->double_data[d];
        if(derivatives[1])
            for(d = 0; d < Npoints; d++){
                int x = d % model->image->width;

                if(model->good[d] >= 0)
                    derivatives[1][model->good[d]] = -x*1.0/model->errors->double_data[d];
            }
        if(derivatives[2])
            for(d = 0; d < Npoints; d++){
                int x = d % model->image->width;
                int y = (d - x) / model->image->width;

                if(model->good[d] >= 0)
                    derivatives[2][model->good[d]] = -y*1.0/model->errors->double_data[d];
            }
    }

    for(d = 0; d < model->Npeaks; d++){
        int idx = model->Nglobals + model->Npeakpars*d;
        double x = params[idx + 0];
        double y = params[idx + 1];
        double A = params[idx + 2];
        /* PSF stamp */
        image_str *stamp = model->stamps[d];
        /* Sub-pixel shift of PSF center relative to stamp position */
        double dx = x - floor(model->peaks[d]->x - model->x0);
        double dy = y - floor(model->peaks[d]->y - model->y0);
        /* Origin of PSF stamp inside the image */
        int x0 = floor(model->peaks[d]->x - model->x0) - 0.5*stamp->width;
        int y0 = floor(model->peaks[d]->y - model->y0) - 0.5*stamp->height;
        int x1;
        int y1;

        /* Simple caching of stamps if sub-pixel adjustments is unchanged */
        if((fabs(A - model->A[d]) >= 0.01*A && model->psf->sz) || dx != model->dx[d] || dy != model->dy[d]){
            if(fabs(A - model->A[d]) >= 0.01*A && model->psf->sz){
                image_delete(model->psfs[d]);
                model->psfs[d] = psf_sampled_image(model->psf, model->peaks[d]->x, model->peaks[d]->y);
                model->A[d] = A;
            }
            psf_image_fill(model->psf, model->psfs[d], model->stamps[d], dx, dy);
            model->dx[d] = dx;
            model->dy[d] = dy;
        }

        if(derivatives && derivatives[idx + 2])
            for(y1 = MAX(0, -y0); y1 < MIN(stamp->height, model->image->height - y0); y1++)
                for(x1 = MAX(0, -x0); x1 < MIN(stamp->width, model->image->width - x0); x1++){
                    double value = PIXEL_DOUBLE(stamp, x1, y1);
                    int idx1 = (y0 + y1)*model->image->width + (x0 + x1);

                    imdata[idx1] += value*A;
                    if(model->good[idx1] >= 0)
                        derivatives[idx + 2][model->good[idx1]] = -value/model->errors->double_data[idx1];
                }
        else
            for(y1 = MAX(0, -y0); y1 < MIN(stamp->height, model->image->height - y0); y1++)
                for(x1 = MAX(0, -x0); x1 < MIN(stamp->width, model->image->width - x0); x1++){
                    double value = PIXEL_DOUBLE(stamp, x1, y1);
                    int idx1 = (y0 + y1)*model->image->width + (x0 + x1);

                    imdata[idx1] += value*A;
                }
    }
}

int fn_mpfit(int Npoints, int Nparams, double *params, double *residuals, double **derivatives, void *data)
{
    model_str *model = (model_str *)data;
    double *imdata = (double *)malloc(sizeof(double)*model->image->width*model->image->height);
    int d;

    model_fill_image(model, params, imdata, derivatives);

    for(d = 0; d < model->image->width*model->image->height; d++){
        if(model->good[d] >= 0){
            residuals[model->good[d]] = (model->image->double_data[d] - imdata[d])/model->errors->double_data[d];
        }
    }

    free(imdata);

    return 0;
}

void model_dump_to_fits(model_str *model, char *filename, char *peaksname)
{
    image_str *image = image_create_double(model->image->width, model->image->height);

    model_fill_image(model, model->params, image->double_data, NULL);
    image_dump_to_fits(image, filename);

    image_delete(image);

    if(peaksname){
        FILE *file = fopen(peaksname, "w");

        int d;
        for(d = 0; d < model->Npeaks; d++)
            fprintf(file, "%g %g %g %g %d\n",
                    model->params[model->Nglobals + model->Npeakpars*d + 0],
                    model->params[model->Nglobals + model->Npeakpars*d + 1],
                    model->params[model->Nglobals + model->Npeakpars*d + 2],
                    model->peaks[d]->chisq,
                    model->peaks[d]->id);
        fclose(file);
    }
}

void model_fit(model_str *model)
{
    int Nparams = model->Nglobals + model->Npeakpars*model->Npeaks;
    struct mp_config_struct conf;

    memset(&conf, 0, sizeof(conf));
    conf.ftol = 1e-16;
    conf.xtol = 1e-16;
    conf.gtol = 1e-16;
    conf.maxiter = 1000;

    memset(&model->result, 0, sizeof(mp_result));
    model->result.xerror = model->parerrs;
    model->result.resid = model->residuals;

    if(debug){
        image_dump_to_fits(model->image, "out.crop.fits");
        image_dump_to_fits(model->errors, "out.croperr0.fits");
        image_dump_to_fits(model->mask, "out.cropmask.fits");
        model_dump_to_fits(model, "out.model.fits", "out.peaks0.txt");
        image_dump_to_fits(model->psfs[0], "out.psf.fits");
    }

    //mpfit((mp_func)fn_mpfit, model->image->width*model->image->height, Nparams, model->params, model->paropts, NULL /* &conf */, model, &model->result);
    mpfit((mp_func)fn_mpfit, model->Ngood, Nparams, model->params, model->paropts, NULL /* &conf */, model, &model->result);

    if(debug){
        model_dump_to_fits(model, "out.result.fits", "out.peaks.txt");
        image_dump_to_fits(model->errors, "out.croperr.fits");
    }

    if(debug){
        dprintf("%g %g: %g -> %g, chisq = %g / %g, (%d/%d, %d peaks), %d/%d runs, status=%d\n",
                model->x0 + 0.5*model->size, model->y0 + 0.5*model->size,
                model->result.orignorm, model->result.bestnorm,
                model->result.bestnorm/(model->result.nfunc - model->result.nfree),
                1.0 - gsl_cdf_chisq_P(model->result.bestnorm, model->result.nfunc - model->result.nfree),
                model->result.nfunc, model->result.nfree, model->Npeaks,
                model->result.niter, model->result.nfev, model->result.status);
    }

    /* if(model->Npeaks > 1) */
    /*     exit(1); */
}

double model_peak_chisq(model_str *model, peak_str *peak, int *nbad_ptr)
{
    double chisq = 0;
    int Npoints = 0;
    int nbad = 0;

    double size = 2*model->psf->fwhm;
    int x;
    int y;

    for(x = MAX(0, floor(peak->x - model->x0 - 0.5*size)); x < MIN(model->image->width, ceil(peak->x - model->x0 + 0.5*size)); x++)
        for(y = MAX(0, floor(peak->y - model->y0 - 0.5*size)); y < MIN(model->image->height, ceil(peak->y - model->y0 + 0.5*size)); y++){
            if(model->good[y*model->image->width + x] < 0)
                continue;

            chisq += pow(model->residuals[model->good[y*model->image->width + x]], 2);
            Npoints ++;

            if(debug && peak->id == peak_to_stop && FALSE)
                printf("%g\n", model->residuals[model->good[y*model->image->width + x]]);

            /* Count the bad pixels near peak top */
            /* FIXME: make it configurable and optional */
            if(x >= floor(peak->x - model->x0 - 1) && x <= ceil(peak->x - model->x0 + 1) &&
               y >= floor(peak->y - model->y0 - 1) && y <= ceil(peak->y - model->y0 + 1) &&
               model->residuals[model->good[y*model->image->width + x]] < -5.0){
                if(debug)
                    dprintf("bad in peak %d: %d %d - %g\n", peak->id, x + model->x0, y + model->y0, model->residuals[model->good[y*model->image->width + x]]);
                nbad ++;
            }
        }

    if(debug && Npoints > model->Npeakpars)
        dprintf("chisq=%g / %g df=%d p-value=%g for peak %d\n", chisq, chisq/(Npoints - model->Npeakpars), Npoints - model->Npeakpars, 1.0-gsl_cdf_chisq_P(chisq, Npoints - model->Npeakpars), peak->id);

    if(nbad_ptr)
        *nbad_ptr = nbad;

    if(Npoints > model->Npeakpars)
        return 1.0 - gsl_cdf_chisq_P(chisq, (Npoints - model->Npeakpars));
    else
        return -1;
}

/* Effectively mask some pixels near the peak top having too small values - Andor sCMOS problem */
void mask_bad_pixels(model_str *model, peak_str *peak)
{
    int x;
    int y;

    for(x = MAX(0, floor(peak->x - model->x0 - 1)); x < MIN(model->image->width, ceil(peak->x - model->x0 + 2)); x++)
        for(y = MAX(0, floor(peak->y - model->y0 - 1)); y < MIN(model->image->height, ceil(peak->y - model->y0 + 2)); y++){
            if(model->residuals[y*model->image->width + x] < -5.0){
                PIXEL_DOUBLE(model->errors, x, y) = psf_satur_level;//*fabs(model->residuals[y*model->image->width + x]);
            }
        }
}

int model_update_peaks(model_str *model, image_str *image, int is_mask_bad)
{
    int Nprocessed = 0; /* Number of either fitted or failed peaks */
    double edgesize = 2.0;
    int d;

    double Pchisq0 = 0.0;

    int is_stop = FALSE;

    if(model->result.status <= 0 || model->result.status > 4)
        /* Fit has not converged properly, let's skip this peak without marking anything good or bad */
        return 0;

    if(debug){
        dprintf("globals:");
        for(d = 0; d < model->Nglobals; d++)
            dprintf(" %g", model->params[d]);
        dprintf("\n");
    }

    Pchisq0 = 1.0 - gsl_cdf_chisq_P(model->result.bestnorm, model->result.nfunc - model->result.nfree);

    for(d = 0; d < model->Npeaks; d++){
        peak_str *peak = model->peaks[d];
        int nbad = 0;
        double chisq = model_peak_chisq(model, peak, &nbad);

        if(peak->state == PEAK_INITIAL){
            if(chisq >= 0 && /* Pchisq0 > 1e-15 && */ /* chisq > 0.01 && */
               peak->x - model->x0 >= edgesize &&
               peak->x - model->x0 < model->image->width - edgesize &&
               peak->y - model->y0 >= edgesize &&
               peak->y - model->y0 < model->image->height - edgesize &&

               model->params[model->Nglobals + model->Npeakpars*d + 0] >= edgesize &&
               model->params[model->Nglobals + model->Npeakpars*d + 0] < model->image->width - edgesize &&
               model->params[model->Nglobals + model->Npeakpars*d + 1] >= edgesize &&
               model->params[model->Nglobals + model->Npeakpars*d + 1] < model->image->height - edgesize){
                /* Peak is _potentially_ measured well, more checks below */

                peak->bg = model->params[0];
                peak->dbg = model->parerrs[0];

                if(debug){
                    dprintf("updating flux for peak %d / %d: %g -> %g\n", peak->id, peak->state, peak->flux, model->params[model->Nglobals + model->Npeakpars*d + 2]);
                }

                peak->x = model->params[model->Nglobals + model->Npeakpars*d + 0] + model->x0;
                peak->dx = model->parerrs[model->Nglobals + model->Npeakpars*d + 0];
                peak->y = model->params[model->Nglobals + model->Npeakpars*d + 1] + model->y0;
                peak->dy = model->parerrs[model->Nglobals + model->Npeakpars*d + 1];
                peak->flux = model->params[model->Nglobals + model->Npeakpars*d + 2];
                peak->dflux = model->parerrs[model->Nglobals + model->Npeakpars*d + 2];

                peak->chisq = chisq;

                /* Bad pixels?.. */
                /* FIXME: make it configurable and optional */
                if(is_mask_bad && nbad <= 2){
                    mask_bad_pixels(model, peak);

                    return -1;
                }

                if(peak->flux == 0 || peak->dflux == 0 || fabs(peak->dflux/peak->flux) > 1.0/model->threshold /* || nbad > 2 */ /* || peak->chisq < 1e-3 */){
                    /* Flux accuracy is below the limit */
                    peak->state = PEAK_FAILED;
                } else if(peak->state != PEAK_MEASURED){
                    peak->state = PEAK_MEASURED;
                }

                Nprocessed ++;
            }

            if(debug && peak->id == peak_to_stop && // peak->state == PEAK_MEASURED &&
               peak->x - model->x0 >= edgesize &&
               peak->x - model->x0 < model->image->width - edgesize &&
               peak->y - model->y0 >= edgesize &&
               peak->y - model->y0 < model->image->height - edgesize){
                is_stop = TRUE;

                image_dump_to_fits(model->psfs[d], "out.psf.fits");
                image_dump_to_fits(model->stamps[d], "out.stamp.fits");
            }

            if(debug)
                dprintf("peak %d : %g %g - %g %g - %g %g - (%g +/- %g) - flags %d - %d\n\n",
                        peak->id, peak->x, peak->y, peak->x - model->x0, peak->y - model->y0,
                        chisq, fabs(peak->dflux/peak->flux), model->params[model->Nglobals + model->Npeakpars*d + 2],
                        model->parerrs[model->Nglobals + model->Npeakpars*d + 2], peak->flags, peak->state);
        }
    }

    if(debug && is_stop){
        dprintf("Stopping at peak %d\n", peak_to_stop);
        exit(1);
    }


    return Nprocessed;
}

void normalize_and_subtract_peaks_from_image(image_str *image, psf_str *psf, struct list_head *peaks, int size, int is_subtract)
{
    image_str *stamp = image_create_double(size, size);
    peak_str *peak;

    /* dprintf("Normalizing fluxes and subtracting the model\n"); */

    foreach(peak, *peaks){
        if(peak->state == PEAK_MEASURED){
            int x0 = floor(peak->x); /* Peak sub-window center */
            int y0 = floor(peak->y);
            double dx = peak->x - x0; /* Sub-pixel adjustments for PSF model */
            double dy = peak->y - y0;
            int x;
            int y;
            double I = 0.0;

            image_str *psfi = psf_sampled_image(psf, peak->x, peak->y);

            peak->fwhm = image_fwhm(psfi);
            peak->ellipticity = image_ellipticity(psfi);

            psf_image_fill(psf, psfi, stamp, dx, dy);

            for(x = 0; x < stamp->width; x++)
                for(y = 0; y < stamp->height; y++){
                    int x1 = x0 - 0.5*stamp->width + x;
                    int y1 = y0 - 0.5*stamp->height + y;

                    double value = PIXEL_DOUBLE(stamp, x, y);

                    I += value;

                    if(is_subtract && !(peak->flags & FLAG_SATURATED) &&
                       x1 >= 0 && x1 < image->width &&
                       y1 >= 0 && y1 < image->height)
                        PIXEL_DOUBLE(image, x1, y1) -= value*peak->flux;
                }

            for(x = x0 - 2; x <= x0 + 2; x++)
                for(y = y0 - 2; y <= y0 + 2; y++)
                    if(x >= 0 && x < image0->width &&
                       y >= 0 && y < image0->height)
                        peak->maxflux = MAX(peak->maxflux, PIXEL_DOUBLE(image0, x, y));

            peak->minflux = peak->maxflux;

            for(x = x0 - 10; x <= x0 + 10; x++)
                for(y = y0 - 10; y <= y0 + 10; y++)
                    if(x >= 0 && x < image0->width &&
                       y >= 0 && y < image0->height)
                        peak->minflux = MIN(peak->minflux, PIXEL_DOUBLE(image0, x, y));

            peak->minflux = image_min_value(psfi);

            if(debug){
                dprintf("peak %d: %g +/- %g -> %g +/- %g (I=%g)\n", peak->id, peak->flux, peak->dflux, peak->flux*I, peak->dflux*I, I);
            }

            /* peak->flux *= I; */
            /* peak->dflux *= I; */

            image_delete(psfi);
        }
    }

    image_delete(stamp);
}

int model_deblend(model_str *model)
{
    image_str *image = image_create_double(model->image->width, model->image->height);
    image_str *res = image_create_double(model->image->width, model->image->height);
    image_str *smooth = image_create_double(model->image->width, model->image->height);
    image_str *kernel = image_create_double(model->image->width, model->image->height);
    peak_str *peak;

    int d;
    int x;
    int y;
    int result = FALSE;

    model_fill_image(model, model->params, image->double_data, NULL);
    psf_image_fill(model->psf, model->psfs[0], kernel, 0, 0);

    for(d = 0; d < image->width*image->height; d++)
        res->double_data[d] = (model->image->double_data[d] - image->double_data[d]);

    /* Smooth the residuals with PSF */
    for(y = 0; y < image->height; y++)
        for(x = 0; x < image->width; x++){
            int dx;
            int dy;
            double value = 0;
            double kernel_sum = 0;

            for(dy = -model->image->height/2; dy < model->image->height/2; dy++)
                for(dx = -model->image->width/2; dx < model->image->width/2; dx++){
                    double weight;

                    if(x + dx < 0 || x + dx >= model->image->width ||
                       y + dy < 0 || y + dy >= model->image->height)
                        continue;

                    weight = PIXEL_DOUBLE(kernel, dx + model->image->width/2, dy + model->image->height/2);

                    if(!(PIXEL(model->mask, x + dx, y + dy) & FLAG_SATURATED)){
                        value += weight*PIXEL_DOUBLE(res, x + dx, y + dy);
                        kernel_sum += weight;
                    }
                }

            PIXEL_DOUBLE(smooth, x, y) = value/kernel_sum;
        }

    if(debug){
        image_dump_to_fits(res, "out.res.fits");
        image_dump_to_fits(smooth, "out.rsmooth.fits");
    }

    /* Find peaks in smoothed residuals */
    find_peaks(image, smooth, model->errors, model->mask, model->threshold, &model->newpeaks);

    if(debug && list_length(&model->newpeaks) > 0)
        dprintf("\t->\t%d new peaks found\n", list_length(&model->newpeaks));

    /* Inject new peaks into the model */
    foreach(peak, model->newpeaks){
        if(peak->x > 2 && peak->x < image->width-2 &&
           peak->y > 2 && peak->y < image->height-2){
            int is_close = FALSE;

            for(d = 0; d < model->Npeaks; d++){
                if(hypot(peak->x + model->x0 - model->peaks[d]->x,
                         peak->y + model->y0 - model->peaks[d]->y) < 1){
                    is_close = TRUE; /* New peak is too close to existing one */
                    break;
                }
            }

            if(!is_close){
                if(debug)
                    dprintf("\tadding peak at %g %g\n", peak->x, peak->y);

                peak->x += model->x0;
                peak->y += model->y0;

                model_add_peak(model, peak);

                result = TRUE;
            }
        }
    }

    image_delete(image);
    image_delete(res);
    image_delete(smooth);
    image_delete(kernel);

    return result;
}

void psf_fit_peaks(image_str *image, image_str *errors, psf_str *psf, image_str *mask, double threshold, int is_deblend, struct list_head *peaks)
{
    struct kdtree *kd = kd_create(2);
    time_str time_start = time_current();
    peak_str *peak;
    int size = 2*floor(0.5*psf->width*psf->pix_step); /* It should be always even */
    int edgesize = 1.0;
    int Npeaks = 0;
    int Nprocessed = 0;
    int maxid = 0;

    Npeaks = list_length(peaks);

    dprintf("Window size is %d x %d pixels\n", size, size);

    /* Fill kd-tree with all peaks */
    foreach(peak, *peaks){
        kd_insert2(kd, peak->x, peak->y, peak);
        maxid = MAX(maxid, peak->id);
    }

    /* Cycle all unfitted peaks, skipping ones near the edges */
    foreach(peak, *peaks)
        if(peak->state == PEAK_INITIAL &&
           peak->x > 2*edgesize && peak->x < image->width - 2*edgesize &&
           peak->y > 2*edgesize && peak->y < image->height - 2*edgesize){
            /* Re-center the sub-window to exclude image edges */
            int xc = MAX(0.5*size + edgesize, MIN(image->width - 0.5*size - edgesize, round(peak->x)));
            int yc = MAX(0.5*size + edgesize, MIN(image->height - 0.5*size - edgesize, round(peak->y)));
            struct kdres *res = kd_nearest_range2(kd, xc, yc, 2.0*size);
            int N = 0;
            int Np = 0;

            image_str *image_cropped = image_crop(image, xc - 0.5*size, yc - 0.5*size, xc + 0.5*size, yc + 0.5*size);
            image_str *errors_cropped = image_crop(errors, xc - 0.5*size, yc - 0.5*size, xc + 0.5*size, yc + 0.5*size);
            image_str *mask_cropped = image_crop(mask, xc - 0.5*size, yc - 0.5*size, xc + 0.5*size, yc + 0.5*size);

            model_str *model = model_create(image_cropped, errors_cropped, mask_cropped, psf, xc - 0.5*size, yc - 0.5*size, size);

            int has_saturated = FALSE;

            model->threshold = threshold;

            if(debug)
                dprintf("\nmodel for peak %d (flags %d): %g %g - %d %d\n", peak->id, peak->flags, peak->x, peak->y, xc, yc);

            while(!kd_res_end(res)){
                peak_str *peak1 = kd_res_item_data(res);

                if((/* peak1->state != PEAK_FAILED && */
                    peak1->x >= xc - 1.0*size && peak1->x < xc + 1.0*size &&
                    peak1->y >= yc - 1.0*size && peak1->y < yc + 1.0*size) &&
                   peak1->flux / peak->flux > 0.01){
                    /* Simple heuristic to exclude faraway faint stars from model */
                    if(peak1->flags & FLAG_SATURATED ||
                       (peak1->x >= xc - 0.7*size && peak1->x < xc + 0.7*size &&
                        peak1->y >= yc - 0.7*size && peak1->y < yc + 0.7*size)){
                        /* Should include the peak in model */
                        model_add_peak(model, peak1);

                        if(peak1->flags & FLAG_SATURATED)
                            has_saturated = TRUE;

                        N ++;
                    }
                }

                kd_res_next(res);
            }

            /* Here we should have at least one unfitted peak */
            model_fit(model);

            /* TODO: Here we should somehow analyze fit results and perform
               the detection and injection of additional peaks into model */
            if(!has_saturated && is_deblend)
                if(model_deblend(model)){
                    if(debug)
                        dprintf("Re-fitting model\n");
                    model_fit(model);
                }

            /* We may need to re-fit the model if it has been somehow updated/modified */
            if((Np = model_update_peaks(model, image, FALSE)) < 0){
                model_fit(model);
                Np = model_update_peaks(model, image, FALSE);
            }

            {
                /* Process deblended peaks if any */
                peak_str *peak1;
                int is_quit = FALSE;

                foreach(peak1, model->newpeaks){
                    if(debug)
                        dprintf("newpeak at %g %g - %d\n", peak1->x, peak1->y, peak1->state);

                    if(peak1->state == PEAK_MEASURED){
                        peak->flags |= FLAG_DEBLENDED;
                        peak1->flags |= FLAG_DEBLENDED;
                        peak1->id = maxid + 1;
                        maxid ++;

                        Npeaks ++;

                        //dprintf("%d - %d\n", Npeaks, list_length(peaks));

                        is_quit = TRUE;

                        kd_insert2(kd, peak->x, peak->y, peak1);
                        del_from_list_in_foreach_and_run(peak1, add_to_list(*peaks, peak1));
                    }
                }

                /* if(is_quit) */
                /*     exit(1); */
            }

            Nprocessed += Np;

            if(isatty(fileno(stderr)))
                dprintf("\r  %d / %d - %d peaks - %g s        \r", Nprocessed, Npeaks, N, 1e-3*time_interval(time_start, time_current()));

            model_delete(model);

            image_delete(image_cropped);
            image_delete(errors_cropped);
            image_delete(mask_cropped);

            kd_res_free(res);
        }

    if(isatty(fileno(stderr)))
        dprintf("\n");

    /* Peak amplitudes should be normalized to fluxes (we assume the integral over PSF to be it) */
    normalize_and_subtract_peaks_from_image(image, psf, peaks, size, TRUE);

    if(debug){
        dump_peaks_to_file(peaks, "out.peaks.result.txt", PEAK_MEASURED);
        dump_peaks_to_file(peaks, "out.peaks.failed.txt", PEAK_FAILED);
        dump_peaks_to_file(peaks, "out.peaks.initial.txt", PEAK_INITIAL);

        image_dump_to_fits(image, "out.residuals.fits");
    }

    {
        int Nmeasured = 0;
        int Nfailed = 0;
        int Nunprocessed = 0;

        foreach(peak, *peaks){
            if(peak->state == PEAK_INITIAL)
                Nunprocessed ++;
            else if(peak->state == PEAK_MEASURED)
                Nmeasured ++;
            else if(peak->state == PEAK_FAILED)
                Nfailed ++;
        }

        dprintf("%d measured, %d failed, %d not processed\n", Nmeasured, Nfailed, Nunprocessed);
    }

    kd_free(kd);
}

image_str *image_errors_estimate(image_str *image)
{
    image_str *res = image_create_double(image->width, image->height);
    double saturation = image_max_value(image);
    int x;
    int y;

    for(y = 0; y < image->height; y++)
        for(x = 0; x < image->width; x++){
            double sum = 0;
            double sum2 = 0;
            int N = 0;
            int dx;
            int dy;

            for(dx = -1; dx <= 1; dx++)
                for(dy = -1; dy <= 1; dy++)
                    if(x + dx >= 0 && x + dx < image->width &&
                       y + dy >= 0 && y + dy < image->height &&
                       PIXEL_DOUBLE(image, x + dx, y + dy) < saturation){
                        sum += PIXEL_DOUBLE(image, x + dx, y + dy);
                        sum2 += PIXEL_DOUBLE(image, x + dx, y + dy)*PIXEL_DOUBLE(image, x + dx, y + dy);
                        N ++;
                    }

            if(N > 2)
                PIXEL_DOUBLE(res, x, y) = sqrt((sum2 - sum*sum/N)/(N - 1));
            else
                PIXEL_DOUBLE(res, x, y) = 0;
        }

    return res;
}

int main(int argc, char **argv)
{
    char *filename = NULL;
    char *psfname = NULL;
    char *outname = NULL;
    char *darkname = NULL;
    char *flatname = NULL;
    char *maskname = NULL;
    char *simplename = NULL;
    char *preprocessedname = NULL;
    char *errorsname = NULL;

    int x0 = 0;
    int y0 = 0;
    int width = 0;
    int height = 0;

    double ra0 = -1;
    double dec0 = 0;

    int is_andor = TRUE;
    int is_psfex = FALSE;
    int is_simple = FALSE;
    int is_gauss = FALSE;
    int is_keep_psf = FALSE;
    int is_keep_simple = FALSE;
    int is_keep_preprocessed = FALSE;
    int is_keep_mask = FALSE;
    int is_keep_errors = FALSE;
    int is_keep_residuals = FALSE;
    int is_keep_bg = FALSE;
    int is_keep_unsharp = FALSE;
    int is_unsharp = FALSE;
    int is_header = TRUE;
    int is_preprocess = FALSE;
    int is_sex = FALSE;
    int is_aper = FALSE;
    int is_deblend = FALSE;

    int is_preprocessed = FALSE;

    double bias = 0.0;
    double gain = 0.0;
    double readnoise = 0.0;
    double threshold = 5.0;

    double saturation = 0.0;
    double welldepth = 20000.0;

    double max_fwhm = 6.0;

    image_str *image = NULL;
    image_str *errors = NULL;
    image_str *mask = NULL;
    image_str *dark = NULL;
    image_str *flat = NULL;
    image_str *bg = NULL;

    psf_str *psf = NULL;

    char *peaksname = NULL;
    int is_radec = TRUE;

    parse_args(argc, argv,
               "x0=%d", &x0,
               "y0=%d", &y0,
               "ra0=%lf", &ra0,
               "dec0=%lf", &dec0,
               "width=%d", &width,
               "height=%d", &height,
               "psf=%s", &psfname,
               "preprocessed=%s", &preprocessedname,
               "errors=%s", &errorsname,

               "-andor", &is_andor,
               "-psfex", &is_psfex,
               "-simple", &is_simple,
               "-gauss", &is_gauss,
               "-keep_psf", &is_keep_psf,
               "-keep_simple", &is_keep_simple,
               "-keep_preprocessed", &is_keep_preprocessed,
               "-keep_mask", &is_keep_mask,
               "-keep_errors", &is_keep_errors,
               "-keep_residuals", &is_keep_residuals,
               "-keep_bg", &is_keep_bg,
               "-keep_unsharp", &is_keep_unsharp,
               "-unsharp", &is_unsharp,
               "-header", &is_header,
               "-preprocess", &is_preprocess,
               "-sex", &is_sex,
               "-aper", &is_aper,
               "-deblend", &is_deblend,
               "-debug", &debug,

               "peak=%d", &peak_to_stop,

               "bias=%lf", &bias,
               "gain=%lf", &gain,
               "readnoise=%lf", &readnoise,
               "threshold=%lf", &threshold,

               "saturation=%lf", &saturation,
               "welldepth=%lf", &welldepth,

               "max_fwhm=%lf", &max_fwhm,

               "psf_var_degree=%d", &psf_var_degree,
               "psf_vignet_size=%d", &psf_vignet_size,
               "psf_size=%d", &psf_size,
               "psf_satur_level=%lf", &psf_satur_level,
               "psf_gain=%lf", &psf_gain,
               "psf_phot_aper=%f", &psf_phot_aper,
               "-psf_subtract_bg", &psf_subtract_bg,

               "dark=%s", &darkname,
               "flat=%s", &flatname,
               "mask=%s", &maskname,

               "peaks=%s", &peaksname,
               "-radec", &is_radec,

               "%s", &filename,
               "%s", &outname,
               NULL);

    if(!filename)
        return EXIT_SUCCESS;

    image = image_create_from_fits(filename);

    if(!image){
        dprintf("Can't load FITS image from %s!\n", filename);
        return EXIT_FAILURE;
    } else if(image->type != IMAGE_DOUBLE){
        /* The code expects IMAGE_DOUBLE */
        image_str *new = image_convert_to_double(image);

        image_delete(image);
        image = new;

        dprintf("Image converted to DOUBLE\n");
    }

    image0 = image_copy(image);

    dprintf("%s: %d x %d, min = %g, max = %g, mean = %g, median = %g, sigma = %g\n", filename, image->width, image->height,
            image_min_value(image), image_max_value(image), image_mean(image), image_median(image), image_sigma(image));

    if(ra0 >= 0){
        double x = 0;
        double y = 0;

        coords_get_x_y(&image->coords, ra0, dec0, &x, &y);

        x0 = x - 0.5*width;
        y0 = y - 0.5*height;

        dprintf("Setting crop center to %d %d\n", (int)x, (int)y);
    }

    /* If no GAIN provided on command line, try to guess it from the header */
    if(is_header && !gain){
        if(image_keyword_find(image, "GAIN")){
            gain = image_keyword_get_double(image, "GAIN");

            dprintf("Setting GAIN to %g from GAIN header keyword\n", gain);
        }

        if(image_keyword_find(image, "RDNOISE")){
            readnoise = image_keyword_get_double(image, "RDNOISE");

            dprintf("Setting READNOISE to %g from RDNOISE header keyword\n", readnoise);
        }

        if(image_keyword_find(image, "SHUTTER")){
            if(image_keyword_get_int(image, "SHUTTER") == 0){
                gain = 0.67;
                readnoise = 1.4;
            } else {
                gain = 1.91;
                readnoise = 2.46;
            }

            dprintf("Setting GAIN to %g and READNOISE to %g from SHUTTER header keyword\n", gain, readnoise);
        }

        if(image_keyword_find(image, "AVERAGED")){
            int avg = MAX(1, image_keyword_get_int(image, "AVERAGED"));

            if(avg > 1){
                gain *= avg;
                readnoise /= sqrt(avg);
                welldepth *= avg;

                dprintf("The image was AVERAGED over %d frames, corrected the parameters correspondingly\n", avg);
            }
        }
    }

    if(!saturation){
        if(is_header && image_keyword_find(image, "SATURATE")){
            saturation = image_keyword_get_double(image, "SATURATE");

            dprintf("Setting SATURATION to %g from SATURATE header keyword\n", saturation);
        } else {
            saturation = floor(welldepth/gain);

            dprintf("Setting SATURATION to %g from %g depth well capacity\n", saturation, welldepth);
        }
    }

    /* Mask */
    if(maskname){
        mask = image_create_from_fits(maskname);

        if(mask && mask->type != IMAGE_UINT16){
            dprintf("MASK should be INT16, and %s is not\n", maskname);
            image_delete_and_null(mask);
        }

        if(mask)
            dprintf("MASK frame read from %s\n", maskname);
    } else {
        mask = image_create(image->width, image->height);
        image_fill(mask, 0);
    }

    if(gain <= 0){
        dprintf("No GAIN provided! Exiting.\n");
        return EXIT_FAILURE;
    }

    /* Dark */
    if(darkname)
        dark = image_create_from_fits(darkname);

    if(dark){
        dprintf("DARK frame read from %s, mean=%g\n", darkname, image_mean(dark));

        if(is_andor){
            //image_linearize(dark, NULL);
            dprintf("ANDOR dual-amplifier correction applied to DARK frame\n");
        }
    } else {
        dark = image_create_double(image->width, image->height);
        image_fill(dark, bias);
        dprintf("No DARK frame, using bias = %g instead\n", bias);
    }

    /* Flat */
    if(flatname)
        flat = image_create_from_fits(flatname);

    if(flat)
        dprintf("FLAT frame read from %s, mean=%g\n", flatname, image_mean(flat));
    else {
        flat = image_create_double(image->width, image->height);
        image_fill(flat, 1.0);
        dprintf("No FLAT frame, using uniform field instead\n");
    }

    /* Pre-processing */
    {
        int d;

        if(is_andor){
            image_str *emask = image_create(image->width, image->height);

            image_linearize(image, emask);
            dprintf("ANDOR dual-amplifier correction applied to the data frame\n");

            saturation = image_keyword_get_double(image, "SATURATE");

            for(d = 0; d < image->width*image->height; d++)
                if(emask->data[d])
                    //mask->data[d] |= FLAG_SATURATED;
                    mask->data[d] |= FLAG_BAD;

            image_delete(emask);
        }

        for(d = 0; d < image->width*image->height; d++)
            if(image->double_data[d] >= saturation ||
               !isfinite(image->double_data[d]) || image->double_data[d] == -0x80000000)
                mask->data[d] |= FLAG_SATURATED;

        if(is_keep_mask){
            maskname = make_string("%s.mask.fits", filename);

            image_dump_to_fits(mask, maskname);
        }

        dprintf("Saturated pixels masked\n");

        errors = image_create_double(image->width, image->height);

        /* Apply dark */
        for(d = 0; d < image->width*image->height; d++){

            image->double_data[d] = image->double_data[d] - dark->double_data[d];
        }

        /* Assign errors */
        for(d = 0; d < image->width*image->height; d++){
            errors->double_data[d] = hypot(readnoise/gain, sqrt(MAX(0, (image->double_data[d] - bias)/gain)));
        }

        /* Apply flat */
        for(d = 0; d < image->width*image->height; d++){
            if(!isfinite(flat->double_data[d]))
                mask->data[d] |= FLAG_SATURATED;
            else {
                image->double_data[d] /= flat->double_data[d];
                errors->double_data[d] /= flat->double_data[d];
            }
        }

        image_keyword_add_double(image, "GAIN", gain, "Median gain");
        image_keyword_add_double(image, "RDNOISE", readnoise, "Readout noise");
        image_keyword_add_double(image, "SATURATE", saturation, "Saturation level");

        if(debug){
            image_dump_to_fits(image, "out.prebg.fits");
            image_dump_to_fits(errors, "out.eprebg.fits");
        }

        /* Estimate background */
        //bg = get_background(image, mask, 128);
        bg = image_background(image, errors, NULL, 256);
        if(0){
            image_str *mask1 = image_create(image->width, image->height);
            int iter;

            for(iter = 0; iter < 3; iter++){
                image_str *bg1 = image_background(image, NULL /* errors */, mask, 128);

                for(d = 0; d < image->width*image->height; d++)
                    if(fabs(image->double_data[d] - bg1->double_data[d]) > 4.0*errors->double_data[d])
                        mask1->data[d] = TRUE;

                image_delete(bg1);
            }

            image_dump_to_fits(mask1, "out.mask1.fits");

            bg = image_background(image, errors, mask1, 128);

            image_delete(mask1);
        }

        if(!psf_subtract_bg)
            for(d = 0; d < image->width*image->height; d++)
                if(!(mask->data[d] & FLAG_SATURATED))
                    image->double_data[d] -= bg->double_data[d];

        dprintf("Image background estimated and subtracted\n");

        /* Make saturated pixels saturated again */
        for(d = 0; d < image->width*image->height; d++)
            if(mask->data[d] & FLAG_SATURATED || image->double_data[d] >= saturation)
                image->double_data[d] = 1.0*saturation + 1;

        if(debug)
            image_dump_to_fits(bg, "out.bg.fits");

        if(is_keep_bg){
            char *bgname = make_string("%s.bg.fits", filename);

            image_dump_to_fits(bg, bgname);
            dprintf("Image background map stored to %s\n", bgname);

            free(bgname);
        }

        if(is_keep_errors || errorsname){
            if(!errorsname)
                errorsname = make_string("%s.errors.fits", filename);

            image_dump_to_fits(errors, errorsname);
            dprintf("RMS image stored to %s\n", errorsname);
        }

        if(is_keep_preprocessed || is_preprocess || preprocessedname){
            if(!preprocessedname)
                preprocessedname = make_string("%s.processed.fits", filename);

            image_dump_to_fits(image, preprocessedname);
            dprintf("Pre-processed image stored to %s\n", preprocessedname);
        }

        is_preprocessed = TRUE;
        dprintf("Image pre-processed\n");
        dprintf("%s: %d x %d, min = %g, max = %g, mean = %g, median = %g, sigma = %g\n", filename, image->width, image->height,
                image_min_value(image), image_max_value(image), image_mean(image), image_median(image), image_sigma(image));

        if(is_preprocess)
            return EXIT_SUCCESS;
    }

    /* Various filenames */
    if(is_psfex && outname)
        psfname = outname;
    else if(!psfname)
        psfname = make_string("%s.psf", filename);

    simplename = make_string("%s.simple", filename);

    /* Fill PSFEx parameters */
    if(gain && !psf_gain)
        psf_gain = gain;

    if(saturation)
        psf_satur_level = saturation;

    dprintf("bias=%g gain=%g readnoise=%g saturation=%g\n", bias, gain, readnoise, saturation);

    if(is_aper){
        /* SExtractor aperture photometry */
        if(!outname || !*outname)
            outname = make_string("%s.txt", filename);
        if(psf_photometry(image, psfname, outname, FALSE)){
            dprintf("SExtractor aperture photometry results stored to %s\n", outname);
            return EXIT_SUCCESS;
        } else
            return EXIT_FAILURE;
    }

    if(file_exists_and_normal(psfname) && !is_psfex){
        dprintf("Loading PSF from PSFEx file %s\n", psfname);
        psf = psf_create(psfname);
    } else {
        char *tmpname = filename;
        char *etmpname = make_temp_filename("/tmp/img_errors_XXXXXX");

        if(psfname && !is_psfex)
            dprintf("Can't find file %s\n", psfname);

        if(is_preprocessed){
            tmpname = make_temp_filename("/tmp/img_processed_XXXXXX");
            image_dump_to_fits(image, tmpname);
        }

        image_dump_to_fits(errors, etmpname);

        dprintf("Extracting PSF from FITS file %s\n", tmpname);

        if(is_keep_psf || is_psfex){
            psf = psf_create_from_fits_and_save(tmpname, etmpname, psfname);
            dprintf("PSF for %s created and stored to %s\n", tmpname, psfname);
        } else {
            psf = psf_create_from_fits(tmpname, etmpname);
            dprintf("PSF for %s created\n", tmpname);
        }

        if(is_preprocessed)
            free(tmpname);
        free(etmpname);

        if(is_preprocessed)
            unlink(tmpname);
        unlink(etmpname);

        if(is_psfex)
            return EXIT_SUCCESS;
    }

    dprintf("PSF: %d x %d, degree = %d x %d, FWHM=%g\n", psf->width, psf->height, psf->degree, psf->zdegree, psf->fwhm);
    dprintf("PSF: x0=%g sx=%g y0=%g sy=%g z0=%g sz=%g\n", psf->x0, psf->sx, psf->y0, psf->sy, psf->z0, psf->sz);

    if(psf->fwhm <= 0 || !psf_is_positive(psf, image)){
        dprintf("Bad PSF FWHM!\n");
        return EXIT_FAILURE;
    } else if(max_fwhm > 0 && psf->fwhm > max_fwhm){
        dprintf("PSF FWHM is too large!\n");
        return EXIT_FAILURE;
    } else if(is_sex){
        /* SExtractor photometry */
        if(!outname || !*outname)
            outname = make_string("%s.txt", filename);
        if(!psfname || !file_exists_and_normal(psfname)){
            dprintf("PSF is necessary for SExtractor PSF photometry\n");
            return EXIT_FAILURE;
        }
        if(psf_photometry(image, psfname, outname, TRUE)){
            dprintf("SExtractor PSF photometry results stored to %s\n", outname);
            return EXIT_SUCCESS;
        } else
            return EXIT_FAILURE;
    } else {
        /* Main processing code */
        image_str *smooth = NULL;
        struct list_head peaks;
        int d;

        if(is_gauss){
            dprintf("Unsharping using Gaussian\n");
            smooth = image_unsharp(image, 1);
        } else {
            dprintf("Unsharping using exact PSF\n");
            smooth = image_unsharp_psf(image, psf);
        }

        if(is_keep_unsharp){
            char *unsharpname = make_string("%s.unsharp.fits", filename);
            char *smoothname = make_string("%s.smooth.fits", filename);
            image_str *smooth0 = image_smooth_psf(image, psf);

            image_dump_to_fits(smooth, unsharpname);
            dprintf("Unsharped image stored to %s\n", unsharpname);

            image_dump_to_fits(smooth0, smoothname);
            dprintf("Smooth image stored to %s\n", smoothname);

            image_delete(smooth0);
            free(smoothname);
            free(unsharpname);
        }

        if(is_unsharp){
            if(!outname)
                outname = make_string("%s.unsharp.fits", filename);

            image_dump_to_fits(smooth, outname);
            dprintf("Unsharped image written to %s\n", outname);
            return EXIT_SUCCESS;
        }

        if(width > 0 && height > 0){
            image = image_crop(image, x0, y0, x0 + width, y0 + height);
            errors = image_crop(errors, x0, y0, x0 + width, y0 + height);
            smooth = image_crop(smooth, x0, y0, x0 + width, y0 + height);
            mask = image_crop(mask, x0, y0, x0 + width, y0 + height);
        }

        if(debug){
            image_dump_to_fits(image, "out.image.fits");
            image_dump_to_fits(errors, "out.errors.fits");
            image_dump_to_fits(smooth, "out.smooth.fits");
            image_dump_to_fits(mask, "out.mask.fits");
        }

        dprintf("Extracting initial peaks\n");

        if(peaksname) {
            load_peaks(peaksname, image, errors, mask, is_radec, &peaks);
            dprintf("%d image peaks loaded from %s\n", list_length(&peaks), peaksname);
        } else {
            find_peaks(image, smooth, errors, mask, threshold, &peaks);
            dprintf("%d image peaks found\n", list_length(&peaks));
        }

        if(!outname || !*outname)
            outname = make_string("%s.txt", filename);

        if(debug){
            dump_peaks_to_file(&peaks, "out.peaks.full.txt", PEAK_INITIAL | PEAK_MEASURED | PEAK_FAILED);
        }

        if(is_simple){
            dump_peaks_to_file(&peaks, outname, PEAK_INITIAL | PEAK_MEASURED | PEAK_FAILED);
            dprintf("Initial peaks stored to %s\n", outname);

            if(debug)
                dump_peaks_to_file(&peaks, "out.peaks.result.txt", PEAK_INITIAL | PEAK_MEASURED | PEAK_FAILED);
        } else {
            int d;

            if(is_keep_simple){
                dump_peaks_to_file(&peaks, simplename, PEAK_INITIAL | PEAK_MEASURED | PEAK_FAILED);

                dprintf("Initial peaks stored to %s\n", simplename);
            }

            if(width > 0 && height > 0){
                psf->x0 -= x0;
                psf->y0 -= y0;
            }

            /* Increase errors for danger zone */
            {
                int d;

                for(d = 0; d < mask->width*mask->height; d++){
                    /* Low weight for saturated pixels */
                    if(mask->data[d] & FLAG_SATURATED){
                        int x = d % image->width;
                        int y = (d - x) / image->width;
                        int size = 0;
                        int x1;
                        int y1;

                        errors->double_data[d] = 100*saturation;

                        for(x1 = x - size; x1 <= x + size; x1++)
                            for(y1 = y - size; y1 <= y + size; y1++){
                                if(x1 >= 0 && x1 < image->width &&
                                   y1 >= 0 && y1 < image->height)
                                    PIXEL_DOUBLE(errors, x1, y1) = MAX(saturation, PIXEL_DOUBLE(errors, x1, y1));
                            }
                    }

                    if(mask->data[d] & FLAG_BAD)
                        errors->double_data[d] = saturation;
                }
            }

            psf_fit_peaks(image, errors, psf, mask, threshold, is_deblend, &peaks);

            {
                peak_str *peak;

                foreach(peak, peaks)
                    coords_get_ra_dec(&image->coords, peak->x, peak->y, &peak->ra, &peak->dec);
            }

            dump_peaks_to_file(&peaks, outname, PEAK_MEASURED);

            dprintf("Final peaks stored to %s\n", outname);

            if(is_keep_residuals){
                char *resname = make_string("%s.residuals.fits", filename);

                image_dump_to_fits(image, resname);
                dprintf("Residual image stored to %s\n", resname);

                free(resname);
            }
        }

        free_list(peaks);

        image_delete(smooth);
    }

    image_delete(image);
    image_delete(errors);
    image_delete(mask);

    psf_delete(psf);

    return EXIT_SUCCESS;
}
