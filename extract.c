#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "utils.h"

#include "extract.h"
#include "psf.h"
#include "mpfit.h"
#include "kdtree.h"

static int debug = FALSE;
static int peak_to_stop = -1;

typedef struct {
    image_str *image;
    image_str *errors;
    psf_str *psf;

    int x0;
    int y0;

    int size;

    double threshold; /* Minimal signal/noise for fitted flux */

    int Nglobals;

    int Npeaks;
    peak_str **peaks;

    double *params;
    double *parerrs;
    mp_par *paropts;

    image_str **psfs;
    image_str *stamp;

    mp_result result;
    double *residuals;
} model_str;

model_str *model_create(image_str *image, image_str *errors, psf_str *psf, int x0, int y0, int size)
{
    model_str *model = (model_str *)malloc(sizeof(model_str));
    int d;

    model->image = image;
    model->errors = errors;
    model->psf = psf;

    model->x0 = x0;
    model->y0 = y0;
    model->size = size;

    model->threshold = 2;

    model->Nglobals = 3;//6;
    model->Npeaks = 0;

    model->peaks = NULL;
    model->params = (double *)calloc(model->Nglobals, sizeof(double));
    model->parerrs = (double *)calloc(model->Nglobals, sizeof(double));
    model->paropts = (mp_par *)calloc(model->Nglobals, sizeof(mp_par));

    model->psfs = NULL;
    model->stamp = image_create_double(size, size);

    model->residuals = calloc(model->image->width*model->image->height, sizeof(double));

    /* Background */
    for(d = 0; d < model->Nglobals; d++){
        model->params[d] = 0;
        model->paropts[d].side = 3; /* Analytical derivative available */
    }
    model->params[0] = image_min_value(model->image);

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

    if(model->psfs){
        int d;

        for(d = 0; d < model->Npeaks; d++)
            image_delete(model->psfs[d]);

        free(model->psfs);
    }

    image_delete(model->stamp);

    free(model);
}

void model_add_peak(model_str *model, peak_str *peak)
{
    /* double delta = 1.0; /\* Positional uncertainty *\/ */
    double delta = 0.2*model->psf->fwhm; /* Positional uncertainty */
    int idx = model->Nglobals + 3*model->Npeaks;

    model->params = realloc(model->params, sizeof(double)*(idx + 3));
    model->parerrs = realloc(model->parerrs, sizeof(double)*(idx + 3));
    model->paropts = realloc(model->paropts, sizeof(mp_par)*(idx + 3));
    memset(model->paropts + idx, 0, 3*sizeof(mp_par));

    model->peaks = realloc(model->peaks, sizeof(peak_str *)*(model->Npeaks + 1));
    model->psfs = realloc(model->psfs, sizeof(image_str*)*(model->Npeaks + 1));

    model->peaks[model->Npeaks] = peak;

    model->psfs[model->Npeaks] = psf_sampled_image(model->psf, peak->x, peak->y);
    psf_image_fill(model->psf, model->psfs[model->Npeaks], model->stamp, 0, 0);

    model->params[idx + 0] = peak->x - model->x0;
    model->params[idx + 1] = peak->y - model->y0;
    if(peak->flux > 0)
        model->params[idx + 2] = peak->flux;
    else
        model->params[idx + 2] = peak->A/image_max_value(model->stamp);

    model->params[0] = peak->bg;

    if(model->params[idx + 2] < 0)
        model->params[idx + 2] = 0;

    model->paropts[idx + 0].step = 0.001;
    model->paropts[idx + 1].step = 0.001;
    model->paropts[idx + 2].relstep = 0.001;

    if(peak->state == PEAK_MEASURED){
        /* For already measured peaks we fix all parameters */
        model->paropts[idx + 0].fixed = 1;
        model->paropts[idx + 1].fixed = 1;
        model->paropts[idx + 2].fixed = 1;
    } else {
        if(model->params[idx + 0] < 0 || model->params[idx + 0] >= model->size ||
           model->params[idx + 1] < 0 || model->params[idx + 1] >= model->size){
            /* For peaks outside the the subimage we fix the positions */
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
        model->paropts[idx + 2].side = 3; /* Analytical derivative available */
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

    if(derivatives){
        if(derivatives[0])
            for(d = 0; d < Npoints; d++)
                /* Derivatives should be normalized inplace, as they are used for residuals only */
                derivatives[0][d] = -1.0/model->errors->double_data[d];
        if(derivatives[1])
            for(d = 0; d < Npoints; d++){
                int x = d % model->image->width;

                derivatives[1][d] = -x*1.0/model->errors->double_data[d];
            }
        if(derivatives[2])
            for(d = 0; d < Npoints; d++){
                int x = d % model->image->width;
                int y = (d - x) / model->image->width;

                derivatives[2][d] = -y*1.0/model->errors->double_data[d];
            }
        /* if(derivatives[3]) */
        /*     for(d = 0; d < Npoints; d++){ */
        /*         int x = d % model->image->width; */
        /*         int y = (d - x) / model->image->width; */

        /*         derivatives[3][d] = -x*y/model->errors->double_data[d]; */
        /*     } */
        /* if(derivatives[4]) */
        /*     for(d = 0; d < Npoints; d++){ */
        /*         int x = d % model->image->width; */

        /*         derivatives[4][d] = -x*x/model->errors->double_data[d]; */
        /*     } */
        /* if(derivatives[5]) */
        /*     for(d = 0; d < Npoints; d++){ */
        /*         int x = d % model->image->width; */
        /*         int y = (d - x) / model->image->width; */

        /*         derivatives[5][d] = -y*y/model->errors->double_data[d]; */
        /*     } */
    }

    for(d = 0; d < model->Npeaks; d++){
        int idx = model->Nglobals + 3*d;
        double x = params[idx + 0];
        double y = params[idx + 1];
        double A = params[idx + 2];
        /* Sub-pixel shift of PSF center relative to stamp position */
        double dx = x - floor(model->peaks[d]->x - model->x0);
        double dy = y - floor(model->peaks[d]->y - model->y0);
        /* Origin of PSF stamp inside the image */
        int x0 = floor(model->peaks[d]->x - model->x0) - 0.5*model->stamp->width;
        int y0 = floor(model->peaks[d]->y - model->y0) - 0.5*model->stamp->height;
        int x1;
        int y1;

        psf_image_fill(model->psf, model->psfs[d], model->stamp, dx, dy);

        for(y1 = MAX(0, -y0); y1 < MIN(model->stamp->height, model->image->height - y0); y1++)
            for(x1 = MAX(0, -x0); x1 < MIN(model->stamp->width, model->image->width - x0); x1++){
                double value = PIXEL_DOUBLE(model->stamp, x1, y1);
                int idx1 = (y0 + y1)*model->image->width + (x0 + x1);

                imdata[idx1] += value*A;
            }

        if(derivatives && derivatives[idx + 2])
            for(y1 = MAX(0, -y0); y1 < MIN(model->stamp->height, model->image->height - y0); y1++)
                for(x1 = MAX(0, -x0); x1 < MIN(model->stamp->width, model->image->width - x0); x1++){
                    double value = PIXEL_DOUBLE(model->stamp, x1, y1);
                    int idx1 = (y0 + y1)*model->image->width + (x0 + x1);

                    derivatives[idx + 2][idx1] = -value/model->errors->double_data[idx];
                }
    }
}

int fn_mpfit(int Npoints, int Nparams, double *params, double *residuals, double **derivatives, void *data)
{
    model_str *model = (model_str *)data;
    int d;

    model_fill_image(model, params, residuals, derivatives);

    for(d = 0; d < Npoints; d++){
        residuals[d] = (model->image->double_data[d] - residuals[d])/model->errors->double_data[d];
    }

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
                    model->params[model->Nglobals + 3*d + 0],
                    model->params[model->Nglobals + 3*d + 1],
                    model->params[model->Nglobals + 3*d + 2],
                    model->peaks[d]->chisq,
                    model->peaks[d]->id);
        fclose(file);
    }
}

void model_fit(model_str *model)
{
    int Nparams = model->Nglobals + 3*model->Npeaks;
    struct mp_config_struct conf;

    memset(&conf, 0, sizeof(conf));
    conf.ftol = 1e-16;
    //conf.xtol = 1e-16;
    //conf.gtol = 1e-16;
    conf.maxiter = 1000;

    memset(&model->result, 0, sizeof(mp_result));
    model->result.xerror = model->parerrs;
    model->result.resid = model->residuals;

    if(debug){
        image_dump_to_fits(model->image, "out.crop.fits");
        image_dump_to_fits(model->errors, "out.croperr.fits");
        model_dump_to_fits(model, "out.model.fits", "out.peaks0.txt");
    }

    mpfit((mp_func)fn_mpfit, model->image->width*model->image->height, Nparams, model->params, model->paropts, &conf, model, &model->result);

    if(debug){
        model_dump_to_fits(model, "out.result.fits", "out.peaks.txt");
        image_dump_to_fits(model->errors, "out.ecrop.fits");
    }

    if(debug){
        dprintf("%g %g: %g -> %g, chisq = %g, (%d/%d, %d peaks), %d/%d runs, status=%d\n",
                model->x0 + 0.5*model->size, model->y0 + 0.5*model->size,
                model->result.orignorm, model->result.bestnorm,
                sqrt(model->result.bestnorm/(model->result.nfunc - model->result.nfree)),
                model->result.nfunc, model->result.nfree, model->Npeaks,
                model->result.niter, model->result.nfev, model->result.status);
    }

    /* if(model->Npeaks > 1) */
    /*     exit(1); */
}

double model_peak_chisq(model_str *model, peak_str *peak)
{
    double chisq = 0;
    int Npoints = 0;
    int Nparams = 3;

    double size = 5*model->psf->fwhm;
    int x;
    int y;

    for(x = MAX(0, round(peak->x - model->x0 - 0.5*size)); x < MIN(model->image->width, round(peak->x - model->x0 + 0.5*size)); x++)
        for(y = MAX(0, round(peak->y - model->y0 - 0.5*size)); y < MIN(model->image->height, round(peak->y - model->y0 + 0.5*size)); y++){
            chisq += pow(model->residuals[y*model->image->width + x], 2);
            Npoints ++;
        }

    if(Npoints > Nparams)
        return sqrt(chisq/(Npoints - Nparams));
    else
        return 0;
}

int model_update_peaks(model_str *model, image_str *image)
{
    int Nprocessed = 0; /* Number of either fitted or failed peaks */
    double edgesize = 1.0;
    int d;

    int is_stop = FALSE;

    if(model->result.status <= 0 || model->result.status > 4)
        /* Fit has not converged properly, let's skip this peak without marking anything good or bad */
        return 0;

    for(d = 0; d < model->Npeaks; d++){
        peak_str *peak = model->peaks[d];

        peak->chisq = model_peak_chisq(model, peak);

        if(peak->state == PEAK_INITIAL){
            peak->bg = model->params[0];
            peak->dbg = model->parerrs[0];

            peak->x = model->params[model->Nglobals + 3*d + 0] + model->x0;
            peak->dx = model->parerrs[model->Nglobals + 3*d + 0];
            peak->y = model->params[model->Nglobals + 3*d + 1] + model->y0;
            peak->dy = model->parerrs[model->Nglobals + 3*d + 1];
            peak->flux = model->params[model->Nglobals + 3*d + 2];
            peak->dflux = model->parerrs[model->Nglobals + 3*d + 2];

            if(peak->chisq > 0 && peak->chisq < 3.0 &&
               peak->x - model->x0 >= edgesize &&
               peak->x - model->x0 < model->image->width - edgesize &&
               peak->y - model->y0 >= edgesize &&
               peak->y - model->y0 < model->image->height - edgesize){
                /* Peak is near the center of model sub-image */
                if(peak->flux == 0 || peak->dflux == 0 || fabs(peak->dflux/peak->flux) > 1.0/model->threshold){
                    /* Flux accuracy is below the limit */
                    peak->state = PEAK_FAILED;
                } else if(peak->state != PEAK_MEASURED){
                    /* Subtract the peak from image */
                    /* int x0 = floor(peak->x); /\* Peak sub-window center *\/ */
                    /* int y0 = floor(peak->y); */
                    /* double dx = peak->x - x0; /\* Sub-pixel adjustments for PSF model *\/ */
                    /* double dy = peak->y - y0; */
                    /* int x; */
                    /* int y; */
                    /* double I = 0; /\* Normalization of the PSF model, to correct the flux *\/ */

                    /* image_delete(model->psfs[d]); */
                    /* model->psfs[d] = psf_sampled_image(model->psf, peak->x, peak->y); */
                    /* psf_image_fill(model->psf, model->psfs[d], model->stamp, dx, dy); */

                    /* for(x = 0; x < model->stamp->width; x++) */
                    /*     for(y = 0; y < model->stamp->height; y++){ */
                    /*         int x1 = x0 - 0.5*model->stamp->width + x; */
                    /*         int y1 = y0 - 0.5*model->stamp->height + y; */

                    /*         double value = PIXEL_DOUBLE(model->stamp, x, y); */

                    /*         I += value; */

                    /*         if(x1 >= 0 && x1 < image->width && */
                    /*            y1 >= 0 && y1 < image->height) */
                    /*             PIXEL_DOUBLE(image, x1, y1) -= 0*value*peak->flux; */
                    /*     } */

                    /* Correct the flux and fluxerr using actual sub-pixel-adjusted PSF normalization */
                    /* peak->flux *= 1./I; */
                    /* peak->dflux *= 1./I; */

                    peak->state = PEAK_MEASURED;
                }

                Nprocessed ++;
            }

            if(debug && peak->id == peak_to_stop &&
               peak->x - model->x0 >= edgesize &&
               peak->x - model->x0 < model->image->width - edgesize &&
               peak->y - model->y0 >= edgesize &&
               peak->y - model->y0 < model->image->height - edgesize){
                is_stop = TRUE;
            }

            if(debug)
                dprintf("peak %d : %g %g - %g %g - %g %g - %d\n\n",
                        peak->id, peak->x, peak->y, peak->x - model->x0, peak->y - model->y0,
                        peak->chisq, fabs(peak->dflux/peak->flux), peak->state);
        }
    }

    if(debug && is_stop){
        dprintf("Stopping at peak %d\n", peak_to_stop);
        exit(1);
    }


    return Nprocessed;
}

void subtract_peaks_from_image(image_str *image, psf_str *psf, struct list_head *peaks, int size)
{
    image_str *stamp = image_create_double(size, size);
    peak_str *peak;

    dprintf("Subtracting the model\n");

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

            psf_image_fill(psf, psfi, stamp, dx, dy);

            for(x = 0; x < stamp->width; x++)
                for(y = 0; y < stamp->height; y++){
                    int x1 = x0 - 0.5*stamp->width + x;
                    int y1 = y0 - 0.5*stamp->height + y;

                    double value = PIXEL_DOUBLE(stamp, x, y);

                    I += value;

                    if(x1 >= 0 && x1 < image->width &&
                       y1 >= 0 && y1 < image->height)
                        PIXEL_DOUBLE(image, x1, y1) -= value*peak->flux;
                }

            peak->flux *= 1.0/I;
            peak->dflux *= 1.0/I;

            image_delete(psfi);
        }
    }

    image_delete(stamp);
}

void psf_fit_peaks(image_str *image, image_str *errors, psf_str *psf, double threshold, struct list_head *peaks)
{
    struct kdtree *kd = kd_create(2);
    peak_str *peak;
    int size = 2*floor(0.5*psf->width*psf->pix_step); /* It should be always even */
    int edgesize = 1.0;
    int Npeaks = 0;
    int Nprocessed = 0;

    Npeaks = list_length(peaks);

    dprintf("Window size is %d x %d pixels\n", size, size);

    /* Fill kd-tree with all peaks */
    foreach(peak, *peaks){
        kd_insert2(kd, peak->x, peak->y, peak);
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

            image_str *image_cropped = image_crop(image, xc - 0.5*size, yc - 0.5*size, xc + 0.5*size, yc + 0.5*size);
            image_str *errors_cropped = image_crop(errors, xc - 0.5*size, yc - 0.5*size, xc + 0.5*size, yc + 0.5*size);

            model_str *model = model_create(image_cropped, errors_cropped, psf, xc - 0.5*size, yc - 0.5*size, size);

            model->threshold = threshold;

            if(debug)
                dprintf("\nmodel for peak %d: %g %g - %d %d\n", peak->id, peak->x, peak->y, xc, yc);

            while(!kd_res_end(res)){
                peak_str *peak1 = kd_res_item_data(res);

                if((peak1->state != PEAK_FAILED &&
                    peak1->x >= xc - 0.7*size + 0*edgesize && peak1->x < xc + 0.7*size - 0*edgesize &&
                    peak1->y >= yc - 0.7*size + 0*edgesize && peak1->y < yc + 0.7*size - 0*edgesize)){
                    /* Should include the peak in model */
                    model_add_peak(model, peak1);

                    N ++;
                }

                kd_res_next(res);
            }

            /* dprintf("\nmodel at: %d %d for peak at %g %g\n", model->x0, model->y0, peak->x, peak->y); */

            /* Here we should have at least one unfitted peak */
            model_fit(model);

            /* TODO: Here we should somehow analyze fit results and perform
               the detection and injection of additional peaks into model */

            Nprocessed += model_update_peaks(model, image);

            if(isatty(fileno(stderr)))
                dprintf("\r  %d / %d - %d peaks        \r", Nprocessed, Npeaks, N);

            model_delete(model);

            image_delete(image_cropped);
            image_delete(errors_cropped);

            kd_res_free(res);
        }

    if(isatty(fileno(stderr)))
        dprintf("\n");

    if(debug){
        dump_peaks_to_file(peaks, "out.peaks.result.txt", PEAK_MEASURED);
        dump_peaks_to_file(peaks, "out.peaks.failed.txt", PEAK_FAILED);
        dump_peaks_to_file(peaks, "out.peaks.initial.txt", PEAK_INITIAL);

        subtract_peaks_from_image(image, psf, peaks, size);
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
    int x0 = 0;
    int y0 = 0;
    int width = 0;
    int height = 0;

    int is_psfex = FALSE;
    int is_simple = FALSE;
    int is_gauss = FALSE;
    int is_keep_psf = FALSE;
    int is_keep_simple = FALSE;
    int is_keep_preprocessed = FALSE;
    int is_keep_mask = FALSE;
    int is_unsharp = FALSE;
    int is_header = FALSE;

    int is_preprocessed = FALSE;

    double bias = 0.0;
    double gain = 0.0;
    double readnoise = 0.0;
    double threshold = 5.0;

    double saturation = 0.0;
    double welldepth = 23000.0;

    double max_fwhm = 6.0;

    image_str *image = NULL;
    image_str *mask = NULL;
    psf_str *psf = NULL;

    parse_args(argc, argv,
               "x0=%d", &x0,
               "y0=%d", &y0,
               "width=%d", &width,
               "height=%d", &height,
               "psf=%s", &psfname,
               "-psfex", &is_psfex,
               "-simple", &is_simple,
               "-gauss", &is_gauss,
               "-keep_psf", &is_keep_psf,
               "-keep_simple", &is_keep_simple,
               "-keep_preprocessed", &is_keep_preprocessed,
               "-keep_mask", &is_keep_mask,
               "-unsharp", &is_unsharp,
               "-header", &is_header,
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

               "dark=%s", &darkname,
               "flat=%s", &flatname,
               "mask=%s", &maskname,

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

    dprintf("%s: %d x %d, min = %g, max = %g, mean = %g, sigma = %g\n", filename, image->width, image->height,
            image_min_value(image), image_max_value(image), image_mean(image), image_sigma(image));

    /* If no GAIN provided on command line, try to guess it from the header */
    if(is_header && !gain){
        if(image_keyword_find(image, "GAIN")){
            gain = image_keyword_get_double(image, "GAIN");

            dprintf("Setting GAIN to %g from GAIN header keyword\n", gain);
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

            gain *= avg;
            readnoise /= sqrt(avg);
            welldepth *= avg;

            dprintf("The image was AVERAGED over %d frames, corrected the parameters correspondingly\n", avg);
        }

        /* TODO: estimate the saturation level here */
    }

    if(!saturation){
        saturation = welldepth/gain;

        dprintf("Setting SATURATION to %g from %g depth well capacity\n", saturation, welldepth);
    }

    /* Fill PSFEx parameters */
    if(gain && !psf_gain)
        psf_gain = gain;

    if(saturation)
        psf_satur_level = saturation;

    dprintf("bias=%g gain=%g readnoise=%g saturation=%g\n", bias, gain, readnoise, saturation);

    /* Various filenames */
    if(is_psfex && outname)
        psfname = outname;
    else if(!psfname)
        psfname = make_string("%s.psf", filename);

    simplename = make_string("%s.simple", filename);
    preprocessedname = make_string("%s.processed.fits", filename);

    /* MASK holding saturation/error flags */
    if(maskname){
        mask = image_create_from_fits(maskname);

        if(mask)
            dprintf("MASK frame read from %s\n", maskname);
    } else {
        mask = image_create(image->width, image->height);
    }

    /* We should mask saturated pixels before dark/flat corrections */
    {
        int d;

        for(d = 0; d < image->width*image->height; d++)
            if(image->double_data[d] >= saturation)
                mask->data[d] = FLAG_SATURATED;

        if(is_keep_mask){
            maskname = make_string("%s.mask.fits", filename);

            image_dump_to_fits(mask, maskname);
        }

        dprintf("Saturated pixels masked\n");
    }

    /* Pre-process the image if DARK and FLAT are provided */
    if(darkname){
        image_str *dark = image_create_from_fits(darkname);
        image_str *flat = flatname ? image_create_from_fits(flatname) : NULL;
        int d;

        if(dark)
            dprintf("DARK frame read from %s\n", darkname);

        if(flat)
            dprintf("FLAT frame read from %s\n", flatname);

        if(dark && flat){
            for(d = 0; d < image->width*image->height; d++)
                image->double_data[d] = (image->double_data[d] - dark->double_data[d])/flat->double_data[d];
        } else if(dark){
            for(d = 0; d < image->width*image->height; d++)
                image->double_data[d] = image->double_data[d] - dark->double_data[d];
        }

        /* Make saturated pixels saturated */
        for(d = 0; d < image->width*image->height; d++)
            if(mask->data[d] & FLAG_SATURATED)
                image->double_data[d] = saturation;

        if(dark)
            image_delete(dark);

        if(flat)
            image_delete(flat);

        if(is_keep_preprocessed)
            image_dump_to_fits(image, preprocessedname);

        is_preprocessed = TRUE;
        dprintf("Image pre-processed\n");
    }

    if(file_exists_and_normal(psfname) && !is_psfex){
        dprintf("Loading PSF from PSFEx file %s\n", psfname);
        psf = psf_create(psfname);
    } else {
        char *tmpname = filename;

        if(psfname && !is_psfex)
            dprintf("Can't find file %s\n", psfname);

        if(is_preprocessed){
            tmpname = make_temp_filename("/tmp/img_processed_XXXXXX");
            image_dump_to_fits(image, tmpname);
        }

        dprintf("Extracting PSF from FITS file %s\n", tmpname);

        if(is_keep_psf || is_psfex){
            psf = psf_create_from_fits_and_save(tmpname, psfname);
            dprintf("PSF for %s created and stored to %s\n", tmpname, psfname);
        } else {
            psf = psf_create_from_fits(tmpname);
            dprintf("PSF for %s created\n", tmpname);
        }

        if(is_preprocessed)
            unlink(tmpname);

        if(is_psfex)
            return EXIT_SUCCESS;
    }

    dprintf("PSF: %d x %d, degree = %d, FWHM=%g\n", psf->width, psf->height, psf->degree, psf->fwhm);
    /* dprintf("PSF: x0=%g sx=%g y0=%g sy=%g\n", psf->x0, psf->sx, psf->y0, psf->sy); */

    if(psf->fwhm <= 0){
        dprintf("Bad PSF FWHM!\n");
        return EXIT_FAILURE;
    } else if(max_fwhm > 0 && psf->fwhm > max_fwhm){
        dprintf("PSF FWHM is too large!\n");
        return EXIT_FAILURE;
    } else {
        /* Main processing code */
        image_str *smooth = NULL;
        image_str *errors = NULL;
        struct list_head peaks;

        if(is_gauss){
            dprintf("Unsharping using Gaussian\n");
            smooth = image_unsharp(image, 1);
        } else {
            dprintf("Unsharping using exact PSF\n");
            smooth = image_unsharp_psf(image, psf);
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
            smooth = image_crop(smooth, x0, y0, x0 + width, y0 + height);
            mask = image_crop(mask, x0, y0, x0 + width, y0 + height);
        }

        if(gain > 0)
            errors = image_errors(image, bias, gain, readnoise);
        else {
            dprintf("No GAIN provided, estimating the noise from image itself\n");
            errors = image_errors_estimate(image);
        }

        dprintf("Extracting initial peaks\n");

        find_peaks(image, smooth, errors, mask, threshold, &peaks);

        dprintf("%d image peaks found\n", list_length(&peaks));

        if(debug){
            image_dump_to_fits(image, "out.image.fits");
            image_dump_to_fits(errors, "out.errors.fits");
            image_dump_to_fits(smooth, "out.smooth.fits");
            image_dump_to_fits(mask, "out.mask.fits");

            dump_peaks_to_file(&peaks, "out.peaks.full.txt", PEAK_INITIAL | PEAK_MEASURED | PEAK_FAILED);
        }

        if(is_simple){
            dump_peaks_to_file(&peaks, outname, PEAK_INITIAL | PEAK_MEASURED | PEAK_FAILED);

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

            /* Low weight for saturated pixels */
            for(d = 0; d < image->width*image->height; d++)
                if(mask->data[d] & FLAG_SATURATED)
                    errors->double_data[d] = saturation;

            psf_fit_peaks(image, errors, psf, threshold, &peaks);

            dump_peaks_to_file(&peaks, outname, PEAK_MEASURED);
        }

        free_list(peaks);

        image_delete(errors);
        image_delete(smooth);
    }

    image_delete(image);
    psf_delete(psf);

    return EXIT_SUCCESS;
}
