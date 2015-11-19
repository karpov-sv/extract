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

    model->Nglobals = 6;
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
    double delta = 2.0; /* Positional uncertainty */
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

    if(peak->state == PEAK_MEASURED){
        model->paropts[idx + 0].fixed = 1;
        model->paropts[idx + 1].fixed = 1;
        model->paropts[idx + 2].fixed = 1;
    } else {
        model->paropts[idx + 0].limited[0] = 1;
        model->paropts[idx + 0].limited[1] = 1;
        model->paropts[idx + 0].limits[0] = peak->x - model->x0 - delta;
        model->paropts[idx + 0].limits[1] = peak->x - model->x0 + delta;

        model->paropts[idx + 1].limited[0] = 1;
        model->paropts[idx + 1].limited[1] = 1;
        model->paropts[idx + 1].limits[0] = peak->y - model->y0 - delta;
        model->paropts[idx + 1].limits[1] = peak->y - model->y0 + delta;

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

        imdata[d] = params[0] + x*params[1] + y*params[2] + x*y*params[3] + x*x*params[4] + y*y*params[5];
    }

    if(derivatives){
        if(derivatives[0])
            for(d = 0; d < Npoints; d++)
                /* Derivatives should be normalized inplace, as they are used for residuals only */
                derivatives[0][d] = -1.0/model->errors->double_data[d];
        if(derivatives[1])
            for(d = 0; d < Npoints; d++){
                int x = d % model->image->width;

                derivatives[1][d] = -x/model->errors->double_data[d];
            }
        if(derivatives[2])
            for(d = 0; d < Npoints; d++){
                int x = d % model->image->width;
                int y = (d - x) / model->image->width;

                derivatives[2][d] = -y/model->errors->double_data[d];
            }
        if(derivatives[3])
            for(d = 0; d < Npoints; d++){
                int x = d % model->image->width;
                int y = (d - x) / model->image->width;

                derivatives[3][d] = -x*y/model->errors->double_data[d];
            }
        if(derivatives[4])
            for(d = 0; d < Npoints; d++){
                int x = d % model->image->width;

                derivatives[4][d] = -x*x/model->errors->double_data[d];
            }
        if(derivatives[5])
            for(d = 0; d < Npoints; d++){
                int x = d % model->image->width;
                int y = (d - x) / model->image->width;

                derivatives[5][d] = -y*y/model->errors->double_data[d];
            }
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
            fprintf(file, "%g %g %g %g\n",
                    model->params[model->Nglobals + 3*d + 0],
                    model->params[model->Nglobals + 3*d + 1],
                    model->params[model->Nglobals + 3*d + 2],
                    model->peaks[d]->chisq);
        fclose(file);
    }
}

void model_fit(model_str *model)
{
    int Nparams = model->Nglobals + 3*model->Npeaks;

    memset(&model->result, 0, sizeof(mp_result));
    model->result.xerror = model->parerrs;
    model->result.resid = model->residuals;

    if(debug){
        image_dump_to_fits(model->image, "out.crop.fits");
        image_dump_to_fits(model->errors, "out.croperr.fits");
        model_dump_to_fits(model, "out.model.fits", "out.peaks.txt");
    }

    mpfit((mp_func)fn_mpfit, model->image->width*model->image->height, Nparams, model->params, model->paropts, NULL, model, &model->result);

    if(debug){
        model_dump_to_fits(model, "out.result.fits", NULL);
    }

    if(debug){
        dprintf("%g %g: %g -> %g, chisq = %g, (%d/%d, %d peaks), %d/%d runs, status=%d\n",
                model->x0 + 0.5*model->size, model->y0 + 0.5*model->size,
                model->result.orignorm, model->result.bestnorm,
                sqrt(model->result.bestnorm/(model->result.nfunc - model->result.nfree)),
                model->result.nfunc, model->result.nfree, model->Npeaks,
                model->result.niter, model->result.nfev, model->result.status);
    }
}

double model_peak_chisq(model_str *model, peak_str *peak)
{
    double chisq = 0;
    int Npoints = 0;
    int Nparams = 3;

    double size = 2*model->psf->fwhm;
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
    double edgesize = 2.0;
    int d;

    if(model->result.status <= 0)
        /* Fit has not converged properly, let's skip this peak without marking anything good or bad */
        return 0;

    for(d = 0; d < model->Npeaks; d++){
        peak_str *peak = model->peaks[d];

        peak->chisq = model_peak_chisq(model, peak);

        if(peak->state != PEAK_FAILED){
            peak->bg = model->params[0];
            peak->dbg = model->parerrs[0];

            peak->x = model->params[model->Nglobals + 3*d + 0] + model->x0;
            peak->dx = model->parerrs[model->Nglobals + 3*d + 0];
            peak->y = model->params[model->Nglobals + 3*d + 1] + model->y0;
            peak->dy = model->parerrs[model->Nglobals + 3*d + 1];
            peak->flux = model->params[model->Nglobals + 3*d + 2];
            peak->dflux = model->parerrs[model->Nglobals + 3*d + 2];

            if(peak->chisq > 0 && peak->chisq < 3 &&
               peak->x - model->x0 >= edgesize &&
               peak->x - model->x0 < model->image->width - edgesize &&
               peak->y - model->y0 >= edgesize &&
               peak->y - model->y0 < model->image->height - edgesize){
                /* Peak is near the center of model sub-image */
                if(peak->flux == 0 || peak->dflux == 0 || fabs(peak->dflux/peak->flux) > 1.0/model->threshold){
                    /* Flux accuracy is below the limit */
                    peak->state = PEAK_FAILED;
                } else {
                    /* Subtract the peak from image */
                    int x0 = floor(peak->x); /* Peak sub-window center */
                    int y0 = floor(peak->y);
                    double dx = peak->x - x0; /* Sub-pixel adjustments for PSF model */
                    double dy = peak->y - y0;
                    int x;
                    int y;
                    double I = 0; /* Normalization of the PSF model, to correct the flux */

                    image_delete(model->psfs[d]);
                    model->psfs[d] = psf_sampled_image(model->psf, peak->x, peak->y);
                    psf_image_fill(model->psf, model->psfs[d], model->stamp, dx, dy);

                    for(x = 0; x < model->stamp->width; x++)
                        for(y = 0; y < model->stamp->height; y++){
                            int x1 = x0 - 0.5*model->stamp->width + x;
                            int y1 = y0 - 0.5*model->stamp->height + y;

                            double value = PIXEL_DOUBLE(model->stamp, x, y);

                            I += value;

                            if(x1 >= 0 && x1 < image->width &&
                               y1 >= 0 && y1 < image->height)
                                PIXEL_DOUBLE(image, x1, y1) -= value*peak->flux;
                        }

                    /* Correct the flux and fluxerr using actual sub-pixel-adjusted PSF normalization */
                    peak->flux *= 1./I;
                    peak->dflux *= 1./I;

                    peak->state = PEAK_MEASURED;
                }

                Nprocessed ++;
            }

            if(debug)
                dprintf("%g %g - %g %g - %d\n\n", peak->x, peak->y, peak->chisq, fabs(peak->dflux/peak->flux), peak->state);


        }
    }

    return Nprocessed;
}

void subtract_peaks_from_image(image_str *image, psf_str *psf, struct list_head *peaks)
{
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
            double I = 0; /* Normalization of the PSF model, to correct the flux */

            image_str *psfi = psf_sampled_image(psf, peak->x, peak->y);
            image_str *stamp = image_create_double(30, 30);

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

            image_delete(stamp);
            image_delete(psfi);
        }
    }
}

void process_image_get_peaks(image_str *image, image_str *errors, psf_str *psf, double threshold, struct list_head *peaks)
{
    image_str *smooth = image_smooth(image, 1);
    struct kdtree *kd = kd_create(2);
    peak_str *peak;
    int size = 2*floor(0.5*psf->width*psf->pix_step); /* It should be always even */
    int edgesize = 2.0;
    int Npeaks = 0;
    int Nprocessed = 0;

    find_peaks(image, smooth, peaks);
    Npeaks = list_length(peaks);
    dprintf("%d initial peaks found\n", Npeaks);

    dprintf("Window size is %d x %d pixels\n", size, size);

    if(debug){
        image_dump_to_fits(image, "out.image.fits");
        image_dump_to_fits(errors, "out.errors.fits");
        image_dump_to_fits(smooth, "out.smooth.fits");

        dump_peaks_to_file_full(peaks, "out.peaks.full.txt");
    }

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
            int xc = MAX(0.5*size + edgesize, MIN(image->width - 0.5*size - edgesize, floor(peak->x)));
            int yc = MAX(0.5*size + edgesize, MIN(image->height - 0.5*size - edgesize, floor(peak->y)));
            struct kdres *res = kd_nearest_range2(kd, xc, yc, 0.75*size);
            int N = 0;

            image_str *image_cropped = image_crop(image, xc - 0.5*size, yc - 0.5*size, xc + 0.5*size, yc + 0.5*size);
            image_str *errors_cropped = image_crop(errors, xc - 0.5*size, yc - 0.5*size, xc + 0.5*size, yc + 0.5*size);

            model_str *model = model_create(image_cropped, errors_cropped, psf, xc - 0.5*size, yc - 0.5*size, size);

            model->threshold = threshold;

            while(!kd_res_end(res)){
                peak_str *peak1 = kd_res_item_data(res);

                if((peak1->state == PEAK_INITIAL &&
                    peak1->x >= xc - 0.5*size + 0*edgesize && peak1->x < xc + 0.5*size - 0*edgesize &&
                    peak1->y >= yc - 0.5*size + 0*edgesize && peak1->y < yc + 0.5*size - 0*edgesize)){
                    /* Should include the peak in model */
                    model_add_peak(model, peak1);

                    N ++;
                }

                kd_res_next(res);
            }

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
        dump_peaks_to_file_measured(peaks, "out.peaks.result.txt");
        dump_peaks_to_file_failed(peaks, "out.peaks.failed.txt");
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

    image_delete(smooth);
}

int main(int argc, char **argv)
{
    char *filename = NULL;
    char *psfname = NULL;
    char *outname = NULL;
    int x0 = 0;
    int y0 = 0;
    int width = 0;
    int height = 0;

    int is_psfex = FALSE;
    int is_simple = FALSE;

    double bias = 100.0;
    double gain = 0.6;
    double readnoise = 1.4;
    double threshold = 5;

    double max_fwhm = 6.0;

    parse_args(argc, argv,
               "x0=%d", &x0,
               "y0=%d", &y0,
               "width=%d", &width,
               "height=%d", &height,
               "psf=%s", &psfname,
               "-psfex", &is_psfex,
               "-simple", &is_simple,
               "-debug", &debug,

               "bias=%lf", &bias,
               "gain=%lf", &gain,
               "readnoise=%lf", &readnoise,
               "threshold=%lf", &threshold,

               "max_fwhm=%lf", &max_fwhm,

               "%s", &filename,
               "%s", &outname,
               NULL);

    if(!filename)
        return EXIT_SUCCESS;

    if(is_psfex){
        psf_str *psf = NULL;

        if(!outname)
            outname = make_string("%s.psf", filename);

        psf = psf_create_from_fits_and_save(filename, outname);

        dprintf("PSF for %s created and stored to %s\n", filename, outname);

        psf_delete(psf);
    } else if(is_simple){
        image_str *image = image_create_from_fits(filename);
        image_str *smooth = NULL;
        struct list_head peaks;

        if(image->type != IMAGE_DOUBLE){
            /* The code expects IMAGE_DOUBLE */
            image_str *new = image_convert_to_double(image);

            image_delete(image);
            image = new;
        }

        if(width > 0 && height > 0){
            image = image_crop(image, x0, y0, x0 + width, y0 + height);
        }

        smooth = image_smooth(image, 1);

        find_peaks(image, smooth, &peaks);

        dump_peaks_to_file_full(&peaks, outname);

        free_list(peaks);
        image_delete(smooth);
        image_delete(image);
    } else {
        image_str *image = image_create_from_fits(filename);
        image_str *errors = NULL;
        psf_str *psf = NULL;
        struct list_head peaks;

        if(!psfname)
            psfname = make_string("%s.psf", filename);

        if(file_exists_and_normal(psfname)) {
            dprintf("Loading PSF from PSFEx file %s\n", psfname);
            psf = psf_create(psfname);
        } else {
            if(psfname)
                dprintf("Can't find file %s\n", psfname);
            dprintf("Extracting PSF from FITS file %s\n", filename);
            psf = psf_create_from_fits(filename);
        }

        if(width > 0 && height > 0){
            image = image_crop(image, x0, y0, x0 + width, y0 + height);

            psf->x0 -= x0;
            psf->y0 -= y0;
        }

        if(image->type != IMAGE_DOUBLE){
            /* The code expects IMAGE_DOUBLE */
            image_str *new = image_convert_to_double(image);

            image_delete(image);
            image = new;
        }

        dprintf("PSF: %d x %d, degree = %d, FWHM=%g\n", psf->width, psf->height, psf->degree, psf->fwhm);
        dprintf("PSF: x0=%g sx=%g y0=%g sy=%g\n", psf->x0, psf->sx, psf->y0, psf->sy);

        dprintf("%s: %d x %d, min = %g, max = %g\n", filename, image->width, image->height,
                image_min_value(image), image_max_value(image));

        if(max_fwhm > 0 && psf->fwhm > max_fwhm){
            dprintf("PSF FWHM is too large!\n");
            return EXIT_FAILURE;
        }

        errors = image_errors(image, bias, gain, readnoise);

        process_image_get_peaks(image, errors, psf, threshold, &peaks);

        dump_peaks_to_file(&peaks, outname);

        free_list(peaks);

        psf_delete(psf);
        image_delete(errors);
        image_delete(image);
    }

    return EXIT_SUCCESS;
}