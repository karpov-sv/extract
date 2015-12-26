#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "utils.h"

#include "psf.h"
#include "image.h"
#include "random.h"

int main(int argc, char **argv)
{
    char *outname = "out.fits";
    char *psfname = "example.psf";
    image_str *image = NULL;
    psf_str *psf = NULL;

    int size = 300;
    int N = 50000;
    double delta = 0.1;
    double bias = 0.0;
    double gain = 1.9;
    double readnoise = 0.0;
    double bg = 1000.0;
    double saturation = 0.0;

    parse_args(argc, argv,
               "size=%d", &size,
               "N=%d", &N,
               "delta=%lf", &delta,
               "bias=%lf", &bias,
               "gain=%lf", &gain,
               "readnoise=%lf", &readnoise,
               "saturation=%lf", &saturation,
               "bg=%lf", &bg,
               "psf=%s", &psfname,
               "%s", &outname,
               NULL);

    dprintf("Loading PSF from %s\n", psfname);
    psf = psf_create(psfname);

    if(!psf){
        dprintf("Can't load PSF from %s!\n", psfname);
        return EXIT_FAILURE;
    } else {
        int psfsize = floor(0.5*psf->width*psf->pix_step);

        dprintf("PSF: %d x %d, degree = %d, FWHM=%g\n", psf->width, psf->height, psf->degree, psf->fwhm);
        dprintf("PSF: x0=%g sx=%g y0=%g sy=%g\n", psf->x0, psf->sx, psf->y0, psf->sy);
        dprintf("PSF window size: %d x %d\n", psfsize, psfsize);
    }

    if(!saturation)
        saturation = 26000/gain;

    image = image_create_double(size, size);
    image_keyword_add_double(image, "BIAS", bias, "Bias used to generate noise");
    image_keyword_add_double(image, "GAIN", gain, "Gain used to generate noise");
    image_keyword_add_double(image, "READNOISE", readnoise, "Read-out noise level");
    image_keyword_add_double(image, "BACKGROUND", bg, "Background level");
    image_keyword_add_double(image, "SATURATION", saturation, "Saturation level");

    if(isatty(fileno(stderr)))
        dprintf("Filling the background...\n");

    image_fill(image, bg);

    if(isatty(fileno(stderr)))
        dprintf("Placing the stars...\n");

    {
        image_str *sampled = psf_sampled_image(psf, psf->x0, psf->y0);
        int psfsize = floor(0.5*psf->width*psf->pix_step);
        int d;

        srandom(time(NULL));

        for(d = 0; d < N; d++){
            double x0 = size*1.0*random()/RAND_MAX;
            double y0 = size*1.0*random()/RAND_MAX;
            double A = pow(10.0, 9.0 - 10.0*pow(1.0*random()/RAND_MAX, delta));
            image_str *kernel = psf_image(psf, sampled, x0-floor(x0), y0-floor(y0), 2*psfsize+1);
            int x1;
            int y1;

            for(x1 = MAX(0, x0-psfsize); x1 < MIN(image->width, x0+psfsize+1); x1++)
                for(y1 = MAX(0, y0-psfsize); y1 < MIN(image->height, y0+psfsize+1); y1++){
                    PIXEL_DOUBLE(image, x1, y1) += A*PIXEL_DOUBLE(kernel, x1 - x0 + psfsize, y1 - y0 + psfsize);;
                }

            if(isatty(fileno(stderr)))
                dprintf("\r %d / %d - %.2lf %.2lf - %g             \r", d, N, x0, y0, A);
        }

        if(isatty(fileno(stderr)))
            dprintf("\n");
    }

    if(isatty(fileno(stderr)))
        dprintf("Generating Poissonian noise...\n");

    {
        image_str *errors = image_errors(image, bias, gain, readnoise);
        int d;

        for(d = 0; d < image->width*image->height; d++){
            image->double_data[d] += random_gauss(errors->double_data[d]);

            image->double_data[d] = MIN(saturation, image->double_data[d]);
        }

        image_dump_to_fits(errors, "out.errors.fits");
    }

    if(isatty(fileno(stderr)))
        dprintf("Saving generated image to %s\n", outname);

    image_dump_to_fits(image, outname);

    return EXIT_SUCCESS;
}
