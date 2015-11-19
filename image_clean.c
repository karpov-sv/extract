#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "utils.h"

#include "image.h"

#define CLEAN_HORIZONTAL
/* #define CLEAN_VERTICAL */

void image_clean_stripes(image_str *image)
{
    int x;
    int y;
    double m = image_median(image);
    double mean = image_mean(image);

    /* dprintf("Image median=%g, mean=%g, mode=%g\n", m, mean, 2.5*m - 1.5*mean); */

    /* We process DOUBLE images only */
    if(image->type != IMAGE_DOUBLE)
        return;

#ifdef CLEAN_VERTICAL
    for(x = 0; x < image->width; x++){
        double *values = (double *)malloc(sizeof(double)*image->height);
        double m1;
        double m2;

        for(y = 0; y < image->height; y++)
            values[y] = PIXEL_DOUBLE(image, x, y);

        m1 = get_median(values, image->height/2);
        m2 = get_median(values + image->height/2, image->height/2);

        for(y = 0; y < image->height/2; y++)
            PIXEL_DOUBLE(image, x, y) -= m1 - m;

        for(y = image->height/2; y < image->height; y++)
            PIXEL_DOUBLE(image, x, y) -= m2 - m;

        free(values);
    }
#endif

#ifdef CLEAN_HORIZONTAL
    for(y = 0; y < image->height; y++){
        double m1 = get_median(image->double_data + y*image->width, image->width);

        for(x = 0; x < image->width; x++)
            PIXEL_DOUBLE(image, x, y) -= m1 - m;
    }
#endif
}
