#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#include "utils.h"

#include "image.h"

#define IS_MASK_NANS FALSE

/* static int N = 40; */
/* static double h[81] = */
/*     { -1.07159639e-19,   4.77863127e-06,   2.12074335e-05, */
/*          5.27535252e-05,   1.03393196e-04,   1.77684646e-04, */
/*          2.80813440e-04,   4.18604169e-04,   5.97493094e-04, */
/*          8.24458160e-04,   1.10690462e-03,   1.45250655e-03, */
/*          1.86900658e-03,   2.36397847e-03,   2.94455879e-03, */
/*          3.61715632e-03,   4.38714879e-03,   5.25857825e-03, */
/*          6.23385676e-03,   7.31349474e-03,   8.49586400e-03, */
/*          9.77700676e-03,   1.11505011e-02,   1.26073911e-02, */
/*          1.41361892e-02,   1.57229540e-02,   1.73514459e-02, */
/*          1.90033601e-02,   2.06586320e-02,   2.22958095e-02, */
/*          2.38924829e-02,   2.54257611e-02,   2.68727807e-02, */
/*          2.82112332e-02,   2.94198955e-02,   3.04791454e-02, */
/*          3.13714492e-02,   3.20818025e-02,   3.25981133e-02, */
/*          3.29115131e-02,   3.30165867e-02,   3.29115131e-02, */
/*          3.25981133e-02,   3.20818025e-02,   3.13714492e-02, */
/*          3.04791454e-02,   2.94198955e-02,   2.82112332e-02, */
/*          2.68727807e-02,   2.54257611e-02,   2.38924829e-02, */
/*          2.22958095e-02,   2.06586320e-02,   1.90033601e-02, */
/*          1.73514459e-02,   1.57229540e-02,   1.41361892e-02, */
/*          1.26073911e-02,   1.11505011e-02,   9.77700676e-03, */
/*          8.49586400e-03,   7.31349474e-03,   6.23385676e-03, */
/*          5.25857825e-03,   4.38714879e-03,   3.61715632e-03, */
/*          2.94455879e-03,   2.36397847e-03,   1.86900658e-03, */
/*          1.45250655e-03,   1.10690462e-03,   8.24458160e-04, */
/*          5.97493094e-04,   4.18604169e-04,   2.80813440e-04, */
/*          1.77684646e-04,   1.03393196e-04,   5.27535252e-05, */
/*          2.12074335e-05,   4.77863127e-06,  -1.07159639e-19}; */

static int N = 1;
static double h[3] = {0.0, 0.0, 0.0};

#define CLEAN_HORIZONTAL
/* #define CLEAN_VERTICAL */

void image_clean_stripes(image_str *image)
{
    int x;
    int y;
    double m = image_median(image);
    double mean = image_mean(image);
    double *values = NULL;
    double *filtered = NULL;

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
    values = calloc(image->height, sizeof(double));
    filtered = calloc(image->height, sizeof(double));

    for(y = 0; y < image->height; y++){
        double m1 = get_median(image->double_data + y*image->width, image->width);

        values[y] = m1;
    }

    for(y = 0; y < image->height; y++){
        int i;

        for(i = -N; i <= N; i++)
            if(y + i >= 0 && y + i < image->height)
                filtered[y + i] += values[y]*h[N + i];
    }

    for(y = 0; y < image->height; y++){
        for(x = 0; x < image->width; x++)
            PIXEL_DOUBLE(image, x, y) -= (values[y] - filtered[y]) - m;
    }

    free(values);
    free(filtered);

#endif
}

image_str *get_calib_image(image_str *image, char *name)
{
    int shutter = image_keyword_get_int(image, "SHUTTER");
    int channel_id = image_keyword_get_int(image, "CHANNEL ID");

    char *filename = make_string("calibrations/shutter_%d_channel_%d_%s.fits", shutter, channel_id, name);

    image_str *result = image_create_from_fits(filename);

    if(!result){
        dprintf("Can't open file %s\n", filename);
        exit(-1);
    } else {
        dprintf("%s read: mean %g min %g max %g\n", filename, image_mean(result), image_min_value(result), image_max_value(result));
    }

    free(filename);

    if(result->type != IMAGE_DOUBLE){
        image_str *new = image_convert_to_double(result);

        image_delete(result);

        result = new;
    }

    return result;
}

void image_linearize(image_str *image, image_str *omask)
{
    int shutter = image_keyword_get_int(image, "SHUTTER");

    image_str *bias = get_calib_image(image, "bias");
    image_str *thresh = get_calib_image(image, "thresh");
    image_str *mask = get_calib_image(image, "mask");

    image_str *scale1_0 = get_calib_image(image, "scale1_0");
    image_str *scale1_1 = get_calib_image(image, "scale1_1");
    image_str *scale1_2 = get_calib_image(image, "scale1_2");
    image_str *scale1_3 = get_calib_image(image, "scale1_3");
    image_str *scale1_4 = get_calib_image(image, "scale1_4");

    image_str *scale2_0 = get_calib_image(image, "scale2_0");
    image_str *scale2_1 = get_calib_image(image, "scale2_1");
    image_str *scale2_2 = get_calib_image(image, "scale2_2");
    image_str *scale2_3 = get_calib_image(image, "scale2_3");
    image_str *scale2_4 = get_calib_image(image, "scale2_4");

    int saturation = image_keyword_get_int(bias, "SATURATE");
    int bias0 = image_keyword_get_int(bias, "BIAS0");

    int d;

    if(!image_keyword_find(bias, "SATURATE") || !image_keyword_find(bias, "BIAS0")){
        dprintf("Incorrect calibration frames for shutter %d and channel %d!\n", shutter, image_keyword_get_int(image, "CHANNEL ID"));
        return;
    }

    for(d = 0; d < image->width*image->height; d++){
        double v = image->double_data[d];
        double lv = log10(v - bias0);

        if(omask && omask->type == IMAGE_UINT16)
            omask->data[d] = mask->double_data[d];
        else if(omask && omask->type == IMAGE_DOUBLE)
            omask->double_data[d] = mask->double_data[d];

        if((mask->double_data[d] > 0 && IS_MASK_NANS) ||
           (mask->double_data[d] > 1 && !IS_MASK_NANS) ||
           v >= saturation){
            image->double_data[d] = saturation;
            continue;
        }

        if(v < thresh->double_data[d])
            v /= scale1_0->double_data[d]*lv*lv*lv*lv + scale1_1->double_data[d]*lv*lv*lv + scale1_2->double_data[d]*lv*lv + scale1_3->double_data[d]*lv + scale1_4->double_data[d];
        else
            v /= scale2_0->double_data[d]*lv*lv*lv*lv + scale2_1->double_data[d]*lv*lv*lv + scale2_2->double_data[d]*lv*lv + scale2_3->double_data[d]*lv + scale2_4->double_data[d];

        v -= bias->double_data[d];

        if(v > saturation || !isfinite(v))
            v = saturation;

        image->double_data[d] = v;
    }

    image_keyword_add_double(image, "SATURATE", saturation, "Saturation level");

    image_delete(bias);
    image_delete(thresh);
    image_delete(mask);

    image_delete(scale1_0);
    image_delete(scale1_1);
    image_delete(scale1_2);
    image_delete(scale1_3);
    image_delete(scale1_4);

    image_delete(scale2_0);
    image_delete(scale2_1);
    image_delete(scale2_2);
    image_delete(scale2_3);
    image_delete(scale2_4);
}
