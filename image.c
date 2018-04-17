#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/mman.h>

#include "time_str.h"
#include "image.h"
#include "utils.h"

image_str *image_create_with_data(int width, int height, u_int16_t *data)
{
    image_str *image = (image_str *)malloc(sizeof(image_str));

    image->type = IMAGE_UINT16;

    image->width  = width;
    image->height = height;

    /* Dummy data, you need to replace it with your values */
    image->time = time_zero();
    image->coords = coords_empty();

    image->Nkeywords = 0;
    image->keywords = NULL;

    image->data = data;

    if(image->data)
        madvise(image->data, image->width*image->height*sizeof(u_int16_t), MADV_SEQUENTIAL);

    return image;
}

image_str *image_create_with_type(int width, int height, enum image_type type)
{
    image_str *image = image_create_with_data(width, height, NULL);

    image->type = type;

    if(image->type == IMAGE_DOUBLE){
        image->double_data = (double *)calloc(width*height, sizeof(double));
        madvise(image->double_data, image->width*image->height*sizeof(double), MADV_SEQUENTIAL);
    } else {
        image->data = (u_int16_t *)calloc(width*height, sizeof(u_int16_t));
        madvise(image->data, image->width*image->height*sizeof(u_int16_t), MADV_SEQUENTIAL);
    }

    return image;
}

image_str *image_create(int width, int height)
{
    return image_create_with_type(width, height, IMAGE_UINT16);
}

image_str *image_create_double(int width, int height)
{
    return image_create_with_type(width, height, IMAGE_DOUBLE);
}

void image_delete(image_str *image)
{
    if(!image)
        return;

    if(image->keywords)
        free(image->keywords);

    if(image->data)
        free(image->data);

    free(image);
}

void image_fill(image_str *image, double val)
{
    int d;

    if(!image)
        return;

    if(image->type == IMAGE_DOUBLE)
        for(d = 0; d < image->width*image->height; d++)
            image->double_data[d] = val;
    else
        for(d = 0; d < image->width*image->height; d++)
            image->data[d] = val;
}

void image_fill_nans(image_str *image, double val)
{
    int d;

    if(!image)
        return;

    if(image->type == IMAGE_DOUBLE)
        for(d = 0; d < image->width*image->height; d++)
            if(!isfinite(image->double_data[d]))
                image->double_data[d] = val;
    else
        for(d = 0; d < image->width*image->height; d++)
            if(!isfinite(image->data[d]))
                image->data[d] = val;
}

void image_clean(image_str *image)
{
    if(!image)
        return;

    if(image->type == IMAGE_DOUBLE)
        memset(image->double_data, 0, image->width*image->height*sizeof(double));
    else
        memset(image->data, 0, image->width*image->height*sizeof(u_int16_t));
}

void image_copy_properties(image_str *from, image_str *to)
{
    to->time = from->time;
    to->coords = from->coords;

    to->keywords = realloc(to->keywords, sizeof(image_keyword_str)*from->Nkeywords);
    to->Nkeywords = from->Nkeywords;
    memcpy(to->keywords, from->keywords, sizeof(image_keyword_str)*from->Nkeywords);
}

image_str *image_copy(image_str *image_in)
{
    image_str *image = NULL;

    if(image_in){
        image = image_create_with_type(image_in->width, image_in->height, image_in->type);

        if(image_in->type == IMAGE_DOUBLE)
            memcpy(image->double_data, image_in->double_data, sizeof(double)*image_in->width*image_in->height);
        else
            memcpy(image->data, image_in->data, sizeof(u_int16_t)*image_in->width*image_in->height);

        image->time = image_in->time;
        image_copy_properties(image_in, image);
    }

    return image;
}

image_str *image_crop(image_str *orig, int x1, int y1, int x2, int y2)
{
    image_str *image = image_create_with_type(x2 - x1, y2 - y1, orig->type);
    int x;
    int y;

    if(orig->type == IMAGE_DOUBLE)
        for(x = x1; x < x2; x++)
            for(y = y1; y < y2; y++)
                PIXEL_DOUBLE(image, x - x1, y - y1) =
                    (x >= 0 && x < orig->width &&
                     y >= 0 && y < orig->height)
                    ? PIXEL_DOUBLE(orig, x, y)
                    : 0;
    else
        for(x = x1; x < x2; x++)
            for(y = y1; y < y2; y++)
                PIXEL(image, x - x1, y - y1) =
                    (x >= 0 && x < orig->width &&
                     y >= 0 && y < orig->height)
                    ? PIXEL(orig, x, y)
                    : 0;

    image_copy_properties(orig, image);

    if(!coords_is_empty(&orig->coords)){
        /* Adjust coordinate transformation */
        image->coords.CRPIX1 -= x1;
        image->coords.CRPIX2 -= y1;
    }

    return image;
}

image_str *image_scale(image_str *orig, int factor)
{
    image_str *image = image_create_with_type(orig->width*factor, orig->height*factor, orig->type);
    int x;
    int y;

    if(orig->type == IMAGE_DOUBLE)
        for(x = 0; x < image->width; x++)
            for(y = 0; y < image->height; y++){
                int xx = x/factor;
                int yy = y/factor;

                PIXEL_DOUBLE(image, x, y) += PIXEL_DOUBLE(orig, xx, yy);
            }
    else
        for(x = 0; x < image->width; x++)
            for(y = 0; y < image->height; y++){
                int xx = x/factor;
                int yy = y/factor;

                PIXEL(image, x, y) += PIXEL(orig, xx, yy);
            }

/*     for(x = 0; x < image->width; x++) */
/*         for(y = 0; y < image->height; y++) */
/*             PIXEL(image, x, y) /= factor*factor; */

    image_copy_properties(orig, image);

    return image;
}

image_str *image_downscale(image_str *orig, int factor)
{
    image_str *image = image_create_with_type(orig->width/factor, orig->height/factor, orig->type);
    int x;
    int y;

    if(orig->type == IMAGE_DOUBLE)
        for(x = 0; x < orig->width; x++)
            for(y = 0; y < orig->height; y++){
                int xx = x/factor;
                int yy = y/factor;

                PIXEL_DOUBLE(image, xx, yy) += PIXEL_DOUBLE(orig, x, y)/factor/factor;
            }
    else
        for(x = 0; x < orig->width; x++)
            for(y = 0; y < orig->height; y++){
                int xx = x/factor;
                int yy = y/factor;

                PIXEL(image, xx, yy) += PIXEL(orig, x, y)/factor/factor;
            }

    image_copy_properties(orig, image);

    return image;
}

image_str *image_convert_to_double(image_str *orig)
{
    image_str *image = image_create_double(orig->width, orig->height);

    if(orig->type == IMAGE_DOUBLE)
        memcpy(image->double_data, orig->double_data, sizeof(double)*orig->width*orig->height);
    else {
        int d;

        for(d = 0; d < orig->width*orig->height; d++){
            image->double_data[d] = orig->data[d];
            if(orig->data[d] == -0x80000000)
                image->double_data[d] = NAN;
        }
    }

    image_copy_properties(orig, image);

    return image;
}

void image_combine(image_str *image, double val, image_str *second, enum image_combine_op op)
{
    int d;

    if(second &&
       ((image->width != second->width) ||
        (image->height != second->height)))
        return;

    for(d = 0; d < image->width*image->height; d++){
        double op2 = val * (second ? (second->type == IMAGE_DOUBLE ? second->double_data[d] : second->data[d]) : 1);
        double op1 = image->type == IMAGE_DOUBLE ? image->double_data[d] : image->data[d];

        switch(op){
        case IMAGE_OP_ADD:
            op1 += op2;
            break;
        case IMAGE_OP_SUB:
            op1 -= op2;
            break;
        case IMAGE_OP_MUL:
            op1 *= op2;
            break;
        case IMAGE_OP_DIV:
            if(op2)
                op1 /= op2;
            else
                op1 = 0;
            break;
        }

        if(image->type == IMAGE_DOUBLE)
            image->double_data[d] = op1;
        else
            image->data[d] = MAX(0, MIN(op1, 0xffff));
    }
}

/* Stripped-down faster version for co-adding two images */
void image_add(image_str *image, image_str *second)
{
    int d;

#define WORKER(type1, type2) \
    for(d = 0; d < image->width*image->height; d++) \
        image->type1##_data[d] += second->type2##_data[d]

    if(image->type == IMAGE_DOUBLE && second->type == IMAGE_DOUBLE)
        WORKER(double, double);
    else if(image->type == IMAGE_DOUBLE && second->type == IMAGE_UINT16)
        WORKER(double, uint16);
    else if(image->type == IMAGE_UINT16 && second->type == IMAGE_DOUBLE)
        WORKER(uint16, double);
    else if(image->type == IMAGE_UINT16 && second->type == IMAGE_UINT16)
        WORKER(uint16, uint16);
#undef WORKER
}

/* FIXME: add support for IMAGE_DOUBLE type */
double image_max_value(image_str *image)
{
    int d;
    double value = 0;

    if(image->type == IMAGE_DOUBLE){
        value = image->double_data[0];

        for(d = 0; d < image->width*image->height; d++)
            if(isfinite(image->double_data[d]))
                value = MAX(value, image->double_data[d]);
    } else {
        value = image->data[0];

        for(d = 0; d < image->width*image->height; d++)
            value = MAX(value, image->data[d]);
    }

    return value;
}

/* FIXME: add support for IMAGE_DOUBLE type */
double image_min_value(image_str *image)
{
    int d;
    double value = 0;

    if(image->type == IMAGE_DOUBLE){
        value = image->double_data[0];

        for(d = 0; d < image->width*image->height; d++)
            value = MIN(value, image->double_data[d]);
    } else {
        value = image->data[0];

        for(d = 0; d < image->width*image->height; d++)
            value = MIN(value, image->data[d]);
    }

    return value;
}

double image_mean(image_str *image)
{
    int d;
    double sum = 0;
    int N = 0;

    if(image->type == IMAGE_DOUBLE)
        for(d = 0; d < image->width*image->height; d++)
            if(isfinite(image->double_data[d])){
                sum += image->double_data[d];
                N ++;
            }
    else
        for(d = 0; d < image->width*image->height; d++)
            if(isfinite(image->data[d])){
                sum += image->data[d];
                N ++;
            }

    return sum*1./N;
}

double image_sigma(image_str *image)
{
    int d;
    double sum = 0;
    double sum2 = 0;
    int N = image->width*image->height;

    if(image->type == IMAGE_DOUBLE)
        for(d = 0; d < image->width*image->height; d++){
            sum += image->double_data[d];
            sum2 += image->double_data[d]*image->double_data[d];
        }
    else
        for(d = 0; d < image->width*image->height; d++){
            sum += image->data[d];
            sum2 += image->data[d]*image->data[d];
        }

    return sqrt((sum2 - sum*sum*1./N)*1./(N - 1));
}

double image_sum(image_str *image)
{
    int d;
    double sum = 0;

    if(image->type == IMAGE_DOUBLE)
        for(d = 0; d < image->width*image->height; d++)
            sum += image->double_data[d];
    else
        for(d = 0; d < image->width*image->height; d++)
            sum += image->data[d];

    return sum;
}

double image_median(image_str *image)
{
    int N = image->width*image->height;
    int *idx = (int *)malloc(sizeof(int)*N);
    int d;
    double value = 0;

    for(d = 0; d < N; d++)
        idx[d] = d;

    if(image->type == IMAGE_DOUBLE){
        value = get_median(image->double_data, image->width*image->height);
    } else {
        int compare_fn(const void *v1, const void *v2)
        {
            int res = image->data[*(int *)v1] - image->data[*(int *)v2];

            return (res > 0 ? 1 : (res < 0 ? -1 : 0));
        }

        qsort(idx, N, sizeof(int), compare_fn);
        value = image->data[idx[(int)floor(N/2)]];

    }

    free(idx);

    return value;
}

image_str *image_smooth(image_str *image, double sigma)
{
    image_str *res = image_create_double(image->width, image->height);
    double saturation = 65535;
    int size = ceil(sigma*4);
    int x;
    int y;

    image_str *kernel = image_create_double(2*size+1, 2*size+1);

    /* Guess saturation level from image keywords */
    if(image_keyword_find(image, "SATURATE"))
        saturation = image_keyword_get_double(image, "SATURATE");

    for(y = -size; y < size+1; y++)
        for(x = -size; x < size+1; x++){
            double arg = (x*x + y*y)/2/sigma/sigma;
            double value = exp(-arg)/2/M_PI/sigma/sigma;

            PIXEL_DOUBLE(kernel, x + size, y + size) = value;
        }

#pragma omp parallel for private(x)
    for(y = 0; y < image->height; y++)
        for(x = 0; x < image->width; x++){
            int x1 = 0;
            int y1 = 0;
            double value = 0;
            double kernel_sum = 0;

            for(x1 = MAX(0, x-size); x1 < MIN(image->width, x+size+1); x1++)
                for(y1 = MAX(0, y-size); y1 < MIN(image->height, y+size+1); y1++){
                    double weight = PIXEL_DOUBLE(kernel, x1 - x + size, y1 - y + size);

                    if(image->type == IMAGE_DOUBLE){
                        if(PIXEL_DOUBLE(image, x1, y1) < saturation){
                            value += PIXEL_DOUBLE(image, x1, y1)*weight;
                            kernel_sum += weight;
                        }
                    } else {
                        if(PIXEL(image, x1, y1) < saturation){
                            value += PIXEL(image, x1, y1)*weight;
                            kernel_sum += weight;
                        }
                    }
                }

            PIXEL_DOUBLE(res, x, y) = value/kernel_sum;
        }

    image_copy_properties(image, res);

    image_delete(kernel);

    return res;
}

image_str *image_unsharp(image_str *image, double sigma)
{
    image_str *smooth = image_smooth(image, sigma);
    int d;

    if(image->type == IMAGE_DOUBLE)
        for(d = 0; d < image->width*image->height; d++)
            smooth->double_data[d] = image->double_data[d] - smooth->double_data[d];
    else
        for(d = 0; d < image->width*image->height; d++)
            smooth->double_data[d] = image->data[d] - smooth->double_data[d];

    return smooth;
}

image_str *image_errors(image_str *image, double bias, double gain, double readnoise)
{
    image_str *errors = image_create_double(image->width, image->height);
    int d;

    if(image->type == IMAGE_DOUBLE)
        for(d = 0; d < image->width*image->height; d++){
            double value = (image->double_data[d] - bias)/gain;

            errors->double_data[d] = hypot(readnoise/gain, sqrt(MAX(0, value)));
        }
    else
        for(d = 0; d < image->width*image->height; d++){
            double value = (image->data[d] - bias)/gain;

            errors->double_data[d] = hypot(readnoise/gain, sqrt(MAX(0, value)));
        }

    image_copy_properties(image, errors);

    return errors;
}
