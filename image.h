/* Image operations */
#ifndef IMAGE_H
#define IMAGE_H

#include <fitsio.h>

#include "time_str.h"
#include "coords.h"

typedef struct {
    char key[FLEN_KEYWORD];
    char value[FLEN_VALUE];
    char unquoted[FLEN_VALUE];
    char comment[FLEN_COMMENT];
    char type; /* FITSIO type */
} image_keyword_str;

typedef struct {
    /* Geometric size */
    int width;
    int height;

    /* Exposition start time */
    time_str time;

    /* Coordinate transformation matrix */
    coords_str coords;

    /* Optional keywords */
    int Nkeywords;
    image_keyword_str *keywords;

    /* Type: 0 = u_int16_t, 1 = double */
    enum image_type {
        IMAGE_UINT16,
        IMAGE_DOUBLE
    } type;

    /* Data */
    union {
        u_int16_t *data; /* 2-byte type */
        u_int16_t *uint16_data; /* 2-byte type */
        double *double_data; /* 8-byte DOUBLE type */
    };
} image_str;

#define PIXEL(image, x, y) ((image)->data[(int)(y)*(image)->width + (int)(x)])
#define PIXEL_UINT16(image, x, y) ((image)->data[(int)(y)*(image)->width + (int)(x)])
#define PIXEL_DOUBLE(image, x, y) ((image)->double_data[(int)(y)*(image)->width + (int)(x)])

image_str *image_create(int , int );
image_str *image_create_double(int , int );
image_str *image_create_with_type(int , int , enum image_type );
image_str *image_create_with_data(int , int , u_int16_t *);
void image_delete(image_str *);

#define image_delete_and_null(image) do {\
        if(image) image_delete(image); (image) = NULL;  \
    } while(0)

void image_fill(image_str *, double );
void image_fill_nans(image_str *, double );
void image_clean(image_str *);

void image_copy_properties(image_str *, image_str *);

image_str *image_copy(image_str *);
image_str *image_crop(image_str *, int , int , int , int );
image_str *image_scale(image_str *, int );
image_str *image_downscale(image_str *, int );

image_str *image_convert_to_double(image_str *);

/* Keyword stuff */
image_keyword_str *image_keyword_add(image_str *, char *, char *, char *);
void image_keyword_add_int(image_str *, char *, int , char *);
void image_keyword_add_int64(image_str *, char *, u_int64_t , char *);
void image_keyword_add_double(image_str *, char *, double , char *);
void image_keyword_add_time(image_str *, char *, time_str , char *);
image_keyword_str *image_keyword_find(image_str *, char *);
char *image_keyword_get_string(image_str *, char *);
int image_keyword_get_int(image_str *, char *);
u_int64_t image_keyword_get_int64(image_str *, char *);
double image_keyword_get_double(image_str *, char *);
double image_keyword_get_sexagesimal(image_str *, char *);
time_str image_keyword_get_time(image_str *, char *);
char *image_keywords_as_hstore(image_str *);
void image_keyword_add_coords(image_str *, coords_str );
coords_str image_keyword_get_coords(image_str *);

/* Input/output routines */
image_str *image_create_from_fits(char *);
void image_dump_to_fits(image_str *, char *);
void image_dump_to_jpeg(image_str *, char *);
void image_convert_to_jpeg(image_str *, unsigned char **, int *);
void image_jpeg_set_scale(int );
void image_jpeg_set_percentile(double, double );

/* FITS-related stuff */
void image_fits_get_info(char *, double *, double *, double *, u_int64_t *, time_str *);
/* Thread-safety */
void cfitsio_lock();
void cfitsio_unlock();

/* Combine image with number and second image, result will be in first one */
enum image_combine_op {
    IMAGE_OP_ADD,
    IMAGE_OP_SUB,
    IMAGE_OP_MUL,
    IMAGE_OP_DIV,
};

void image_combine(image_str *, double , image_str *, enum image_combine_op );
void image_add(image_str *, image_str *);

/* Some image statistics */
double image_max_value(image_str *);
double image_min_value(image_str *);
double image_mean(image_str *);
double image_sigma(image_str *);
double image_sum(image_str *);
double image_median(image_str *);

/* ANDOR-specific cleanups */
void image_clean_stripes(image_str *);
void image_linearize(image_str *, image_str *);

/* Image smoothing - Gaussian */
image_str *image_smooth(image_str *, double );
image_str *image_unsharp(image_str *, double );

/* Errors image from data image using bias, gain and readout noise */
image_str *image_errors(image_str *, double , double , double );

/* Image background estimation */
image_str *image_background(image_str *, image_str *, image_str *, int );

#endif /* IMAGE_H */
