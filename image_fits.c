#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

#include <fitsio2.h>

#include "utils.h"
#include "time_str.h"
#include "image.h"

static int include_length = 1;
static char *include_list[] = {"*"};

static int exclude_length = 33;
static char *exclude_list[] = {"SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2",
                               "EXTEND", "COMMENT", "BZERO", "BSCALE",
                               "ZIMAGE", "ZBITPIX", "ZNAXIS", "ZNAXIS1", "ZNAXIS2",
                               "ZTILE1", "ZTILE2", "ZQUANTIZ", "ZDITHER0", "ZCMPTYPE",
                               "ZNAME1", "ZVAL1", "ZNAME2", "ZVAL2", "XTENSION", /* 24 */
                               "PCOUNT", "GCOUNT", "TFIELDS", "TTYPE1", "TFORM1",
                               "TTYPE2", "TFORM2", "TTYPE3", "TFORM3"};

/* Mutex for library protection - will use if it is not reentrant */
static pthread_mutex_t cfitsio_mutex = PTHREAD_MUTEX_INITIALIZER;

void cfitsio_lock()
{
    if(!fits_is_reentrant())
        pthread_mutex_lock(&cfitsio_mutex);
}

void cfitsio_unlock()
{
    if(!fits_is_reentrant())
        pthread_mutex_unlock(&cfitsio_mutex);
}

static void image_keywords_from_fits(image_str *image, fitsfile *fits)
{
    char card[FLEN_CARD];
    int status = 0;

    fits_read_record(fits, 0, card, &status);

    while(TRUE){
        char key[FLEN_KEYWORD];
        char value[FLEN_VALUE];
        char comment[FLEN_COMMENT];
        int length = 0;
        image_keyword_str *kw = NULL;

        status = 0;
        fits_find_nextkey(fits, include_list, include_length, exclude_list, exclude_length, card, &status);

        if(status == KEY_NO_EXIST)
            break;

        fits_get_keyname(card, key, &length, &status);
        fits_parse_value(card, value, comment, &status);

        if(key[0] != '_'){
            kw = image_keyword_add(image, key, value, comment);

            fits_get_keytype(value, &kw->type, &status);
        }
    }
}

static void image_keywords_to_fits(image_str *image, fitsfile *fits)
{
    int status = 0;
    int d;

    for(d = 0; d < image->Nkeywords; d++){
        image_keyword_str *kw = &image->keywords[d];

        switch(kw->type){
        case 'I':
            {
                int value = 0;

                sscanf(kw->unquoted, "%d", &value);
                fits_update_key(fits, TINT, kw->key, &value, kw->comment, &status);
            break;
            }
        case 'F':
            {
                double value = 0;

                sscanf(kw->unquoted, "%lf", &value);
                fits_update_key(fits, TDOUBLE, kw->key, &value, kw->comment, &status);
            break;
            }
        default:
            status = 0;
            fits_update_key(fits, TSTRING, kw->key, kw->unquoted, kw->comment, &status);
            break;
        }
    }
}

image_str *image_create_from_fits(char *filename)
{
    fitsfile *fits;
    char buf[FLEN_CARD];
    image_str *image = NULL;

    int status = 0;  /* Error code for CFITSIO library */
    int Ndims = 0;
    int bitpix = 0;
    long Npixels = 0;
    long firstpix[2] = {1, 1};
    long fits_dims[2] = {1, 1};

    if(!filename || !file_exists_and_normal(filename))
        return NULL;

    fits_open_file(&fits, filename, READONLY, &status);

    do {
        /* read dimensions */
        fits_get_img_dim(fits, &Ndims, &status);
        fits_get_img_size(fits, 2, fits_dims, &status);
        fits_get_img_type(fits, &bitpix, &status);

        if(Ndims < 2){
            int Nhdus = 0;
            int hdu = 0;
            int type = 0;

            fits_get_num_hdus(fits, &Nhdus, &status);
            fits_get_hdu_num(fits, &hdu);

            if(hdu < Nhdus){
                /* dprintf("Moving to HDU %d of %d\n", hdu + 1, Nhdus); */
                fits_movabs_hdu(fits, hdu + 1, &type, &status);
            } else {
                /* dprintf("Error: FITS images with less than 2 dimensions are not supported\n"); */
                /* fits_close_file(fits, &status); */
                /* return NULL; */
                break;
            }
        }
    } while(Ndims < 2);

    if(status || Ndims > 2){
        if(Ndims > 2)
            dprintf("Error: FITS images with more than 2 dimensions are not supported\n");
        fits_close_file(fits, &status);

        return image;
    }

    Npixels = fits_dims[0]*fits_dims[1];

    if(bitpix == DOUBLE_IMG || bitpix == FLOAT_IMG){
        double *data = (double *) malloc(Npixels * sizeof(double));

        fits_read_pix(fits, TDOUBLE, firstpix, Npixels, NULL, data, NULL, &status);
        image = image_create_with_data(fits_dims[0], fits_dims[1], (u_int16_t *)data);
        image->type = IMAGE_DOUBLE;
    } else {
        u_int16_t *data = (u_int16_t *) malloc(Npixels * sizeof(u_int16_t));

        fits_read_pix(fits, TUSHORT, firstpix, Npixels, NULL, data, NULL, &status);
        image = image_create_with_data(fits_dims[0], fits_dims[1], data);
    }

    /* Read all keywords we may need */
    image_keywords_from_fits(image, fits);
    /* Rewind the header */
    fits_read_record(fits, 0, buf, &status);

    image->coords = image_keyword_get_coords(image);
    if(image->coords.ra0 == 0 && image->coords.dec0 == 0){
        image->coords.ra0 = 15.0*image_keyword_get_sexagesimal(image, "RA");
        image->coords.dec0 = image_keyword_get_sexagesimal(image, "DEC");
    }
    image->time = time_str_from_date_time(image_keyword_get_string(image, "TIME"));

    fits_close_file(fits, &status);

    return image;
}

static char *make_ra_string(double ra)
{
    double ra_hr = ra/15;
    int ra1 = (int)ra_hr;
    int ra2 = (int)(60*(ra_hr - ra1));
    double ra3 = 3600*(ra_hr - ra1) - 60*ra2;

    return make_string("%02d:%02d:%05.2f", ra1, ra2, ra3);
}

static char *make_dec_string(double dec)
{
    int dec1 = (int)fabs(dec);
    int dec2 = (int)((fabs(dec) - dec1)*60);
    double dec3 = ((fabs(dec) - dec1)*60 - dec2)*60;

    return make_string("%03d:%02d:%05.2f", dec1*(int)sign(dec), dec2, dec3);
}

void image_dump_to_fits(image_str *image, char *filename)
{
    fitsfile *fits;
    int status;
    long naxis = 2;
    long naxes[2] = {image->width, image->height};
    char *filename_exclam = make_string("!%s", filename);

    if(!filename)
        return;

    cfitsio_lock();

    status = 0;
    fits_create_file(&fits, filename_exclam, &status);
    free(filename_exclam);

    /* Creating 3-dimensional FITS file of necessary type. */
    if(image->type == IMAGE_DOUBLE){
        fits_create_img(fits, DOUBLE_IMG, naxis, naxes, &status);
        fits_write_2d_dbl(fits, 0, image->width, image->width, image->height,
                          image->double_data, &status);
    } else {
        fits_create_img(fits, USHORT_IMG, naxis, naxes, &status);
        fits_write_2d_usht(fits, 0, image->width, image->width, image->height,
                           image->data, &status);
    }

    if(!coords_is_empty(&image->coords))
        image_keyword_add_coords(image, image->coords);

    /* Write keywords */
    image_keywords_to_fits(image, fits);

    /* Write header */
    {
        double ra = image->coords.ra0;
        double dec = image->coords.dec0;
        char *ra_string = NULL;
        char *dec_string = NULL;
        char *time_string = time_str_get_date_time(image->time);
        char *date_obs_string = make_string("%04d-%02d-%02d", image->time.year, image->time.month, image->time.day);
        char *time_obs_string = make_string("%02d:%02d:%04.2lf", image->time.hour, image->time.minute, image->time.second + 1e-6*image->time.microsecond);
        double jd = time_str_get_JD(image->time);

        if(!coords_is_empty(&image->coords))
            coords_get_ra_dec(&image->coords, 0.5*image->width, 0.5*image->height, &ra, &dec);

        ra_string = make_ra_string(ra);
        dec_string = make_dec_string(dec);

        fits_update_key(fits, TSTRING, "RA", ra_string, "Rough center coordinates", &status);
        fits_update_key(fits, TSTRING, "DEC", dec_string, "Rough center coordinates", &status);
        fits_update_key(fits, TSTRING, "TIME", time_string, "Date and time", &status);
        fits_update_key(fits, TSTRING, "DATE-OBS", date_obs_string, NULL, &status);
        fits_update_key(fits, TSTRING, "TIME-OBS", time_obs_string, NULL, &status);

        fits_update_key(fits, TDOUBLE, "JD", &jd, "Julianic date, millisecond precision", &status);

        free(ra_string);
        free(dec_string);
        free(time_string);
        free(date_obs_string);
        free(time_obs_string);
    }

    fits_close_file(fits, &status);

    /* Print out error messages if any */
    fits_report_error(stderr, status);

    cfitsio_unlock();
}

void image_fits_get_info(char *filename, double *ra_ptr, double *dec_ptr, double *exp_ptr, u_int64_t *number_ptr, time_str *time_ptr)
{
    fitsfile *fits;
    image_str *image = image_create(0, 0);
    int status = 0;  /* Error code for CFITSIO library */

    fits_open_file(&fits, filename, READONLY, &status);

    /* Read all keywords we may need */
    image_keywords_from_fits(image, fits);

    if(exp_ptr)
        *exp_ptr = image_keyword_get_double(image, "EXPOSURE");
    if(ra_ptr)
        *ra_ptr = 15.0*image_keyword_get_sexagesimal(image, "RA");
    if(dec_ptr)
        *dec_ptr = image_keyword_get_sexagesimal(image, "DEC");
    if(time_ptr)
        *time_ptr = time_str_from_date_time(image_keyword_get_string(image, "TIME"));
    if(number_ptr)
        *number_ptr = image_keyword_get_int64(image, "FRAMENUMBER");

    fits_close_file(fits, &status);
    image_delete(image);
}

double image_fits_get_ra(char *filename)
{
    double ra = 0;

    image_fits_get_info(filename, &ra, NULL, NULL, NULL, NULL);

    return ra;
}

double image_fits_get_dec(char *filename)
{
    double dec = 0;

    image_fits_get_info(filename, NULL, &dec, NULL, NULL, NULL);

    return dec;
}
