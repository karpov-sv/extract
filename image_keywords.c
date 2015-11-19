#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include <fitsio2.h>

#include "utils.h"
#include "image.h"

image_keyword_str *image_keyword_add(image_str *image, char *key, char *value, char *comment)
{
    image_keyword_str *kw = image_keyword_find(image, key);
    int status = 0;

    if(!kw){
        image->keywords = realloc(image->keywords, sizeof(image_keyword_str)*(image->Nkeywords + 1));
        kw = &image->keywords[image->Nkeywords];
        image->Nkeywords ++;

        if(key)
            strncpy(kw->key, key, FLEN_KEYWORD);
        else
            *kw->key = '\0';
    }

    if(value)
        strncpy(kw->value, value, FLEN_VALUE);
    else
        *kw->value = '\0';

    kw->value[FLEN_VALUE - 1] = '\0';

    ffc2s(kw->value, kw->unquoted, &status);
    if(comment)
        strncpy(kw->comment, comment, FLEN_COMMENT);
    else
        *kw->comment = '\0';

    kw->type = 'C'; /* String */

    return kw;
}

void image_keyword_add_int(image_str *image, char *key, int value, char *comment)
{
    char *string = make_string("%d", value);
    image_keyword_str *kw = image_keyword_add(image, key, string, comment);

    kw->type = 'I'; /* Int */

    free(string);
}

void image_keyword_add_int64(image_str *image, char *key, u_int64_t value, char *comment)
{
    char *string = make_string("%lld", value);
    image_keyword_str *kw = image_keyword_add(image, key, string, comment);

    kw->type = 'I'; /* Int */

    free(string);
}

void image_keyword_add_double(image_str *image, char *key, double value, char *comment)
{
    char *string = make_string("%.12g", value);
    image_keyword_str *kw = image_keyword_add(image, key, string, comment);

    kw->type = 'F'; /* Float */

    free(string);
}

void image_keyword_add_time(image_str *image, char *key, time_str value, char *comment)
{
    char *string = time_str_get_date_time(value);

    image_keyword_add(image, key, string, comment);

    free(string);
}

image_keyword_str *image_keyword_find(image_str *image, char *key)
{
    image_keyword_str *kw = NULL;
    int d;

    for(d = 0; d < image->Nkeywords; d++)
        if(!strcmp(image->keywords[d].key, key)){
            kw = &image->keywords[d];
            break;
        }

    return kw;
}

char *image_keyword_get_string(image_str *image, char *key)
{
    char *value = NULL;
    image_keyword_str *kw = image_keyword_find(image, key);

    if(kw)
        value = kw->unquoted;

    return value;
}

int image_keyword_get_int(image_str *image, char *key)
{
    image_keyword_str *kw = image_keyword_find(image, key);
    int value = 0;

    /* TODO: check keyword type?.. */
    if(kw && kw->value)
        sscanf(kw->unquoted, "%d", &value);

    return value;
}

u_int64_t image_keyword_get_int64(image_str *image, char *key)
{
    image_keyword_str *kw = image_keyword_find(image, key);
    long long int value = 0;

    if(kw && kw->value)
        sscanf(kw->unquoted, "%lld", &value);

    return value;
}

double image_keyword_get_double(image_str *image, char *key)
{
    image_keyword_str *kw = image_keyword_find(image, key);
    double value = 0;

    if(kw && kw->value)
        sscanf(kw->unquoted, "%lf", &value);

    return value;
}

double image_keyword_get_sexagesimal(image_str *image, char *key)
{
    image_keyword_str *kw = image_keyword_find(image, key);
    double t1 = 0;
    double t2 = 0;
    double t3 = 0;

    if(kw && kw->value)
        sscanf(kw->unquoted, "%lf %lf %lf", &t1, &t2, &t3);

    return t1 + sign(t1)*(t2*1./60 + t3*1./3600);
}

char *image_keywords_as_hstore(image_str *image)
{
    char *hstring = NULL;
    char *result = NULL;
    int d;

    /* FIXME: move to somewhere else */
    if(!coords_is_empty(&image->coords))
        image_keyword_add_coords(image, image->coords);

    for(d = 0; d < image->Nkeywords; d++){
        if(hstring)
            add_to_string(&hstring, ", ");
        add_to_string(&hstring, "\"%s\" => \"%s\"", image->keywords[d].key, image->keywords[d].unquoted);
    }

    result = make_string("'%s'::hstore", hstring);

    free(hstring);

    return result;
}

void image_keyword_add_coords(image_str *image, coords_str coords)
{
    image_keyword_add(image, "CTYPE1", "'RA---TAN'", NULL);
    image_keyword_add_double(image, "CRPIX1", image->coords.CRPIX1 + 1, NULL);
    image_keyword_add_double(image, "CD1_1", image->coords.CD11, NULL);
    image_keyword_add_double(image, "CD1_2", image->coords.CD12, NULL);
    image_keyword_add_double(image, "CRVAL1", image->coords.ra0, NULL);
    image_keyword_add(image, "CUNIT1", "'deg'", NULL);

    image_keyword_add(image, "CTYPE2", "'DEC--TAN'", NULL);
    image_keyword_add_double(image, "CRPIX2", image->coords.CRPIX2 + 1, NULL);
    image_keyword_add_double(image, "CD2_1", image->coords.CD21, NULL);
    image_keyword_add_double(image, "CD2_2", image->coords.CD22, NULL);
    image_keyword_add_double(image, "CRVAL2", image->coords.dec0, NULL);
    image_keyword_add(image, "CUNIT1", "'deg'", NULL);

    image_keyword_add_double(image, "LONPOLE", 180.0, "Native longitude of celestial pole");
    image_keyword_add_double(image, "LATPOLE", -90.0, "Native latitude of celestial pole");
    image_keyword_add_double(image, "EQUINOX", 2000.0, "Equinox of equatorial coordinates");

    /* Write SIP polynomial coefficients, if any */
    if(coords.A_ORDER){
        int i;
        int j;

        image_keyword_add(image, "CTYPE1", "'RA---TAN-SIP'", NULL);
        image_keyword_add(image, "CTYPE2", "'DEC--TAN-SIP'", NULL);

        image_keyword_add_int(image, "A_ORDER", image->coords.A_ORDER, "Polynomial order, axis 1");
        image_keyword_add_int(image, "B_ORDER", image->coords.B_ORDER, "Polynomial order, axis 2");
        image_keyword_add_int(image, "AP_ORDER", image->coords.AP_ORDER, "Inv polynomial order, axis 1");
        image_keyword_add_int(image, "BP_ORDER", image->coords.BP_ORDER, "Inv polynomial order, axis 2");

        for(i = 0; i <= MAX_SIP_ORDER; i++)
            for(j = 0; j <= MAX_SIP_ORDER; j++){
                if(i + j <= coords.A_ORDER ){
                    char *name = make_string("A_%d_%d", i, j);
                    image_keyword_add_double(image, name, coords.A[i + j*MAX_SIP_ORDER], NULL);
                    free(name);
                }

                if(i + j <= coords.AP_ORDER ){
                    char *name = make_string("AP_%d_%d", i, j);
                    image_keyword_add_double(image, name, coords.AP[i + j*MAX_SIP_ORDER], NULL);
                    free(name);
                }

                if(i + j <= coords.B_ORDER ){
                    char *name = make_string("B_%d_%d", i, j);
                    image_keyword_add_double(image, name, coords.B[i + j*MAX_SIP_ORDER], NULL);
                    free(name);
                }

                if(i + j <= coords.BP_ORDER ){
                    char *name = make_string("BP_%d_%d", i, j);
                    image_keyword_add_double(image, name, coords.BP[i + j*MAX_SIP_ORDER], NULL);
                    free(name);
                }
            }
    }

    image_keyword_add_double(image, "PIXSCALE", coords_get_pixscale(&image->coords), "Estimated pixel scale, degrees");

    image_keyword_add_int(image, "MAG_FILTER", image->coords.filter, "Filter (0 - Clear, 1 - B, 2 - V, 3 - R)");
    /* image_keyword_add_double(image, "MAG_0", image->coords.mag0, "Cat = MAG_0 + MAG_SCALE*Instr"); */
    /* image_keyword_add_double(image, "MAG_0_ERR", image->coords.mag0_err, "MAG_0 stderr"); */
    /* image_keyword_add_double(image, "MAG_SCALE", image->coords.mag_scale, "Cat = MAG_0 + MAG_SCALE*Instr"); */
    /* image_keyword_add_double(image, "MAG_SCALE_ERR", image->coords.mag_scale_err, "MAG_SCALE stderr"); */
    /* image_keyword_add_double(image, "MAG_COVAR", image->coords.mag_covar, "Covariance of MAG_SCALE anf MAG_0"); */
    image_keyword_add_double(image, "MAG_SIGMA", image->coords.mag_sigma, "Standard deviation of Cat - Instr");
}

coords_str image_keyword_get_coords(image_str *image)
{
    coords_str coords = coords_empty();

    coords.CD11 = image_keyword_get_double(image, "CD1_1");
    coords.CD12 = image_keyword_get_double(image, "CD1_2");
    coords.ra0 = image_keyword_get_double(image, "CRVAL1");

    coords.CD21 = image_keyword_get_double(image, "CD2_1");
    coords.CD22 = image_keyword_get_double(image, "CD2_2");
    coords.dec0 = image_keyword_get_double(image, "CRVAL2");

    coords.CRPIX1 = image_keyword_get_double(image, "CRPIX1") - 1;
    coords.CRPIX2 = image_keyword_get_double(image, "CRPIX2") - 1;

    coords.A_ORDER = image_keyword_get_int(image, "A_ORDER");
    coords.B_ORDER = image_keyword_get_int(image, "B_ORDER");
    coords.AP_ORDER = image_keyword_get_int(image, "AP_ORDER");
    coords.BP_ORDER = image_keyword_get_int(image, "BP_ORDER");

    if(coords.A_ORDER && coords.A_ORDER <= MAX_SIP_ORDER){
        /* Read SIP polynomial coefficients */
        int i;
        int j;

        for(i = 0; i <= MAX_SIP_ORDER; i++)
            for(j = 0; j <= MAX_SIP_ORDER; j++){
                if(i + j <= coords.A_ORDER ){
                    char *name = make_string("A_%d_%d", i, j);
                    coords.A[i + j*MAX_SIP_ORDER] = image_keyword_get_double(image, name);
                    free(name);
                }

                if(i + j <= coords.AP_ORDER ){
                    char *name = make_string("AP_%d_%d", i, j);
                    coords.AP[i + j*MAX_SIP_ORDER] = image_keyword_get_double(image, name);
                    free(name);
                }

                if(i + j <= coords.B_ORDER ){
                    char *name = make_string("B_%d_%d", i, j);
                    coords.B[i + j*MAX_SIP_ORDER] = image_keyword_get_double(image, name);
                    free(name);
                }

                if(i + j <= coords.BP_ORDER ){
                    char *name = make_string("BP_%d_%d", i, j);
                    coords.BP[i + j*MAX_SIP_ORDER] = image_keyword_get_double(image, name);
                    free(name);
                }
            }
    }

    coords.filter = image_keyword_get_int(image, "MAG_FILTER");
    /* coords.mag0 = image_keyword_get_double(image, "MAG_0"); */
    /* coords.mag0_err = image_keyword_get_double(image, "MAG_0_ERR"); */
    /* coords.mag_scale = image_keyword_get_double(image, "MAG_SCALE"); */
    /* coords.mag_scale_err = image_keyword_get_double(image, "MAG_SCALE_ERR"); */
    /* coords.mag_covar = image_keyword_get_double(image, "MAG_COVAR"); */
    coords.mag_sigma = image_keyword_get_double(image, "MAG_SIGMA");

    return coords;
}
