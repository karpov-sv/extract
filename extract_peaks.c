#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "utils.h"

#include "extract.h"
#include "kdtree.h"

/* Delta arrays - possible steps from given pixel */
static int dx[] = {0, 0,  0, 1, 1,  1, -1, -1, -1};
static int dy[] = {0, 1, -1, 1, 0, -1,  1,  0, -1};
static int dN = 9;

#define IS_PIXEL_VALID(image, x, y) ((x) >= 0 && (x) < (image)->width && (y) >= 0 && (y) < (image)->height)

/* Single image pixel */
typedef struct pixel_str {
    LIST_HEAD(struct pixel_str);

    int x;
    int y;
    double flux;
    double flux_err;
    double smooth;

    int is_edge;
} pixel_str;

/* Region on an image */
typedef struct region_str {
    LIST_HEAD(struct region_str);

    int npixels;
    int npixels_bg;

    int is_edge;
    int is_noise;

    double min_flux;
    double max_flux;
    double mean_flux;

    double bg_mean;
    double bg_sigma;

    pixel_str *saddle;
    pixel_str *max;
    double prominence;

    struct list_head pixels;
} region_str;

static inline int get_max_n_idx(image_str *image, image_str *mask, int x, int y)
{
    int d;
    int max_value = 0;
    int max_id = 0;

    for(d = 0; d < dN; d++){
        if(IS_PIXEL_VALID(image, x + dx[d], y + dy[d])){
            int value = PIXEL_DOUBLE(image, x + dx[d], y + dy[d]);

            if(d == 0 || ((!mask || !PIXEL(mask, x + dx[d], y + dy[d])) && value > max_value)){
                max_value = value;
                max_id = d;
            }
        }
    }

    return max_id;
}

static inline void get_reachable_maximum_xy(image_str *image, image_str *mask, int x0, int y0, int *x_ptr, int *y_ptr, int maxsteps)
{
    int x = x0;
    int y = y0;
    int id = 0;
    int Nsteps = 0;

    while((id = get_max_n_idx(image, mask, x, y)) && (!maxsteps || Nsteps < maxsteps)){
        x += dx[id];
        y += dy[id];
        Nsteps ++;
    }

    *x_ptr = x;
    *y_ptr = y;
}

static peak_str *peak_from_region(region_str *region, image_str *mask)
{
    peak_str *peak = (peak_str *)calloc(1, sizeof(peak_str));

    int N = 0;
    double sum_x = 0;
    double sum_y = 0;
    double sum_x2 = 0;
    double sum_y2 = 0;
    double sum_xy = 0;
    double sum_weight = 0;
    double sum_flux = 0;
    double sum_flux_err = 0;
    pixel_str *pixel;

    double X2 = 0;
    double Y2 = 0;
    double XY = 0;

    double bg = region->bg_mean + region->bg_sigma;

    int flags = 0;

    foreach(pixel, region->pixels){
        double value = pixel->smooth;

        if(value >= 1e-10){
            double x = pixel->x;
            double y = pixel->y;

            sum_x += value*x;
            sum_x2 += value*x*x;
            sum_y += value*y;
            sum_y2 += value*y*y;
            sum_xy += value*x*y;
            sum_weight += value;

            sum_flux += value;
            sum_flux_err += pixel->flux_err*pixel->flux_err;

            flags |= PIXEL(mask, pixel->x, pixel->y);
            if(x == 0 || x == mask->width-1 || y == 0 || y == mask->height-1)
                flags |= FLAG_TRUNCATED;

            N ++;
        }
    }

    peak->flags = flags;

    if(sum_weight <= 0){
        peak->x = pixel->x;
        peak->y = pixel->y;
        peak->dx = 0.5;
        peak->dy = 0.5;
        peak->bg = bg;
        peak->excess = 0;

        return peak;
    }

    /* Position and flux */
    peak->x = sum_x/sum_weight;
    peak->y = sum_y/sum_weight;

    peak->flux = sum_flux; /* In principle it may be negative */
    peak->dflux = sqrt(sum_flux_err);

    peak->bg = bg;

    /* Rough estimate of peak significance */
    peak->excess = peak->flux/peak->dflux;
    /* if(region->bg_sigma > 0) */
    /*     peak->excess = peak->flux/region->bg_sigma/sqrt(region->npixels); */
    /* else */
    /*     peak->excess = 0; */

    /* if(peak->excess > 5.0) */
    /*     dprintf("%g - %g %g %d\n", peak->excess, peak->flux, region->bg_sigma, region->npixels); */

    peak->A = region->max_flux - region->min_flux; /* Conservative estimate for peak amplitude, always positive */

    /* Second-order moments */
    X2 = sum_x2/sum_weight - peak->x*peak->x;
    Y2 = sum_y2/sum_weight - peak->y*peak->y;
    XY = sum_xy/sum_weight - peak->x*peak->y;

    /* Handling singular cases */
    if(X2*Y2 - XY*XY < 1./12/12){
        X2 += 1./12;
        Y2 += 1./12;
    }

    /* Formulae from the SExtractor manual, p.29 */
    peak->theta = 0.5*atan2(2*XY, X2 - Y2);
    peak->a = sqrt(fabs(0.5*(X2 + Y2) + sqrt(0.25*(X2 - Y2)*(X2 - Y2) + XY*XY)));
    peak->b = sqrt(fabs(0.5*(X2 + Y2) - sqrt(0.25*(X2 - Y2)*(X2 - Y2) + XY*XY)));

    peak->bg = bg;

    peak->state = PEAK_INITIAL;

    return peak;
}

static inline void add_pixel_to_region(int x, int y, region_str *region, image_str *image, image_str *smooth, image_str *errors, region_str **map)
{
    pixel_str *pixel = (pixel_str *)malloc(sizeof(pixel_str));

    map[x + y*image->width] = region;

    pixel->x = x;
    pixel->y = y;

    pixel->flux = PIXEL_DOUBLE(image, x, y);
    pixel->flux_err = PIXEL_DOUBLE(errors, x, y);

    pixel->smooth = PIXEL_DOUBLE(smooth, x, y);

    pixel->is_edge = FALSE;

    add_to_list(region->pixels, pixel);
}

static void update_region_stats(region_str *region, image_str *image, region_str **map)
{
    pixel_str *pixel = list_first_item(region->pixels);
    double bg_sum = 0;
    double bg_sum2 = 0;

    region->npixels = 0;
    region->npixels_bg = 0;

    region->min_flux = pixel->flux;
    region->max_flux = pixel->flux;
    region->mean_flux = 0;

    region->bg_mean = 0;
    region->bg_sigma = 0;

    region->max = NULL;
    region->saddle = NULL;

    foreach(pixel, region->pixels){
        int i;

        pixel->is_edge = FALSE;

        region->min_flux = MIN(region->min_flux, pixel->flux);
        region->max_flux = MAX(region->max_flux, pixel->flux);
        region->mean_flux += pixel->flux;
        region->npixels ++;

        if(!region->max || region->max->flux < pixel->flux)
            region->max = pixel;

        for(i = 1; i < dN; i++){
            int x1 = pixel->x + dx[i];
            int y1 = pixel->y + dy[i];

            if(IS_PIXEL_VALID(image, x1, y1) && map[x1 + y1*image->width] && map[x1 + y1*image->width] != region){
                pixel->is_edge = TRUE;
                bg_sum += pixel->flux;
                bg_sum2 += pixel->flux*pixel->flux;
                region->npixels_bg ++;

                if(!region->saddle || region->saddle->flux < pixel->flux)
                    region->saddle = pixel;

                break;
            }
        }
    }

    region->mean_flux /= region->npixels;

    if(region->npixels_bg > 0){
        region->bg_mean = bg_sum/region->npixels_bg;
    }

    if(region->npixels_bg > 1){
        region->bg_sigma = sqrt((bg_sum2 - bg_sum*bg_sum/region->npixels_bg)/(region->npixels_bg - 1));
    }

    {
        double *v = (double *)malloc(sizeof(double)*region->npixels_bg);
        int i = 0;

        foreach(pixel, region->pixels){
            if(pixel->is_edge)
                v[i++] = pixel->flux;
        }

        region->bg_mean = get_median_mad(v, region->npixels_bg, &region->bg_sigma);
        region->bg_sigma *= 1.4826;

        free(v);
    }

    if(region->saddle){
        region->prominence = (region->max_flux - region->saddle->flux);
    }
}

static void dump_region_a(region_str *region, char *filename)
{
    FILE *file = fopen(filename, "a");
    pixel_str *pixel = NULL;

    foreach(pixel, region->pixels)
        if(pixel != region->saddle)
            fprintf(file, "%d %d %g %g %d\n", pixel->x, pixel->y, pixel->flux, pixel->flux - region->bg_mean, pixel->is_edge);
        else
            fprintf(file, "%d %d %g %g %d\n", pixel->x, pixel->y, pixel->flux, pixel->flux - region->bg_mean, 2);

    fclose(file);
}

void measure_peak(peak_str *peak, image_str *image)
{
    double flux = 0;
    double flux_err = 0;
    double r0 = 2;
    int x;
    int y;

    for(x = MAX(floor(peak->x - r0), 0); x <= MIN(ceil(peak->x + r0), image->width - 1); x++)
        for(y = MAX(floor(peak->y - r0), 0); y <= MIN(ceil(peak->y + r0), image->height - 1); y++)
            if(hypot(x - peak->x, y - peak->y) < r0)
                flux += PIXEL_DOUBLE(image, x, y);

    if(flux > 0)
        peak->flux = flux;
    else
        peak->excess = 0;
}

void find_peaks(image_str *image, image_str *smooth, image_str *errors, image_str *mask, double threshold, struct list_head *peaks)
{
    region_str **map = (region_str **)calloc(image->width*image->height, sizeof(region_str *));
    struct list_head regions;
    region_str *region = NULL;
    int x = 0;
    int y = 0;

    int Ntotal = 0;
    int Ngood = 0;

    init_list(*peaks);
    init_list(regions);

    for(y = 0; y < smooth->height; y++)
        for(x = 0; x < smooth->width; x++){
            int d;

            if(!map[x + y*image->width] && PIXEL_DOUBLE(smooth, x, y) > 1.5*PIXEL_DOUBLE(errors, x, y)){
                int x0;
                int y0;

                get_reachable_maximum_xy(smooth, NULL, x, y, &x0, &y0, 0);

                if(mask && ((PIXEL(mask, x, y) & FLAG_BAD) || (PIXEL(mask, x0, y0) & FLAG_BAD)))
                    /* Skip BAD pixels */
                    continue;

                region = map[x0 + y0*image->width];

                if(!region){
                    region = calloc(1, sizeof(region_str));
                    init_list(region->pixels);
                    add_to_list(regions, region);

                    add_pixel_to_region(x0, y0, region, image, smooth, errors, map);
                }

                if(x != x0 || y != y0)
                    add_pixel_to_region(x, y, region, image, smooth, errors, map);
            }
        }

    dprintf("%d initial regions\n", list_length(&regions));

    foreach(region, regions)
        update_region_stats(region, image, map);

    foreach(region, regions){
        int is_merged = FALSE;

        /* do { */
        /*     is_merged = FALSE; */

        /*     if(region->max->is_edge){ */
        /*         int x0 = region->saddle->x; */
        /*         int y0 = region->saddle->y; */
        /*         int i; */
        /*         int id = 0; */
        /*         double max = 0; */

        /*         for(i = 1; i < dN; i++){ */
        /*             if(IS_PIXEL_VALID(image, x0+dx[i], y0+dy[i]) /\* && map[x0+dx[i] + (y0+dy[i])*image->width] != region *\/){ */
        /*                 if(!id || PIXEL_DOUBLE(smooth, x0+dx[i], y0+dy[i]) > max){ */
        /*                     id = i; */
        /*                     max = PIXEL_DOUBLE(smooth, x0+dx[i], y0+dy[i]); */
        /*                 } */
        /*             } */
        /*         } */

        /*         if(id){ */
        /*             region_str *r = map[x0+dx[id] + (y0+dy[id])*image->width]; */

        /*             if(r && r != region){ */
        /*                 pixel_str *pixel = NULL; */

        /*                 foreach(pixel, r->pixels) */
        /*                     add_pixel_to_region(pixel->x, pixel->y, region, image, smooth, errors, map); */

        /*                 del_from_list(r); */
        /*                 free_list(r->pixels); */
        /*                 free(r); */

        /*                 is_merged = TRUE; */

        /*                 update_region_stats(region, image, map); */
        /*             } */
        /*         } */
        /*     } */
        /* } while(is_merged); */

        while(region->saddle){
            int x0 = region->saddle->x;
            int y0 = region->saddle->y;
            int i;

            is_merged = FALSE;

            for(i = 1; i < dN; i++)
                if(IS_PIXEL_VALID(image, x0+dx[i], y0+dy[i]) &&
                   map[x0+dx[i] + (y0+dy[i])*image->width] &&
                   map[x0+dx[i] + (y0+dy[i])*image->width] != region){
                    region_str *r = map[x0+dx[i] + (y0+dy[i])*image->width];

                    //dprintf("%d %d - %g - %i - %g\n", x0, y0, region->prominence, i, r->prominence);

                    /* if(r->max_flux - r->bg_mean < region->max_flux - region->bg_mean) */
                    /*     dprintf("%d %d - %d %d - %g\n", region->max->x, region->max->y, r->max->x, r->max->y, */
                    /*             (r->max_flux - r->saddle->flux)/(region->max_flux - region->bg_mean)); */

                    if((region->max_flux - region->saddle->flux)/(region->max_flux - 0*region->bg_mean) < 0.5 ||
                       (r->max_flux - r->saddle->flux)/(r->max_flux - 0*r->min_flux) < 0.5){
                        pixel_str *pixel = NULL;

                        /* dprintf("merging %d %d into %d %d - %g\n", r->max->x, r->max->y, region->max->x, region->max->y, */
                        /*         (r->max_flux - r->saddle->flux)/(region->max_flux - region->min_flux)); */
                        foreach(pixel, r->pixels)
                            add_pixel_to_region(pixel->x, pixel->y, region, image, smooth, errors, map);

                        del_from_list(r);
                        free_list(r->pixels);
                        free(r);

                        is_merged = TRUE;

                        update_region_stats(region, image, map);
                    } else {
                        /* dprintf("not merging %d %d into %d %d - %g\n", r->max->x, r->max->y, region->max->x, region->max->y, */
                        /*         (r->max_flux - r->saddle->flux)/(region->max_flux - region->min_flux)); */
                    }

                    /* dprintf("%d %d: min %g max %g N %d Nbg %d\n", */
                    /*         region->max->x, region->max->y, region->min_flux, region->max_flux, */
                    /*         region->npixels, region->npixels_bg); */
                }

            if(!is_merged)
                break;
        }
    }

    dprintf("%d regions after merging\n", list_length(&regions));

    system("rm -f out.all.txt");
    system("touch out.all.txt");

    foreach(region, regions){
        peak_str *peak = NULL;

        peak = peak_from_region(region, mask);

        measure_peak(peak, smooth);

        if(peak->excess > threshold || TRUE){
            dump_region_a(region, "out.all.txt");
        }

        if(peak->excess > threshold && /* region->npixels - region->npixels_bg > 10 && */
           /* region->max_flux - region->bg_mean > 5.0*region->bg_sigma && */ TRUE){
            peak->id = Ntotal++;
            add_to_list(*peaks, peak);
        } else
            free(peak);

        free_list(region->pixels);
    }

    dprintf("%d peaks\n", list_length(peaks));

    free_list(regions);
}

void dump_peaks_to_file(struct list_head *peaks, char *filename, int state)
{
    FILE *file = (!filename || filename[0] == '-') ? stdout : fopen(filename, "w");
    peak_str *peak = NULL;

    foreach(peak, *peaks){
        if(peak->state & state)
            fprintf(file, "%g %g %g %g %d %g %g %g %g %g %g %d\n", peak->x, peak->y, peak->flux, peak->dflux, peak->flags, peak->dx, peak->dy, peak->a, peak->b, peak->theta, peak->excess, peak->id);
    }

    if(file != stdout)
        fclose(file);
}
