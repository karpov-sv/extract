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

static int dx5[] = {0, 0,  0, 1, -1};
static int dy5[] = {0, 1, -1, 0, 0};
static int dN5 = 5;

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

region_str *extract_region(image_str *image, image_str *smooth, image_str *mask, int x0, int y0, double threshold)
{
    region_str *region = (region_str *)calloc(1, sizeof(region_str));
    int stack_length = 100;
    int stack_pos = 0;
    int *stack_x = (int *)malloc(stack_length*sizeof(int));
    int *stack_y = (int *)malloc(stack_length*sizeof(int));
    int x = x0;
    int y = y0;

    double bg_sum = 0;
    double bg_sum2 = 0;

    pixel_str *pixel = NULL;

    /* dprintf("\n%d %d - %g\n", x0, y0, threshold); */

    init_list(region->pixels);

    region->npixels = 0;
    region->npixels_bg = 0;

    region->is_edge = FALSE;
    region->is_noise = FALSE;

    region->min_flux = PIXEL_DOUBLE(image, x, y);
    region->max_flux = PIXEL_DOUBLE(image, x, y);
    region->mean_flux = 0;

    region->bg_mean = 0;
    region->bg_sigma = 0;

    inline void claim_pixel(x, y)
    {
        pixel_str *pixel = (pixel_str *)malloc(sizeof(pixel_str));

        PIXEL(mask, x, y) = TRUE;

        pixel->x = x;
        pixel->y = y;

        pixel->flux = PIXEL_DOUBLE(image, x, y);
        pixel->flux_err = 0;//

        pixel->is_edge = FALSE;

        region->min_flux = MIN(region->min_flux, pixel->flux);
        region->max_flux = MAX(region->max_flux, pixel->flux);
        region->mean_flux += pixel->flux;

        add_to_list(region->pixels, pixel);

        region->npixels ++;
    }

    claim_pixel(x, y);

    /* Recursive descent */
    while(TRUE){
        int i = 0;
        double value0 = PIXEL_DOUBLE(smooth, x, y);
        int Nadded = 0;

        /* Find and claim all destinations leading down from this point */
        for(i = 1; i < dN5; i++){
            int x1 = x + dx5[i];
            int y1 = y + dy5[i];

            if(IS_PIXEL_VALID(smooth, x1, y1) && !PIXEL(mask, x1, y1)){
                if((PIXEL_DOUBLE(smooth, x1, y1) <= value0 + 1e-10) ||
                   (threshold > 0 && PIXEL_DOUBLE(smooth, x1, y1) >= threshold + 1e-10)){
                    if(stack_pos == stack_length){
                        stack_length *= 2;
                        stack_x = realloc(stack_x, sizeof(int)*stack_length);
                        stack_y = realloc(stack_y, sizeof(int)*stack_length);
                    }

                    claim_pixel(x1, y1);

                    stack_x[stack_pos] = x1;
                    stack_y[stack_pos] = y1;

                    stack_pos ++;

                    Nadded ++;
                }
            } else
                region->is_edge = TRUE;
        }

        /* Go to the stack upper point, if any */
        if(stack_pos > 0){
            stack_pos --;
            x = stack_x[stack_pos];
            y = stack_y[stack_pos];
        } else
            /* Processed all points */
            break;
    }

    foreach(pixel, region->pixels){
        int i;

        pixel->is_edge = FALSE;

        for(i = 1; i < dN; i++){
            int x1 = pixel->x + dx[i];
            int y1 = pixel->y + dy[i];

            if(IS_PIXEL_VALID(smooth, x1, y1) && !PIXEL(mask, x1, y1)){
                pixel->is_edge = TRUE;
                bg_sum += pixel->flux;
                bg_sum2 += pixel->flux*pixel->flux;
                region->npixels_bg ++;

                break;
            }
        }
    }

    region->mean_flux /= region->npixels;

    region->bg_mean = bg_sum/region->npixels_bg;
    region->bg_sigma = sqrt((bg_sum2 - bg_sum*bg_sum/region->npixels_bg)/(region->npixels_bg - 1));

    foreach(pixel, region->pixels){
        PIXEL(mask, pixel->x, pixel->y) = FALSE;
    }

    free(stack_x);
    free(stack_y);

    /* dprintf("%d pixels, %d bg, %g %g\n", region->npixels, region->npixels_bg, region->bg_mean, region->bg_sigma); */

    if(!threshold && region->npixels_bg > 10){
        threshold = region->bg_mean + 3.0*region->bg_sigma;

        if(isfinite(threshold) && threshold > 0){
            free_list(region->pixels);
            free(region);

            region = extract_region(image, smooth, mask, x0, y0, threshold);
        }
    }

    return region;
}

void mask_region(image_str *mask, region_str *region)
{
    pixel_str *pixel = NULL;

    foreach(pixel, region->pixels)
        if(!pixel->is_edge){
            PIXEL(mask, pixel->x, pixel->y) = TRUE;
        }
}

void dump_region(region_str *region, char *filename)
{
    FILE *file = fopen(filename, "w");
    pixel_str *pixel = NULL;

    foreach(pixel, region->pixels)
        fprintf(file, "%d %d %g %g %d\n", pixel->x, pixel->y, pixel->flux, pixel->flux - region->bg_mean, pixel->is_edge);

    fclose(file);
}

static peak_str *peak_from_region(region_str *region)
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
    pixel_str *pixel;

    double X2 = 0;
    double Y2 = 0;
    double XY = 0;

    double bg = region->bg_mean + region->bg_sigma;

    foreach(pixel, region->pixels){
        /* double value = pixel->flux - bg; */
        double value = pixel->smooth;
        //value = MAX(0, value);

        if(value >= 1e-10){
            sum_x += value*pixel->x;
            sum_x2 += value*pixel->x*pixel->x;
            sum_y += value*pixel->y;
            sum_y2 += value*pixel->y*pixel->y;
            sum_xy += value*pixel->x*pixel->y;
            sum_weight += value;

            sum_flux += value;//pixel->flux - bg;

            N ++;
        }
    }

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

    peak->bg = bg;

    /* Rough estimate of peak significance */
    peak->excess = peak->flux/region->bg_sigma/sqrt(region->npixels);

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

static inline void add_pixel_to_region(int x, int y, region_str *region, image_str *image, image_str *smooth, region_str **map)
{
    pixel_str *pixel = (pixel_str *)malloc(sizeof(pixel_str));

    map[x + y*image->width] = region;

    pixel->x = x;
    pixel->y = y;

    pixel->flux = PIXEL_DOUBLE(image, x, y);
    pixel->flux_err = 0;//

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

            if(IS_PIXEL_VALID(image, x1, y1) && map[x1 + y1*image->width] != region){
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

void find_peaks(image_str *image, image_str *smooth, struct list_head *peaks)
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

            if(!map[x + y*image->width] && PIXEL(smooth, x, y) > 0){
                int x0;
                int y0;

                get_reachable_maximum_xy(smooth, NULL, x, y, &x0, &y0, 0);

                region = map[x0 + y0*image->width];

                if(!region){
                    region = calloc(1, sizeof(region_str));
                    init_list(region->pixels);
                    add_to_list(regions, region);

                    add_pixel_to_region(x0, y0, region, image, smooth, map);
                }

                if(x != x0 || y != y0)
                    add_pixel_to_region(x, y, region, image, smooth, map);
            }
        }

    dprintf("%d initial regions\n", list_length(&regions));

    foreach(region, regions)
        update_region_stats(region, image, map);

    foreach(region, regions){
        peak_str *peak = NULL;

        peak = peak_from_region(region);

        if(peak->excess > 3.0 && region->npixels - region->npixels_bg > 1 &&
           /* region->max_flux - region->bg_mean > 5.0*region->bg_sigma && */ TRUE){
            add_to_list(*peaks, peak);
        } else
            free(peak);

        free_list(region->pixels);
    }

    dprintf("%d peaks\n", list_length(peaks));

    free_list(regions);
}

void dump_peaks_to_file(struct list_head *peaks, char *filename)
{
    FILE *file = (!filename || filename[0] == '-') ? stdout : fopen(filename, "w");
    peak_str *peak = NULL;

    foreach(peak, *peaks){
        if(peak->state == PEAK_MEASURED)
            fprintf(file, "%g %g %g %g %g %g %g %g\n", peak->x, peak->y, peak->flux, peak->bg, peak->dx, peak->dy, peak->dflux, peak->dbg);
    }

    if(file != stdout)
        fclose(file);
}

void dump_peaks_to_file_measured(struct list_head *peaks, char *filename)
{
    FILE *file = (!filename || filename[0] == '-') ? stdout : fopen(filename, "w");
    peak_str *peak = NULL;

    foreach(peak, *peaks){
        if(peak->state == PEAK_MEASURED)
            fprintf(file, "%g %g %g %g %g %g %g %g %g %g %g %g\n", peak->x, peak->y, peak->flux, peak->bg, peak->dx, peak->dy, peak->dflux, peak->dbg, peak->a, peak->b, peak->theta, peak->excess);
    }

    if(file != stdout)
        fclose(file);
}

void dump_peaks_to_file_failed(struct list_head *peaks, char *filename)
{
    FILE *file = (!filename || filename[0] == '-') ? stdout : fopen(filename, "w");
    peak_str *peak = NULL;

    foreach(peak, *peaks){
        if(peak->state == PEAK_FAILED)
            fprintf(file, "%g %g %g %g %g %g %g %g %g %g %g %g\n", peak->x, peak->y, peak->flux, peak->bg, peak->dx, peak->dy, peak->dflux, peak->dbg, peak->a, peak->b, peak->theta, peak->excess);
    }

    if(file != stdout)
        fclose(file);
}

void dump_peaks_to_file_full(struct list_head *peaks, char *filename)
{
    FILE *file = (!filename || filename[0] == '-') ? stdout : fopen(filename, "w");
    peak_str *peak = NULL;

    foreach(peak, *peaks){
        fprintf(file, "%g %g %g %g %g %g %g %g %g %g %g %g\n", peak->x, peak->y, peak->flux, peak->bg, peak->dx, peak->dy, peak->dflux, peak->dbg, peak->a, peak->b, peak->theta, peak->excess);
    }

    if(file != stdout)
        fclose(file);
}
