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
/* static int dx[] = {0, 0,  0, 1, -1}; */
/* static int dy[] = {0, 1, -1, 0, 0}; */
/* static int dN = 5; */

#define IS_PIXEL_VALID(image, x, y) ((x) >= 0 && (x) < (image)->width && (y) >= 0 && (y) < (image)->height)

/* Single image pixel */
typedef struct pixel_str {
    LIST_HEAD(struct pixel_str);

    int x;
    int y;
    double flux;
    double flux_err;
    double smooth;
    double excess;

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

    int min_x;
    int max_x;
    int min_y;
    int max_y;

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

    /* if(list_length(&region->pixels) > 100) */
    /*     dprintf("len = %d\n", list_length(&region->pixels)); */

    foreach(pixel, region->pixels){
        double value = pixel->smooth;

        if(value >= 1e-10){
            double x = pixel->x;
            double y = pixel->y;

            sum_flux += value;
            sum_flux_err += pixel->flux_err*pixel->flux_err;

            if(!(PIXEL(mask, pixel->x, pixel->y) & FLAG_SATURATED)){
                sum_x += value*x;
                sum_x2 += value*x*x;
                sum_y += value*y;
                sum_y2 += value*y*y;
                sum_xy += value*x*y;
                sum_weight += value;
            }

            flags |= PIXEL(mask, pixel->x, pixel->y);
            //if(x == 0 || x == mask->width-1 || y == 0 || y == mask->height-1)
            if(x <= 1 || x >= mask->width-2 || y <= 1 || y >= mask->height-2)
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

    pixel->excess = PIXEL_DOUBLE(smooth, x, y)/PIXEL_DOUBLE(errors, x, y);

    pixel->is_edge = FALSE;

    //dprintf("%d %d - %g\n", pixel->x, pixel->y, pixel->flux);

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

    region->min_x = pixel->x;
    region->max_x = pixel->x;
    region->min_y = pixel->y;
    region->max_y = pixel->y;

    foreach(pixel, region->pixels){
        int i;

        pixel->is_edge = FALSE;

        region->min_flux = MIN(region->min_flux, pixel->flux);
        region->max_flux = MAX(region->max_flux, pixel->flux);
        region->mean_flux += pixel->flux;
        region->npixels ++;

        region->min_x = MIN(region->min_x, pixel->x);
        region->max_x = MAX(region->min_x, pixel->x);
        region->min_y = MIN(region->min_y, pixel->y);
        region->max_y = MAX(region->min_y, pixel->y);

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

    if(region->npixels_bg > 4){
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
    //FILE *file = fopen(filename, "a");
    FILE *file = fopen(filename, "w");
    pixel_str *pixel = NULL;

    foreachback(pixel, region->pixels)
        if(pixel != region->saddle)
            fprintf(file, "%d %d %g %g %g %d\n", pixel->x, pixel->y, pixel->flux, pixel->flux - region->bg_mean, pixel->excess, pixel->is_edge);
        else
            fprintf(file, "%d %d %g %g %g %d\n", pixel->x, pixel->y, pixel->flux, pixel->flux - region->bg_mean, pixel->excess, 2);

    fclose(file);
}

void measure_peak(peak_str *peak, image_str *image, image_str *errors, image_str *mask)
{
    double flux = 0;
    double flux_err = 0;
    double flux0 = 0;
    double flux_err0 = 0;
    double r0 = 3;
    int flags = 0;
    int N = 0;
    int x;
    int y;

    for(x = MAX(floor(peak->x - r0), 0); x <= MIN(ceil(peak->x + r0), image->width - 1); x++)
        for(y = MAX(floor(peak->y - r0), 0); y <= MIN(ceil(peak->y + r0), image->height - 1); y++)
            if(hypot(x - peak->x, y - peak->y) < r0){
                flux += PIXEL_DOUBLE(image, x, y);
                flux_err += PIXEL_DOUBLE(errors, x, y)*PIXEL_DOUBLE(errors, x, y);
                flags |= PIXEL(mask, x, y);
                N ++;
            }

    /* for(x = MAX(floor(peak->x - r0), 0); x <= MIN(ceil(peak->x + r0), image->width - 1); x++) */
    /*     for(y = MAX(floor(peak->y - r0), 0); y <= MIN(ceil(peak->y + r0), image->height - 1); y++) */
    /*         if(hypot(x - peak->x, y - peak->y) < r0){ */
    /*             PIXEL_DOUBLE(image, x, y) = flux/N; */
    /*         } */

    if(flux > 0){
        peak->flux = flux;
        peak->dflux = sqrt(flux_err);
        peak->flags |= flags;

        peak->excess = peak->flux/peak->dflux;
    } else
        peak->excess = 0;
}

void dump_peak(peak_str *peak, image_str *image, image_str *errors, image_str *mask, char *filename)
{
    FILE *file = fopen(filename, "w");
    double r0 = 3;
    int x;
    int y;

    for(x = MAX(floor(peak->x - r0), 0); x <= MIN(ceil(peak->x + r0), image->width - 1); x++)
        for(y = MAX(floor(peak->y - r0), 0); y <= MIN(ceil(peak->y + r0), image->height - 1); y++)
            if(hypot(x - peak->x, y - peak->y) < 2*r0){
                fprintf(file, "%d %d %g %g %g %d\n", x, y, PIXEL_DOUBLE(image, x, y), PIXEL_DOUBLE(errors, x, y), PIXEL_DOUBLE(image, x, y)/PIXEL_DOUBLE(errors, x, y), hypot(x - peak->x, y - peak->y) < r0);
            }

    fclose(file);
}

static inline double sinc(double x)
{
    if(fabs(x) < 1e-4)
        return 1;
    else
        return sin(M_PI*x)/M_PI/x;
}

void measure_peak_new(peak_str *peak, image_str *image, image_str *mask)
{
    double flux = 0;
    //double flux_err = 0;
    double r0 = 3.5;
    int flags = 0;
    int N = 0;
    double x;
    double y;

    inline double interpolate_image(double x, double y)
    {
        double sum = 0;
        int ex = 4;
        int x1;
        int y1;

        for(x1 = MAX(0, round(x) - ex); x1 <= MIN(image->width - 1, round(x) + ex); x1++)
            for(y1 = MAX(0, round(y) - ex); y1 <= MIN(image->height - 1, round(y) + ex); y1++)
                sum += PIXEL_DOUBLE(image, x1, y1)*sinc(x - x1)*sinc(y - y1);

        if(!isfinite(sum)){
            dprintf("%g %g - %g\n", x, y, sum);

            exit(1);
        }

        return sum;
    }

    /* printf("\n"); */

    for(x = MAX(floor(peak->x - r0 - 1), 0); x <= MIN(ceil(peak->x + r0 + 1), image->width - 1); x++)
        for(y = MAX(floor(peak->y - r0 - 1), 0); y <= MIN(ceil(peak->y + r0 + 1), image->height - 1); y++){
            int nsteps = 3;
            double step = 1.0/nsteps;
            double dx;
            double dy;

            for(dy = -0.5; dy < 0.5-step; dy += step)
                for(dx = -0.5; dx < 0.5-step; dx += step){
                    if(hypot(x - peak->x, y - peak->y) < r0){
                        flux += interpolate_image(x + dx + 0.5*step, y + dy + 0.5*step)*step*step;
                    }
                }

            /* Flags */
            if(hypot(x - peak->x, y - peak->y) <= 3.0*r0){
                //flux += PIXEL_DOUBLE(image, x, y);
                flags |= PIXEL(mask, x, y);
                /* printf("%g %g - %d - %d\n", x, y, PIXEL(mask, x, y), flags); */
                N ++;
            }

        }


    if(flux > 0){
        peak->flux = flux;
        peak->flags |= flags;
        /* printf("%g %g %g %d\n", peak->x, peak->y, peak->flux, peak->flags); */
    } else
        peak->excess = 0;
}

int compare_peaks_fn(const void *v1, const void *v2)
{
    int res = (*(peak_str**)v1)->flux - (*(peak_str**)v2)->flux;

    return (res > 0 ? 1 : (res < 0 ? -1 : 0));
}

void sort_peaks(struct list_head *list)
{
    int N = list_length(list);
    peak_str **peaks = (peak_str **)malloc(sizeof(peak_str *)*N);
    peak_str *peak;
    int d = 0;

    foreach(peak, *list){
        peaks[d++] = peak;
        del_from_list_in_foreach_and_run(peak, {});
    }

    qsort(peaks, N, sizeof(peak_str *), compare_peaks_fn);

    for(d = 0; d < N; d++)
        add_to_list(*list, peaks[d]);
}

void find_peaks(image_str *image, image_str *smooth, image_str *errors, image_str *mask, double threshold, struct list_head *peaks)
{
    region_str **map = (region_str **)calloc(image->width*image->height, sizeof(region_str *));
    image_str *weighted = image_create_double(image->width, image->height);
    struct list_head regions;
    region_str *region = NULL;
    int x = 0;
    int y = 0;

    int Ntotal = 0;
    int Ngood = 0;

    init_list(*peaks);
    init_list(regions);

    for(x = 0; x < image->width*image->height; x++)
        weighted->double_data[x] = smooth->double_data[x]/errors->double_data[x];

    if(weighted->width > 50)
        image_dump_to_fits(weighted, "out.weighted.fits");
    else
        image_dump_to_fits(weighted, "out.weighted.small.fits");

    for(y = 0; y < smooth->height; y++)
        for(x = 0; x < smooth->width; x++){
            int d;

            if(!map[x + y*image->width] && PIXEL_DOUBLE(smooth, x, y) > 2.0*PIXEL_DOUBLE(errors, x, y)){
                int x0;
                int y0;

                get_reachable_maximum_xy(weighted, mask /* NULL */, x, y, &x0, &y0, 0);

                if(mask && ((PIXEL(mask, x, y) & FLAG_BAD) || (PIXEL(mask, x, y) & FLAG_SATURATED)))
                    /* Skip BAD pixels */
                    continue;

                region = map[x0 + y0*image->width];

                if(!region){
                    //dprintf("\n\n");
                    region = calloc(1, sizeof(region_str));
                    init_list(region->pixels);
                    add_to_list(regions, region);

                    add_pixel_to_region(x0, y0, region, image, smooth, errors, map);
                }

                if(x != x0 || y != y0)
                    add_pixel_to_region(x, y, region, image, smooth, errors, map);
            }
        }

    //dprintf("%d initial regions\n", list_length(&regions));

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

    foreach(region, regions){
        if(list_length(&region->pixels) > 1000){
            free_list(region->pixels);
            del_from_list_in_foreach_and_run(region, free(region));
        }
    }

    foreach(region, regions){
        double fill = 1.0*region->npixels/(region->max_x + 1 - region->min_x)/(region->max_y + 1 - region->min_y);

        /* dprintf("%d %g\n", region->npixels, fill); */

        /* if(fill < 0.25){ */
        /*     free_list(region->pixels); */
        /*     del_from_list_in_foreach_and_run(region, free(region)); */
        /* } */
    }

    //dprintf("%d regions after merging\n", list_length(&regions));

    /* system("rm -f out.all.txt"); */
    /* system("touch out.all.txt"); */

    foreach(region, regions){
        peak_str *peak = NULL;

        peak = peak_from_region(region, mask);

        measure_peak(peak, image, errors, mask);
        //measure_peak_new(peak, image, mask);

        /* if(peak->excess > threshold || TRUE){ */
        /*     dump_region_a(region, "out.all.txt"); */
        /* } */

        if(peak->excess > threshold &&
           (region->npixels > 2 || (!(peak->flags & FLAG_SATURATED) && region->npixels > 2)) &&
           /* region->npixels - region->npixels_bg > 10 && */
           /* region->max_flux - region->bg_mean > 5.0*region->bg_sigma && */ TRUE){
            peak->id = Ntotal++;
            add_to_list(*peaks, peak);

            /* if(peak->id == 4){ */
            /*     dump_peak(peak, image, errors, mask, "out.all2.txt"); */
            /*     dump_region_a(region, "out.all.txt"); */
            /*     //exit(1); */
            /* } */
        } else
            free(peak);

        free_list(region->pixels);
    }

    sort_peaks(peaks);

    //dprintf("%d peaks\n", list_length(peaks));

    image_delete(weighted);

    free_list(regions);
}

void load_peaks(char *filename, image_str *image, image_str *errors, image_str *mask, int is_radec, struct list_head *peaks)
{
    FILE *file = fopen(filename, "r");
    int Ntotal = 0;

    init_list(*peaks);

    while(!feof(file)){
        double x;
        double y;

        fscanf(file, "%lf %lf\n", &x, &y);

        if(is_radec)
            coords_get_x_y(&image->coords, x, y, &x, &y);

        if(x >= 0 && x < image->width &&
           y >= 0 && y < image->height){
            peak_str *peak = calloc(1, sizeof(peak_str));

            peak->x = x;
            peak->y = y;
            peak->state = PEAK_INITIAL;
            peak->id = Ntotal++;

            measure_peak(peak, image, errors, mask);
            //measure_peak_new(peak, image, mask);

            add_to_list(*peaks, peak);
        }
    }

    sort_peaks(peaks);

    fclose(file);
}

void dump_peaks_to_file(struct list_head *peaks, char *filename, int state)
{
    FILE *file = (!filename || filename[0] == '-') ? stdout : fopen(filename, "w");
    peak_str *peak = NULL;

    foreach(peak, *peaks){
        if(peak->state & state)
            fprintf(file, "%g %g %g %g %d %g %g %g %g %g %d %g %g\n", peak->x, peak->y, peak->flux, peak->dflux, peak->flags, peak->a, peak->b, peak->theta, peak->bg, peak->chisq, peak->id, peak->ra, peak->dec);
    }

    if(file != stdout)
        fclose(file);
}
