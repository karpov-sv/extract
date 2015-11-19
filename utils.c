/* Misc functions - strings etc */

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <signal.h>
#include <math.h>
#include <ctype.h>
#include <ftw.h>
#include <sys/statvfs.h>
#include <fcntl.h>

#include "utils.h"

/* #include "popen_noshell.h" */

/* Parse the command-line arguments
 * Accepts an argc, argv and the "name=format",&value pairs.
 * Finally, NULL should be the last argument.
 * All arguments used become a null strings.
 */
void parse_args(int argc, char **argv, ...)
{
    va_list l;
    char *pattern;
    int *already_parsed = NULL;
    int print_help = FALSE;
    char **strings = malloc(sizeof(char *)*argc);
    int nstrings = argc;
    int d;

    /* Fill strings array to process */
    for(d = 0; d < argc; d++){
        strings[d] = argv[d];
    }

    /* Add strings from config file, if any */
    {
        char *filename = make_string("%s.opts", argv[0]);

        if(file_exists_and_normal(filename)){
            FILE *file = fopen(filename, "r");

            /* dprintf("Reading options from %s\n", filename); */

            while(!feof(file)){
                char *buffer = calloc(1024, sizeof(char));

                fgets(buffer, 1024, file);

                if(buffer){
                    int pos = strlen(buffer) - 1;
                    char *string = buffer;

                    /* Strip trailing whitespaces */
                    while(pos >= 0 && isspace(*(buffer + pos)))
                        *(buffer + pos) = '\0';

                    /* Strip leading ones */
                    while(isspace(*string))
                        string ++;

                    if(*string){
                        /* dprintf("-> %s\n", string); */
                        strings = realloc(strings, sizeof(char*)*(nstrings + 1));
                        strings[nstrings] = string;
                        nstrings ++;
                    } else
                        free(buffer);
                }
            }

            fclose(file);
        }

        free(filename);
    }

    /* Process them */

    va_start(l, argv);

    for(d = 1; d < argc; d++)
        if(!strcmp(argv[d], "-help"))
            print_help = TRUE;

    if(print_help)
        printf("Usage: %s <options>\n"
               "where options are:\n", argv[0]);

    already_parsed = (int*)calloc(nstrings, sizeof(int));

    while((pattern = va_arg(l, char*))){
        void *pointer = va_arg(l, void*);
        char *pos = strchr(pattern, '=');
        int pattern_used = FALSE;

        if(print_help)
            printf("\t%s\n", pattern);
        else if(pos){
            for(d = 1; d < nstrings; d++)
                if(!strncmp(strings[d], pattern, pos - pattern + 1)){

                    /* dprintf(":: %s == %s\n", strings[d], pattern); */

                    if(!pattern_used){
                        if(strstr(pattern, "%s")){
                            /* Special handling for strings... */
                            char *tmp = make_string("%s", strings[d] + (pos - pattern + 1));

                            *(char**)pointer = tmp;
                        } else
                            sscanf(strings[d], pattern, pointer);

                        pattern_used = TRUE;
                    }

                    already_parsed[d] = TRUE;

                    //break;
                }
        } else if(!strstr(pattern, "%") && pattern[0] == '-'){
            /* Boolean argument - set variable to 'TRUE' if found */
            /* Must be integer!!! */
            for(d = 1; d < nstrings; d++)
                if((strings[d][0] == '-' || strings[d][0] == '+') && !strncmp(strings[d] + 1, pattern + 1, pos - pattern)){
                    if(!pattern_used){
                        *(int*)pointer = (strings[d][0] == '-');

                        pattern_used = TRUE;
                    }

                    /* dprintf(":: %s == %s\n", strings[d], pattern); */
                    already_parsed[d] = TRUE;

                    //break;
                }
        } else if(!strcmp(pattern, "%s")){
            /* Parse first unused element to a string */
            for(d = 1; d < nstrings; d++)
                if(!already_parsed[d]){
                    if(!pattern_used){
                        /* FIXME: this will cause a memleak later */
                        char *tmp = make_string("%s", strings[d]);

                        /* dprintf(":: %s == %s\n", strings[d], pattern); */

                        *(char**)pointer = tmp;
                        pattern_used = TRUE;
                    }

                    already_parsed[d] = TRUE;

                    break;
                }
        }
    }

    /* Make all used args an empty strings */
    for(d = 1; d < nstrings; d++)
        if(already_parsed[d])
            strings[d][0] = '\0';

    va_end(l);

    if(already_parsed)
        free(already_parsed);

    if(strings)
        free(strings);

    if(print_help)
        exit(EXIT_SUCCESS);
}

/* Exit with error message */
void exit_with_error(const char *template, ...)
{
    va_list ap;
    char *buffer = NULL;

    va_start(ap, template);
    vasprintf((char**)&buffer, template, ap);
    va_end(ap);

    perror(buffer);

    exit(EXIT_FAILURE);
}

/* Print an error message and return -1 */
int return_with_error(const char *template, ...)
{
    va_list ap;
    char *buffer = NULL;

    va_start(ap, template);
    vasprintf((char**)&buffer, template, ap);
    va_end(ap);

    perror(buffer);

    return -1;
}

/* Make string by printf-like template
 *
 * vasprintf is a GNU extension, so it maybe
 * need to be rewritten in platform-independent way ?
 */
char *make_string(const char *template, ...)
{
    va_list ap;
    char *buffer = NULL;

    va_start(ap, template);
    vasprintf((char**)&buffer, template, ap);
    va_end(ap);

    return buffer;
}

/* Append to existing string */
void add_to_string(char **string, const char *template, ...)
{
    va_list ap;
    char *buffer = NULL;
    int length = 0;

    va_start(ap, template);
    vasprintf((char**)&buffer, template, ap);
    va_end(ap);

    length = (*string) ? strlen(*string) : 0;

    *string = realloc(*string, length + strlen(buffer) + 1);
    memcpy((char*)*string + length, buffer, strlen(buffer) + 1);

    if(buffer)
        free(buffer);
}

char *make_hex_string(char *data, int length)
{
    char *str = NULL;
    char *chars = NULL;
    int d;

    for(d = 0; d < length; d++){
        unsigned char value = data[d];

        add_to_string(&str, "%02X ", value);

        if(data[d] > 31)
            add_to_string(&chars, "%c", value);
        else
            add_to_string(&chars, ".");
    }

    if(chars){
        add_to_string(&str, "| %s", chars);
        free(chars);
    }

    return str;
}

char *make_long_hex_string(char *data, int length)
{
    char *str = NULL;
    char *chars = NULL;
    int d;

    for(d = 0; d < length; d++){
        unsigned char value = data[d];

        if(! (d % 16))
            add_to_string(&str, "%04X : ", d);
        add_to_string(&str, "%02X ", value);
        if(value > 31)
            add_to_string(&chars, "%c", value);
        else
            add_to_string(&chars, ".");
        if(! ((d - 15) % 16) || d == length - 1){
            add_to_string(&str, "| %s \n", chars);
            free_and_null(chars);
        }
    }

    return str;
}

char *make_sexagesimal_string(double x)
{
    double val = fabs(x);
    int deg = floor(val);
    int min = floor((val - deg)*60.0);
    double sec = 60.*((val - deg)*60 - min);

    if(x < 0)
        deg = -deg;

    return make_string("%02d %02d %05.2lf", deg, min, sec);
}

char *make_binary_string(int x)
{
    char *str = NULL;
    int z;

    for (z = 128; z > 0; z >>= 1){
        add_to_string(&str, ((x & z) == z) ? "1" : "0");
    }

    return str;
}

/* Check whether the file exists and is "normal" */
int file_exists_and_normal(char *filename)
{
    struct stat status;
    int result = !stat(filename, &status);

    if (result &&
        !S_ISREG(status.st_mode) &&
        !S_ISLNK(status.st_mode))
        result = FALSE;

    if(result &&
       access(filename, R_OK) != 0)
        /* Access denied */
        result = FALSE;

    return result;
}

long long int get_file_size(char *filename)
{
    struct stat status;
    int result = !stat(filename, &status);
    long long int size = -1;

    if(result)
        size = status.st_size;

    return size;
}

/* mkstemp wrapper - creates temp file and returns its path */
char *make_temp_filename(const char *pattern)
{
    char *filename = make_string("%s", pattern);
    int fd = mkstemp(filename);

    close(fd);

    return filename;
}

/* mkdtemp wrapper - creates temp dir and returns its path */
char *make_temp_dirname(const char *pattern)
{
    char *dirname = make_string("%s", pattern);

    mkdtemp(dirname);

    return dirname;
}

/* Parse an hh:mm:ss string */
double parse_angular_string(char *string)
{
    char *tmp = make_string("%s", string);
    char *pos = tmp;
    double result = 0;
    double sign = 1;
    double hour = 0;
    double min = 0;
    double sec = 0;
    double msec = 0;

    /* replace all delimiters with spaces */
    while(*pos){
        if(*pos == '.' || *pos == ':' || *pos == 'h' ||
           *pos == 'd' || *pos == 'm' || *pos == 's' ||
           *pos == '\'' || *pos == '"' /* || *pos == 0x9c */)
            *pos = ' ';
        pos ++;
    }

    sscanf(tmp, "%lf %lf %lf %lf",
           &hour, &min, &sec, &msec);

    if(hour < 0){
        hour = -hour;
        sign = -1;
    }

    result = sign*(hour + (min + (sec + msec*1./60)*1./60)*1./60);

    free(tmp);

    return result;
}

inline double mod(double x, double y)
{
    int k = floor(x / y);

    return x - k * y;
}

inline double sign(double x)
{
    return (x > 0 ? 1 : (x < 0 ? -1 : 0));
}

double get_mean(double *x, int N, double *sd_ptr)
{
    double sum = 0;
    double sum2 = 0;
    int d;

    for(d = 0; d < N; d++){
        sum += x[d];
        sum2 += x[d]*x[d];
    }

    if(sd_ptr)
        *sd_ptr = sqrt((sum2 - sum*sum/N)/(N - 1));

    return sum/N;
}

unsigned char crc8(unsigned char *data, int length)
{
    unsigned char crc = 0xFF;
    unsigned int i;

    while (length--){
        crc ^= *data++;

        for (i = 0; i < 8; i++)
            crc = crc & 0x80 ? (crc << 1) ^ 0x31 : crc << 1;
    }

    return crc;
}

int system_run(const char *template, ...)
{
    va_list ap;
    char *buffer;
    int result = 0;

    va_start(ap, template);
    vasprintf((char**)&buffer, template, ap);
    va_end(ap);

    result = fast_system(buffer, TRUE);

    free(buffer);

    return result;
}

int system_run_silently(const char *template, ...)
{
    va_list ap;
    char *buffer;
    int result = 0;

    va_start(ap, template);
    vasprintf((char**)&buffer, template, ap);
    va_end(ap);

    result = fast_system(buffer, TRUE);

    free(buffer);

    return result;
}

int unlink_cb(const char *fpath, const struct stat *sb, int typeflag, struct FTW *ftwbuf)
{
    int rv = remove(fpath);

    if (rv)
        perror(fpath);

    return rv;
}

int remove_dir(char *path)
{
    return nftw(path, unlink_cb, 64, FTW_DEPTH | FTW_PHYS);
}

u_int64_t free_disk_space(char *path)
{
    struct statvfs s;

    statvfs(path, &s);

    return s.f_bavail*s.f_frsize;
}

int fast_system(char *command, int verbose)
{
/* #ifndef __linux__ */
    return system(command);
/* #else */
/*     FILE *fp; */
/*     struct popen_noshell_pass_to_pclose pclose_arg; */
/*     char buf[256]; */
/*     int status; */


/*     fp = popen_noshell_compat(command, "r", &pclose_arg); */

/*     if(!fp){ */
/*         dprintf("Command execution failed: %s\n", command); */

/*         return -1; */
/*     } */

/*     while (fgets(buf, sizeof(buf)-1, fp)) { */
/*         if(verbose) */
/*             printf("%s", buf); /\* Dump all subprocess output to simulate system() behaviour *\/ */
/*     } */

/*     status = pclose_noshell(&pclose_arg); */

/*     return status; */
/* #endif */
}

/*
 *  This Quickselect routine is based on the algorithm described in
 *  "Numerical recipes in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 *  This code by Nicolas Devillard - 1998. Public domain.
 */

#define ELEM_SWAP(a,b) { register double t=(a);(a)=(b);(b)=t; }

double quick_select(double arr[], int n)
{
    int low, high ;
    int median;
    int middle, ll, hh;

    low = 0 ; high = n-1 ; median = (low + high) / 2;
    for (;;) {
        if (high <= low) /* One element only */
            return arr[median] ;

        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]) ;
            return arr[median] ;
        }

        /* Find median of low, middle and high items; swap into position low */
        middle = (low + high) / 2;
        if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]) ;
        if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]) ;
        if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]) ;

        /* Swap low item (now in position middle) into position (low+1) */
        ELEM_SWAP(arr[middle], arr[low+1]) ;

        /* Nibble from each end towards middle, swapping items when stuck */
        ll = low + 1;
        hh = high;
        for (;;) {
            do ll++; while (arr[low] > arr[ll]) ;
            do hh--; while (arr[hh]  > arr[low]) ;

            if (hh < ll)
                break;

            ELEM_SWAP(arr[ll], arr[hh]) ;
        }

        /* Swap middle item (in position low) back into correct position */
        ELEM_SWAP(arr[low], arr[hh]) ;

        /* Re-set active partition */
        if (hh <= median)
            low = ll;
        if (hh >= median)
            high = hh - 1;
    }
}

#undef ELEM_SWAP

double get_median(double *x, int N)
{
    double *arr = (double *)malloc(sizeof(double)*N);
    double value = x[0];

    memcpy(arr, x, N*sizeof(double));

    value = quick_select(arr, N);

    free(arr);

    return value;
}

double get_median_mad(double *x, int N, double *mad_ptr)
{
    double *arr = (double *)malloc(sizeof(double)*N);
    double value = x[0];

    memcpy(arr, x, N*sizeof(double));

    value = quick_select(arr, N);

    if(mad_ptr){
        int i;

        for(i = 0; i < N; i++)
            arr[i] = fabs(arr[i] - value);

        *mad_ptr = quick_select(arr, N);
    }

    free(arr);

    return value;
}

int copy_file(char *from_name, char *to_name)
{
    char buf[BUFSIZ];
    size_t size;

    int source = open(from_name, O_RDONLY, 0);
    int dest = open(to_name, O_WRONLY | O_CREAT | O_TRUNC, 0644);

    while ((size = read(source, buf, BUFSIZ)) > 0) {
        write(dest, buf, size);
    }

    close(source);
    close(dest);

    return 0;
}
