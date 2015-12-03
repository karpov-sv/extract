#ifndef UTILS_H
#define UTILS_H

/* Misc functions - strings etc */
#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>

/* Optional debug messages */
#ifdef DEBUG
#define dprintf(x...) fprintf(stderr, x)
#else
#define dprintf(x...) {;}
#endif

#ifndef	FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE  1
#endif

/* FIXME: possible side-effects in macros, as these calls x or y twice. Safer way:
   #define MAX(a,b) \
       ({ typeof (a) _a = (a); \
           typeof (b) _b = (b); \
         _a > _b ? _a : _b; })
*/
#ifndef MIN
//#define MIN(x, y) (((x) > (y)) ? y : x)
#define MIN(a, b) \
    ({ typeof (a) _a = (a);     \
        typeof (b) _b = (b);    \
        _a < _b ? _a : _b; })
#endif /* MIN */
#ifndef MAX
//#define MAX(x, y) (((x) > (y)) ? x : y)
#define MAX(a, b) \
    ({ typeof (a) _a = (a);     \
        typeof (b) _b = (b);    \
        _a > _b ? _a : _b; })
#endif /* MAX */

/* Parse the command-line arguments
 * Accepts an argc, argv and the "name=format",&value pairs.
 * Finally, NULL should be the last argument
 */
void parse_args(int , char **, ...);

/* Exit with error message */
void exit_with_error(const char *, ...) __attribute__ ((format (printf, 1, 2)));
/* Print an error message and return -1 */
int return_with_error(const char *, ...) __attribute__ ((format (printf, 1, 2)));

/* String by printf-like template */
char *make_string(const char *, ...) __attribute__ ((format (printf, 1, 2)));

/* Append to existing string */
void add_to_string(char **, const char *, ...);
/* Hexadecimal string */
char *make_hex_string(char *, int );
/* Memory dump, pretty-printed */
char *make_long_hex_string(char *, int );
/* DD MM SS.SS */
char *make_sexagesimal_string(double );
/* Binary */
char *make_binary_string(int );

/* Check whether the file exists and is "normal" */
int file_exists_and_normal(char *);
long long int get_file_size(char *);

/* mkstemp wrapper - creates temp file and returns its path */
char *make_temp_filename(const char *);
/* mkdtemp wrapper */
char *make_temp_dirname(const char *);

/* Parse an hh:mm:ss string */
double parse_angular_string(char *);

/* Modulo */
inline double mod(double , double );

/* Sign */
inline double sign(double );

/* Rough estimation of the median */
double get_median(double *, int );
double get_median_mad(double *, int , double *);
double get_quantile(double *, int , double );

/* Mean and variance of the array */
double get_mean(double *, int , double *);

/* Create a lambda function.  Note: unlike lambdas in functional
   languages, this lambda does not capture the containing
   environment.  Thus, if you access the enclosing environment, you
   must ensure that the lifetime of this lambda is bound by the
   lifetime of the enclosing environment (i.e., until the enclosing
   function returns).  This means that if you access local
   variables, bad things will happen.  If you don't access local
   variables, you're fine.  */
#define lambda(l_ret_type, l_arguments, l_body)                 \
    ({                                                          \
        l_ret_type l_anonymous_functions_name l_arguments       \
            l_body                                              \
            &l_anonymous_functions_name;                        \
    })

/* CRC8 */
unsigned char crc8(unsigned char *, int );

/* Run system() with printf() arguments */
int system_run(const char *, ...);
int system_run_silently(const char *, ...);

/* free() the pointer and set its value to NULL */
#define free_and_null(ptr) do {\
        if(ptr) free(ptr); (ptr) = NULL;  \
    } while(0)


int remove_dir(char *);

u_int64_t free_disk_space(char *);

#ifdef __CYGWIN__
/* Cygwin does not have some long double functions */
#define fmodl(x, y) fmod((x), (y))
#define floorl(x) floor((x))
#define powl(x, y) pow((x), (y))
#endif

int fast_system(char *, int );

int copy_file(char *, char *);

#endif /* UTILS_H */
