#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "random.h"

/* This code is extracted from Emmanuel Bertin' SkyMaker source */
/* FIXME: replace with something faster, more accurate and with less restrictive license */

inline double random_gauss(double sigma)
{
    double x,y,z, r;

    while((z = pow(x = random_double() - 0.5, 2.0) + pow(y = random_double() - 0.5, 2.0)) > 0.25);
    while((r = random_double()) <= 0.0);

    return sigma*sqrt(-2.0*log(r)/z)*x;
}

/* inline double random_gauss(double s) */
/* { */
/*     double m = 0; */
/*     const double p0 = 0.322232431088;     const double q0 = 0.099348462606; */
/*     const double p1 = 1.0;                const double q1 = 0.588581570495; */
/*     const double p2 = 0.342242088547;     const double q2 = 0.531103462366; */
/*     const double p3 = 0.204231210245e-1;  const double q3 = 0.103537752850; */
/*     const double p4 = 0.453642210148e-4;  const double q4 = 0.385607006340e-2; */
/*     double u, t, p, q, z; */

/*     u   = random_double(); */
/*     if (u < 0.5) */
/*         t = sqrt(-2.0 * log(u)); */
/*     else */
/*         t = sqrt(-2.0 * log(1.0 - u)); */
/*     p   = p0 + t * (p1 + t * (p2 + t * (p3 + t * p4))); */
/*     q   = q0 + t * (q1 + t * (q2 + t * (q3 + t * q4))); */
/*     if (u < 0.5) */
/*         z = (p / q) - t; */
/*     else */
/*         z = t - (p / q); */
/*     return (m + s * z); */
/* } */

inline int random_int()
{
    return (int)rand();
}

inline double random_double()
{
    return (double)rand()/RAND_MAX;
}

void random_initialize(int seed)
{
    if(seed)
        srand((unsigned int)seed);
    else
        srand((unsigned int)time(NULL));
}

static double gammln(double xx)
{
    double x,tmp,ser;
    static double cof[6]={76.18009173, -86.50532033, 24.01409822,
                          -1.231739516, 0.120858003e-2, -0.536382e-5};
    int j;

    tmp = (x = xx-1.0) + 5.5;
    tmp -= (x + 0.5)*log(tmp);
    ser = 1.0;
    for (j = 0; j < 6; j++)
        ser += cof[j]/(x += 1.0);

    return log(2.50662827465*ser) - tmp;
}

inline double random_poisson(double xm)
{
    double gammln();
    double sq, alxm, g, oldm, em, t, y;

    sq = alxm = g = 0.0;
    oldm = -1.0;

    if (xm < 12.0){
        if (xm != oldm){
            oldm = xm;
            g = exp(-xm);
        }

        em = -1.0;
        t = 1.0;

        do {
            em += 1.0;
            t *= random_double();
        } while (t > g);
    } else {
        if (xm != oldm){
            oldm = xm;
            sq = sqrt(2.0*xm);
            alxm = log(xm);
            g = xm*alxm - gammln(xm + 1.0);
        } do{
            do{
                y = tan(M_PI*random_double());
                em = sq*y + xm;
            } while (em < 0.0);

            em = floor(em);
            t = 0.9*(1.0 + y*y)*exp(em*alxm - gammln(em + 1.0) - g);
        } while (random_double() > t);
    }

    return em;
}
