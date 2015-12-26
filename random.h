#ifndef RANDOM_H
#define RANDOM_H

#ifndef	RAND_MAX
#define	RAND_MAX	0x7fffffffL	/* default dynamic range of rand() */
#endif

double random_double();
double random_gauss(double );
double random_poisson(double );

int random_int();

void random_initialize(int );

#endif /* RANDOM_H */
