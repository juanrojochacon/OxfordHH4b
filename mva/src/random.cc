// random.cc - nh 22/10/14

#include "random.h"

#include <iostream>

#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

static gsl_rng *random_generator = 0;

static void rng_check()
{
	if (random_generator == 0)
	{
		std::cout << "RNG Warning: generator not initialised"<<std::endl;
		std::cout << "Initialsing RNG with default seed"<<std::endl;
		rng_init(12424645645243);
	}
}

// Initialise RNG
void rng_init(unsigned long int const& s)
{
	if (random_generator != 0)
		gsl_rng_free(random_generator);

	random_generator = gsl_rng_alloc(gsl_rng_ranlux);
	gsl_rng_set(random_generator, s);

	std::cout << "RNG Initialised. Seed = " <<s<<"."<<std::endl;
}

//return an integer from 0 to n-1
unsigned long int rng_uniform(unsigned long int const& n)
{
	rng_check();
	return gsl_rng_uniform_int(random_generator, n);
}

// return a double from [0,1)
double rng_uniform()
{
	rng_check();
	return gsl_rng_uniform(random_generator);
}

// return a double from [a,b)
double rng_uniform(double const& a, double const& b)
{
	rng_check();
	return (b-a)*gsl_rng_uniform(random_generator) + a;
}

// Return box-muller gaussian random number
double rng_gaussian(double const& sigma)
{
	rng_check();
	return gsl_ran_gaussian(random_generator, sigma);
}
