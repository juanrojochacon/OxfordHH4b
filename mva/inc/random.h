// random.h nh 22/10/14
// Random number generator interface

#pragma once

// Initialise RNG
void rng_init(unsigned long int const& s);

//return an integer from 0 to n-1
unsigned long int rng_uniform(unsigned long int const& n);

// return a double from [0,1)
double rng_uniform();

// return a double from [a,b)
double rng_uniform(double const& a, double const& b);

// Return box-muller gaussian random number
double rng_gaussian(double const& sigma);