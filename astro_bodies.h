#ifndef _ASTRO_BODIES_H_
#define _ASTRO_BODIES_H_

#include "bhtree.h"

#define LYR_TO_M (9.461e15)
#define SOLAR_MASS_KG (1.98855e30)

typedef struct
{
	double radius_lyr;
	double min_star_mass;
	double max_star_mass;
	double density_sigma;
	double density_mu;
	double tangential_velocity;
} galaxy_parameters_t;

void galaxy_create(bhtree_t* tree, int nstars, galaxy_parameters_t params);

#endif
