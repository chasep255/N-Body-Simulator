#include "astro_bodies.h"
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <stdio.h>

static double frand()
{
	const double RMAX = 1.0 / RAND_MAX;
	return rand() * RMAX;
}

double sfrand()
{
	return rand() % 2 ? frand() : -frand();
}

double gaussian(double mu, double sigma)
{
	const double epsilon = DBL_MIN;
	const double two_pi = 2.0 * 3.14159265358979323846;

	static double z0, z1;
	static int generate;
	generate = !generate;

	if (!generate)
	   return z1 * sigma + mu;

	double u1, u2;
	do
	{
		u1 = frand();
		u2 = frand();
	}
	while (u1 <= epsilon);

	z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
	return z0 * sigma + mu;
}


void galaxy_create(bhtree_t* tree, int nstars, galaxy_parameters_t params)
{
	double mass_sigma = (params.max_star_mass - params.min_star_mass) / 2.0;
	double rad_meters = params.radius_lyr * LYR_TO_M;
	double rad_sigma = params.density_sigma * rad_meters;
	double z_sigma = 0.05 * rad_sigma;
	
	for(int i = 0; i < nstars; i++)
	{
		double radius = fabs(gaussian(0, rad_sigma));
		double theta = 2 * 3.141592653589 * frand();
		double z = gaussian(0, z_sigma);
		double mass = fabs(params.min_star_mass - gaussian(params.min_star_mass, mass_sigma * fmod(theta, 3.14159265))) + params.min_star_mass;
		
		radius = fmin(2 * rad_meters, radius);
		z = z < 0.0 ? fmax(-3.0 * z_sigma, z) : fmin(3.0 * z_sigma, z);
		
		double sin_theta = sin(theta);
		double cos_theta = cos(theta);
		
		double v = atan(2 * radius / rad_meters) * params.tangential_velocity;
		
		double x = cos_theta * radius;
		double y = sin_theta * radius;
		double vx = -sin_theta * v;
		double vy = cos_theta * v;
		double vz = 0;
		
		object_t obj = object((vec3d_t){x, y, z}, (vec3d_t){vx, vy, vz}, mass);
		bhtree_push(tree, obj);
	}
}
