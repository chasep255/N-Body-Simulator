#ifndef _PHYSICS_H_
#define _PHYSICS_H_

#include <math.h>
#include "vec3d.h"

typedef struct
{
	double mass;
	union
	{
		struct
		{
			double x, y, z;
		};
		vec3d_t position;
	};
	
	union
	{
		struct
		{
			double vx, vy, vz;
		}__attribute__((packed));
		vec3d_t velocity;
	};
} object_t;

object_t object(vec3d_t x, vec3d_t v, double m);
object_t object_step(object_t o, vec3d_t a, double dt);

void physics_set_gravitational_constant(double g);
double physics_get_gravitational_constant();

void physics_set_light_speed(double c);
double physics_get_light_speed();

void physics_set_softening(double e);
double physics_get_softening();

double physics_gamma(double v);
vec3d_t physics_gravity_acceleration_vector(object_t a, object_t b);
#endif
