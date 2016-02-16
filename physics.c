#include "physics.h"

static double c = 1.0;
static double over_c2 = 1.0;
static double g = 1.0;
static double e = 1.0;

object_t object(vec3d_t x, vec3d_t v, double m)
{
	return (object_t){.position = x, .velocity = v, .mass = m};
}

vec3d_t physics_gravity_acceleration_vector(object_t a, object_t b)
{
	vec3d_t r = vec3d_sub(b.position, a.position);
	double r2 = vec3d_mag2(r);
	return vec3d_scaler_mul(r, g * b.mass / (r2 * sqrt(r2 + e)));
}

object_t object_step(object_t o, vec3d_t a, double dt)
{
	vec3d_t v1 = o.velocity;
	
	o.velocity.x += a.x * dt;
	o.velocity.y += a.y * dt;
	o.velocity.z += a.z * dt;
	
	dt *= 0.5;
	o.position.x += (o.velocity.x + v1.x) * dt;
	o.position.y += (o.velocity.y + v1.y) * dt;
	o.position.z += (o.velocity.z + v1.z) * dt;
	return o;
}


void physics_set_gravitational_constant(double _g)
{
	g = _g;
}

double physics_get_gravitational_constant()
{
	return g;
}

void physics_set_light_speed(double _c)
{
	over_c2 = 1.0 / (_c * _c);
	c = _c;
}

double physics_get_light_speed()
{
	return c;
}

void physics_set_softening(double _e)
{
	e = _e;
}

double physics_get_softening()
{
	return e;
}
