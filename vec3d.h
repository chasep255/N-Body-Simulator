#ifndef _VEC3D_H_
#define _VEC3D_H_

#include <math.h>

typedef struct
{
	double x, y, z;
}__attribute__((packed)) vec3d_t;

static inline vec3d_t vec3d(double x, double y, double z)
{
	return (vec3d_t){.x = x, .y = y, .z = z};
}

static inline double vec3d_dist2(vec3d_t a, vec3d_t b)
{
	double dx, dy, dz;
	
	dx = a.x - b.x;
	dy = a.y - b.y;
	dz = a.z - b.z;
	
	return dx * dx + dy * dy + dz * dz;
}

static inline double vec3d_dist(vec3d_t a, vec3d_t b)
{	
	return sqrt(vec3d_dist2(a, b));
}

static inline vec3d_t vec3d_scaler_mul(vec3d_t a, double b)
{
	return (vec3d_t){.x = a.x * b, .y = a.y * b, .z = a.z * b};
}

static inline double vec3d_mag2(vec3d_t a)
{
	return a.x * a.x + a.y * a.y + a.z * a.z;
}

static inline double vec3d_mag(vec3d_t a)
{
	return sqrt(vec3d_mag2(a));
}

static inline vec3d_t vec3d_add(vec3d_t a, vec3d_t b)
{
	return (vec3d_t){.x = a.x + b.x, .y = a.y + b.y, .z = a.z + b.z};
}

static inline vec3d_t vec3d_sub(vec3d_t a, vec3d_t b)
{
	return (vec3d_t){.x = a.x - b.x, .y = a.y - b.y, .z = a.z - b.z};
}

static inline vec3d_t vec3d_negate(vec3d_t a)
{
	return (vec3d_t){.x = -a.x, .y = -a.y, .z = -a.z};
}

static inline vec3d_t vec3d_normalize(vec3d_t a)
{
	return vec3d_scaler_mul(a, 1.0 / vec3d_mag(a));
}

#endif
