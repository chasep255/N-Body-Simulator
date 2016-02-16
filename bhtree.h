#ifndef _BHTREE_H_
#define _BHTREE_H_

#include "physics.h"
#include <stddef.h>

typedef struct bhtree_node
{
	char is_object;
	double total_mass;
	vec3d_t center_of_mass;
	
	struct bhtree_node* fne;
	struct bhtree_node* fnw;
	struct bhtree_node* fsw;
	struct bhtree_node* fse;
	struct bhtree_node* bne;
	struct bhtree_node* bnw;
	struct bhtree_node* bsw;
	struct bhtree_node* bse;
} bhtree_node_t;

typedef struct bhtree_object
{
	char is_object;
	union
	{	
		struct
		{
			double total_mass;
			vec3d_t center_of_mass;
			vec3d_t velocity;
		};
		object_t obj;
	};
} bhtree_object_t;

typedef struct
{
	bhtree_node_t* root;
	double root_size;
	size_t count;
	
	bhtree_node_t* container_pool;
	bhtree_object_t* object_pool;
	size_t cpool_used;
	size_t opool_used;
} bhtree_t;

bhtree_t bhtree(size_t max_n, double initial_size);
bhtree_node_t* bhtree_push(bhtree_t* tree, object_t obj);
void bhtree_update_com(bhtree_t* tree);
void bhtree_free(bhtree_t* tree);
vec3d_t bhtree_calc_gravity_vector(bhtree_t* tree, bhtree_object_t* node, double threshold_angle);
void bhtree_time_step(bhtree_t* tree, double dt, double threshold_angle);

#endif
