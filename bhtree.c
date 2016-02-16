#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <x86intrin.h>
#include "bhtree.h"

bhtree_node_t* bhtree_push_recursive(bhtree_t* tree, bhtree_node_t* restrict root, bhtree_node_t* restrict node, vec3d_t corner, double dim);
static void bhtree_update_com_recursive(bhtree_node_t* node);
static vec3d_t bhtree_calc_gravity_vector_recursive(bhtree_node_t* current_node, double dim);
static void bhtree_populate_list_recursive(bhtree_node_t* node, bhtree_node_t** list, size_t* size);

static __thread struct
{
	double threshold_tangent;
	bhtree_object_t* object_node;
} calc_grav_params;

static __thread char duplicate_insert_flag;

bhtree_node_t* bhtree_push_recursive(bhtree_t* tree, bhtree_node_t* restrict root, bhtree_node_t* restrict node, vec3d_t corner, double dim)
{
	if(root == 0)
		return node;
	
	double dim_2 = dim * 0.5;
	vec3d_t center = vec3d_add(corner, vec3d(dim_2, dim_2, dim_2));
	
	if(root->is_object)
	{
		assert(root->center_of_mass.x >= corner.x);
		assert(root->center_of_mass.y >= corner.y);
		assert(root->center_of_mass.z >= corner.z);
		assert(root->center_of_mass.x <= corner.x + dim);
		assert(root->center_of_mass.y <= corner.y + dim);
		assert(root->center_of_mass.z <= corner.z + dim);
		
		//probably won't happen but would crash program if not handled
		if(root->center_of_mass.x == node->center_of_mass.x &&
				root->center_of_mass.y == node->center_of_mass.y &&
				root->center_of_mass.z == node->center_of_mass.z)
		{
			duplicate_insert_flag = 1;
			root->total_mass += node->total_mass;
			return root;
		}
		
		bhtree_node_t* container = tree->container_pool + tree->cpool_used++;
		memset(container, 0, sizeof(bhtree_node_t));
		
		container = bhtree_push_recursive(tree, container, root, corner, dim);
		container = bhtree_push_recursive(tree, container, node, corner, dim);
		
		return container;
	}

	bhtree_node_t** quad;
	vec3d_t coff;
	vec3d_t diff = vec3d_sub(node->center_of_mass, center);
	
	if(diff.x < 0)
	{
		if(diff.y < 0)
		{
			if(diff.z < 0)
			{
				quad = &(root->bsw);
				coff = vec3d(0, 0, 0);
			}
			else
			{
				quad = &(root->fsw);
				coff = vec3d(0, 0, dim_2);
			}
		}
		else
		{
			if(diff.z < 0)
			{
				quad = &(root->bnw);
				coff = vec3d(0, dim_2, 0);
			}
			else
			{
				quad = &(root->fnw);
				coff = vec3d(0, dim_2, dim_2);
			}
		}
	}
	else
	{
		if(diff.y < 0)
		{
			if(diff.z < 0)
			{
				quad = &(root->bse);
				coff = vec3d(dim_2, 0, 0);
			}
			else
			{
				quad = &(root->fse);
				coff = vec3d(dim_2, 0, dim_2);
			}
		}
		else
		{
			if(diff.z < 0)
			{
				quad = &(root->bne);
				coff = vec3d(dim_2, dim_2, 0);
			}
			else
			{
				quad = &(root->fne);
				coff = vec3d(dim_2, dim_2, dim_2);
			}
		}
	}
	
	vec3d_t next_corner = vec3d_add(corner, coff);
	*quad = bhtree_push_recursive(tree, *quad, node, next_corner, dim_2);
	return root;
}

bhtree_t bhtree(size_t max_n, double initial_size)
{
	bhtree_t tr;
	tr.root = 0;
	tr.count = 0;
	tr.root_size = initial_size;
	tr.object_pool = malloc(sizeof(bhtree_object_t) * max_n);
	tr.container_pool = malloc(sizeof(bhtree_node_t) * max_n);
	tr.opool_used = 0;
	tr.cpool_used = 0;
	return tr;
}

bhtree_node_t* bhtree_push(bhtree_t* tree, object_t obj)
{
	double dim_2 = tree->root_size * 0.5;
	if(fabs(obj.x) >= dim_2 || fabs(obj.y) >= dim_2 || fabs(obj.z) >= dim_2)
	{
		return 0;
	}
	
	bhtree_object_t* node = tree->object_pool + tree->opool_used++;
	memset(node, 0, sizeof(bhtree_object_t));
	node->center_of_mass = obj.position;
	node->total_mass = obj.mass;
	node->velocity = obj.velocity;
	node->is_object = 1;
	duplicate_insert_flag = 0;
	tree->root = bhtree_push_recursive(tree, tree->root, (bhtree_node_t*)node, vec3d(-dim_2, -dim_2, -dim_2), tree->root_size);
	if(!duplicate_insert_flag) 
		tree->count++;
	return (bhtree_node_t*)node;
}

void bhtree_update_com_recursive(bhtree_node_t* node)
{
	if(node == 0 || node->is_object)
		return;
	
	vec3d_t com = {0};
	double m = 0;
	
	#define _update_quad(q)\
	if(q)\
	{\
		bhtree_update_com_recursive(q);\
		com = vec3d_add(com, vec3d_scaler_mul(q->center_of_mass, q->total_mass));\
		m += q->total_mass;\
	}
	
	_update_quad(node->bne);
	_update_quad(node->bnw);
	_update_quad(node->bse);
	_update_quad(node->bsw);
	_update_quad(node->fne);
	_update_quad(node->fnw);
	_update_quad(node->fse);
	_update_quad(node->fsw);
	
	#undef _update_quad
	
	node->total_mass = m;
	node->center_of_mass = m ? vec3d_scaler_mul(com, 1.0 / m) : vec3d(0, 0, 0);
}

void bhtree_update_com(bhtree_t* tree)
{
	bhtree_update_com_recursive(tree->root);
}

void bhtree_free(bhtree_t* tree)
{
	free(tree->object_pool);
	free(tree->container_pool);
	memset(tree, 0, sizeof(bhtree_t));
}

vec3d_t bhtree_calc_gravity_vector_avx_recursive(bhtree_node_t* node, double dim)
{
	
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wmaybe-uninitialized" 
	const __m256d ox = _mm256_set1_pd(calc_grav_params.object_node->center_of_mass.x);
	const __m256d oy = _mm256_set1_pd(calc_grav_params.object_node->center_of_mass.y);
	const __m256d oz = _mm256_set1_pd(calc_grav_params.object_node->center_of_mass.z);
	const __m256d g = _mm256_set1_pd(physics_get_gravitational_constant());
	const __m256d threshold = _mm256_set1_pd(calc_grav_params.threshold_tangent);
	const __m256d epislon = _mm256_set1_pd(physics_get_softening());
	
	__m256d bx, by, bz, bm, bdim;
	__m256d fx, fy, fz, fm, fdim;
	
	double dim_2 = dim * 0.5;
	
	#define _load(q, i, _x, _y, _z, _d, _m)\
	if(q)\
	{\
		_x[i] = q->center_of_mass.x;\
		_y[i] = q->center_of_mass.y;\
		_z[i] = q->center_of_mass.z;\
		_m[i] = q->total_mass;\
		_d[i] = q->is_object ? 0 : dim_2;\
	}\
	else\
	{\
		_x[i] = 0;\
		_y[i] = 0;\
		_z[i] = 0;\
		_m[i] = 0;\
		_d[i] = 0;\
	}
	
	_load(node->bne, 0, bx, by, bz, bdim, bm);
	_load(node->bnw, 1, bx, by, bz, bdim, bm);
	_load(node->bse, 2, bx, by, bz, bdim, bm);
	_load(node->bsw, 3, bx, by, bz, bdim, bm);
	_load(node->fne, 0, fx, fy, fz, fdim, fm);
	_load(node->fnw, 1, fx, fy, fz, fdim, fm);
	_load(node->fse, 2, fx, fy, fz, fdim, fm);
	_load(node->fsw, 3, fx, fy, fz, fdim, fm);
	
	#undef _load
	
	__m256d bdx = _mm256_sub_pd(bx, ox);
	__m256d bdy = _mm256_sub_pd(by, oy);
	__m256d bdz = _mm256_sub_pd(bz, oz);
	
	__m256d bdx2 = _mm256_mul_pd(bdx, bdx);
	__m256d bdy2 = _mm256_mul_pd(bdy, bdy);
	__m256d bdz2 = _mm256_mul_pd(bdz, bdz);
	
	__m256d bd2 = _mm256_add_pd(bdx2, _mm256_add_pd(bdy2, bdz2));
	__m256d bd = _mm256_sqrt_pd(bd2);
	
	__m256d fdx = _mm256_sub_pd(fx, ox);
	__m256d fdy = _mm256_sub_pd(fy, oy);
	__m256d fdz = _mm256_sub_pd(fz, oz);
	
	__m256d fdx2 = _mm256_mul_pd(fdx, fdx);
	__m256d fdy2 = _mm256_mul_pd(fdy, fdy);
	__m256d fdz2 = _mm256_mul_pd(fdz, fdz);
	
	__m256d fd2 = _mm256_add_pd(fdx2, _mm256_add_pd(fdy2, fdz2));
	__m256d fd = _mm256_sqrt_pd(fd2);
	
	__m256d ftan = _mm256_div_pd(fdim, fd);
	__m256d btan = _mm256_div_pd(bdim, bd);
	
	__m256d bcmp = _mm256_cmp_pd(btan, threshold, _CMP_LE_OQ);
	__m256d fcmp = _mm256_cmp_pd(ftan, threshold, _CMP_LE_OQ);
	
	if(!_mm256_testz_pd(fcmp, _mm256_set1_pd(0.0)))
		goto fskip;
	
	__m256d fgm = _mm256_mul_pd(g, fm);
	__m256d fr2sqrtr2e = _mm256_mul_pd(fd2, _mm256_sqrt_pd(_mm256_add_pd(fd2, epislon)));
	__m256d fa = _mm256_div_pd(fgm, fr2sqrtr2e);
	
	__m256d fax = _mm256_mul_pd(fa, fdx);
	__m256d fay = _mm256_mul_pd(fa, fdy);
	__m256d faz = _mm256_mul_pd(fa, fdz);
	
	fskip:
	
	if(!_mm256_testz_pd(bcmp, _mm256_set1_pd(0.0)))
		goto bskip;
	
	__m256d bgm = _mm256_mul_pd(g, bm);
	__m256d br2sqrtr2e = _mm256_mul_pd(bd2, _mm256_sqrt_pd(_mm256_add_pd(bd2, epislon)));
	__m256d ba = _mm256_div_pd(bgm, br2sqrtr2e);
	
	__m256d bax = _mm256_mul_pd(ba, bdx);
	__m256d bay = _mm256_mul_pd(ba, bdy);
	__m256d baz = _mm256_mul_pd(ba, bdz);
	
	
	vec3d_t a = {0};
	
	#define _check_quad(q, i, _ax, _ay, _az, _tan, _cmp)\
	if(q && (void*)q != (void*)calc_grav_params.object_node)\
	{\
		if(_cmp[i])\
			a = vec3d_add(a, vec3d(_ax[i], _ay[i], _az[i]));\
		else\
			a = vec3d_add(a, bhtree_calc_gravity_vector_avx_recursive(q, dim_2));\
	}
	
	bskip:
	
	_check_quad(node->bne, 0, bax, bay, baz, btan, bcmp);
	_check_quad(node->bnw, 1, bax, bay, baz, btan, bcmp);
	_check_quad(node->bse, 2, bax, bay, baz, btan, bcmp);
	_check_quad(node->bsw, 3, bax, bay, baz, btan, bcmp);
	_check_quad(node->fne, 0, fax, fay, faz, ftan, fcmp);
	_check_quad(node->fnw, 1, fax, fay, faz, ftan, fcmp);
	_check_quad(node->fse, 2, fax, fay, faz, ftan, fcmp);
	_check_quad(node->fsw, 3, fax, fay, faz, ftan, fcmp);
	
	#undef _check_quad
	
	return a;
	#pragma GCC diagnostic pop
}

vec3d_t bhtree_calc_gravity_vector_recursive(bhtree_node_t* current_node, double dim)
{
	const vec3d_t ZERO = {0};
	
	double dist = vec3d_dist(calc_grav_params.object_node->center_of_mass, current_node->center_of_mass);
	double tangent = dim / dist;
	if(tangent < calc_grav_params.threshold_tangent)
	{
		return physics_gravity_acceleration_vector(calc_grav_params.object_node->obj, 
				object(current_node->center_of_mass, ZERO, current_node->total_mass));
	}
	
	double dim_2 = dim * 0.5;
	
	#define _check_quad(q, n)\
	vec3d_t n = {0};\
	if(q && (void*)q != (void*)calc_grav_params.object_node)\
	{\
		if(q->is_object)\
			n = physics_gravity_acceleration_vector(calc_grav_params.object_node->obj,\
				((bhtree_object_t*)q)->obj);\
		else\
			n = bhtree_calc_gravity_vector_recursive(q, dim_2);\
	}
	
	_check_quad(current_node->bne, a1);
	_check_quad(current_node->bnw, a2);
	_check_quad(current_node->bse, a3);
	_check_quad(current_node->bsw, a4);
	_check_quad(current_node->fne, a5);
	_check_quad(current_node->fnw, a6);
	_check_quad(current_node->fse, a7);
	_check_quad(current_node->fsw, a8);
	
	#undef _check_quad
	
	a1 = vec3d_add(a1, a2);
	a3 = vec3d_add(a3, a4);
	a5 = vec3d_add(a5, a6);
	a7 = vec3d_add(a7, a8);
	
	a1 = vec3d_add(a1, a3);
	a5 = vec3d_add(a5, a7);
	
	return vec3d_add(a1, a5);
}

vec3d_t bhtree_calc_gravity_vector(bhtree_t* tree, bhtree_object_t* node, double threshold_angle)
{
	calc_grav_params.threshold_tangent = tan(threshold_angle);
	calc_grav_params.object_node = node;
	return bhtree_calc_gravity_vector_recursive(tree->root, tree->root_size);
}

void bhtree_populate_list_recursive(bhtree_node_t* node, bhtree_node_t** list, size_t* size)
{
	if(node == 0)
		return;
	else if(node->is_object)
		list[(*size)++] = node;
	else
	{
		bhtree_populate_list_recursive(node->bne, list, size);
		bhtree_populate_list_recursive(node->bnw, list, size);
		bhtree_populate_list_recursive(node->bse, list, size);
		bhtree_populate_list_recursive(node->bsw, list, size);
		bhtree_populate_list_recursive(node->fne, list, size);
		bhtree_populate_list_recursive(node->fnw, list, size);
		bhtree_populate_list_recursive(node->fse, list, size);
		bhtree_populate_list_recursive(node->fsw, list, size);
	}
}

void bhtree_time_step(bhtree_t* tree, double dt, double threshold_angle)
{	
	bhtree_update_com(tree);
	object_t* objs = malloc(sizeof(object_t) * tree->count);
	bhtree_node_t** nodes =  malloc(sizeof(bhtree_node_t*) * tree->count);
	size_t size = 0;
	bhtree_populate_list_recursive(tree->root, nodes, &size);
	
	size = 0;
	#pragma omp parallel for schedule(dynamic)
	for(size_t i = 0; i < tree->count; i++)
	{
		vec3d_t a = bhtree_calc_gravity_vector(tree, (bhtree_object_t*)nodes[i], threshold_angle);
		object_t obj = object_step(((bhtree_object_t*)nodes[i])->obj, a, dt);
		#pragma omp critical
			objs[size++] = obj;
	}
	
	double total_mass = 0.0;
	double x = 0.0, y = 0.0, z = 0.0;
	
	#pragma omp parallel for reduction(+ : total_mass, x, y, z)
	for(size_t i = 0; i < tree->count; i++)
	{
		total_mass += objs[i].mass;
		x += objs[i].mass * objs[i].position.x;
		y += objs[i].mass * objs[i].position.y;
		z += objs[i].mass * objs[i].position.z;
	}

	x /= total_mass;
	y /= total_mass;
	z /= total_mass;
	
	double max_side = 0;
	#pragma omp parallel for reduction(max : max_side)
	for(size_t i = 0; i < tree->count; i++)
	{
		objs[i].x -= x;
		objs[i].y -= y;
		objs[i].z -= z;
		
//		max_side = fmax(max_side, fabs(objs[i].x));
//		max_side = fmax(max_side, fabs(objs[i].y));
//		max_side = fmax(max_side, fabs(objs[i].z));
	}
	
	size_t old_count = tree->count;
	tree->cpool_used = 0;
	tree->opool_used = 0;
	//tree->root_size = max_side * 2.01;
	tree->count = 0;
	tree->root = 0;
	
	for(size_t i = 0; i < old_count; i++)
		bhtree_push(tree, objs[i]);
	
	free(nodes);
	free(objs);
}
