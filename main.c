#include "physics.h"
#include "bhtree.h"
#include "astro_bodies.h"
#include <math.h>
#include <GL/glut.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <pthread.h>

double max_mass = 0.0;
double min_mass = 0;
double angle = 0.0;
bhtree_t tree;
void draw_tree(bhtree_node_t* n, double dim, vec3d_t center)
{
	if(n == NULL)
		return;
//	glColor3f(1, 1, 1);
//	glPointSize(3);
//	glPushMatrix();
//	glTranslated(center.x, center.y, center.z);
//	glutWireCube(dim);
//	glPopMatrix();
	if(n->is_object)
	{
		double m = n->total_mass - min_mass;
		glColor3f((max_mass - m) / max_mass, 0.75, m / max_mass);
		glBegin(GL_POINTS);
		{
			glVertex3d(n->center_of_mass.x, n->center_of_mass.y, n->center_of_mass.z);
		}
		glEnd();
	}
	else
	{
		double dim_2 = dim / 2;
		double dim_4 = dim / 4;
		draw_tree(n->bne, dim_2, vec3d_add(center, vec3d(dim_4, dim_4, -dim_4)));
		draw_tree(n->bnw, dim_2, vec3d_add(center, vec3d(-dim_4, dim_4, -dim_4)));
		draw_tree(n->bse, dim_2, vec3d_add(center, vec3d(dim_4, -dim_4, -dim_4)));
		draw_tree(n->bsw, dim_2, vec3d_add(center, vec3d(-dim_4, -dim_4, -dim_4)));
		draw_tree(n->fne, dim_2, vec3d_add(center, vec3d(dim_4, dim_4, dim_4)));
		draw_tree(n->fnw, dim_2, vec3d_add(center, vec3d(-dim_4, dim_4, dim_4)));
		draw_tree(n->fse, dim_2, vec3d_add(center, vec3d(dim_4, -dim_4, dim_4)));
		draw_tree(n->fsw, dim_2, vec3d_add(center, vec3d(-dim_4, -dim_4, dim_4)));
	}
}

double zoom;
void display()
{
	glClearColor(0, 0, 0, 0);
	glClear(GL_COLOR_BUFFER_BIT);
	
	glLoadIdentity();
	
	gluLookAt(tree.root->center_of_mass.x, tree.root->center_of_mass.y, tree.root->center_of_mass.z - zoom, tree.root->center_of_mass.x, tree.root->center_of_mass.y, tree.root->center_of_mass.z, 0, 1, 0);
	glScalef(0.1, 0.1, 0.1);
	glRotated(angle, 1, 0, 0);
	
	glColor4f(1, 1, 1, 0.25);
	draw_tree(tree.root, tree.root_size, vec3d(0, 0, 0));
	
	glBegin(GL_POINTS);
	{
		glColor3f(1.0, 0.0, 0.0);
		glVertex3d(tree.root->center_of_mass.x, tree.root->center_of_mass.y, tree.root->center_of_mass.z);
	}
	glEnd();
	glutSwapBuffers();
}

void reshape(int w, int h)
{
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60, (double)w / (double)h, 1.0, zoom * 1e9);
	glMatrixMode(GL_MODELVIEW);
}

void keyboard(unsigned char key, int x, int y)
{
	if(key == '-')
	{
		zoom *= 1.5;
		glutPostRedisplay();
	}
	else if(key == '=')
	{
		zoom /= 1.5;
		glutPostRedisplay();
	}
	else if(key == 'w')
	{
		angle += 2;
		glutPostRedisplay();
	}
	else if(key == 's')
	{
		angle -= 2;
		glutPostRedisplay();
	}
}
int update = 1;
void* time_step(void* args)
{
	static double tlast = 0.0;
	static long years = 0;
	
	double start = omp_get_wtime();
	bhtree_time_step(&tree, 60.0 * 60.0 * 24.0 * 365 * 500, 25 * 3.14159 / 180.0);
	double end = omp_get_wtime();
	
	years += 500;
	printf("dt = %lf\t fps = %lf \t years = %ld \t stars = %zd\n", end - start, 1.0 / (omp_get_wtime() - tlast), years, tree.count);
	tlast = omp_get_wtime();
	update = 1;
	pthread_exit(0);
}

void idle()
{		
	if(update)
	{
		update = 0;
		glutPostRedisplay();
		pthread_t thread;
		pthread_create(&thread, 0, time_step, 0);
		pthread_detach(thread);
	}
}

int main(int argc, char** argv)
{
	omp_set_num_threads(4);
	
	srand(time(0));
	glutInit(&argc, argv);
	glutInitDisplayMode(GL_DOUBLE);
	
	glutInitWindowSize(1000, 1000);
	glutCreateWindow("N-Body");
	
	glutReshapeFunc(reshape);
	glutDisplayFunc(display);
	glutIdleFunc(idle);
	glutKeyboardFunc(keyboard);
	
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glPointSize(1.0);
	
	physics_set_gravitational_constant(6.67408e-11);
	physics_set_softening(1e9);
	physics_set_light_speed(3.0e8);
	
	tree = bhtree(10000, 100000 * LYR_TO_M);
	galaxy_parameters_t params;
	params.max_star_mass = 1000 * SOLAR_MASS_KG;
	params.min_star_mass = SOLAR_MASS_KG;
	params.tangential_velocity = 1e5;
	params.radius_lyr = 100;
	params.density_sigma = 0.3;
	galaxy_create(&tree, 10000, params);
	
	zoom = params.radius_lyr * LYR_TO_M * 3.0;
	max_mass = params.max_star_mass;
	min_mass = params.min_star_mass;
	glutMainLoop();
	bhtree_free(&tree);
}
