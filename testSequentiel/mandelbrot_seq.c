/**
 * @file mandelbrot_seq.c
 * @brief Simple Mandelbrot fractal renderer.
 * @author FG
 * @date November 2014
 * @version 0.1
 */

#define _GNU_SOURCE
#include "gfx.h"
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <semaphore.h>

// Size of the window
#define WIDTH 1280
#define HEIGHT 960

// Coordinates and size of the window within the Mandelbrot plane.
// For more examples of coordinates, check out http://www.cuug.ab.ca/dewara/mandelbrot/images.html
struct params_st {
	double xc;	// x position in the complex plane
	double yc;	// y position in the complex plane
	double size;	// size of the window in the complex plane
	long max_iter;	// maximum number of iterations
	double dcol;	// color increment (> 0)
};
typedef struct params_st params_t;

// Definition of a color map
struct colormap_st {
	uint32 *map;
	int length;
};
typedef struct colormap_st colormap_t;

struct Param_master {
    int workers;
    int Nbblocks;
    params_t *p;
    SURFACE *surface;
    colormap_t *colmap;
};
typedef struct Param_master Param_master;

struct Param_worker {
    int NbBloc;
    params_t *p;
    SURFACE *surface;
    colormap_t *colmap;
};
typedef struct Param_worker Param_worker;

int  CountBlock= 0;

pthread_mutex_t lock;
sem_t mutex;

int getIndex(int blocks) {
    sem_wait(&mutex);
    int tmp = CountBlock;
    if (CountBlock < blocks) {
	printf("%d -- ", CountBlock);
	CountBlock++;
    sem_post(&mutex);
        return tmp;
    }
    sem_post(&mutex);
    return tmp;
}
/**
 * Create the colormap.
 * This tool can be used to visualize gradients: http://angrytools.com/gradient/
 */
void create_colormap(colormap_t *colmap) {
	int shade_count = 256;
	int step_count = 3;
	colmap->length = shade_count * step_count;
	colmap->map = (uint32 *) malloc(sizeof(uint32) * colmap->length);
	uint32 steps[] = { COLOR(6,88,189), COLOR(6,214,100), COLOR(255,99,133) };
	int c = 0;
	int j;
	for (j = 0; j < step_count; j++) {
		double dr = ((double)COLOR_GET_R(steps[(j+1) % step_count]) - (double)COLOR_GET_R(steps[j]))/(double)shade_count;
		double dg = ((double)COLOR_GET_G(steps[(j+1) % step_count]) - (double)COLOR_GET_G(steps[j]))/(double)shade_count;
		double db = ((double)COLOR_GET_B(steps[(j+1) % step_count]) - (double)COLOR_GET_B(steps[j]))/(double)shade_count;
		int i;
		for (i = 0; i < shade_count; i++) {
			uint8 r = (uint8)((double)COLOR_GET_R(steps[j]) + dr * (double)i);
			uint8 g = (uint8)((double)COLOR_GET_G(steps[j]) + dg * (double)i);
			uint8 b = (uint8)((double)COLOR_GET_B(steps[j]) + db * (double)i);
			colmap->map[c++] = COLOR(r,g,b);
		}
	}
}

/**
 * Free the memory used by the color map.
 */
void free_colormap(colormap_t *colmap) {
	free(colmap->map);
	colmap->map = 0;
}

/**
 * Render the Mandelbrot set.
 */
void *mandelbrot(void *arg) {

    // création du objet type param_worker qui contiendra les infos contenue dans le param master
    Param_worker *data_mand = (Param_worker*) arg;
    double x1 = data_mand->p->xc - data_mand->p->size;
    double x2 = data_mand->p->xc + data_mand->p->size;
    double y1 = data_mand->p->yc - data_mand->p->size;
    double y2 = data_mand->p->yc + data_mand->p->size;

	double dx = (x2 - x1) / WIDTH;
	double dy = (y2 - y1) / HEIGHT;

    int blocID = getIndex(data_mand->NbBloc);

    while ((blocID < data_mand->NbBloc)) {
	double y = y1;

    int i;
    double ratio = 1.0 / data_mand->NbBloc;

    double tailleBloc=WIDTH/data_mand->NbBloc;
	for ( i = 0; i < HEIGHT; i++) {
		double x =blocID*ratio*dx*WIDTH+ x1;
        int j;
        int a;
        if(blocID==data_mand->NbBloc-1)
        {
            //Si on a atteind le dernier bloc a calculer on stocke WIDTH dans a
            a=WIDTH;
        }
        else
        {
            //Sinon a est égal blocID*tailleBloc+tailleBloc
            a=blocID*tailleBloc+tailleBloc;
        }
        // on initialise j à la valeur de départ du bloc
		for ( j = blocID*tailleBloc; j <a ; j++) {
			double zx = 0;
			double zy = 0;
			uint32 color = COLOR(0,0,0);
            long depth;
			for ( depth = 0; depth < data_mand->p->max_iter; depth++) {
				double zx_new = (zx*zx) - (zy*zy) + x;
				double zy_new = 2.0*zx*zy + y;
				zx = zx_new;
				zy = zy_new;
				// Did the pixel diverge (go to infinity)?
				if ((zx*zx + zy*zy) > 4.0) {
					color = data_mand->colmap->map[((int)((double)depth*data_mand->p->dcol)) % data_mand->colmap->length];
					break;
				}
			}
			gfx_setpix(data_mand->surface, j, i, color);
			x += dx;
		}
		y += dy;

		// Every 32 lines: present surface to screen and check keyboard
		if(i%32==0)
        {
            gfx_present(data_mand->surface);
			if (gfx_is_esc_pressed()) {
				break;
			}
        }

	}
        blocID = getIndex(data_mand->NbBloc);
    }
}

void* master_func(void *arg) {

    Param_master *data = (Param_master*) arg;

    CountBlock = 0;

    Param_worker data_worker;
    data_worker.p = data->p;// on remplis le parametre Params_t par celui initialiser plus haut (p)
    data_worker.colmap = data->colmap; // on remplis le parametre colmap par celui initialiser plus haut
    data_worker.surface = data->surface; // on remplis le parametre surface par celui initialiser plus haut
    data_worker.NbBloc=data->Nbblocks;  // on remplis le parametre NbBblocks avec le nombre de block

    // allocation de la memoire pour les threads worker
    pthread_t *thread_worker = (pthread_t*) malloc(data->workers * sizeof (pthread_t));
	int i;
    for(i=0;i<data->workers;i++)
    {
        // Creation des threads workers
        pthread_create(&thread_worker[i],NULL,mandelbrot,&data_worker);
    }
    for(i=0;i<data->workers;i++)
    {
        // join des threads workers
        pthread_join(thread_worker[i],NULL);
    }
    //Liberation de la mémoir
    free(thread_worker);

        FILE *stream ;
    if((stream = freopen("file.txt", "w", stdout)) == NULL)
      exit(-1);
    int temps =SDL_GetTicks();
    printf("Il faut %d ms pour calculer le total",temps);
    stream = freopen("CON", "w", stdout);
    printf("Il faut %d ms pour calculer le total",temps);
}
/**
 * Program's entry point.
 * @param argc number of arguments.
 * @param arguments (as an array of strings).
 * @return status code.
 */
int main(int argc, char **argv) {
	colormap_t colmap; // instenciation d'un colormap
    create_colormap(&colmap); // remplissage des attribut du colormap

    freopen("CON", "w", stdout);
    freopen("CON", "w", stderr);

	    // création d'une nouvelle surface avec en parametre la taille définit dans les defines.
	SURFACE *surface = gfx_init("Mandelbrot", WIDTH, HEIGHT);
	if (surface == NULL) {
		fprintf(stderr, "Failed initializing video mode!\n");
		return EXIT_FAILURE;
	}
    // Mandelbrot parameters
	params_t p = {
		0.2929859127507,
		0.6117848324958,
		1.0E-11,
		8000,
		0.9 };

	// Classic coordinates
	params_t r = {
		-0.65,
		-0.0,
		1.2,
		1500,
		10 };



	// Longer computation
	params_t q = {
		-0.17476469999956,
		-1.0713151000007,
		5.095053e-13,
		8000,
		0.35 };



    sem_init(&mutex, 0, 1); // initialisation du semaphore
    pthread_t thread_master; // thread master

    // création d'un param master qui contiendra le nombre de bloc, le nombre de worker
	Param_master *param_master = malloc(sizeof (Param_master));

	param_master->workers= atoi(argv[1]);  // définit le nombre de thread worker
	param_master->Nbblocks=atoi(argv[2]);  // définit le nombre de bloc
	param_master->p = &r;
    param_master->colmap = &colmap;
    param_master->surface = surface;

    pthread_create(&thread_master,NULL,master_func,param_master);
    pthread_join(thread_master,NULL);


    gfx_present(surface);
    system("PAUSE");
	gfx_close();

	free_colormap(&colmap);
    free(param_master);
	return EXIT_SUCCESS;
}
