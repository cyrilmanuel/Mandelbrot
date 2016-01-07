#define _GNU_SOURCE
#include "gfx.h"
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

// Size of the window
#define WIDTH 1280
#define HEIGHT 960

int leave = 0; // Global int for exit master thread //

// Coordinates and size of the window within the Mandelbrot plane.
// For more examples of coordinates, check out http://www.cuug.ab.ca/dewara/mandelbrot/images.html

/**
 *   Struct of the parameters of the map
 *   @xc x position in the complex plane
 *   @yc y position in the complex plane
 *   @size size of the window in the complex plane
 *   @max_iter maximum number of iterations
 *   @dcol color increment (> 0)
 */
struct Params_st {
    double xc;
    double yc;
    double size;
    long max_iter;
    double dcol;
};
typedef struct Params_st Params_t;

/**
 *   Struct of the color_map
 *   @map 32 bit value of map
 *   @length int length of the map
 */
struct Colormap_st {
    uint32 *map;
    int length;
};
typedef struct Colormap_st Colormap_t;

/**
 *   Struct of paramateres for a thread worker
 *   @surface surface of working
 *   @colmap color_map
 *   @width_min Witdh_min for the working thread
 *   @width_max Witdh_max for the working thread
 *   @p p pointers on Params_t
 */
struct Param_worker {
    int blocks;
    double ratio;
    Params_t *p;
    SURFACE *surface;
    Colormap_t *colmap;
};
typedef struct Param_worker Param_worker;

typedef struct Bloc {
    int bloc_id;
    int width_min;
    int width_max;
} Bloc;

/** Global !! */

Bloc* blocs = NULL;
int shared_index = 0;

pthread_mutex_t lock;

int getIndex(int blocks) {


    pthread_mutex_lock(&lock);
    int tmp = shared_index;
    if (shared_index < blocks) {
	printf("%d -- ", shared_index);
	shared_index++;
	pthread_mutex_unlock(&lock);
        return tmp;
    }

    pthread_mutex_unlock(&lock);
    return tmp;
}

/*********************************************************/

/**
 *   Struct wich will be pass in parameter in the main thread
 *   @workers number of workers
 *   @blocks numbers of blocks
 *   @width_img width of the img (WIDTH)
 *   @width_blocks widths for each blocks
 */
struct Param_master {
    int workers;
    int blocks;
    int width_img;
    Params_t *p;
    SURFACE *surface;
    Colormap_t *colmap;
};
typedef struct Param_master Param_master;

/**
 * Create the colormap.
 * This tool can be used to visualize gradients: http://angrytools.com/gradient/
 */
void create_colormap(Colormap_t *colmap) {
    int shade_count = 256;
    int step_count = 3;
    colmap->length = shade_count * step_count;
    colmap->map = (uint32 *) malloc(sizeof (uint32) * colmap->length);
    uint32 steps[] = {COLOR(6, 88, 189), COLOR(6, 214, 100), COLOR(255, 99, 133)};
    int c = 0;
    int j;
    for ( j = 0; j < step_count; j++) {
        double dr = ((double) COLOR_GET_R(steps[(j + 1) % step_count]) - (double) COLOR_GET_R(steps[j])) / (double) shade_count;
        double dg = ((double) COLOR_GET_G(steps[(j + 1) % step_count]) - (double) COLOR_GET_G(steps[j])) / (double) shade_count;
        double db = ((double) COLOR_GET_B(steps[(j + 1) % step_count]) - (double) COLOR_GET_B(steps[j])) / (double) shade_count;
        int i;
        for ( i = 0; i < shade_count; i++) {
            uint8 r = (uint8) ((double) COLOR_GET_R(steps[j]) + dr * (double) i);
            uint8 g = (uint8) ((double) COLOR_GET_G(steps[j]) + dg * (double) i);
            uint8 b = (uint8) ((double) COLOR_GET_B(steps[j]) + db * (double) i);
            colmap->map[c++] = COLOR(r, g, b);
        }
    }
}

/**
 * Free the memory used by the color map.
 */
void free_colormap(Colormap_t *colmap) {
    free(colmap->map);
    colmap->map = 0;
}

/**
 * Exit when ESC is pressed
 */
void* exit_esc(void* arg) {
    leave = 0;
    while (!gfx_is_esc_pressed()) {
        usleep(100000); // Check every 0.1 sec. (10 Hz)
    }
    leave = 1;
    pthread_exit(EXIT_SUCCESS);
}

/**
 * Refresh at 25hz
 */
void* present_25hz(void* arg) {
    while (1) {
        SURFACE *surface = (SURFACE*) arg;
        usleep(1); // Check every 0.4 sec. (25 Hz)
        gfx_present(surface);
    }
}

/**
 * Render the Mandelbrot set.
 */
void* mandelbrot(void *arg) {

    // Struct will be pass in argument in the function
    Param_worker *data_mand = (Param_worker*) arg;

    // Calculate of a the width min and max for x
    double x1 = data_mand->p->xc - data_mand->p->size;
    double x2 = data_mand->p->xc + data_mand->p->size;

    // Calculate of a the width min and max for y
    double y1 = data_mand->p->yc - data_mand->p->size;
    double y2 = data_mand->p->yc + data_mand->p->size;

    // Calculate the dx and dy for Mandelbrot
    double dx = (x2 - x1) / WIDTH;
    double dy = (y2 - y1) / HEIGHT;

    // Define y1 in block_data

    Bloc currentBloc;
    int current_element = getIndex(data_mand->blocks);

    while ((current_element < data_mand->blocks) && leave!=1) {
	double y = y1;

        currentBloc = blocs[current_element];


        // Ratio of the numbers of workers

        /**
         *	Loop who starts the calculation of mandelbrot set
         */
         int h;
        for ( h = 0; h < HEIGHT; h++) {

            // Defined by the number of blocks and the id of thread
            double x = currentBloc.bloc_id * dx * data_mand->ratio * WIDTH + x1;
            int r;
            for ( r = currentBloc.width_min; r < currentBloc.width_max; r++) {

                double zx = 0;
                double zy = 0;
                uint32 color = COLOR(0, 0, 0);
                long depth;
                for ( depth = 0; depth < data_mand->p->max_iter; depth++) {
                    double zx_new = (zx * zx) - (zy * zy) + x;
                    double zy_new = 2.0 * zx * zy + y;
                    zx = zx_new;
                    zy = zy_new;
                    // Did the pixel diverge (go to infinity)?
                    if ((zx * zx + zy * zy) > 4.0) {
                        // Set the color for the mandelbrot function
                        color = data_mand->colmap->map[((int) ((double) depth * data_mand->p->dcol)) % data_mand->colmap->length];
                        break;
                    }
                }
                gfx_setpix(data_mand->surface, r, h, color); // Set the surface and color of SDL window
                x += dx;
                //if (leave) return 0; // Leave the programm
            }
            y += dy;
        }

        current_element = getIndex(data_mand->blocks);
    }

    return 0;
}

/**
 *	Master function which will create the workers threads and launch the threads
 */
void* master_func(void *arg) {

    // Struct will be pass in argument in the function
    Param_master *data = (Param_master*) arg;

    blocs = (Bloc*) malloc(data->blocks * sizeof (Bloc));
    shared_index = 0;
    double ratio = 1 / (double) data->blocks;

    // Tab creation of data_worker base on a struct
    Param_worker data_worker;
    data_worker.p = data->p; // we send his Params_st
    data_worker.colmap = data->colmap; // same with colmap
    data_worker.surface = data->surface; // and the surface
    data_worker.ratio = ratio;
    data_worker.blocks = data->blocks;

    int add = 0;
    int i = 0;

    // Calcul width for each blocks
    int width = data->width_img / data->blocks;

    // We intitialize the first and the last block apart
    blocs[0].width_min = add; // with the start of his block
    blocs[0].width_max = width; // and the end
    blocs[0].bloc_id = 0;

    // We initialize others blocks
    for (i = 1; i < (data->blocks) - 1; i++) {
        blocs[i].width_min = blocs[i - 1].width_min + width;
        blocs[i].width_max = blocs[i - 1].width_max + width;
        blocs[i].bloc_id = i;
    }

    // We initialize the last block and give him everything left
    blocs[i].width_min = blocs[i - 1].width_min + width;
    blocs[i].width_max = data->width_img - blocs[i].width_max;
    blocs[i].bloc_id = i;

    // Memory for workers' thread
    pthread_t* thread_worker = (pthread_t*) malloc(data->workers * sizeof (pthread_t));



    // Creation of thread workers
    for ( i = 0; i < data->workers; i++) {
        if (pthread_create(&thread_worker[i], NULL, mandelbrot, &data_worker) != 0) {
            fprintf(stderr, "Erreur dans la creation du thread\n");
            return NULL;
        }
    }
    // Laucnh the thread workers
    for ( i = 0; i < data->workers; i++) {
        if (pthread_join(thread_worker[i], NULL) != 0) {
            perror("Join error");
            return 0;
        }
    }

    free(thread_worker);

    return 0;
}

/*
 *
 */
int main(int argc, char** argv) {
    // We set Param_master for initializing some int
    Param_master *param_master = malloc(sizeof (Param_master));

    pthread_t escape; // Thread ESC
    pthread_t present; // Thread that will refresh
    pthread_t master; // Master thread

    param_master->width_img = WIDTH; // Define the width of the img
    param_master->workers =  1; // Define the number of workers
    param_master->blocks = 50; // Define the number of blocks

    // Creation of the thread which is checking for the ESC
    if (pthread_create(&escape, NULL, exit_esc, NULL) != 0) {
        fprintf(stderr, "Erreur dans la creation du thread\n");
        return EXIT_FAILURE;
    }

    Colormap_t colmap; // Set colormap_t
    create_colormap(&colmap); // Create colormap_t

    // Rendering surface
    SURFACE *surface = gfx_init("Mandelbrot", WIDTH, HEIGHT);

    // If surface=NULL then we print "Failed initializing video mode"
    if (surface == NULL) {
        fprintf(stderr, "Failed initializing video mode!\n");
        return EXIT_FAILURE;
    }

    // Creation of the thread which defines the Hz display
    if (pthread_create(&present, NULL, present_25hz, surface) != 0) {
        fprintf(stderr, "Erreur dans la creation du thread\n");
        return EXIT_FAILURE;
    }
    // Mandelbrot computation parameters
    Params_t p = {
        0.2929859127507,
        0.6117848324958,
        1.0E-11,
        8000,
        0.9
    };

    // Longer computation
    Params_t q = {
            -0.17476469999956,
            -1.0713151000007,
            5.095053e-13,
            8000,
            0.35 };
    // Classic coordinates
    Params_t r = {
            -0.65,
            -0.0,
            1.2,
            150,
            10 };

    param_master->p = &q; // We take Params_st
    param_master->colmap = &colmap; // We take colmap
    param_master->surface = surface; // We take surface

    // Creation of thread master
    if (pthread_create(&master, NULL, master_func, param_master) != 0) {
        fprintf(stderr, "Erreur dans la creation du thread\n");
        return EXIT_FAILURE;
    }

    // Execute thread master
    if (pthread_join(master, NULL) != 0) {
        perror("Join error");
        return 0;
    }


    // Free memory
    free(blocs);
    while(leave != 1){
    	usleep(10);
    }
    gfx_present(surface);
    system("PAUSE");
    gfx_close();
    free_colormap(&colmap);
    return 0;
}
