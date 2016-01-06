/**
 * @file gfx.c
 * @brief Minimalistic SDL wrapper providing pixel drawing functions.
 * @author F.G
 * @date October, 2014
 * @version 0.1
 */

#include "gfx.h"

static int DEPTH = 32;

/**
 * Compute pixel address for coordinates (x,y).
 * @param surface for which to compute the position.
 * @param x coordinate.
 * @param y coordinate.
 * @return pixel address.
 */
static uint32 *pix_addr(SURFACE *surface, int x, int y) {
	return (uint32 *) surface->pixels + y * surface->pitch / (DEPTH / 8) + x;
}

/**
 * Set pixel color at coordinates (x,y).
 * @param surface for which to draw the pixel.
 * @param x coordinate.
 * @param y coordinate.
 * @param col pixel color (RGB).
 */
void gfx_setpix(SURFACE *surface, int x, int y, uint32 col) {
	uint32 *pix = pix_addr(surface, x, y);
	*pix = SDL_MapRGB(surface->format, COLOR_GET_R(col), COLOR_GET_G(col), COLOR_GET_B(col));
}

/**
 * Get pixel color at coordinates (x,y).
 * @param surface for which to read the pixel.
 * @param x coordinate.
 * @param y coordinate.
 * @return pixel color (RGB).
 */
uint32 gfx_getpix(SURFACE *surface, int x, int y) {
	uint32 *pix = pix_addr(surface, x, y);
	uint8 r, g, b;
	SDL_GetRGB(*pix, surface->format, &r, &g, &b);
	return COLOR(r, g, b);
}

/**
 * Initialize a video mode.
 * @param title of the window.
 * @param width of the window.
 * @param height of the window.
 * @return pointer to the initialized surface or 0 if the call failed.
 */
SURFACE *gfx_init(char *title, int width, int height) {
	if (SDL_Init(SDL_INIT_VIDEO) == -1) {
		printf("Unable to initialize SDL!\n");
		return NULL;
	}

	// HACK: software surfaces shouldn't require locking
	SURFACE *surface = SDL_SetVideoMode(width, height, DEPTH, SDL_SWSURFACE);
	if (surface == NULL) {
		printf("Unable to initialize SDL video mode!\n");
		SDL_Quit();
		return NULL;
	}

	SDL_WM_SetCaption(title, 0);

	return surface;
}

/**
 * Check whether the ESC key was pressed or windows close button was clicked.
 * @return true if a key was pressed, false otherwise.
 */
bool gfx_is_esc_pressed() {
	SDL_Event event;
	if (SDL_PollEvent(&event)) {
		// Exit when window's close button is clicked
		if (event.type == SDL_QUIT) {
			return true;
		}
		// Exit when ESC is pressed
		else if (event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_ESCAPE) {
			return true;
		}
	}
	return false;
}

/**
 * Present the surface.
 * @param surface to present.
 */
void gfx_present(SURFACE *surface) {
	SDL_Flip(surface);
}

/**
 * Close the graphic mode.
 */
void gfx_close() {
	SDL_Quit();
}
