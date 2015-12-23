/**
 * @file gfx.h
 * @brief Minimalistic SDL wrapper providing pixel drawing functions.
 * @author Florent Gluck, florent.gluck@hesge.ch
 * @date October, 2014
 * @version 0.1
 */

#ifndef __GFX__H__
#define __GFX__H__

#include <unistd.h>
#include <time.h>
#include <stdbool.h>
#include <SDL/SDL.h>

// Macros to pack and unpack color components.
#define COLOR(r,g,b) ((uint32)b | ((uint32)g<<8) | ((uint32)r<<16))
#define COLOR_GET_B(col) (col & 0xFF)
#define COLOR_GET_G(col) ((col & 0xFF00)>>8)
#define COLOR_GET_R(col) ((col & 0xFF0000)>>16)

typedef Uint32 uint32;
typedef Uint8 uint8;
typedef SDL_Surface SURFACE;

extern SURFACE *gfx_init(char *title, int width, int height);
extern void gfx_setpix(SURFACE *surface, int x, int y, uint32 col);
extern uint32 gfx_getpix(SURFACE *surface, int x, int y);
extern bool gfx_lock(SURFACE *surface);
extern void gfx_unlock(SURFACE *surface);
extern void gfx_present(SURFACE *surface);
extern bool gfx_is_esc_pressed();
extern void gfx_close();

#endif
