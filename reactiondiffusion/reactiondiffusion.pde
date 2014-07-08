/*
	Copyright (C) 2014 Stephen M. Cameron
	Author: Stephen M. Cameron

	This file is part of reactiondiffusion.

	Reactiondiffusion is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 2 of the License, or
	(at your option) any later version.

	Reactiondiffusion is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with reactiondiffusion; if not, write to the Free Software
	Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/* This program is based on the work of Karl Sims.
 * See: http://www.karlsims.com/rd.html
 */

int xdim = 200;
int ydim = 100;
int current = 0;

float[][][] a;
float[][][] b;

float diffusion_rate_a = 1.0;
float diffusion_rate_b = 0.5;
float feed_rate_a = 0.055;
float kill_rate_b = 0.062; 
float delta_t = 1.00;

/* convolution matrix */
float [][] conv;

void add_spot(int c, int x, int y)
{
	int i, j, tx, ty;

	for (i = x - 3; i <= x + 3; i++) {
		tx = i;
		if (tx < 0) {
			tx = 0;
		} else if (tx > xdim - 1) {
			tx = xdim - 1;
		}
		for (j = y - 3; j <= y + 3; j++) {
			ty = j;
			if (ty < 0) {
				ty = 0;
			} else if (ty > ydim - 1) {
				ty = ydim - 1;
			}
			a[c][tx][ty] = 1.0;
			b[c][tx][ty] = 1.0;
		}
	}
}

void initialize(int c)
{
	int x, y;

	for (x = 0; x < xdim; x++) {
		for (y = 0; y < ydim; y++) {
			a[c][x][y] = 1.0;
			b[c][x][y] = 0.0;
		}
	}

	add_spot(c, xdim / 2, ydim / 2);
	add_spot(c, xdim / 3, ydim / 2);
	add_spot(c, 2 * xdim / 3, ydim / 3);
}

void setup()
{

	size(xdim, ydim);
	background(0);

	a = new float[2][xdim][ydim];
	b = new float[2][xdim][ydim];

	conv = new float[3][3];

	conv[0][0] = 0.05;
	conv[1][0] = 0.2;
	conv[2][0] = 0.05;
	conv[0][1] = 0.2;
	conv[1][1] = -1.0;
	conv[2][1] = 0.2;
	conv[0][2] = 0.05;
	conv[1][2] = 0.2;
	conv[2][2] = 0.05;

	initialize(current);
	initialize(current + 1);
}

void draw_current(int current)
{
	float mcolor;
	int x, y;
	int c;

	noSmooth();
	background(0);
	for (x = 0; x < xdim; x++) {
		for (y = 0; y < ydim; y++) {
			mcolor = 255.0 * (a[current][x][y] / (1.0 - b[current][x][y]));
			c = int(mcolor);
			stroke(c);
			point(x, y);
		}
	}
}

float convolve(int from, int x, int y, int use_a)
{
	int ix, iy, tx, ty, sx, sy;
	float sum;

	sum = 0.0;
	for (ix = x - 1; ix <= x + 1; ix++) {
		sx = ix - x + 1;
		for (iy = y - 1; iy <= y + 1; iy++) {
			sy = iy - y + 1;
			tx = ix;
			if (tx < 0) {
				tx = 0;
			} else if (tx > xdim - 1) {
				tx = xdim - 1;
			}
			ty = iy;
			if (ty < 0) {
				ty = 0;
			} else if (ty > ydim -1) {
				ty = ydim - 1;
			}

			/* better way of selecting a or b with java? Pass by ref? */
			if (use_a != 0) {
				sum +=  a[from][tx][ty] * conv[sx][sy];
			} else {
				sum +=  b[from][tx][ty] * conv[sx][sy];
			}
		}
	}
	return sum;
}

void iterate_simulation()
{
	int new_one;
	int x, y;
	float cva, cvb, newa, newb, olda, oldb;

	new_one = (current + 1) % 2;

	for (x = 0; x < xdim; x++) {
		for (y = 0; y < ydim; y++) {
		
			olda = a[current][x][y];
			oldb = b[current][x][y];
	
			cva = convolve(current, x, y, 1);
			cvb = convolve(current, x, y, 0);

			newa = olda + (diffusion_rate_a * cva -
					olda * oldb * oldb + feed_rate_a * (1 - olda)) * delta_t;
			newb = oldb + (diffusion_rate_b * cvb +
					olda * oldb * oldb -
					(kill_rate_b + feed_rate_a) * oldb) * delta_t;

			a[new_one][x][y] = newa;
			b[new_one][x][y] = newb;
		}
	}
	current = new_one;
}

void draw()
{
	draw_current(current);
	iterate_simulation();
}

