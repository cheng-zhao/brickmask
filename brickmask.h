/*******************************************************************************
* 
* brickmask: assign bit codes defined on Legacy Survey brick pixels
*            to a catalogue with sky coordinates.
*
* Github repository:
*       https://github.com/cheng-zhao/brickmask
*
* Copyright (c) 2020 Cheng Zhao <zhaocheng03@gmail.com>
* 
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
* 
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*
*******************************************************************************/

#ifndef _BRICKMASK_H_
#define _BRICKMASK_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fitsio.h>
#include "define.h"
#ifdef OMP
#include <omp.h>
#endif

/********** read_data.c **********/

int read_ascii(const char *, const char, DATA **, size_t *, const int);

int read_fits(const char *, DATA **, size_t *, const int);

int read_list(const char *, char ***, double **, double **, double **,
    double **, size_t *, const int);

int read_wcs_header(fitsfile *, WCS *);

/********** read_data.c **********/

/********** baselib.c **********/

size_t safe_strcpy(char *, const char *, const size_t);

int wcs_world2pix(const WCS *, const double, const double, double *, double *);

/********** baselib.c **********/

/********** find_brick.c **********/

long find_brick(const double, const double, const double *, const double *,
    const double *, const double *, const size_t);

int compare_pos(const double, const double, const double, const double,
    const double, const double);

int compare_idx(const void *, const void *);

/********** find_brick.c **********/

/********** save_res.c **********/

int write_ascii(const char *, const DATA *, const size_t, const int);

int write_fits(const char *, const char *, const DATA *, const long, const int);

/********** save_res.c **********/

#endif
