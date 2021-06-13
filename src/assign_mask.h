/*******************************************************************************
* assign_mask.h: this file is part of the brickmask program.

* brickmask: assign bit codes defined on Legacy Survey brick pixels
             to a catalogue with sky coordinates.

* Github repository:
        https://github.com/cheng-zhao/brickmask

* Copyright (c) 2020 -- 2021 Cheng Zhao <zhaocheng03@gmail.com>
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

*******************************************************************************/

#ifndef __ASSIGN_MASK_H__
#define __ASSIGN_MASK_H__

#include "load_conf.h"
#include "get_brick.h"
#include "data_io.h"
#include <stddef.h>

/*============================================================================*\
                          Data structures for maskbits
\*============================================================================*/

/* Data structure for the WCS header for coordinate conversion. */
typedef struct {
  double ang[8];        /* factors for computing the celestial angles */
  double r[2];          /* reference pixel coordinates                */
  double m[2][2];       /* linear transformation matrix               */
  double idetm;         /* inversed determinant of `m`                */
} WCS;

/* Data structure for the maskbits. */
typedef struct {
  long size;                    /* total number of maskbits             */
  int dtype;                    /* data type of maskbit values          */
  long dim[2];                  /* dimensions of maskbits               */
  uint64_t mnull;               /* bit code for objects outside bricks  */
  unsigned char *bit;           /* maskbit values                       */
  WCS *wcs;                     /* WCS parameters                       */
} MASK;


/*============================================================================*\
                        Interface for assigning maskbits
\*============================================================================*/

/******************************************************************************
Function `assign_mask`:
  Assign maskbits to the data catalogue.
Arguments:
  * `brick`:    structure for bricks;
  * `data`:     structure for the data catalogue;
  * `verbose`:  indicate whether to show detailed standard outputs.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int assign_mask(const BRICK *brick, DATA *data, const bool verbose);

#endif
