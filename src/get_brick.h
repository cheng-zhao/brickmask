/*******************************************************************************
* get_brick.h: this file is part of the brickmask program.

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

#ifndef __GET_BRICK_H__
#define __GET_BRICK_H__

#include "load_conf.h"
#include "data_io.h"
#include <stddef.h>
#include <stdint.h>

/*============================================================================*\
                           Data structure for bricks
\*============================================================================*/

/* Data structure for bricks. */
typedef struct {
  size_t n;             /* number of bricks                             */
  double *ra1;          /* minimum right ascension of the bricks        */
  double *ra2;          /* maximum right ascension of the bricks        */
  double *dec1;         /* minimum declination of the bricks            */
  double *dec2;         /* maximum declination of the bricks            */
  char **name;          /* name of the bricks                           */
  int nsp;              /* number of subsamples                         */
  int *subid;           /* IDs of subsamples                            */
  size_t *nmask;        /* number of maskbit files for each subsample   */
  char ***fmask;        /* names of maskbit files                       */
  uint64_t mnull;       /* bit code for objects outside maskbit bricks  */
#ifdef MPI
  int nlen;             /* length of the brick names                    */
  size_t *mlen;         /* total length of filenames for each subsample */
#endif
} BRICK;

/*============================================================================*\
                        Interfaces for setting up bricks
\*============================================================================*/

/******************************************************************************
Function `get_brick`:
  Get brick information from files.
Arguments:
  * `conf`:     structure for storing configurations.
Return:
  Address of the structure for bricks on success; NULL on error.
******************************************************************************/
BRICK *get_brick(const CONF *conf);

/******************************************************************************
Function `sort_data`:
  Sort the input data sample based on the brick IDs.
Arguments:
  * `brick`:    structure for bricks;
  * `data`:     structure for the data catalogue;
  * `verbose`:  indicate whether to show detailed outputs.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int sort_data(BRICK *brick, DATA *data, const bool verbose);

/******************************************************************************
Function `brick_destroy`:
  Deconstruct the structure for bricks.
Arguments:
  * `data`:     structure for bricks.
******************************************************************************/
void brick_destroy(BRICK *brick);

#endif
