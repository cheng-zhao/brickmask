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

#include <stddef.h>
#include <stdint.h>
#include "load_conf.h"

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
  unsigned char *photsys;       /* photometric region (N or S)          */
  char *bpath;          /* path of the bricks                           */
  uint64_t mnull;       /* bit code for objects outside maskbit bricks  */
  uint64_t enull;       /* bit code for objects outside nexp bricks     */
#ifdef MPI
  int nlen;             /* length of the brick names                    */
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
Function `brick_destroy`:
  Deconstruct the structure for bricks.
Arguments:
  * `data`:     structure for bricks.
******************************************************************************/
void brick_destroy(BRICK *brick);

#endif
