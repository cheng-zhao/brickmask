/*******************************************************************************
* data_io.h: this file is part of the brickmask program.

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

#ifndef __DATA_IO_H__
#define __DATA_IO_H__

#include "load_conf.h"
#include <stddef.h>
#include <stdint.h>

/*============================================================================*\
                         Data structures for the inputs
\*============================================================================*/

/* Type of the input catalogue. */
typedef enum {
  BRICKMASK_FFMT_ASCII = 0,
  BRICKMASK_FFMT_FITS = 1
} BRICKMASK_ffmt_t;

/* Data structure for the input catalogue. */
typedef struct {
  BRICKMASK_ffmt_t fmt; /* format of the input data catalogue           */
  int mtype;            /* data type of maskbits                        */
  int etype[3];         /* data type of nexp bits                       */
  size_t n;             /* number of objects read from file             */
  size_t nmax;          /* number of objects allocated for the data     */
  size_t csize;         /* size of the content read from file           */
  size_t cmax;          /* number of bytes allocated for the content    */
  size_t nbrick;        /* total number of bricks containing the data   */
  double *ra;           /* right ascension                              */
  double *dec;          /* declination                                  */
  size_t *idx;          /* index of the data before sorting             */
  size_t *cidx;         /* ASCII: index for the rest of the columns     */
  size_t *iidx;         /* index range for different input catalogues   */
  long *id;             /* brick ID, signed type for sorting comparison */
  uint64_t *mask;       /* maskbit value                                */
  uint64_t *nexp[3];    /* nexp values                                  */
  void *content;        /* ASCII: address for the rest of the columns
                           FITS:  properties of output columns          */
} DATA;

/* Data structure for information of FITS columns. */
typedef struct {
  int n;                /* column number                                */
  long i;               /* starting position of the column              */
  long w;               /* width (number of bytes) of the column        */
} FITS_COL_t;

/*============================================================================*\
               Interfaces for reading and saving data catalogues
\*============================================================================*/

/******************************************************************************
Function `read_data`:
  Read data from the input catalogue.
Arguments:
  * `conf`:     structure for storing configurations.
Return:
  Address of the structure for the input catalogue on success; NULL on error.
******************************************************************************/
DATA *read_data(const CONF *conf);

/******************************************************************************
Function `save_data`:
  Save data to the output catalogue.
Arguments:
  * `conf`:     structure for storing configurations;
  * `data`:     structure for the data catalogue.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int save_data(const CONF *conf, DATA *data);

/******************************************************************************
Function `data_destroy`:
  Deconstruct the structure for the input data.
Arguments:
  * `data`:     structure for the input data catalogue.
******************************************************************************/
void data_destroy(DATA *data);

#endif
