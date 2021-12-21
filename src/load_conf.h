/*******************************************************************************
* load_conf.h: this file is part of the brickmask program.

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

#ifndef __LOAD_CONF_H__
#define __LOAD_CONF_H__

#include <stdbool.h>

/*============================================================================*\
                   Data structure for storing configurations
\*============================================================================*/

typedef struct {
  char *fconf;          /* Name of the configuration file. */
  char *flist;          /* BRICK_LIST           */
  char **fmask;         /* MASKBIT_FILES        */
  int mnull;            /* MASKBIT_NULL         */
  int nsub;             /* Number of subsamples. */
  int *subid;           /* SUBSAMPLE_ID         */
  char *input;          /* INPUT                */
  int ftype;            /* FILE_TYPE            */
  char comment;         /* ASCII_COMMENT        */
  char **cname;         /* COORD_COLUMN         */
  int cnum[2];          /* Column number of (RA,Dec) for ASCII input. */
  char *output;         /* OUTPUT               */
  char **ocol;          /* OUTPUT_COLUMN        */
  int ncol;             /* Number of output columns. */
  int *onum;            /* Column numbers to be saved to the output. */
  char *mcol;           /* MASKBIT_COLUMN       */
  int ovwrite;          /* OVERWRITE            */
  bool verbose;         /* VERBOSE              */
} CONF;


/*============================================================================*\
                      Interface for loading configurations
\*============================================================================*/

/******************************************************************************
Function `load_conf`:
  Read, check, and print configurations.
Arguments:
  * `argc`:     number of arguments passed via command line;
  * `argv`:     array of command line arguments.
Return:
  The structure for storing configurations.
******************************************************************************/
CONF *load_conf(const int argc, char *const *argv);

/******************************************************************************
Function `conf_destroy`:
  Release memory allocated for the configurations.
Arguments:
  * `conf`:     the structure for storing configurations.
******************************************************************************/
void conf_destroy(CONF *conf);

#endif
