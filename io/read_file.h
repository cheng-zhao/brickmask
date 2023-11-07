/*******************************************************************************
* read_file.h: this file is part of the brickmask program.

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

#ifndef __READ_FILE_H__
#define __READ_FILE_H__

#include "load_conf.h"
#include "data_io.h"
#include "get_brick.h"
#include "assign_mask.h"

/*============================================================================*\
                       Interfaces for reading input files
\*============================================================================*/

/******************************************************************************
Function `read_ascii`:
  Read data from the input ASCII catalogue.
Arguments:
  * `fname`:    filename of the input catalogue;
  * `conf`:     structure for storing configurations;
  * `data`:     structure for the input data catalog.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int read_ascii(const char *fname, const CONF *conf, DATA *data);

/******************************************************************************
Function `read_fname`:
  Read filenames from a text file.
Arguments:
  * `fname`:    name of the file storing filenames;
  * `output`:   string array for storing filenames;
  * `num`:      number of filenames read from file.
Return:
  Number of bytes read on success; zero on error.
******************************************************************************/
size_t read_fname(const char *fname, char ***output, size_t *num);

/******************************************************************************
Function `read_brick`:
  Read the brick name and range of (RA, Dec) from a brick list file.
Arguments:
  * `fname`:    the filename of the brick list;
  * `brick`:    structure for storing information of bricks.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int read_brick(const char *fname, BRICK *brick);

/******************************************************************************
Function `read_fits`:
  Read data from the input FITS catalogue.
Arguments:
  * `fname`:    filename of the input catalogue;
  * `conf`:     structure for storing configurations;
  * `data`:     structure for the input data catalogue.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int read_fits(const char *fname, const CONF *conf, DATA *data);

/******************************************************************************
Function `read_mask`:
  Read a maskbit file.
Arguments:
  * `fname`:    name of a masbit file;
  * `mask`:     structure for maskbits.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int read_mask(const char *fname, MASK *mask);

#endif
