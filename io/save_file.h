/*******************************************************************************
* save_res.h: this file is part of the brickmask program.

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

#ifndef __SAVE_RES_H__
#define __SAVE_RES_H__

#include "data_io.h"

/*============================================================================*\
                    Interfaces for saving output catalogues
\*============================================================================*/

/******************************************************************************
Function `save_ascii`:
  Write the data catalogue to an ASCII file.
Arguments:
  * `conf`:     structure for storing configurations;
  * `data`:     structure for the the data catalogue.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int save_ascii(const CONF *conf, const DATA *data);

/******************************************************************************
Function `save_fits`:
  Write the data catalogue to a FITS file.
Arguments:
  * `conf`:     structure for storing configurations;
  * `data`:     structure for the the data catalogue.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int save_fits(const CONF *conf, DATA *data);

#endif
