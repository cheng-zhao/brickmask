/*******************************************************************************
* mpi_schedule.h: this file is part of the brickmask program.

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

#ifdef MPI
#ifndef __MPI_SCHEDULE_H__
#define __MPI_SCHEDULE_H__

#include "get_brick.h"
#include "data_io.h"

/*============================================================================*\
                 Function for assigning maskbits with MPI tasks
\*============================================================================*/

/******************************************************************************
Function `mpi_init_worker`:
  Initialise MPI workers with brick lists and the input data.
Arguments:
  * `brick`:    structure for bricks;
  * `data`:     structure for the data catalogue;
  * `verbose`:  indicate whether to show detailed standard outputs.
******************************************************************************/
void mpi_init_worker(BRICK **brick, DATA **data, const bool verbose);

/******************************************************************************
Function `mpi_gather_data`:
  Gather data from different MPI tasks.
Arguments:
  * `brick`:    structure for bricks;
  * `data`:     structure for the data catalogue.
******************************************************************************/
void mpi_gather_data(BRICK *brick, DATA *data);

#endif
#endif
