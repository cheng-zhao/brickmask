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

#include "brickmask.h"

/******************************************************************************
Function `find_brick`:
  Find the index of the brick that contains the current coordinate, using
  binary search.

Arguments:
  * `ra`:       RA of the current coordinate;
  * `dec`:      Dec of the current coordinate;
  * `ra1`:      lower limit of RA for the bricks;
  * `ra2`:      upper limit of RA for the bricks;
  * `dec1`:     lower limit of Dec for the bricks;
  * `dec2`:     upper limit of Dec for the bricks;
  * `n`:        the number of bricks.
Return:
  If the brick is found, return the index; else return -1.
******************************************************************************/
long find_brick(const double ra, const double dec, const double *ra1,
    const double *ra2, const double *dec1, const double *dec2,
    const size_t n) {
  long i, l, u;
  l = 0;
  u = n - 1;

  while (l <= u) {
    i = (l + u) >> 1;
    if (compare_pos(ra, dec, ra1[i], ra2[i], dec1[i], dec2[i]) > 0) l = i + 1;
    else if (compare_pos(ra, dec, ra1[i], ra2[i], dec1[i], dec2[i]) < 0)
      u = i - 1;
    else return i;
  }

  return -1;
}


/******************************************************************************
Function `compare_pos`:
  Compare a coordinate with a range (first Dec, and then RA).

Arguments:
  * `ra`:       RA of the current coordinate;
  * `dec`:      Dec of the current coordinate;
  * `ra1`:      lower limit of the RA range;
  * `ra2`:      upper limit of the RA range;
  * `dec1`:     lower limit of the Dec range;
  * `dec2`:     upper limit of the Dec range.
Return:
  +1 if the coordinate is larger than the range, -1 if smaller, and 0 if the
  coordinate is inside the range.
******************************************************************************/
int compare_pos(const double ra, const double dec, const double ra1,
    const double ra2, const double dec1, const double dec2) {
  if (dec < dec1) return -1;
  if (dec >= dec2) return 1;
  if (ra < ra1) return -1;
  if (ra >= ra2) return 1;
  return 0;
}


/******************************************************************************
Function `compare_idx`:
  Compare the indices of bricks, to be used for sorting the input data.

Arguments:
  * `a`:        a pointer to one data record;
  * `b`:        a pointer to another data record.
Return:
  +1 if the index of a is larger than b, -1 of smaller, and 0 if the two
  indices are identical.
******************************************************************************/
int compare_idx(const void *a, const void *b) {
  if (((DATA *) a)->idx > ((DATA *) b)->idx) return 1;
  if (((DATA *) a)->idx < ((DATA *) b)->idx) return -1;
  return 0;
}

