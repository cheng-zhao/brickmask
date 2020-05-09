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
Function `safe_strcpy`:
  Copy the first `num`-1 characters from `src` to `dest`, with null
  termination.

Arguments:
  * `dest`:     pointer to the destination array for the resulting string;
  * `src`:      the source string to be copied;
  * `num`:      maximum number of bytes to be written in `dest`.
Return:
  Number of characters that would have been written if `num` had been
  sufficiently large, not counting the terminating null character
  (see also `snprintf`).
******************************************************************************/
size_t safe_strcpy(char *dest, const char *src, const size_t num) {
  size_t i = 0;

  while (i < num - 1 && src[i] != '\0') {
    dest[i] = src[i];
    i++;
  }
  dest[i] = '\0';
  while (src[i] != '\0') i++;
  return i;
}


/******************************************************************************
Function `wcs_world2pix`:
  Convert WCS world coordinates to pixels.

Arguments:
  * `wcs`:      the pre- read/computed WCS keywords.
  * `ra`:       the input RA;
  * `dec`:      the input Dec;
  * `x`:        the output x in pixels;
  * `y`:        the output y in pixels;
Return:
  Non-zero integer if there is problem.
******************************************************************************/
int wcs_world2pix(const WCS *wcs, const double ra, const double dec,
    double *x, double *y) {
  int i, j;
  double alpha[2], theta, phi, b[3], l[3], xx, yy;

  alpha[0] = ra * M_PI / 180.0;
  alpha[1] = dec * M_PI / 180.0;

  if (dec == 90 || dec == -90) l[0] = l[1] = 0;
  else {
    l[0] = (ra == 90 || ra == 270) ? 0 : cos(alpha[1]) * cos(alpha[0]);
    l[1] = (ra == 0 || ra == 180) ? 0 : cos(alpha[1]) * sin(alpha[0]);
  }
  l[2] = (dec == 0) ? 0 : sin(alpha[1]);

  b[0] = b[1] = b[2] = 0;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++) b[j] += wcs->r[j][i] * l[i];

  for (i = 0; i < 3; i++) {
    if (b[i] > 1) b[i] = 1;
    if (b[i] < -1) b[i] = -1;
  }

  theta = asin(b[2]);
  phi = atan2(b[1], b[0]);
  if (theta == 0) {
    P_ERR("failed to convert coordinate (%g, %g)\n", ra, dec);
    return WCS_ERROR;
  }

  theta = 180.0 / (tan(theta) * M_PI);  /* R_theta */
  xx = theta * sin(phi);
  yy = -theta * cos(phi);
  *x = wcs->cd[1][1] * xx - wcs->cd[0][1] * yy;
  *y = -wcs->cd[1][0] * xx + wcs->cd[0][0] * yy;
  *x = wcs->crpix[0] + (*x) / wcs->cdtmp - 1;
  *y = wcs->crpix[1] + (*y) / wcs->cdtmp - 1;

  return 0;
}

