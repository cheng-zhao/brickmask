#!/usr/bin/env python3
#
# Copyright (c) 2020 - 2021 Cheng Zhao <zhaocheng03@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

import numpy as np
from scipy import in1d
import healpy as hp
import pymangle
import copy

def eBOSS_ELG_mask(ra, dec, mskdir, mskbit=None):
  '''
  Assign extra masks to eBOSS ELGs.

  Parameters
  ----------
  ra: array_like
        Right Ascension of the catalogue.
  dec: array_like
        Declination of the catalogue.
  mskdir: str
        Directory containing the following mask files:
        * ELG_centerpost.ply
        * ELG_TDSSFES_62arcsec.pix.snap.balk.ply
        * ebosselg_badphot.26Aug2019.ply
        These files are available at:
        https://data.sdss.org/sas/dr16/eboss/lss/catalogs/DR16/ELGmasks
  mskbit: array_like, optional
        Maskbit codes assigned by the `brickmask` program.
        See https://github.com/cheng-zhao/brickmask

  Returns
  -------
  mask: array
        Maskbit codes with extra masks.

  References
  ----------
  * https://doi.org/10.1093/mnras/staa3336 (https://arxiv.org/abs/2007.09007)
  * https://doi.org/10.1093/mnras/stab510  (https://arxiv.org/abs/2007.08997)
  '''
  # Validate inputs
  ra = np.atleast_1d(ra)
  dec = np.atleast_1d(dec)
  if mskbit is None:
    mask = np.zeros_like(ra, dtype='int16')
  else:
    mask = copy.copy(np.atleast_1d(mskbit)).astype('int16')

  if len(ra) != len(dec) or len(ra) != len(mask):
    raise ValueError('ra, dec, and mskbit (if set) must have the same length')

  # Add 2**8 for discrepancy between mskbit and anymask
  mask_pixel = [2981667,3464728,3514005,3645255,4546075,4685432,5867869, \
        5933353,6031493,6072514,6080368,6092477,6301369,6408277,6834661, \
        2907700,3583785,3587880,4067035,4669088,6007074,6186688,6190785, \
        6199270,6371066,6547876,6551972,6645991,6711673,6735965,6744444, \
        6744445,6748540,6752636,6769023,6773119,6781133]
  theta = np.radians(90 - dec)
  phi = np.radians(360 - ra)
  pix = hp.pixelfunc.ang2pix(1024, theta, phi, nest=False, lonlat=True)
  bit = in1d(pix, mask_pixel)
  mask += bit * 2**8

  # Add bit codes for polygon masks
  plys = ['ELG_centerpost.ply', 'ELG_TDSSFES_62arcsec.pix.snap.balk.ply', \
         'ebosselg_badphot.26Aug2019.ply']
  ply_codes = [9, 10, 11]       # Bit code for the polygon masks

  for ply, code in zip(plys, ply_codes):
    iply = '{}/{}'.format(mskdir, ply)
    m = pymangle.Mangle(iply)
    bit = m.polyid(ra, dec) != -1
    mask += bit * 2**code

  return mask
