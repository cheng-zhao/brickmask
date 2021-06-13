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

from astropy.io import fits
import numpy as np
import os

dr = 'dr8'
caps = ['north', 'south']

legacy_dir = '/global/project/projectdirs/cosmo/data/legacysurvey'
bricks_fmt = '{}/{}/survey-bricks-{}-{}.fits.gz'
maskbit_fmt = '{}/{}/coadd/{}/{}/legacysurvey-{}-maskbits.fits.fz'

output_fmt = 'legacysurvey_maskbits_{}_{}.txt'

ilist = f'/global/project/projectdirs/cosmo/data/legacysurvey/{dr}/{{}}/survey-bricks-{dr}-{{}}.fits.gz'
olist = f'legacysurvey_maskbits_{dr}_{{}}.txt'

for cap in caps:
  ifile = '{}/{}'.format(legacy_dir, bricks_fmt.format(dr,cap,dr,cap))
  ofile = output_fmt.format(dr, cap)

  if not os.path.isfile(ifile):
    print(f'Error: cannot access brick list: {ifile}')
    exit(1)

  names = np.unique(fits.open(ifile)[1].data['brickname'])
  fnames = ['{}/{}\n'.format(legacy_dir, maskbit_fmt.format(dr,cap,n[:3],n,n)) for n in names]

  with open(ofile, 'w') as f:
    f.writelines(fnames)

