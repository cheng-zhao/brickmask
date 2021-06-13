/*******************************************************************************
* bit_code.c: this file is part of the brickmask program.

* brickmask: assign bit codes defined on Legacy Survey brick pixels
             to a catalogue with sky coordinates.

* Github repository:
        https://github.com/cheng-zhao/brickmask

* Copyright (c) 2020 -- 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

/* Macros for the template functions. */
#ifdef BRICKMASK_MASKBIT_DTYPE

/* Macros for generating function names. */
#ifndef CONCAT_FNAME
  #define CONCAT_FNAME(a,b)             a##_##b
#endif

#ifndef MASKBIT_FUNC
  #define MASKBIT_FUNC(a,b)             CONCAT_FNAME(a,b)
#endif

/*============================================================================*\
                  Function for constructing cut-sky catalogue
\*============================================================================*/

/******************************************************************************
Function `assign_bitcode_<BRICKMASK_MASKBIT_DTYPE>`:
  Assign maskbit codes to objects.
Arguments:
  * `mask`:     structure for maskbits;
  * `data`:     structure for the data catalogue;
  * `imin`:     starting index of the data to be processed;
  * `imax`:     ending index of the data to be processed;
  * `subid`:    ID of the subsample.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int MASKBIT_FUNC(assign_bitcode, BRICKMASK_MASKBIT_DTYPE)
    (const MASK *mask, DATA *data, const size_t imin, const size_t imax,
    const unsigned char subid) {
  const BRICKMASK_MASKBIT_DTYPE *bits = (BRICKMASK_MASKBIT_DTYPE *) mask->bit;
  for (size_t i = imin; i < imax; i++) {
    double x, y;
    world2pix(mask->wcs, data->ra[i], data->dec[i], &x, &y);
    long rx = round(x);
    long ry = round(y);
    if (rx < 0 || rx >= mask->dim[0] || ry < 0 || ry >= mask->dim[1]) {
      P_ERR("invalid pixel value (%ld, %ld) for coordinate (" OFMT_DBL ", "
          OFMT_DBL ")\n", rx, ry, data->ra[i], data->dec[i]);
      return BRICKMASK_ERR_MASK;
    }
    uint64_t bit = bits[rx + ry * mask->dim[0]];
#ifdef EBOSS
    if (!(EBOSS_MASK_VALID(bit))) continue;
    if (EBOSS_XYBUG_VALID(bit)) data->mask[i] += bit - EBOSS_XYBUG_BIT;
    else data->mask[i] += bit;

    long ix = (long) x;
    long iy = (long) y;
    if (ix < 0 || ix >= mask->dim[0] || iy < 0 || iy >= mask->dim[1]) {
      P_ERR("invalid pixel value (%ld, %ld) for coordinate (" OFMT_DBL ", "
          OFMT_DBL ")\n", ix, iy, data->ra[i], data->dec[i]);
      return BRICKMASK_ERR_MASK;
    }
    bit = bits[ix + iy * mask->dim[0]];
    if (EBOSS_XYBUG_VALID(bit)) data->mask[i] += EBOSS_XYBUG_BIT;
    /* Assign subsample ID if necessary. */
    if (data->subid) data->subid[i] = subid;
#else
    data->mask[i] += bit;
    /* Assign subsample ID if necessary. */
    if (!(bit & mask->mnull) && data->subid) data->subid[i] = subid;
#endif
  }

  return 0;
}

#undef BRICKMASK_MASKBIT_DTYPE

#endif
