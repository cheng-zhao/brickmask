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
  * `ra`:       right ascension of the data catalogue;
  * `dec`:      declination of the data catalogue;
  * `code`:     array for maskbits;
  * `ndata`:    number of data points to be processed.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int MASKBIT_FUNC(assign_bitcode, BRICKMASK_MASKBIT_DTYPE)
    (const MASK *mask, const double *ra, const double *dec, uint64_t *code,
     const size_t ndata) {
  const BRICKMASK_MASKBIT_DTYPE *bits = (BRICKMASK_MASKBIT_DTYPE *) mask->bit;
  for (size_t i = 0; i < ndata; i++) {
    double x, y;
    world2pix(mask->wcs, ra[i], dec[i], &x, &y);
    long rx = round(x);
    long ry = round(y);
    if (rx < 0 || rx >= mask->dim[0] || ry < 0 || ry >= mask->dim[1]) {
      P_ERR("invalid pixel value (%ld, %ld) for coordinate (" OFMT_DBL ", "
          OFMT_DBL ")\n", rx, ry, ra[i], dec[i]);
      return BRICKMASK_ERR_MASK;
    }
    uint64_t bit = bits[rx + ry * mask->dim[0]];
    code[i] += bit;
  }

  return 0;
}

#undef BRICKMASK_MASKBIT_DTYPE

#endif
