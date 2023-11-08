/*******************************************************************************
* assign_mask.c: this file is part of the brickmask program.

* brickmask: assign bit codes defined on Legacy Survey brick pixels
             to a catalogue with sky coordinates.

* Github repository:
        https://github.com/cheng-zhao/brickmask

* Copyright (c) 2020 -- 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "assign_mask.h"
#include "read_file.h"
#include <fitsio.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#ifdef MPI
#include <mpi.h>
#define BRICKMASK_QUIT(status)  MPI_Abort(MPI_COMM_WORLD, status); exit(status);
#else
#define BRICKMASK_QUIT(status)  return status;
#endif

/*============================================================================*\
                       Functions for processing maskbits
\*============================================================================*/

/******************************************************************************
Function `mask_destroy`:
  Deconstruct the structure for maskbits.
Arguments:
  * `mask`:     structure for maskbits.
******************************************************************************/
static void mask_destroy(MASK *mask) {
  if (!mask) return;
  if (mask->bit) free(mask->bit);
  if (mask->wcs) free(mask->wcs);
  free(mask);
}

/******************************************************************************
Function `mask_init`:
  Initialise the structure for maskbits.
Arguments:
  * `mnull`:    bit code for objects outside bricks;
  * `enull`:    bit code for objects without nexp.
Return:
  Address of the structure for maskbits on success; NULL on error.
******************************************************************************/
static inline MASK *mask_init(const uint64_t mnull, const uint64_t enull) {
  MASK *mask = calloc(1, sizeof(MASK));
  if (!mask) {
    P_ERR("failed to allocate memory for maskbits\n");
    return NULL;
  }
  WCS *wcs = calloc(1, sizeof(WCS));
  if (!wcs) {
    P_ERR("failed to allocate memory for WCS parameters\n");
    free(mask);
    return NULL;
  }

  mask->wcs = wcs;
  mask->bit = NULL;
  mask->mnull = mnull;
  mask->enull = enull;

  return mask;
}

/******************************************************************************
Function `world2pix`:
  Convert world coordinates (RA, Dec) to pixel indices with the 'TAN' scheme.
  Ref: https://doi.org/10.1051/0004-6361:20021327
Arguments:
  * `wcs`:      structure for WCS parameters;
  * `ra`:       the input right ascension;
  * `dec`:      the input declination;
  * `x`:        the output pixel index along the x direction;
  * `y`:        the output pixel index along the y direction.
******************************************************************************/
static inline void world2pix(const WCS *wcs, const double ra, const double dec,
    double *x, double *y) {
  double r = ra * DEGREE_2_RAD;
  double d = dec * DEGREE_2_RAD;
  double sina = sin(r);
  double cosa = cos(r);
  double sind = sin(d);
  double cosd = cos(d);

  double fac1 = cosa * cosd;
  double fac2 = sina * cosd;

  double theta = sind * wcs->ang[0] + fac1 * wcs->ang[1] + fac2 * wcs->ang[2];
  double phi1 = sind * wcs->ang[3] + fac1 * wcs->ang[4] + fac2 * wcs->ang[5];
  double phi2 = fac1 * wcs->ang[6] + fac2 * wcs->ang[7];

  theta = (theta >= 1 || theta <= -1) ? 0 :
      sqrt(1 - theta * theta) / theta * RAD_2_DEGREE;
  double fac = (phi1 == 0 && phi2 == 0) ? 0 :
      theta / sqrt(phi1 * phi1 + phi2 * phi2);
  double xx = fac * phi2;
  double yy = -fac * phi1;

  *x = (xx * wcs->m[1][1] - yy * wcs->m[0][1]) * wcs->idetm + wcs->r[0] - 1;
  *y = (-xx * wcs->m[1][0] + yy * wcs->m[0][0]) * wcs->idetm + wcs->r[1] - 1;
}

/******************************************************************************
Function `num_digit`:
  Compute the number of digits of an unsigned integer.
Arguments:
  * `num`:      the integer to be examined.
Return:
  Number of digits.
******************************************************************************/
static inline int num_digit(size_t num) {
  int n;
  for (n = 1; num > 9; n++) num /= 10;
  return n;
}


/*============================================================================*\
                   Template functions for assigning maskbits
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

#ifdef BRICKMASK_MASKBIT_DTYPE
  #undef BRICKMASK_MASKBIT_DTYPE
#endif

#define BRICKMASK_MASKBIT_DTYPE uint8_t
#include "bit_code.c"

#define BRICKMASK_MASKBIT_DTYPE uint16_t
#include "bit_code.c"

#define BRICKMASK_MASKBIT_DTYPE uint32_t
#include "bit_code.c"

#define BRICKMASK_MASKBIT_DTYPE uint64_t
#include "bit_code.c"


/*============================================================================*\
                        Function for assigning maskbits
\*============================================================================*/

/******************************************************************************
Function `assign_mask_from_file`:
  Read maskbits from file and assign them to the data catalogue.
Arguments:
  * `fname`:    filename of the maskbits file;
  * `mask`:     structure for maskbits;
  * `ra`:       right ascension of the data catalogue;
  * `dec`:      declination of the data catalogue;
  * `code`:     array for maskbits;
  * `vnull`:    bit code for objects outside the mask bricks;
  * `dtype`:    data type of bit codes;
  * `ndata`:    number of data points to be processed.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int assign_mask_from_file(const char *fname, MASK *mask,
    const double *ra, const double *dec, uint64_t *code, const uint64_t vnull,
    int *dtype, const size_t ndata) {
  /* Check if the maskbit file exists. */
  if (access(fname, R_OK)) {
    for (size_t i = 0; i < ndata; i++) code[i] = vnull;
    return 0;
  }
  /* Read maskbits. */
  if (read_mask(fname, mask, dtype)) return BRICKMASK_ERR_MASK;

  /* Choose the maskbit code assignment function given the data type. */
  int (*assign_bitcode_func) (const MASK *, const double *, const double *,
      uint64_t *, size_t) = NULL;
  switch (*dtype) {
    case TBYTE:  assign_bitcode_func = assign_bitcode_uint8_t;  break;
    case TSHORT: assign_bitcode_func = assign_bitcode_uint16_t; break;
    case TINT:   assign_bitcode_func = assign_bitcode_uint32_t; break;
    case TLONG:  assign_bitcode_func = assign_bitcode_uint64_t; break;
    default:
      P_ERR("unexpected data type for maskbits: %d\n", *dtype);
      return BRICKMASK_ERR_MASK;
  }

  if (assign_bitcode_func(mask, ra, dec, code, ndata))
    return BRICKMASK_ERR_MASK;
  return 0;
}

/*============================================================================*\
                        Interface for assigning maskbits
\*============================================================================*/

/******************************************************************************
Function `assign_mask`:
  Assign maskbits to the data catalogue.
Arguments:
  * `brick`:    structure for bricks;
  * `data`:     structure for the data catalogue;
  * `verbose`:  indicate whether to show detailed standard outputs.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int assign_mask(const BRICK *brick, DATA *data, const bool verbose) {
#ifdef MPI
  int size, rank;
  size = rank = 0;
  if (MPI_Comm_size(MPI_COMM_WORLD, &size) ||
      MPI_Comm_rank(MPI_COMM_WORLD, &rank))
    MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);

  if (rank == BRICKMASK_MPI_ROOT) {
#endif
    printf("Assigning maskbits to the data ...");
    if (verbose) printf("\n");
    fflush(stdout);
#ifdef MPI
  }
#endif

  if (!brick || ! data) {
    P_ERR("the bricks or input data catalogue is not initialised\n");
    BRICKMASK_QUIT(BRICKMASK_ERR_INIT);
  }

  /* Initialise progress printing. */
  size_t cnt = 0;
  size_t next = data->nbrick / BRICKMASK_PROGRESS_NUM;
  if (!next) next = 1;
  size_t step = next;
  int ndg = num_digit(data->nbrick);
  int wcol = ndg * 2 + 3;               /* column width of the progress */

#ifdef MPI
  if (verbose && rank == BRICKMASK_MPI_ROOT) {
    printf("  Bricks processed by root task: %*zu / %zu",
        ndg, cnt, data->nbrick);
    fflush(stdout);
  }
#else
  if (verbose) {
    printf("  Bricks processed: %*zu / %zu", ndg, cnt, data->nbrick);
    fflush(stdout);
  }
#endif

  if (!brick->n || !data->n) {
#ifndef MPI
    P_WRN("no brick or data is available\n");
    printf(FMT_DONE);
#endif
    return 0;
  }

  /* Initialise maskbits. */
  MASK *mask = mask_init(brick->mnull, brick->enull);
  if (!mask) {
    BRICKMASK_QUIT(BRICKMASK_ERR_MEMORY);
  }

  /* Read and assign maskbits. */
  char fname[BRICKMASK_MAX_FILENAME];
  size_t imin, imax;
  imin = 0;
  while (imin < data->n) {
    /* Get the index range of data sharing the same brick. */
    for (imax = imin + 1; imax < data->n; imax++) {
      if (data->id[imax] != data->id[imin]) break;
    }
    size_t bid = data->id[imin];        /* ID of the corresponding brick. */

    /* Determine filename of the maskbits brick. */
    char photsys = brick->photsys[bid];
    if (photsys != 'N' && photsys != 'S') {     /* no photometric data */
      for (size_t i = imin; i < imax; i++) {
        data->mask[i] = mask->mnull;
        data->nexp[0][i] = data->nexp[1][i] = data->nexp[2][i] = mask->enull;
        imin = imax;
        continue;
      }
    }
    const char *ns = (photsys == 'N') ? "north" : "south";
    int nchar = snprintf(fname, BRICKMASK_MAX_FILENAME,
        "%s/%s/coadd/%.3s/%s/legacysurvey-%s-maskbits.fits.fz",
        brick->bpath, ns, brick->name[bid], brick->name[bid], brick->name[bid]);
    if (nchar >= BRICKMASK_MAX_FILENAME) {
      P_ERR("the brickmask filename is too long\n"
          "Please enlarge `BRICKMASK_MAX_FILENAME` in `src/define.h`\n");
      mask_destroy(mask);
      BRICKMASK_QUIT(BRICKMASK_ERR_FILE);
    }

    if (assign_mask_from_file(fname, mask, data->ra + imin, data->dec + imin,
        data->mask + imin, mask->mnull, &(mask->dtype), imax - imin)) {
      mask_destroy(mask);
      BRICKMASK_QUIT(BRICKMASK_ERR_MASK);
    }

    /* Read NOBS_G. */
    const char nobs[3] = {'g', 'r', 'z'};
    for (int k = 0; k < 3; k++) {
      nchar = snprintf(fname, BRICKMASK_MAX_FILENAME,
          "%s/%s/coadd/%.3s/%s/legacysurvey-%s-nexp-%c.fits.fz",
          brick->bpath, ns, brick->name[bid], brick->name[bid],
          brick->name[bid], nobs[k]);
      if (nchar >= BRICKMASK_MAX_FILENAME) {
        P_ERR("the brickmask filename is too long\n"
            "Please enlarge `BRICKMASK_MAX_FILENAME` in `src/define.h`\n");
        mask_destroy(mask);
        BRICKMASK_QUIT(BRICKMASK_ERR_FILE);
      }
      if (assign_mask_from_file(fname, mask, data->ra + imin, data->dec + imin,
          data->nexp[k] + imin, mask->enull, mask->etype + k, imax - imin)) {
        mask_destroy(mask);
        BRICKMASK_QUIT(BRICKMASK_ERR_MASK);
      }
    }

    imin = imax;

    /* Print the reading progress. */
#ifdef MPI
    if (verbose && rank == BRICKMASK_MPI_ROOT)
#else
    if (verbose)
#endif
    {
      if (++cnt >= next) {
        next += step;
        printf("\x1B[%dD%*zu / %zu", wcol, ndg, cnt, data->nbrick);
        fflush(stdout);
      }
    }
  }

#ifdef MPI
  if (verbose && rank == BRICKMASK_MPI_ROOT)
#else
  if (verbose)
#endif
    printf("\x1B[%dD%*zu / %zu\n", wcol, ndg, cnt, data->nbrick);

  if (data->mtype < mask->dtype) data->mtype = mask->dtype;
  for (int k = 0; k < 3; k++)
    if (data->etype[k] < mask->etype[k]) data->etype[k] = mask->etype[k];

  mask_destroy(mask);
#ifdef MPI
  if (rank == BRICKMASK_MPI_ROOT)
#endif
  printf(FMT_DONE);
  return 0;
}
