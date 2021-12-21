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
  * `mnull`:    bit code for objects outside bricks.
Return:
  Address of the structure for maskbits on success; NULL on error.
******************************************************************************/
static inline MASK *mask_init(uint32_t mnull) {
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

  return mask;
}

/******************************************************************************
Function `get_maskbit_fname`:
  Find maskbit files containing a given brick name.
Arguments:
  * `brick`:    structure for bricks;
  * `bname`:    name of the brick to be found;
  * `fname`:    pointers to maskbit filenames that are found;
  * `nsp`:      number of subsamples containing the brick name.
******************************************************************************/
static void get_maskbit_fname(const BRICK *brick, const char *bname,
    char **fname, unsigned char *subid, int *nsp) {
  int n = 0;
  for (int i = 0; i < brick->nsp; i++) {
    for (size_t j = 0; j < brick->nmask[i]; j++) {
      char *fcheck = brick->fmask[i][j];
      if (!(*fcheck)) continue;         /* the file has been visited */
      if (strstr(fcheck, bname)) {      /* the file is found */
        if (brick->subid) subid[n] = brick->subid[i];
        fname[n++] = fcheck;
        break;                          /* stop searching this subsample */
      }
    }
  }
  *nsp = n;
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

  theta = (theta >= 1) ? 0 : sqrt(1 - theta * theta) / theta * RAD_2_DEGREE;
  double fac = theta / sqrt(phi1 * phi1 + phi2 * phi2);
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

  /* Allocate memory for pointers to maskbit filenames, and subsample IDs. */
  char **fname = malloc(brick->nsp * sizeof(char *));
  if (!fname) {
    P_ERR("failed to allocate memory for maskbit file pointers\n");
    BRICKMASK_QUIT(BRICKMASK_ERR_MEMORY);
  }
  unsigned char *subid = calloc(brick->nsp, sizeof(unsigned char));
  if (!subid) {
    P_ERR("failed to allocate memory for subsample IDs\n");
    free(fname);
    BRICKMASK_QUIT(BRICKMASK_ERR_MEMORY);
  }
  int nsp = 0;

  /* Initialise maskbits. */
  MASK *mask = mask_init(brick->mnull);
  if (!mask) {
    free(fname); free(subid);
    BRICKMASK_QUIT(BRICKMASK_ERR_MEMORY);
  }

  /* Read and assign maskbits. */
  bool has_null = false;
  size_t imin, imax;
  imin = 0;
  while (imin < data->n) {
    /* Get the index range of data sharing the same brick. */
    for (imax = imin + 1; imax < data->n; imax++) {
      if (data->id[imax] != data->id[imin]) break;
    }
    size_t bid = data->id[imin];        /* ID of the corresponding brick. */

    /* Erase maskbit filenames for the previous brick. */
    for (int i = 0; i < nsp; i++) fname[i][0] = '\0';
    /* Get maskbit filenames corresponding to this brick. */
    get_maskbit_fname(brick, brick->name[bid], fname, subid, &nsp);
    if (!nsp) {                 /* no maskbit file for this object */
      has_null = true;
      for (size_t i = imin; i < imax; i++) data->mask[i] = mask->mnull;
      imin = imax;
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
      continue;
    }

    for (int i = 0; i < nsp; i++) {
      /* Read maskbits for each subsample. */
      if (read_mask(fname[i], mask)) {
        free(fname); free(subid); mask_destroy(mask);
        BRICKMASK_QUIT(BRICKMASK_ERR_MASK);
      }

      /* Choose the maskbit code assigning function given the data type. */
      int (*assign_bitcode_func) (const MASK *, DATA *, const size_t,
          const size_t, const uint8_t) = NULL;
      switch (mask->dtype) {
        case TBYTE:  assign_bitcode_func = assign_bitcode_uint8_t;  break;
        case TSHORT: assign_bitcode_func = assign_bitcode_uint16_t; break;
        case TINT:   assign_bitcode_func = assign_bitcode_uint32_t; break;
        case TLONG:  assign_bitcode_func = assign_bitcode_uint64_t; break;
        default:
          P_ERR("unexpected data type for maskbits: %d\n", mask->dtype);
          free(fname); free(subid); mask_destroy(mask);
          BRICKMASK_QUIT(BRICKMASK_ERR_MASK);
      }

      /* Assign maskbits. */
      assign_bitcode_func(mask, data, imin, imax, subid[i]);
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

  if (!has_null) data->mtype = mask->dtype;
  else if (data->mtype < mask->dtype) data->mtype = mask->dtype;

  free(fname);
  free(subid);
  mask_destroy(mask);
#ifdef MPI
  if (rank == BRICKMASK_MPI_ROOT)
#endif
  printf(FMT_DONE);
  return 0;
}
