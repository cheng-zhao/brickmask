/*******************************************************************************
* sort_data.c: this file is part of the brickmask program.

* brickmask: assign bit codes defined on Legacy Survey brick pixels
             to a catalogue with sky coordinates.

* Github repository:
        https://github.com/cheng-zhao/brickmask

* Copyright (c) 2020 -- 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <fitsio.h>
#include "define.h"
#include "get_brick.h"
#include "data_io.h"

/*============================================================================*\
                    Functions for finding bricks of the data
\*============================================================================*/

/******************************************************************************
Function `compare_pos`:
  Compare a coordinate with a range (first Dec, and then RA).
Arguments:
  * `ra1`:      lower limit of the RA range;
  * `ra2`:      upper limit of the RA range;
  * `dec1`:     lower limit of the Dec range;
  * `dec2`:     upper limit of the Dec range;
  * `ra`:       RA to be examined;
  * `dec`:      Dec to be examined.
Return:
  Zero if the coordinate is inside the range; +1/-1 for above/below the range.
******************************************************************************/
static inline int compare_pos(const double ra1, const double ra2,
    const double dec1, const double dec2, const double ra, const double dec) {
  if (dec < dec1) return -1;
  if (dec >= dec2) return 1;
  if (ra < ra1) return -1;
  if (ra >= ra2) return 1;
  return 0;
}

/******************************************************************************
Function `find_brick`:
  Find the index of the brick that contains a coordinate with binary search.
Arguments:
  * `brick`:    structure for bricks;
  * `ra`:       RA to be examined;
  * `dec`:      Dec to be examined.
Return:
  Index of the brick on success; -1 on error.
******************************************************************************/
static long find_brick(const BRICK *brick, const double ra,
    const double dec) {
  size_t l, u;
  l = 0;
  u = brick->n - 1;

  while (l <= u) {
    size_t i = (l + u) >> 1;
    if (compare_pos(brick->ra1[i], brick->ra2[i], brick->dec1[i],
        brick->dec2[i], ra, dec) > 0) l = i + 1;
    else if (compare_pos(brick->ra1[i], brick->ra2[i], brick->dec1[i],
        brick->dec2[i], ra, dec) < 0) u = i - 1;
    else return i;
  }
  return -1;
}

/******************************************************************************
Function `get_brick_id`:
  Get brick IDs for the input data sample.
Arguments:
  * `brick`:    structure for bricks;
  * `data`:     structure for the data catalogue.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static inline int get_brick_id(BRICK *brick, DATA *data) {
  for (size_t i = 0; i < data->n; i++) {
    if (data->ra[i] == 360) data->ra[i] -= BRICKMASK_TOL;
    if (data->dec[i] == 90) data->dec[i] -= BRICKMASK_TOL;
    data->id[i] = find_brick(brick, data->ra[i], data->dec[i]);
    if (data->id[i] < 0) {
      P_ERR("cannot find the brick for coordinate (%g, %g)\n",
          data->ra[i], data->dec[i]);
      return BRICKMASK_ERR_BRICK;
    }
  }
  free(brick->ra1);
  free(brick->ra2);
  free(brick->dec1);
  free(brick->dec2);
  brick->ra1 = brick->ra2 = brick->dec1 = brick->dec2 = NULL;
  return 0;
}

/*============================================================================*\
                  Definitions for sorting the data by brick ID
\*============================================================================*/

/* Structure for storing the properties of a data object. */
typedef struct {
  double ra;
  double dec;
  size_t idx;
  long id;
} DATA_ELEMENT;

/* Data type for the two arrays. */
#define TIMSORT_DTYPE1                  long
#define TIMSORT_DTYPE2                  DATA
/* Data type of the variable for binding elements from both arrays. */
#define TIMSORT_BIND_DTYPE              DATA_ELEMENT
/* Macro for changing the starting indices of the arrays. */
#define TIMSORT_OFFSET(x,y,s)                                           \
  (x) += (s); (y)->ra += (s); (y)->dec += (s); (y)->idx += (s);
/* Macro for reseting the starting indices of the arrays. */
#define TIMSORT_RESET(x,y,s)                                            \
  (y)->ra -= (s); (y)->dec -= (s); (y)->idx -= (s);
/* Macro for comparing array elements with indices i and j. */
#define TIMSORT_CMP_IDX(x,y,i,j)        ((x)[(i)] - (x)[(j)])
/* Macro for comparing binded value with the arrays with index i. */
#define TIMSORT_CMP_BIND(x,y,i,b)       ((x)[(i)] - (b).id)
/* Macro for assigning values with index i to those with index j. */
#define TIMSORT_ASSIGN_IDX(x,y,i,j)                                     \
  (x)[(j)] = (x)[(i)];           (y)->ra[(j)] = (y)->ra[(i)];           \
  (y)->dec[(j)] = (y)->dec[(i)]; (y)->idx[(j)] = (y)->idx[(i)];
/* Macro for assigning array elements with index i, to the binded variable. */
#define TIMSORT_ASSIGN_BIND(x,y,i,b)                                    \
  (b).id = (x)[(i)];       (b).ra = (y)->ra[(i)];                       \
  (b).dec = (y)->dec[(i)]; (b).idx = (y)->idx[(i)];
/* Macro for assigning binded value to array elements with index i. */
#define TIMSORT_GET_BIND(x,y,i,b)                                       \
  (x)[(i)] = (b).id;       (y)->ra[(i)] = (b).ra;                       \
  (y)->dec[(i)] = (b).dec; (y)->idx[(i)] = (b).idx;
/* Macro for swapping array elements with indices i and j. */
#define TIMSORT_SWAP(x,y,i,j) {                                         \
  DATA_ELEMENT _tmp;                                                    \
  TIMSORT_ASSIGN_BIND(x,y,i,_tmp);                                      \
  (x)[(i)] = (x)[(j)];           (y)->ra[(i)] = (y)->ra[(j)];           \
  (y)->dec[(i)] = (y)->dec[(j)]; (y)->idx[(i)] = (y)->idx[(j)];         \
  TIMSORT_GET_BIND(x,y,j,_tmp);                                         \
}

#include "timsort.c"


/*============================================================================*\
                 Functions for restoring the order of the data
\*============================================================================*/

/******************************************************************************
Function `reorder_mask_reduce`:
  Reorder maskbits, and reduce the length of the data type if applicable.
Arguments:
  * `code`:     array for maskbits;
  * `dtype`:    datatype of maskbits;
  * `idx`:      original index of the data points;
  * `ndata`:    number of data points.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static inline int reorder_mask_reduce(uint64_t **code, const int dtype,
    const size_t *idx, const size_t ndata) {
  void *mask;
  switch (dtype) {
    case TBYTE:
      if (!(mask = malloc(ndata * sizeof(unsigned char)))) {
        P_ERR("failed to allocate memory for saving maskbits\n");
        return BRICKMASK_ERR_MEMORY;
      }
      for (size_t i = 0; i < ndata; i++)
        ((unsigned char *) mask)[idx[i]] = (*code)[i];
      break;
    case TSHORT:
      if (!(mask = malloc(ndata * sizeof(uint16_t)))) {
        P_ERR("failed to allocate memory for saving maskbits\n");
        return BRICKMASK_ERR_MEMORY;
      }
      for (size_t i = 0; i < ndata; i++)
        ((uint16_t *) mask)[idx[i]] = (*code)[i];
      break;
    case TINT:
      if (!(mask = malloc(ndata * sizeof(uint32_t)))) {
        P_ERR("failed to allocate memory for saving maskbits\n");
        return BRICKMASK_ERR_MEMORY;
      }
      for (size_t i = 0; i < ndata; i++)
        ((uint32_t *) mask)[idx[i]] = (*code)[i];
      break;
    case TLONG:
      if (!(mask = malloc(ndata * sizeof(uint64_t)))) {
        P_ERR("failed to allocate memory for saving maskbits\n");
        return BRICKMASK_ERR_MEMORY;
      }
      for (size_t i = 0; i < ndata; i++)
        ((uint64_t *) mask)[idx[i]] = (*code)[i];
      break;
    default:
      P_ERR("unexpected data type for maskbits: %d\n", dtype);
      return BRICKMASK_ERR_UNKNOWN;
  }

  free(*code);
  *code = mask;
  return 0;
}

/******************************************************************************
Function `reorder_mask`:
  Restore the original maskbits before data sorting.
Arguments:
  * `code`:     array for maskbits;
  * `idx`:      original index of the data points;
  * `ndata`:    number of data points.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static inline int reorder_mask(uint64_t **code, const size_t *idx,
    const size_t ndata) {
  uint64_t *mask = malloc(ndata * sizeof(uint64_t));
  if (!mask) {
    P_ERR("failed to allocate memory for saving maskbits\n");
    return BRICKMASK_ERR_MEMORY;
  }

  for (size_t i = 0; i < ndata; i++) mask[idx[i]] = (*code)[i];
  free(*code);
  *code = mask;
  return 0;
}


/*============================================================================*\
                         Interface for sorting the data
\*============================================================================*/

/******************************************************************************
Function `sort_data`:
  Sort the input data sample based on the brick IDs.
Arguments:
  * `brick`:    structure for bricks;
  * `data`:     structure for the data catalogue;
  * `verbose`:  indicate whether to show detailed outputs.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int sort_data(BRICK *brick, DATA *data, const bool verbose) {
  printf("Sorting the input data based on brick IDs ...");
  if (!brick || !data) {
    P_ERR("the bricks or input data is not initialised\n");
    return BRICKMASK_ERR_INIT;
  }
  if (verbose) printf("\n");
  fflush(stdout);

  /* Get brick ID and sort the data. */
  if (get_brick_id(brick, data)) return BRICKMASK_ERR_BRICK;
  tim_sort(data->id, data, data->n);

  /* Count the total number of bricks for the data. */
  long prev = -1;
  for (size_t i = 0; i < data->n; i++) {
    if (data->id[i] != prev) {
      data->nbrick++;
      prev = data->id[i];
    }
  }
  if (verbose) printf("  %zu bricks contain data points\n", data->nbrick);

  printf(FMT_DONE);
  return 0;
}

/******************************************************************************
Function `reorder_data`:
  Restore the original order of the input data sample.
Arguments:
  * `data`:     structure for the data catalogue.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int reorder_data(DATA *data) {
  printf("Restoring the original order of the input data ...");
  if (!data) {
    P_ERR("the input data is not initialised\n");
    return BRICKMASK_ERR_INIT;
  }
  fflush(stdout);

  switch (data->fmt) {
    case BRICKMASK_FFMT_ASCII:
      if (reorder_mask(&(data->mask), data->idx, data->n))
        return BRICKMASK_ERR_MEMORY;
      if (reorder_mask(data->nexp, data->idx, data->n))
        return BRICKMASK_ERR_MEMORY;
      if (reorder_mask(data->nexp + 1, data->idx, data->n))
        return BRICKMASK_ERR_MEMORY;
      if (reorder_mask(data->nexp + 2, data->idx, data->n))
        return BRICKMASK_ERR_MEMORY;
      break;
    case BRICKMASK_FFMT_FITS:
      if (reorder_mask_reduce(&(data->mask), data->mtype, data->idx, data->n))
        return BRICKMASK_ERR_MEMORY;
      if (reorder_mask_reduce(data->nexp, data->etype[0], data->idx, data->n))
        return BRICKMASK_ERR_MEMORY;
      if (reorder_mask_reduce(data->nexp + 1, data->etype[1], data->idx,
          data->n)) return BRICKMASK_ERR_MEMORY;
      if (reorder_mask_reduce(data->nexp + 2, data->etype[2], data->idx,
          data->n)) return BRICKMASK_ERR_MEMORY;
      break;
  }

  printf(FMT_DONE);
  return 0;
}
