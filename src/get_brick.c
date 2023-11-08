/*******************************************************************************
* get_brick.c: this file is part of the brickmask program.

* brickmask: assign bit codes defined on Legacy Survey brick pixels
             to a catalogue with sky coordinates.

* Github repository:
        https://github.com/cheng-zhao/brickmask

* Copyright (c) 2020 -- 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "get_brick.h"
#include "read_file.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>

/*============================================================================*\
                        Functions for setting up bricks
\*============================================================================*/

/******************************************************************************
Function `brick_init`:
  Initialise the structure for bricks.
Arguments:
  * `conf`:     structure for storing configurations.
Return:
  Address of the structure for bricks on success; NULL on error.
******************************************************************************/
static BRICK *brick_init(const CONF *conf) {
  BRICK *brick = calloc(1, sizeof(BRICK));
  if (!brick) {
    P_ERR("failed to allocate memory for bricks\n");
    return NULL;
  }

  brick->ra1 = brick->ra2 = brick->dec1 = brick->dec2 = NULL;
  brick->name = NULL;
  brick->photsys = NULL;
  brick->mnull = conf->mnull;
  brick->enull = conf->enull;

  size_t len = strlen(conf->bpath) + 1;
  if (!(brick->bpath = calloc(len, sizeof(char)))) {
    P_ERR("failed to allocate memory for bricks\n");
    brick_destroy(brick);
    return NULL;
  }
  memcpy(brick->bpath, conf->bpath, len);

  return brick;
}

/******************************************************************************
Function `check_brick`:
  Adjust edges of bricks to eliminate gaps, and check orders of coordinates.
Arguments:
  * `brick`:    structure for bricks.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static inline int check_brick(BRICK *brick) {
#ifdef MPI
  if (brick->n > BRICKMASK_MAX_BRICK ||
      brick->n * (brick->nlen + 1) > BRICKMASK_MAX_BRICK) {
    P_ERR("there are too many bricks: %zu\n", brick->n);
    return BRICKMASK_ERR_BRICK;
  }
#else
  if (brick->n > LONG_MAX) {
    P_ERR("there are too many bricks: %zu\n", brick->n);
    return BRICKMASK_ERR_BRICK;
  }
#endif
  for (size_t i = 0; i < brick->n; i++) {
    brick->ra1[i] = round(brick->ra1[i] / BRICKMASK_TOL) * BRICKMASK_TOL;
    brick->ra2[i] = round(brick->ra2[i] / BRICKMASK_TOL) * BRICKMASK_TOL;
    brick->dec1[i] = round(brick->dec1[i] / BRICKMASK_TOL) * BRICKMASK_TOL;
    brick->dec2[i] = round(brick->dec2[i] / BRICKMASK_TOL) * BRICKMASK_TOL;
  }

  /* Both RA and Dec should be in ascending order. */
  for (size_t i = 1; i < brick->n; i++) {
    if (brick->dec1[i - 1] > brick->dec1[i]) {
      P_ERR("invalid declination order in the brick list file\n");
      return BRICKMASK_ERR_BRICK;
    }
    else if (brick->dec1[i - 1] == brick->dec1[i] &&
        brick->ra1[i - 1] > brick->ra1[i]) {
      P_ERR("invalid right ascension order in the brick list file\n");
      return BRICKMASK_ERR_BRICK;
    }
  }
  return 0;
}

/******************************************************************************
Function `get_brick`:
  Get brick information from files.
Arguments:
  * `conf`:     structure for storing configurations.
Return:
  Address of the structure for bricks on success; NULL on error.
******************************************************************************/
BRICK *get_brick(const CONF *conf) {
  printf("Getting information of bricks ...");
  if (!conf) {
    P_ERR("configuration parameters are not loaded\n");
    return NULL;
  }
  if (conf->verbose) printf("\n");
  fflush(stdout);

  /* Initialise the structure for bricks. */
  BRICK *brick = brick_init(conf);
  if (!brick) return NULL;

  /* Read coordinate ranges and brick names from file. */
  if (read_brick(conf->flist, brick) || check_brick(brick)) {
    brick_destroy(brick);
    return NULL;
  }
  if (conf->verbose)
    printf("  Brick information is loaded from file: `%s'\n", conf->flist);

  printf(FMT_DONE);
  return brick;
}

/******************************************************************************
Function `brick_destroy`:
  Deconstruct the structure for bricks.
Arguments:
  * `data`:     structure for bricks.
******************************************************************************/
void brick_destroy(BRICK *brick) {
  if (!brick) return;
  if (brick->ra1) free(brick->ra1);
  if (brick->ra2) free(brick->ra2);
  if (brick->dec1) free(brick->dec1);
  if (brick->dec2) free(brick->dec2);
  if (brick->name) {
    if (*(brick->name)) free(*(brick->name));
    free(brick->name);
  }
  if (brick->photsys) free(brick->photsys);
  if (brick->bpath) free(brick->bpath);
  free(brick);
}
