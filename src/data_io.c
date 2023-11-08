/*******************************************************************************
* data_io.c: this file is part of the brickmask program.

* brickmask: assign bit codes defined on Legacy Survey brick pixels
             to a catalogue with sky coordinates.

* Github repository:
        https://github.com/cheng-zhao/brickmask

* Copyright (c) 2020 -- 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "data_io.h"
#include "read_file.h"
#include "save_file.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <fitsio.h>

/*============================================================================*\
                   Functions for reading the input catalogue
\*============================================================================*/

/******************************************************************************
Function `data_init`:
  Initialise the structure for the input data.
Arguments:
  * `conf`:     structure for storing configurations.
Return:
  Address of the structure for the input catalogue on success; NULL on error.
******************************************************************************/
static DATA *data_init(const CONF *conf) {
  DATA *data = calloc(1, sizeof(DATA));
  if (!data) {
    P_ERR("failed to allocate memory for the input data catalog\n");
    return NULL;
  }

  data->fmt = conf->ftype;
  data->ra = data->dec = NULL;
  data->idx = data->cidx = data->iidx = NULL;
  data->id = NULL;
  data->mask = NULL;
  data->content = NULL;

  if (!(data->iidx = calloc(conf->ncat + 1, sizeof(size_t)))) {
    P_ERR("failed to allocate memory for the input data catalog\n");
    data_destroy(data);
    return NULL;
  }
  if (data->fmt == BRICKMASK_FFMT_ASCII) {
    data->nmax = BRICKMASK_DATA_INIT_NUM;
    data->cmax = BRICKMASK_CONTENT_INIT_SIZE;
    if (!(data->ra = malloc(data->nmax * sizeof(double))) ||
        !(data->dec = malloc(data->nmax * sizeof(double))) ||
        !(data->cidx = malloc(data->nmax * sizeof(size_t))) ||
        !(data->content = malloc(data->cmax))) {
      P_ERR("failed to allocate memory for the input data catalog\n");
      data_destroy(data);
      return NULL;
    }
  }

  /* Check the data type of masks required by `MASKBIT_NULL` and `NEXP_NULL`. */
  uint64_t mnull = conf->mnull;
  if (mnull < UINT8_MAX) data->mtype = TBYTE;
  else if (mnull < UINT16_MAX) data->mtype = TSHORT;
  else if (mnull < UINT32_MAX) data->mtype = TINT;
  else data->mtype = TLONG;

  uint64_t enull = conf->enull;
  if (enull < UINT8_MAX)
    data->etype[0] = data->etype[1] = data->etype[2] = TBYTE;
  else if (enull < UINT16_MAX)
    data->etype[0] = data->etype[1] = data->etype[2] = TSHORT;
  else if (enull < UINT32_MAX)
    data->etype[0] = data->etype[1] = data->etype[2] = TINT;
  else
    data->etype[0] = data->etype[1] = data->etype[2] = TLONG;

  return data;
}

/******************************************************************************
Function `read_data`:
  Read data from the input catalogue.
Arguments:
  * `conf`:     structure for storing configurations.
Return:
  Address of the structure for the input catalogue on success; NULL on error.
******************************************************************************/
DATA *read_data(const CONF *conf) {
  printf("Reading objects from the input catalogs ...");
  if (!conf) {
    P_ERR("configuration parameters are not loaded\n");
    return NULL;
  }
  if (conf->verbose)
    printf("\n  Input catalogs obtained from: `%s'\n", conf->ilist);
  fflush(stdout);

  /* Initialise the structure for the input catalogue. */
  DATA *data = data_init(conf);
  if (!data) return NULL;
  size_t ndata = 0;

  /* Read the input catalogue. */
  switch (data->fmt) {
    case BRICKMASK_FFMT_ASCII:
      for (int i = 0; i < conf->ncat; i++) {
        if (read_ascii(conf->input[i], conf, data)) {
          data_destroy(data);
          return NULL;
        }
        data->iidx[i + 1] = data->n;
        if (conf->verbose) {
          printf("  %zu objects read from `%s'\n", data->n - ndata,
              conf->input[i]);
          ndata = data->n;
        }
      }
      break;
    case BRICKMASK_FFMT_FITS:
      for (int i = 0; i < conf->ncat; i++) {
        if (read_fits(conf->input[i], conf, data)) {
          data_destroy(data);
          return NULL;
        }
        data->iidx[i + 1] = data->n;
        if (conf->verbose) {
          printf("  %zu objects read from `%s'\n", data->n - ndata,
              conf->input[i]);
          ndata = data->n;
        }
      }
      break;
  }

  if (!data->n) {
    P_ERR("no valid object read from the input catalogs\n");
    data_destroy(data);
    return NULL;
  }

  /* Reduce the memory cost if applicable. */
  if (data->fmt == BRICKMASK_FFMT_ASCII) {
    double *dtmp = realloc(data->ra, data->n * sizeof(double));
    if (dtmp) data->ra = dtmp;
    dtmp = realloc(data->dec, data->n * sizeof(double));
    if (dtmp) data->dec = dtmp;
    size_t *stmp = realloc(data->cidx, data->n * sizeof(size_t));
    if (stmp) data->cidx = stmp;
    char *ctmp = realloc(data->content, data->csize);
    if (ctmp) data->content = ctmp;
  }

#ifdef MPI
  if (data->n > BRICKMASK_MAX_DATA) {
    P_ERR("too many objects in the catalogs: %zu\n", data->n);
    data_destroy(data);
    return NULL;
  }
#endif

  /* Allocate memory for the rest of the properties. */
  if (!(data->idx = malloc(data->n * sizeof(size_t))) ||
      !(data->id = malloc(data->n * sizeof(long))) ||
      !(data->mask = calloc(data->n, sizeof(uint64_t))) ||
      !(data->nexp[0] = calloc(data->n, sizeof(uint64_t))) ||
      !(data->nexp[1] = calloc(data->n, sizeof(uint64_t))) ||
      !(data->nexp[2] = calloc(data->n, sizeof(uint64_t)))) {
    P_ERR("failed to allocate memory for additional columns of the data\n");
    data_destroy(data);
    return NULL;
  }

  for (size_t i = 0; i < data->n; i++) data->idx[i] = i;

  if (conf->verbose)
    printf("  %zu objects read in total from %d files\n", data->n, conf->ncat);
  printf(FMT_DONE);
  return data;
}


/*============================================================================*\
                    Function for saving the output catalogue
\*============================================================================*/

/******************************************************************************
Function `save_data`:
  Save data to the output catalogue.
Arguments:
  * `conf`:     structure for storing configurations;
  * `data`:     structure for the data catalogue.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int save_data(const CONF *conf, DATA *data) {
  printf("Saving objects with maskbits to the output catalog ...");
  if (!conf) {
    P_ERR("configuration parameters are not loaded\n");
    return BRICKMASK_ERR_INIT;
  }
  if (!data) {
    P_ERR("the data catalogue is not initialised\n");
    return BRICKMASK_ERR_INIT;
  }
  if (conf->verbose)
    printf("\n  Output catalogs obtained from: `%s'\n", conf->olist);
  fflush(stdout);

  switch (data->fmt) {
    case BRICKMASK_FFMT_ASCII:
      for (int i = 0; i < conf->ncat; i++) {
        if (save_ascii(conf, data, i)) return BRICKMASK_ERR_SAVE;
        if (conf->verbose) {
          printf("  %zu objects saved to `%s'\n",
              data->iidx[i + 1] - data->iidx[i], conf->output[i]);
        }
      }
      break;
    case BRICKMASK_FFMT_FITS:
      for (int i = 0; i < conf->ncat; i++) {
        if (save_fits(conf, data, i)) return BRICKMASK_ERR_SAVE;
        if (conf->verbose) {
          printf("  %zu objects saved to `%s'\n",
              data->iidx[i + 1] - data->iidx[i], conf->output[i]);
        }
      }
      break;
  }

  data_destroy(data);

  if (conf->verbose) printf("  %d files saved\n", conf->ncat);
  printf(FMT_DONE);
  return 0;
}

/*============================================================================*\
                              Function for cleanup
\*============================================================================*/

/******************************************************************************
Function `data_destroy`:
  Deconstruct the structure for the input data.
Arguments:
  * `data`:     structure for the input data catalogue.
******************************************************************************/
void data_destroy(DATA *data) {
  if (!data) return;
  if (data->ra) free(data->ra);
  if (data->dec) free(data->dec);
  if (data->idx) free(data->idx);
  if (data->cidx) free(data->cidx);
  if (data->iidx) free(data->iidx);
  if (data->id) free(data->id);
  if (data->mask) free(data->mask);
  for (int i = 0; i < 3; i++) {
    if (data->nexp[i]) free(data->nexp[i]);
  }
  if (data->content) free(data->content);
  free(data);
}
