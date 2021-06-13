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
  data->idx = NULL;
  data->id = NULL;
  data->mask = NULL;
  data->subid = NULL;
  data->content = NULL;

  /* Check the data type of masks required by `MASKBIT_NULL`. */
  uint64_t mnull = conf->mnull;
  if (mnull < UINT8_MAX) data->mtype = TBYTE;
  else if (mnull < UINT16_MAX) data->mtype = TSHORT;
  else if (mnull < UINT32_MAX) data->mtype = TINT;
  else data->mtype = TLONG;

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
  printf("Reading objects from the input catalog ...");
  if (!conf) {
    P_ERR("configuration parameters are not loaded\n");
    return NULL;
  }
  if (conf->verbose) printf("\n  Filename: `%s'\n", conf->input);
  fflush(stdout);

  /* Initialise the structure for the input catalogue. */
  DATA *data = data_init(conf);
  if (!data) return NULL;

  /* Read the input catalogue. */
  switch (data->fmt) {
    case BRICKMASK_FFMT_ASCII:
      if (read_ascii(conf, data)) {
        data_destroy(data);
        return NULL;
      }
      break;
    case BRICKMASK_FFMT_FITS:
      if (read_fits(conf, data)) {
        data_destroy(data);
        return NULL;
      }
      break;
  }

#ifdef MPI
  if (data->n > BRICKMASK_MAX_DATA) {
    P_ERR("too many objects in the catalog: %zu\n", data->n);
    data_destroy(data);
    return NULL;
  }
#endif

  /* Allocate memory for the rest of the properties. */
  if (!(data->id = malloc(data->n * sizeof(long))) ||
      !(data->mask = calloc(data->n, sizeof(uint64_t))) ||
      (conf->subid && !(data->subid = calloc(data->n, sizeof(uint8_t))))) {
    P_ERR("failed to allocate memory for additional columns of the data\n");
    data_destroy(data);
    return NULL;
  }

  if (conf->verbose) printf("  %zu objects are read from the file\n", data->n);
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
  if (conf->verbose) printf("\n  Filename: `%s'\n", conf->output);
  fflush(stdout);

  switch (data->fmt) {
    case BRICKMASK_FFMT_ASCII:
      if (save_ascii(conf, data)) return BRICKMASK_ERR_SAVE;
      break;
    case BRICKMASK_FFMT_FITS:
      if (save_fits(conf, data)) return BRICKMASK_ERR_SAVE;
      break;
  }

  data_destroy(data);
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
  if (data->id) free(data->id);
  if (data->mask) free(data->mask);
  if (data->subid) free(data->subid);
  if (data->content) free(data->content);
  free(data);
}
