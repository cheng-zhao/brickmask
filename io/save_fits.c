/*******************************************************************************
* save_fits.c: this file is part of the brickmask program.

* brickmask: assign bit codes defined on Legacy Survey brick pixels
             to a catalogue with sky coordinates.

* Github repository:
        https://github.com/cheng-zhao/brickmask

* Copyright (c) 2020 -- 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "read_file.h"
#include <fitsio.h>
#include <stdio.h>
#include <string.h>

#define FITS_ABORT_SINGLE {                                     \
  P_ERR("cfitsio error: ");                                     \
  fits_report_error(stderr, status);                            \
  status = 0;                                                   \
  if (fp) fits_close_file(fp, &status);                         \
  return BRICKMASK_ERR_FILE;                                    \
}

#define FITS_ABORT {                                            \
  P_ERR("cfitsio error: ");                                     \
  fits_report_error(stderr, status);                            \
  status = 0;                                                   \
  if (ifp) { fits_close_file(ifp, &status); status = 0; }       \
  if (fp) fits_close_file(fp, &status);                         \
  return BRICKMASK_ERR_FILE;                                    \
}

#define FITS_ABORT_MEM {                                        \
  P_ERR("cfitsio error: ");                                     \
  fits_report_error(stderr, status);                            \
  status = 0;                                                   \
  if (ifp) { fits_close_file(ifp, &status); status = 0; }       \
  if (fp) fits_close_file(fp, &status);                         \
  if (memptr) free(memptr);                                     \
  if (chunk) free(chunk);                                       \
  return BRICKMASK_ERR_FILE;                                    \
}

/*============================================================================*\
                        Functions for saving FITS files
\*============================================================================*/

/******************************************************************************
Function `reorder_mask`:
  Reorder maskbits, and reduce the length of the data type if applicable.
Arguments:
  * `data`:     structure for the the data catalogue.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static inline int reorder_mask(DATA *data) {
  void *mask;
  switch (data->mtype) {
    case TBYTE:
      if (!(mask = malloc(data->n * sizeof(unsigned char)))) {
        P_ERR("failed to allocate memory for saving maskbits\n");
        return BRICKMASK_ERR_MEMORY;
      }
      for (size_t i = 0; i < data->n; i++)
        ((unsigned char *) mask)[data->idx[i]] = data->mask[i];
      break;
    case TSHORT:
      if (!(mask = malloc(data->n * sizeof(uint16_t)))) {
        P_ERR("failed to allocate memory for saving maskbits\n");
        return BRICKMASK_ERR_MEMORY;
      }
      for (size_t i = 0; i < data->n; i++)
        ((uint16_t *) mask)[data->idx[i]] = data->mask[i];
      break;
    case TINT:
      if (!(mask = malloc(data->n * sizeof(uint32_t)))) {
        P_ERR("failed to allocate memory for saving maskbits\n");
        return BRICKMASK_ERR_MEMORY;
      }
      for (size_t i = 0; i < data->n; i++)
        ((uint32_t *) mask)[data->idx[i]] = data->mask[i];
      break;
    case TLONG:
      if (!(mask = malloc(data->n * sizeof(uint64_t)))) {
        P_ERR("failed to allocate memory for saving maskbits\n");
        return BRICKMASK_ERR_MEMORY;
      }
      for (size_t i = 0; i < data->n; i++)
        ((uint64_t *) mask)[data->idx[i]] = data->mask[i];
      break;
    default:
      P_ERR("unexpected data type for maskbits: %d\n", data->mtype);
      return BRICKMASK_ERR_UNKNOWN;
  }

  free(data->mask);
  data->mask = mask;
  return 0;
}

/******************************************************************************
Function `reorder_subid`:
  Restore the original subsample ID before data sorting.
Arguments:
  * `data`:     structure for the the data catalogue.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static inline int reorder_subid(DATA *data) {
  if (!data->subid) return 0;           /* subsample ID is not required */
  unsigned char *subid = malloc(data->n * sizeof(unsigned char));
  if (!subid) {
    P_ERR("failed to allocate memory for saving subsample IDs\n");
    return BRICKMASK_ERR_MEMORY;
  }

  for (size_t i = 0; i < data->n; i++) subid[data->idx[i]] = data->subid[i];
  free(data->subid);
  data->subid = subid;
  return 0;
}

/******************************************************************************
Function `force_output`:
  Generate the filename for force overwriting.
Arguments:
  * `fname`:    the original filename.
Return:
  Address of the string for force overwriting.
******************************************************************************/
static inline char *force_output(const char *fname) {
  size_t len = strlen(fname);   /* strings from libcfg are surely terminated */
  char *output = malloc((len + 2) * sizeof(char));
  if (!output) {
    P_ERR("failed to allocate memory for the output filename\n");
    return NULL;
  }
  memcpy(output + 1, fname, (len + 1) * sizeof(char));
  *output = '!';
  return output;
}


/*============================================================================*\
                  Template function for saving a FITS catalog
\*============================================================================*/

/******************************************************************************
Function `fits_write_<BRICKMASK_MASKBIT_DTYPE><FITS_WRITE_SUBID_NAME>
    <FITS_WRITE_ALLCOL_NAME>`:
  Save specific columns of a FITS file to another FITS file, and append
  additional columns.
Arguments:
  * `fname`:    filename for the output FITS file;
  * `conf`:     structure for storing configurations;
  * `data`:     structure for the data catalogue.
Return:
  Zero on success; non-zero on error.
******************************************************************************/

#ifdef BRICKMASK_WFITS_MTYPE
  #undef BRICKMASK_FITS_MTYPE
#endif
#ifdef BRICKMASK_WFITS_SUBID
  #undef BRICKMASK_FITS_SUBID
#endif
#ifdef BRICKMASK_WFITS_OVERWRITE
  #undef BRICKMASK_WFITS_OVERWRITE
#endif
#ifdef BRICKMASK_WFITS_ALLCOL
  #undef BRICKMASK_FITS_ALLCOL
#endif

#define BRICKMASK_WFITS_MTYPE           TBYTE
#define BRICKMASK_WFITS_SUBID           0
#define BRICKMASK_WFITS_OVERWRITE       0
#define BRICKMASK_WFITS_ALLCOL          0
#include "fits_write.c"

#define BRICKMASK_WFITS_MTYPE           TBYTE
#define BRICKMASK_WFITS_SUBID           0
#define BRICKMASK_WFITS_OVERWRITE       0
#define BRICKMASK_WFITS_ALLCOL          1
#include "fits_write.c"

#define BRICKMASK_WFITS_MTYPE           TBYTE
#define BRICKMASK_WFITS_SUBID           0
#define BRICKMASK_WFITS_OVERWRITE       1
#define BRICKMASK_WFITS_ALLCOL          0
#include "fits_write.c"

#define BRICKMASK_WFITS_MTYPE           TBYTE
#define BRICKMASK_WFITS_SUBID           0
#define BRICKMASK_WFITS_OVERWRITE       1
#define BRICKMASK_WFITS_ALLCOL          1
#include "fits_write.c"

#define BRICKMASK_WFITS_MTYPE           TBYTE
#define BRICKMASK_WFITS_SUBID           1
#define BRICKMASK_WFITS_OVERWRITE       0
#define BRICKMASK_WFITS_ALLCOL          0
#include "fits_write.c"

#define BRICKMASK_WFITS_MTYPE           TBYTE
#define BRICKMASK_WFITS_SUBID           1
#define BRICKMASK_WFITS_OVERWRITE       0
#define BRICKMASK_WFITS_ALLCOL          1
#include "fits_write.c"

#define BRICKMASK_WFITS_MTYPE           TBYTE
#define BRICKMASK_WFITS_SUBID           1
#define BRICKMASK_WFITS_OVERWRITE       1
#define BRICKMASK_WFITS_ALLCOL          0
#include "fits_write.c"

#define BRICKMASK_WFITS_MTYPE           TBYTE
#define BRICKMASK_WFITS_SUBID           1
#define BRICKMASK_WFITS_OVERWRITE       1
#define BRICKMASK_WFITS_ALLCOL          1
#include "fits_write.c"


#define BRICKMASK_WFITS_MTYPE           TSHORT
#define BRICKMASK_WFITS_SUBID           0
#define BRICKMASK_WFITS_OVERWRITE       0
#define BRICKMASK_WFITS_ALLCOL          0
#include "fits_write.c"

#define BRICKMASK_WFITS_MTYPE           TSHORT
#define BRICKMASK_WFITS_SUBID           0
#define BRICKMASK_WFITS_OVERWRITE       0
#define BRICKMASK_WFITS_ALLCOL          1
#include "fits_write.c"

#define BRICKMASK_WFITS_MTYPE           TSHORT
#define BRICKMASK_WFITS_SUBID           0
#define BRICKMASK_WFITS_OVERWRITE       1
#define BRICKMASK_WFITS_ALLCOL          0
#include "fits_write.c"

#define BRICKMASK_WFITS_MTYPE           TSHORT
#define BRICKMASK_WFITS_SUBID           0
#define BRICKMASK_WFITS_OVERWRITE       1
#define BRICKMASK_WFITS_ALLCOL          1
#include "fits_write.c"

#define BRICKMASK_WFITS_MTYPE           TSHORT
#define BRICKMASK_WFITS_SUBID           1
#define BRICKMASK_WFITS_OVERWRITE       0
#define BRICKMASK_WFITS_ALLCOL          0
#include "fits_write.c"

#define BRICKMASK_WFITS_MTYPE           TSHORT
#define BRICKMASK_WFITS_SUBID           1
#define BRICKMASK_WFITS_OVERWRITE       0
#define BRICKMASK_WFITS_ALLCOL          1
#include "fits_write.c"

#define BRICKMASK_WFITS_MTYPE           TSHORT
#define BRICKMASK_WFITS_SUBID           1
#define BRICKMASK_WFITS_OVERWRITE       1
#define BRICKMASK_WFITS_ALLCOL          0
#include "fits_write.c"

#define BRICKMASK_WFITS_MTYPE           TSHORT
#define BRICKMASK_WFITS_SUBID           1
#define BRICKMASK_WFITS_OVERWRITE       1
#define BRICKMASK_WFITS_ALLCOL          1
#include "fits_write.c"


#define BRICKMASK_WFITS_MTYPE           TINT
#define BRICKMASK_WFITS_SUBID           0
#define BRICKMASK_WFITS_OVERWRITE       0
#define BRICKMASK_WFITS_ALLCOL          0
#include "fits_write.c"

#define BRICKMASK_WFITS_MTYPE           TINT
#define BRICKMASK_WFITS_SUBID           0
#define BRICKMASK_WFITS_OVERWRITE       0
#define BRICKMASK_WFITS_ALLCOL          1
#include "fits_write.c"

#define BRICKMASK_WFITS_MTYPE           TINT
#define BRICKMASK_WFITS_SUBID           0
#define BRICKMASK_WFITS_OVERWRITE       1
#define BRICKMASK_WFITS_ALLCOL          0
#include "fits_write.c"

#define BRICKMASK_WFITS_MTYPE           TINT
#define BRICKMASK_WFITS_SUBID           0
#define BRICKMASK_WFITS_OVERWRITE       1
#define BRICKMASK_WFITS_ALLCOL          1
#include "fits_write.c"

#define BRICKMASK_WFITS_MTYPE           TINT
#define BRICKMASK_WFITS_SUBID           1
#define BRICKMASK_WFITS_OVERWRITE       0
#define BRICKMASK_WFITS_ALLCOL          0
#include "fits_write.c"

#define BRICKMASK_WFITS_MTYPE           TINT
#define BRICKMASK_WFITS_SUBID           1
#define BRICKMASK_WFITS_OVERWRITE       0
#define BRICKMASK_WFITS_ALLCOL          1
#include "fits_write.c"

#define BRICKMASK_WFITS_MTYPE           TINT
#define BRICKMASK_WFITS_SUBID           1
#define BRICKMASK_WFITS_OVERWRITE       1
#define BRICKMASK_WFITS_ALLCOL          0
#include "fits_write.c"

#define BRICKMASK_WFITS_MTYPE           TINT
#define BRICKMASK_WFITS_SUBID           1
#define BRICKMASK_WFITS_OVERWRITE       1
#define BRICKMASK_WFITS_ALLCOL          1
#include "fits_write.c"


#define BRICKMASK_WFITS_MTYPE           TLONG
#define BRICKMASK_WFITS_SUBID           0
#define BRICKMASK_WFITS_OVERWRITE       0
#define BRICKMASK_WFITS_ALLCOL          0
#include "fits_write.c"

#define BRICKMASK_WFITS_MTYPE           TLONG
#define BRICKMASK_WFITS_SUBID           0
#define BRICKMASK_WFITS_OVERWRITE       0
#define BRICKMASK_WFITS_ALLCOL          1
#include "fits_write.c"

#define BRICKMASK_WFITS_MTYPE           TLONG
#define BRICKMASK_WFITS_SUBID           0
#define BRICKMASK_WFITS_OVERWRITE       1
#define BRICKMASK_WFITS_ALLCOL          0
#include "fits_write.c"

#define BRICKMASK_WFITS_MTYPE           TLONG
#define BRICKMASK_WFITS_SUBID           0
#define BRICKMASK_WFITS_OVERWRITE       1
#define BRICKMASK_WFITS_ALLCOL          1
#include "fits_write.c"

#define BRICKMASK_WFITS_MTYPE           TLONG
#define BRICKMASK_WFITS_SUBID           1
#define BRICKMASK_WFITS_OVERWRITE       0
#define BRICKMASK_WFITS_ALLCOL          0
#include "fits_write.c"

#define BRICKMASK_WFITS_MTYPE           TLONG
#define BRICKMASK_WFITS_SUBID           1
#define BRICKMASK_WFITS_OVERWRITE       0
#define BRICKMASK_WFITS_ALLCOL          1
#include "fits_write.c"

#define BRICKMASK_WFITS_MTYPE           TLONG
#define BRICKMASK_WFITS_SUBID           1
#define BRICKMASK_WFITS_OVERWRITE       1
#define BRICKMASK_WFITS_ALLCOL          0
#include "fits_write.c"

#define BRICKMASK_WFITS_MTYPE           TLONG
#define BRICKMASK_WFITS_SUBID           1
#define BRICKMASK_WFITS_OVERWRITE       1
#define BRICKMASK_WFITS_ALLCOL          1
#include "fits_write.c"


/*============================================================================*\
                 Interface for saving the FITS-format catalogue
\*============================================================================*/

/******************************************************************************
Function `save_fits`:
  Write the data catalogue to a FITS file.
Arguments:
  * `conf`:     structure for storing configurations;
  * `data`:     structure for the the data catalogue.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int save_fits(const CONF *conf, DATA *data) {
  if (!conf) {
    P_ERR("configuration parameters are not loaded\n");
    return BRICKMASK_ERR_INIT;
  }
  if (!data) {
    P_ERR("the data catalog is not read\n");
    return BRICKMASK_ERR_INIT;
  }

  /* Restore the orders of masks and subsample IDs. */
  if (reorder_mask(data) || reorder_subid(data)) return BRICKMASK_ERR_MEMORY;

  /* Generate the filename for force overwriting. */
  char *output = force_output(conf->output);
  if (!output) return BRICKMASK_ERR_MEMORY;

  /* Choose the function for saving the FITS catalogue. */
  int (*save_fits_func) (const char *, const CONF *, const DATA *) = NULL;
  switch (data->mtype) {
    case TBYTE:
      if (strcmp(conf->input, conf->output)) {
        if (data->subid) {
          save_fits_func = (conf->ncol) ? fits_save_uint8_t_subid :
              fits_save_uint8_t_subid_all;
        }
        else {
          save_fits_func = (conf->ncol) ? fits_save_uint8_t :
              fits_save_uint8_t_all;
        }
      }
      else {
        if (data->subid) {
          save_fits_func = (conf->ncol) ? fits_save_uint8_t_subid_overwrite :
              fits_save_uint8_t_subid_overwrite_all;
        }
        else {
          save_fits_func = (conf->ncol) ? fits_save_uint8_t_overwrite :
              fits_save_uint8_t_overwrite_all;
        }
      }
      break;
    case TSHORT:
      if (strcmp(conf->input, conf->output)) {
        if (data->subid) {
          save_fits_func = (conf->ncol) ? fits_save_uint16_t_subid :
              fits_save_uint16_t_subid_all;
        }
        else {
          save_fits_func = (conf->ncol) ? fits_save_uint16_t :
              fits_save_uint16_t_all;
        }
      }
      else {
        if (data->subid) {
          save_fits_func = (conf->ncol) ? fits_save_uint16_t_subid_overwrite :
              fits_save_uint16_t_subid_overwrite_all;
        }
        else {
          save_fits_func = (conf->ncol) ? fits_save_uint16_t_overwrite :
              fits_save_uint16_t_overwrite_all;
        }
      }
      break;
    case TINT:
      if (strcmp(conf->input, conf->output)) {
        if (data->subid) {
          save_fits_func = (conf->ncol) ? fits_save_uint32_t_subid :
              fits_save_uint32_t_subid_all;
        }
        else {
          save_fits_func = (conf->ncol) ? fits_save_uint32_t :
              fits_save_uint32_t_all;
        }
      }
      else {
        if (data->subid) {
          save_fits_func = (conf->ncol) ? fits_save_uint32_t_subid_overwrite :
              fits_save_uint32_t_subid_overwrite_all;
        }
        else {
          save_fits_func = (conf->ncol) ? fits_save_uint32_t_overwrite :
              fits_save_uint32_t_overwrite_all;
        }
      }
      break;
    case TLONG:
      if (strcmp(conf->input, conf->output)) {
        if (data->subid) {
          save_fits_func = (conf->ncol) ? fits_save_uint64_t_subid :
              fits_save_uint64_t_subid_all;
        }
        else {
          save_fits_func = (conf->ncol) ? fits_save_uint64_t :
              fits_save_uint64_t_all;
        }
      }
      else {
        if (data->subid) {
          save_fits_func = (conf->ncol) ? fits_save_uint64_t_subid_overwrite :
              fits_save_uint64_t_subid_overwrite_all;
        }
        else {
          save_fits_func = (conf->ncol) ? fits_save_uint64_t_overwrite :
              fits_save_uint64_t_overwrite_all;
        }
      }
      break;
    default:
      P_ERR("unexpected data type for maskbits: %d\n", data->mtype);
      free(output);
      return BRICKMASK_ERR_UNKNOWN;
  }

  /* Save the FITS catalogue. */
  if (save_fits_func(output, conf, data)) {
    free(output);
    return BRICKMASK_ERR_SAVE;
  }

  free(output);
  return 0;
}
