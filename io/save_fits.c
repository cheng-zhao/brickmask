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
                        Function for saving FITS files
\*============================================================================*/

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
  * `data`:     structure for the the data catalogue;
  * `idx`:      index of the output catalogue.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int save_fits(const CONF *conf, const DATA *data, const int idx) {
  if (!conf) {
    P_ERR("configuration parameters are not loaded\n");
    return BRICKMASK_ERR_INIT;
  }
  if (!data) {
    P_ERR("the data catalog is not read\n");
    return BRICKMASK_ERR_INIT;
  }

  /* Generate the filename for force overwriting. */
  char *output = force_output(conf->output[idx]);
  if (!output) return BRICKMASK_ERR_MEMORY;

  /* Choose the function for saving the FITS catalogue. */
  int (*save_fits_func) (const char *, const CONF *, const DATA *, const int) =
      NULL;
  switch (data->mtype) {
    case TBYTE:
      if (strcmp(conf->input[idx], conf->output[idx])) {
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
      if (strcmp(conf->input[idx], conf->output[idx])) {
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
      if (strcmp(conf->input[idx], conf->output[idx])) {
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
      if (strcmp(conf->input[idx], conf->output[idx])) {
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
  if (save_fits_func(output, conf, data, idx)) {
    free(output);
    return BRICKMASK_ERR_SAVE;
  }

  free(output);
  return 0;
}
