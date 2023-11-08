/*******************************************************************************
* read_fits.c: this file is part of the brickmask program.

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
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include <math.h>

/* Validate macros for data types. */
#if TBYTE >= TSHORT
  #error "Datatype error: TBYTE must be smaller than TSHORT"
#endif
#if TSHORT >= TINT
  #error "Datatype error: TSHORT must be smaller than TINT"
#endif
#if TINT >= TLONG
  #error "Datatype error: TINT must be smaller than TLONG"
#endif
#if CHAR_BIT != 8
  #error "this program works only on architectures with 8 bits per byte"
#endif

#define FITS_ABORT {                                            \
  P_ERR("cfitsio error: ");                                     \
  fits_report_error(stderr, status);                            \
  status = 0;                                                   \
  if (fp) fits_close_file(fp, &status);                         \
  return BRICKMASK_ERR_FILE;                                    \
}

/*============================================================================*\
                     Functions for processing FITS columns
\*============================================================================*/

/******************************************************************************
Function `get_fits_col`:
  Get information of columns in the FITS file.
Arguments:
  * `conf`:     structure for storing configurations;
  * `data`:     structure for the input data catalogue;
  * `fp`:       pointer to the opened FITS file.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static inline int get_fits_col(const CONF *conf, DATA *data, fitsfile *fp) {
  if (!conf->ncol) return 0;    /* copy the whole file directly */
  FITS_COL_t *col = malloc(conf->ncol * sizeof(FITS_COL_t));
  if (!col) {
    P_ERR("failed to allocate memory for the information of FITS columns\n");
    return BRICKMASK_ERR_MEMORY;
  }
  data->content = col;

  /* Get column numbers. */
  int status = 0;
  for (int i = 0; i < conf->ncol; i++) {
    if (fits_get_colnum(fp, BRICKMASK_FITS_CASESEN, conf->ocol[i],
        &col[i].n, &status)) FITS_ABORT;
  }

  /* Check if there are duplicate columns. */
  for (int i = 1; i < conf->ncol; i++) {
    for (int j = 0; j < i; j++) {
      if (col[j].n == col[i].n) {
        P_ERR("FITS columns `%s' and `%s' are essentially identical\n",
            conf->ocol[j], conf->ocol[i]);
        fits_close_file(fp, &status);
        return BRICKMASK_ERR_CFG;
      }
    }
  }

  /* Get the total number of columns. */
  int nc = 0;
  if (fits_get_num_cols(fp, &nc, &status)) FITS_ABORT;

  /* Get column positions and widths. */
  long idx = 1;         /* indices of CFITSIO start from 1 */
  for (int i = 1; i <= nc; i++) {
    long w = 0;         /* width (in bytes) of this column */
    if (fits_get_coltype(fp, i, NULL, NULL, &w, &status)) FITS_ABORT;
    for (int j = 0; j < conf->ncol; j++) {
      if (col[j].n == i) {
        col[j].i = idx;
        col[j].w = w;
        break;
      }
    }
    idx += w;
  }

  return 0;
}

/******************************************************************************
Function `get_fits_coord`:
  Read coordinates from the FITS file.
Arguments:
  * `conf`:     structure for storing configurations;
  * `data`:     structure for the input data catalogue;
  * `ndata`:    number of rows to be read;
  * `fp`:       pointer to the opened FITS file.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int get_fits_coord(const CONF *conf, DATA *data, const size_t ndata,
    fitsfile *fp) {
  /* Get columns for RA and Dec. */
  int status = 0;
  int col[2];
  const char *cname[2] = {"RA", "Dec"};
  for (int i = 0; i < 2; i++) {
    char colname[FLEN_VALUE];
    memset(colname, 0, FLEN_VALUE);
    if (fits_get_colname(fp, BRICKMASK_FITS_CASESEN, conf->cname[i], colname,
        col + i, &status)) FITS_ABORT;

    /* Check if the column is specified by number. */
    if (strncmp(colname, conf->cname[i], FLEN_VALUE)) {
      P_WRN("the FITS column name for %s is: `%s'\n", cname[i], colname);
    }
  }

  /* Get the optimal number of rows to read at one time. */
  long nstep = 0;
  if (fits_get_rowsize(fp, &nstep, &status)) FITS_ABORT;

  /* Read the file and retrieve coordinates. */
  long nread = 1;
  long nrest = ndata;
  int anynul;
  while (nrest) {
    long nrow = (nstep < nrest) ? nstep : nrest;
    if (fits_read_col_dbl(fp, col[0], nread, 1, nrow, 0,
        data->ra + (data->n + nread - 1), &anynul, &status)) FITS_ABORT;
    if (fits_read_col_dbl(fp, col[1], nread, 1, nrow, 0,
        data->dec + (data->n + nread - 1), &anynul, &status)) FITS_ABORT;
    nread += nrow;
    nrest -= nrow;
  }

  data->n += ndata;
  return 0;
}

/*============================================================================*\
                       Functions for reading FITS tables
\*============================================================================*/

/******************************************************************************
Function `read_brick`:
  Read the brick name and range of (RA, Dec) from a brick list file.
Arguments:
  * `fname`:    the filename of the brick list;
  * `brick`:    structure for storing information of bricks.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int read_brick(const char *fname, BRICK *brick) {
  int status = 0;
  fitsfile *fp = NULL;
  long n;

  /* Open FITS file and read the number of bricks. */
  if (fits_open_data(&fp, fname, READONLY, &status)) FITS_ABORT;
  if (fits_get_num_rows(fp, &n, &status)) FITS_ABORT;
  if (!(brick->n = n)) {
    P_ERR("no brick found in file: `%s'\n", fname);
    fits_close_file(fp, &status);
    return BRICKMASK_ERR_FILE;
  }

  /* Allocate memory. */
  if (!(brick->ra1 = malloc(n * sizeof(double))) ||
      !(brick->ra2 = malloc(n * sizeof(double))) ||
      !(brick->dec1 = malloc(n * sizeof(double))) ||
      !(brick->dec2 = malloc(n * sizeof(double))) ||
      !(brick->photsys = malloc(n * sizeof(unsigned char))) ||
      !(brick->name = malloc(n * sizeof(char *)))) {
    P_ERR("failed to allocate memory for the information of bricks\n");
    return BRICKMASK_ERR_MEMORY;
  }
  for (size_t i = 0; i < brick->n; i++) brick->name[i] = NULL;

  /* Get the length of brick names. */
  int col, nlen;
  if (fits_get_colnum(fp, BRICKMASK_FITS_CASESEN, BRICKMASK_FITS_BRICKNAME,
      &col, &status)) FITS_ABORT;
  if (fits_get_col_display_width(fp, col, &nlen, &status)) FITS_ABORT;

  /* Allocate memory for brick names, with null-termination. */
  if (!(brick->name[0] = calloc(n * (nlen + 1), sizeof(char)))) {
    P_ERR("failed to allocate memory for the names of bricks\n");
    return BRICKMASK_ERR_MEMORY;
  }
  for (size_t i = 1; i < brick->n; i++)
    brick->name[i] = brick->name[0] + i * (nlen + 1);

  /* Get columns of coordinate ranges and photsys. */
  int ccol[5];
  if (fits_get_colnum(fp, BRICKMASK_FITS_CASESEN, BRICKMASK_FITS_RAMIN,
      ccol, &status)) FITS_ABORT;
  if (fits_get_colnum(fp, BRICKMASK_FITS_CASESEN, BRICKMASK_FITS_RAMAX,
      ccol + 1, &status)) FITS_ABORT;
  if (fits_get_colnum(fp, BRICKMASK_FITS_CASESEN, BRICKMASK_FITS_DECMIN,
      ccol + 2, &status)) FITS_ABORT;
  if (fits_get_colnum(fp, BRICKMASK_FITS_CASESEN, BRICKMASK_FITS_DECMAX,
      ccol + 3, &status)) FITS_ABORT;
  if (fits_get_colnum(fp, BRICKMASK_FITS_CASESEN, BRICKMASK_FITS_PHOTSYS,
      ccol + 4, &status)) FITS_ABORT;

  /* Get the optimal number of rows to read at one time. */
  long nstep = 0;
  if (fits_get_rowsize(fp, &nstep, &status)) FITS_ABORT;

  /* Read the file and retrieve coordinates. */
  long nread = 1;
  long nrest = n;
  int anynul = 0;
  while (nrest) {
    long nrow = (nstep < nrest) ? nstep : nrest;
    /* Read brick names. */
    if (fits_read_col_str(fp, col, nread, 1, nrow, "", brick->name + nread - 1,
        &anynul, &status)) FITS_ABORT;
    /* Read ranges of sky coordinates. */
    if (fits_read_col_dbl(fp, ccol[0], nread, 1, nrow, 0,
        brick->ra1 + nread - 1, &anynul, &status)) FITS_ABORT;
    if (fits_read_col_dbl(fp, ccol[1], nread, 1, nrow, 0,
        brick->ra2 + nread - 1, &anynul, &status)) FITS_ABORT;
    if (fits_read_col_dbl(fp, ccol[2], nread, 1, nrow, 0,
        brick->dec1 + nread - 1, &anynul, &status)) FITS_ABORT;
    if (fits_read_col_dbl(fp, ccol[3], nread, 1, nrow, 0,
        brick->dec2 + nread - 1, &anynul, &status)) FITS_ABORT;
    /* Read photsys. */
    if (fits_read_col_byt(fp, ccol[4], nread, 1, nrow, '\0',
        brick->photsys + nread - 1, &anynul, &status)) FITS_ABORT;
    nread += nrow;
    nrest -= nrow;
  }

  /* Finished reading the file. */
  if (fits_close_file(fp, &status)) FITS_ABORT;

#ifdef MPI
  brick->nlen = nlen;
#endif
  return 0;
}

/******************************************************************************
Function `read_fits`:
  Read data from the input FITS catalogue.
Arguments:
  * `fname`:    filename of the input catalogue;
  * `conf`:     structure for storing configurations;
  * `data`:     structure for the input data catalogue.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int read_fits(const char *fname, const CONF *conf, DATA *data) {
  if (!conf) {
    P_ERR("configuration parameters are not loaded\n");
    return BRICKMASK_ERR_INIT;
  }
  if (!data) {
    P_ERR("structure for the input data is not initialised\n");
    return BRICKMASK_ERR_INIT;
  }

  /* Open the file for reading. */
  int status = 0;
  fitsfile *fp = NULL;
  if (fits_open_data(&fp, fname, READONLY, &status)) FITS_ABORT;

  /* Check if the maskbit and subsample ID columns are alread in the file. */
  int col = 0;
  if (fits_get_colnum(fp, BRICKMASK_FITS_CASESEN, conf->mcol,
      &col, &status) != COL_NOT_FOUND) {
    P_ERR("the maskbit column (%s) exists in the input catalog.",
        conf->mcol);
    status = 0;
    fits_close_file(fp, &status);
    return BRICKMASK_ERR_FILE;
  }
  status = 0;
  fits_clear_errmsg();

  /* Get the number of objects. */
  long ndata = 0;
  if (fits_get_num_rows(fp, &ndata, &status)) FITS_ABORT;
  if (!ndata) {
    if (fits_close_file(fp, &status)) FITS_ABORT;
    return 0;
  }

  /* Allocate memory. */
  double *tmp[2];
  if (!(tmp[0] = realloc(data->ra, (data->n + ndata) * sizeof(double))) ||
      !(tmp[1] = realloc(data->dec, (data->n + ndata) * sizeof(double)))) {
    P_ERR("failed to allocate memory for the input data catalog\n");
    if (tmp[0]) free(tmp[0]);
    fits_close_file(fp, &status);
    return BRICKMASK_ERR_MEMORY;
  }
  data->ra = tmp[0];
  data->dec = tmp[1];

  /* Get properties of output columns. */
  if (!data->content && get_fits_col(conf, data, fp)) return BRICKMASK_ERR_FILE;

  /* Read coordinates. */
  if (get_fits_coord(conf, data, ndata, fp)) return BRICKMASK_ERR_FILE;

  /* Close file. */
  if (fits_close_file(fp, &status)) FITS_ABORT;
  return 0;
}


/*============================================================================*\
                      Functions for reading maskbit images
\*============================================================================*/

/******************************************************************************
Function `read_wcs_header`:
  Read and preprocess WCS keywords from the FITS header.
Arguments:
  * `fp`:       pointer to the opened FITS file;
  * `wcs`:      structure for WCS parameters.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static inline int read_wcs_header(fitsfile *fp, WCS *wcs) {
  int status = 0;
  double a[2];
  if (fits_read_key_dbl(fp, "CRVAL1", a, NULL, &status)) FITS_ABORT;
  if (fits_read_key_dbl(fp, "CRVAL2", a + 1, NULL, &status)) FITS_ABORT;
  if (fits_read_key_dbl(fp, "CRPIX1", wcs->r, NULL, &status)) FITS_ABORT;
  if (fits_read_key_dbl(fp, "CRPIX2", wcs->r + 1, NULL, &status)) FITS_ABORT;
  if (fits_read_key_dbl(fp, "CD1_1", wcs->m[0], NULL, &status)) FITS_ABORT;
  if (fits_read_key_dbl(fp, "CD1_2", wcs->m[0] + 1, NULL, &status)) FITS_ABORT;
  if (fits_read_key_dbl(fp, "CD2_1", wcs->m[1], NULL, &status)) FITS_ABORT;
  if (fits_read_key_dbl(fp, "CD2_2", wcs->m[1] + 1, NULL, &status)) FITS_ABORT;

  a[0] *= DEGREE_2_RAD;
  a[1] *= DEGREE_2_RAD;
  double sina = sin(a[0]);
  double cosa = cos(a[0]);
  double sind = sin(a[1]);
  double cosd = cos(a[1]);
  wcs->ang[0] = sind;
  wcs->ang[1] = cosa * cosd;
  wcs->ang[2] = sina * cosd;
  wcs->ang[3] = -cosd;
  wcs->ang[4] = cosa * sind;
  wcs->ang[5] = sina * sind;
  wcs->ang[6] = -sina;
  wcs->ang[7] = cosa;
  wcs->idetm = wcs->m[0][0] * wcs->m[1][1] - wcs->m[0][1] * wcs->m[1][0];

  if (!wcs->idetm) {
    P_ERR("the translation matrix is not invertable:\n  " OFMT_DBL "  "
        OFMT_DBL "\n  " OFMT_DBL "  " OFMT_DBL "\n",
        wcs->m[0][0], wcs->m[0][1], wcs->m[1][0], wcs->m[1][1]);
    fits_close_file(fp, &status);
    return BRICKMASK_ERR_MASK;
  }
  wcs->idetm = 1 / wcs->idetm;

  return 0;
}

/******************************************************************************
Function `read_mask`:
  Read a maskbit file.
Arguments:
  * `fname`:    name of a masbit file;
  * `mask`:     structure for maskbits;
  * `dtype`:    data type of bit codes.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int read_mask(const char *fname, MASK *mask, int *dtype) {
  if (!fname) {
    P_ERR("the maskbit filename is not available\n");
    return BRICKMASK_ERR_INIT;
  }
  if (!mask) {
    P_ERR("the structure for maskbits is not initialised\n");
    return BRICKMASK_ERR_INIT;
  }

  /* Open the maskbit file and check the file type. */
  fitsfile *fp = NULL;
  int status, hdutype, naxis, bitpix;
  status = hdutype = naxis = bitpix = 0;

  if (fits_open_data(&fp, fname, READONLY, &status)) FITS_ABORT;

  if (fits_get_hdu_type(fp, &hdutype, &status)) FITS_ABORT;
  if (hdutype != IMAGE_HDU) {
    P_ERR("the first HDU of the maskbit file must be IMAGE_HDU: `%s'\n", fname);
    fits_close_file(fp, &status);
    return BRICKMASK_ERR_MASK;
  }

  /* Check image dimensions and data type. */
  if (fits_get_img_param(fp, 2, &bitpix, &naxis, mask->dim, &status)) FITS_ABORT;
  if (naxis != 2) {
    P_ERR("image dimension of the maskbit file must be 2: `%s'\n", fname);
    fits_close_file(fp, &status);
    return BRICKMASK_ERR_MASK;
  }
  if (mask->dim[0] <= 0 || mask->dim[1] <= 0) {
    P_ERR("invalid image dimension (%ld, %ld) of maskbit file: `%s'\n",
        mask->dim[0], mask->dim[1], fname);
    fits_close_file(fp, &status);
    return BRICKMASK_ERR_MASK;
  }

  /* Check WCS keywords. */
  char value[FLEN_VALUE];
  if (fits_read_key_str(fp, "CTYPE1", value, NULL, &status)) FITS_ABORT;
  if (strncmp(value, "RA---TAN", 9)) {
    P_ERR("unsupported WCS header: CTYPE1 = %s (only 'RA---TAN' is allowed)\n",
        value);
    fits_close_file(fp, &status);
    return BRICKMASK_ERR_MASK;
  }
  if (fits_read_key_str(fp, "CTYPE2", value, NULL, &status)) FITS_ABORT;
  if (strncmp(value, "DEC--TAN", 9)) {
    P_ERR("unsupported WCS header: CTYPE2 = %s (only 'DEC--TAN' is allowed)\n",
        value);
    fits_close_file(fp, &status);
    return BRICKMASK_ERR_MASK;
  }

  /* Allocate memory for the maskbits. */
  if (mask->dim[0] * mask->dim[1] > mask->size) {
    mask->size = mask->dim[0] * mask->dim[1];
    int nbyte = bitpix / CHAR_BIT;
    unsigned char *tmp = realloc(mask->bit, mask->size * nbyte);
    if (!tmp) {
      P_ERR("failed to allocate memory for maskbits\n");
      fits_close_file(fp, &status);
      return BRICKMASK_ERR_MEMORY;
    }
    mask->bit = tmp;
  }

  switch (bitpix) {
    case BYTE_IMG:     *dtype = TBYTE;  break;
    case SHORT_IMG:    *dtype = TSHORT; break;
    case LONG_IMG:     *dtype = TINT;   break;
    case LONGLONG_IMG: *dtype = TLONG;  break;
    default:
      P_ERR("invalid data type (%d) of the maskbit image: `%s'\n",
          bitpix, fname);
      fits_close_file(fp, &status);
      return BRICKMASK_ERR_MASK;
  }

  /* Read and preprocess WCS keywords. */
  if (read_wcs_header(fp, mask->wcs)) return BRICKMASK_ERR_MASK;

  /* Read maskbits. */
#ifdef FAST_FITS_IMG
  /* Low-level image access. */
  if (fits_read_tblbytes(fp, 1, 1, mask->size, mask->bit, &status)) FITS_ABORT;
#else
  /* Flexible image access (for compressed images). */
  if (fits_set_bscale(fp, 1, 0, &status)) FITS_ABORT;
  if (fits_read_img(fp, *dtype, 1, mask->size, 0, mask->bit, 0, &status))
    FITS_ABORT;
#endif

  if (fits_close_file(fp, &status)) FITS_ABORT;
  return 0;
}
