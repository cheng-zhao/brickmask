/*******************************************************************************
* fits_write.c: this file is part of the brickmask program.

* brickmask: assign bit codes defined on Legacy Survey brick pixels
             to a catalogue with sky coordinates.

* Github repository:
        https://github.com/cheng-zhao/brickmask

* Copyright (c) 2020 -- 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

/* Macros for the template functions. */
#if defined(BRICKMASK_WFITS_OVERWRITE) && defined(BRICKMASK_WFITS_ALLCOL)

/* Macros for generating function names. */
#ifndef CONCAT_FNAME
  #define CONCAT_FNAME(a,b,c)           a##b##c
#endif

#ifndef FITS_WRITE_FUNC
  #define FITS_WRITE_FUNC(a,b,c)        CONCAT_FNAME(a,b,c)
#endif


/*============================================================================*\
                             Definition validation
\*============================================================================*/

#ifdef FITS_WRITE_OVERWRITE_NAME
  #undef FITS_WRITE_OVERWRITE_NAME
#endif
#ifdef FITS_WRITE_ALLCOL_NAME
  #undef FITS_WRITE_ALLCOL_NAME
#endif

#if     BRICKMASK_WFITS_OVERWRITE == 1
  #define FITS_WRITE_OVERWRITE_NAME     _overwrite
#elif   BRICKMASK_WFITS_OVERWRITE == 0
  #define FITS_WRITE_OVERWRITE_NAME
#endif

#if     BRICKMASK_WFITS_ALLCOL == 1
  #define FITS_WRITE_ALLCOL_NAME        _all
#elif   BRICKMASK_WFITS_ALLCOL == 0
  #define FITS_WRITE_ALLCOL_NAME
#else
  #error "unexpected definition of `BRICKMASK_WFITS_ALLCOL`"
#endif


/*============================================================================*\
                          Shortcuts for error handling
\*============================================================================*/

#ifdef FITS_WRITE_ABORT
  #undef FITS_WRITE_ABORT
#endif

#if BRICKMASK_WFITS_OVERWRITE == 1
  #define FITS_WRITE_ABORT {                                            \
    P_ERR("cfitsio error: "); fits_report_error(stderr, status);        \
    status = 0; if (fp) { fits_close_file(fp, &status); status = 0; }   \
    if (chunk) free(chunk);                                             \
    if (tab) free(tab);                                                 \
    return BRICKMASK_ERR_FILE;                                          \
  }
#else
  #define FITS_WRITE_ABORT {                                            \
    P_ERR("cfitsio error: "); fits_report_error(stderr, status);        \
    status = 0; if (fp) { fits_close_file(fp, &status); status = 0; }   \
    if (ofp) fits_close_file(ofp, &status);                             \
    if (chunk) free(chunk);                                             \
    if (tab) free(tab);                                                 \
    return BRICKMASK_ERR_FILE;                                          \
  }
#endif


/*============================================================================*\
                       Function for saving a FITS catalog
\*============================================================================*/

/******************************************************************************
Function `fits_save_<BRICKMASK_MASKBIT_DTYPE><FITS_WRITE_SUBID_NAME>
    <FITS_WRITE_OVERWRITE_NAME><FITS_WRITE_ALLCOL_NAME>`:
  Save specific columns of a FITS file to another FITS file, and append
  additional columns.
Arguments:
  * `fname`:    filename for the output FITS file;
  * `conf`:     structure for storing configurations;
  * `data`:     structure for the data catalogue;
  * `icat`:      index of the output catalogue.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int FITS_WRITE_FUNC(fits_save, FITS_WRITE_OVERWRITE_NAME,
    FITS_WRITE_ALLCOL_NAME)
    (const char *fname, const CONF *conf, const DATA *data, const int icat) {
  int status = 0;
  fitsfile *fp = NULL;
#if BRICKMASK_WFITS_OVERWRITE == 1
#else
  fitsfile *ofp = NULL;
#endif
  unsigned char *chunk = NULL;  /* array for reading FITS table */
  unsigned char *tab = NULL;    /* contents of the FITS table to be written */

  /* Define the datatypes of the maskbits columns. */
  int dtypes[4];
  char *cols[4];
  unsigned char *codes[4];
  dtypes[0] = data->mtype;
  cols[0] = conf->mcol;
  codes[0] = (unsigned char *) data->mask;
  for (int i = 0; i < 3; i++) {
    dtypes[i + 1] = data->etype[i];
    cols[i + 1] = conf->ecol[i];
    codes[i + 1] = (unsigned char *) data->nexp[i];
  }

  size_t sizes[4];
  char *tforms[4];
  for (int i = 0; i < 4; i++) {
    switch (dtypes[i]) {
      case TBYTE:
        sizes[i] = sizeof(uint8_t);
        tforms[i] = "B";
        break;
      case TSHORT:
        sizes[i] = sizeof(uint16_t);
        tforms[i] = "I";
        break;
      case TINT:
        sizes[i] = sizeof(uint32_t);
        tforms[i] = "J";
        break;
      case TLONG:
        sizes[i] = sizeof(uint64_t);
        tforms[i] = "K";
        break;
      default:
        P_ERR("unexpected data type for maskbits: %d\n", dtypes[i]);
        return BRICKMASK_ERR_UNKNOWN;
    }
  }

  /* Open the input file to read columns. */
#if BRICKMASK_WFITS_OVERWRITE == 1
  if (fits_open_data(&fp, conf->input[icat], READWRITE, &status))
    FITS_WRITE_ABORT;
#else
  if (fits_open_data(&fp, conf->input[icat], READONLY, &status))
    FITS_WRITE_ABORT;
#endif

  /* Get the total number of columns and rows. */
  int nc = 0;
  if (fits_get_num_cols(fp, &nc, &status)) FITS_WRITE_ABORT;
  long nr = 0;
  if (fits_get_num_rows(fp, &nr, &status)) FITS_WRITE_ABORT;

  /* Compute the total length (in bytes) of input columns. */
  long iwidth = 0;
  for (int i = 1; i <= nc; i++) {
    long w = 0;         /* width (in bytes) of this column */
    if (fits_get_coltype(fp, i, NULL, NULL, &w, &status)) FITS_WRITE_ABORT;
    iwidth += w;
  }

  /* Compute the total length (in bytes) of output columns. */
#if BRICKMASK_WFITS_ALLCOL == 1
  long owidth = iwidth;
#else
  long owidth = 0;
  FITS_COL_t *col = data->content;
  for (int i = 0; i < conf->ncol; i++) owidth += col[i].w;
#endif
  for (int i = 0; i < 4; i++) owidth += sizes[i];

#if BRICKMASK_WFITS_OVERWRITE == 0
  /* Create the output file for receiving columns. */
  if (fits_create_file(&ofp, fname, &status)) FITS_WRITE_ABORT;
  /* Copy the header. */
  if (fits_copy_hdutab(fp, ofp, 1, 0, &status)) FITS_WRITE_ABORT;

  #if BRICKMASK_WFITS_ALLCOL == 0
  /* Duplicate columns to be saved. */
  for (int i = 0; i < conf->ncol; i++) {
    if (fits_copy_col(ofp, ofp, col[i].n, nc + i + 1, 1, &status))
      FITS_WRITE_ABORT;
  }
  /* Remove original columns. */
  for (int i = 0; i < nc; i++) {
    if (fits_delete_col(ofp, 1, &status)) FITS_WRITE_ABORT;
  }
  nc = conf->ncol;
  #endif

  /* Append maskbit columns. */
  for (int i = 0; i < 4; i++) {
    if (fits_insert_col(ofp, nc + i + 1, cols[i], tforms[i], &status))
      FITS_WRITE_ABORT;
  }
#endif

  /* Set the number of rows to be read/written at once */
  long nstep;
  if (fits_get_rowsize(fp, &nstep, &status)) FITS_WRITE_ABORT;
  if (nstep < BRICKMASK_FILE_CHUNK / iwidth)
    nstep = BRICKMASK_FILE_CHUNK / iwidth;
  const long nchunk = nstep * iwidth;
#if BRICKMASK_WFITS_OVERWRITE == 1
  const long ntab = nr * owidth;
#else
  const long ntab = nstep * owidth;
#endif

  /* Allocate memory for reading data by chunks. */
  if (!(chunk = malloc(nchunk))) {
    P_ERR("failed to allocate memory for reading FITS columns\n");
    fits_close_file(fp, &status); status = 0;
#if BRICKMASK_WFITS_OVERWRITE == 0
    fits_close_file(ofp, &status);
#endif
    return BRICKMASK_ERR_MEMORY;
  }
  /* Allocate memory for the output FITS table. */
  if (!(tab = malloc(ntab))) {
    P_ERR("failed to allocate memory for writing FITS columns\n");
    fits_close_file(fp, &status); status = 0;
#if BRICKMASK_WFITS_OVERWRITE == 0
    fits_close_file(ofp, &status);
#endif
    if (chunk) free(chunk);
    return BRICKMASK_ERR_MEMORY;
  }

  /* Copy data and append column in chunks. */
  long nread = 1;
  long nrest = nr;
#if BRICKMASK_WFITS_OVERWRITE == 1
  size_t idx = 0;
#endif
  while (nrest) {
    long nrow = (nstep < nrest) ? nstep : nrest;

    /* Read the chunk at once. */
    if (fits_read_tblbytes(fp, nread, 1, nrow * iwidth, chunk, &status))
      FITS_WRITE_ABORT;

    /* Construct the FITS table to be written. */
#if BRICKMASK_WFITS_OVERWRITE == 0
    size_t idx = 0;
#endif
    for (long i = 0; i < nrow; i++) {
      /* Copy columns. */
      unsigned char *ichunk = chunk + i * iwidth;
#if BRICKMASK_WFITS_ALLCOL == 1
      memcpy(tab + idx, ichunk, iwidth);
      idx += iwidth;
#else
      for (int j = 0; j < conf->ncol; j++) {
        memcpy(tab + idx, ichunk + col[j].i - 1, col[j].w);
        idx += col[j].w;
      }
#endif

      /* Append maskbit value with big endian. */
      const long didx = i + nread - 1;
      for (int j = 0; j < 4; j++) {
#ifdef WITH_BIG_ENDIAN
        memcpy(tab + idx, codes[j] + (data->iidx[icat] + didx) * sizes[j],
            sizes[j]);
        idx += sizes[j];
#else
        unsigned char *msk = codes[j] + (data->iidx[icat] + didx) * sizes[j];
        if (dtypes[j] >= TLONG) {
          tab[idx++] = msk[7];
          tab[idx++] = msk[6];
          tab[idx++] = msk[5];
          tab[idx++] = msk[4];
        }
        if (dtypes[j] >= TINT) {
          tab[idx++] = msk[3];
          tab[idx++] = msk[2];
        }
        if (dtypes[j] >= TSHORT) {
          tab[idx++] = msk[1];
        }
        tab[idx++] = msk[0];
#endif
      }
    }

#if BRICKMASK_WFITS_OVERWRITE == 0
    /* Write the FITS table. */
    if (fits_write_tblbytes(ofp, nread, 1, nrow * owidth, tab, &status))
      FITS_WRITE_ABORT;
#endif

    nread += nrow;
    nrest -= nrow;
  }

  free(chunk);
  chunk = NULL;

#if BRICKMASK_WFITS_OVERWRITE == 1
  /* Delete the existing FITS table. */
  if (fits_delete_rows(fp, 1, nr, &status)) FITS_WRITE_ABORT;
  #if BRICKMASK_WFITS_ALLCOL == 0
  /* Duplicate columns to be saved. */
  for (int i = 0; i < conf->ncol; i++) {
    if (fits_copy_col(fp, fp, col[i].n, nc + i + 1, 1, &status))
      FITS_WRITE_ABORT;
  }
  /* Remove original columns. */
  for (int i = 0; i < nc; i++) {
    if (fits_delete_col(fp, 1, &status)) FITS_WRITE_ABORT;
  }
  nc = conf->ncol;
  #endif

  /* Append maskbit columns. */
  for (int i = 0; i < 4; i++) {
    if (fits_insert_col(fp, nc + i + 1, cols[i], tforms[i], &status))
      FITS_WRITE_ABORT;
  }

  /* Write the FITS table. */
  if (fits_write_tblbytes(fp, 1, 1, ntab, tab, &status)) FITS_WRITE_ABORT;
#endif

  free(tab);
  tab = NULL;
  if (fits_close_file(fp, &status)) FITS_WRITE_ABORT;
#if BRICKMASK_WFITS_OVERWRITE == 0
  /* Close the file directly. */
  if (fits_close_file(ofp, &status)) {
    P_ERR("cfitsio error: "); fits_report_error(stderr, status);
    return BRICKMASK_ERR_FILE;
  }
#endif

  return 0;
}


#undef BRICKMASK_WFITS_OVERWRITE
#undef BRICKMASK_WFITS_ALLCOL

#undef FITS_WRITE_OVERWRITE_NAME
#undef FITS_WRITE_ALLCOL_NAME

#undef FITS_WRITE_ABORT

#endif
