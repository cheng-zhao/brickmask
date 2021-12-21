/*******************************************************************************
* fits_write.c: this file is part of the brickmask program.

* brickmask: assign bit codes defined on Legacy Survey brick pixels
             to a catalogue with sky coordinates.

* Github repository:
        https://github.com/cheng-zhao/brickmask

* Copyright (c) 2020 -- 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

/* Macros for the template functions. */
#if defined(BRICKMASK_WFITS_MTYPE) && defined(BRICKMASK_WFITS_SUBID) && \
  defined(BRICKMASK_WFITS_ALLCOL)

/* Macros for generating function names. */
#ifndef CONCAT_FNAME
  #define CONCAT_FNAME(a,b,c,d)         a##_##b##c##d
#endif

#ifndef FITS_WRITE_FUNC
  #define FITS_WRITE_FUNC(a,b,c,d)      CONCAT_FNAME(a,b,c,d)
#endif


/*============================================================================*\
                             Definition validation
\*============================================================================*/

#ifdef BRICKMASK_MASKBIT_DTYPE
  #undef BRICKMASK_MASKBIT_DTYPE
#endif
#ifdef BRICKMASK_MASKBIT_TFORM
  #undef BRICKMASK_MASKBIT_TFORM
#endif
#ifdef FITS_WRITE_SUBID_NAME
  #undef FITS_WRITE_SUBID_NAME
#endif
#ifdef FITS_WRITE_ALLCOL_NAME
  #undef FITS_WRITE_ALLCOL_NAME
#endif

#if     BRICKMASK_WFITS_MTYPE == TBYTE
  #define BRICKMASK_MASKBIT_DTYPE       uint8_t
  #define BRICKMASK_MASKBIT_TFORM       "B"
#elif   BRICKMASK_WFITS_MTYPE == TSHORT
  #define BRICKMASK_MASKBIT_DTYPE       uint16_t
  #define BRICKMASK_MASKBIT_TFORM       "I"
#elif   BRICKMASK_WFITS_MTYPE == TINT
  #define BRICKMASK_MASKBIT_DTYPE       uint32_t
  #define BRICKMASK_MASKBIT_TFORM       "J"
#elif   BRICKMASK_WFITS_MTYPE == TLONG
  #define BRICKMASK_MASKBIT_DTYPE       uint64_t
  #define BRICKMASK_MASKBIT_TFORM       "K"
#else
  #error "unexpected definition of `BRICKMASK_WFITS_MTYPE`"
#endif

#if     BRICKMASK_WFITS_SUBID == 1
  #define FITS_WRITE_SUBID_NAME         _subid
#elif   BRICKMASK_WFITS_SUBID == 0
  #define FITS_WRITE_SUBID_NAME
#else
  #error "unexpected definition of `BRICKMASK_WFITS_SUBID`"
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

#ifndef FITS_WRITE_ABORT
#define FITS_WRITE_ABORT {                                      \
  P_ERR("cfitsio error: ");                                     \
  fits_report_error(stderr, status);                            \
  status = 0;                                                   \
  if (ifp) { fits_close_file(ifp, &status); status = 0; }       \
  if (fp) fits_close_file(fp, &status);                         \
  if (memptr) free(memptr);                                     \
  if (chunk) free(chunk);                                       \
  return BRICKMASK_ERR_FILE;                                    \
}
#endif


/*============================================================================*\
                       Function for saving a FITS catalog
\*============================================================================*/

/******************************************************************************
Function `fits_save_<BRICKMASK_MASKBIT_DTYPE><FITS_WRITE_SUBID_NAME>
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
static int FITS_WRITE_FUNC(fits_save, BRICKMASK_MASKBIT_DTYPE,
    FITS_WRITE_SUBID_NAME, FITS_WRITE_ALLCOL_NAME)
    (const char *fname, const CONF *conf, const DATA *data) {
  int status = 0;
  fitsfile *fp, *ifp;
  fp = ifp = NULL;
  void *memptr = NULL;          /* pointer for a FITS file in memory */
  unsigned char *chunk = NULL;  /* array for copying FITS data       */

  /* Open the input file to read columns. */
  if (fits_open_data(&ifp, conf->input, READONLY, &status)) FITS_WRITE_ABORT;

  /* Get the total number of columns. */
  int nc = 0;
  if (fits_get_num_cols(ifp, &nc, &status)) FITS_WRITE_ABORT;

  /* Compute the total length (in bytes) of the columns to be read. */
  long iwidth = 0;
#if BRICKMASK_WFITS_ALLCOL == 1
  for (int i = 1; i <= nc; i++) {
    long w = 0;         /* width (in bytes) of this column */
    if (fits_get_coltype(ifp, i, NULL, NULL, &w, &status)) FITS_WRITE_ABORT;
    iwidth += w;
  }
#else
  FITS_COL_t *col = data->content;
  for (int i = 0; i < conf->ncol; i++) iwidth += col[i].w;
#endif

  /* Compute the total length (in bytes) of output columns. */
  long owidth = iwidth + sizeof(BRICKMASK_MASKBIT_DTYPE);
#if BRICKMASK_WFITS_SUBID == 1
  owidth++;
#endif

  /* Create the output file for receiving columns. */
  if (!strcmp(conf->input, conf->output)) {     /* create file in memory */
    size_t memsize = data->n * owidth;  /* initial memory size for the file */
    if (memsize < BRICKMASK_FITS_MIN_MEMSIZE)
      memsize = BRICKMASK_FITS_MIN_MEMSIZE;
    size_t deltasize = memsize / 5;     /* increase 20% of the initial size */
    if (deltasize < BRICKMASK_FITS_MIN_MEMSIZE)
      deltasize = BRICKMASK_FITS_MIN_MEMSIZE;

    if (!(memptr = calloc(memsize, sizeof(char)))) {
      P_ERR("failed to allocate memory for the output file\n");
      fits_close_file(ifp, &status);
      return BRICKMASK_ERR_MEMORY;
    }
    if (fits_create_memfile(&fp, &memptr, &memsize, deltasize, realloc,
        &status)) FITS_WRITE_ABORT;
  }
  else if (fits_create_file(&fp, fname, &status)) FITS_WRITE_ABORT;

  /* Copy the header. */
  if (fits_copy_hdutab(ifp, fp, 1, 0, &status)) FITS_WRITE_ABORT;

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

  /* Append maskbit and subsample ID columns. */
  if (fits_insert_col(fp, nc + 1, conf->mcol,
      BRICKMASK_MASKBIT_TFORM, &status)) FITS_WRITE_ABORT;
#if BRICKMASK_WFITS_SUBID == 1
  if (fits_insert_col(fp, nc + 2, BRICKMASK_FITS_SUBID, "B", &status))
    FITS_WRITE_ABORT;
#endif

  /* Insert all rows. */
  long nr = 0;
  if (fits_get_num_rows(ifp, &nr, &status)) FITS_WRITE_ABORT;
  if (fits_insert_rows(fp, 0, nr, &status)) FITS_WRITE_ABORT;

  /* Get the optimal row number for saving data. */
  long nstep;
  if (fits_get_rowsize(fp, &nstep, &status)) FITS_WRITE_ABORT;

  /* Allocate memory for copying data by chunks. */
  if (!(chunk = malloc(owidth * nstep * sizeof(unsigned char)))) {
    P_ERR("failed to allocate memory for copying FITS columns\n");
    fits_close_file(ifp, &status); status = 0;
    fits_close_file(fp, &status); if (memptr) free(memptr);
    return BRICKMASK_ERR_MEMORY;
  }

  /* Copy data and append column in chunks. */
  long nread = 1;
  long nrest = nr;
  while (nrest) {
    long nrow = (nstep < nrest) ? nstep : nrest;

    /* Read columns. */
    long idx = 0;
    for (long i = nread; i < nread + nrow; i++) {
#if BRICKMASK_WFITS_ALLCOL == 1
      /* Read all columns at once. */
      if (fits_read_tblbytes(ifp, i, 1, iwidth, chunk + idx, &status))
        FITS_WRITE_ABORT;
      idx += iwidth;
#else
      /* Read specific columns. */
      for (int j = 0; j < conf->ncol; j++) {
        if (fits_read_tblbytes(ifp, i, col[j].i, col[j].w, chunk + idx,
            &status)) FITS_WRITE_ABORT;
        idx += col[j].w;
      }
#endif
      /* Append maskbit value without considering endianness. */
      memcpy(chunk + idx, (BRICKMASK_MASKBIT_DTYPE *) (data->mask) + i - 1,
          sizeof(BRICKMASK_MASKBIT_DTYPE));
      idx += sizeof(BRICKMASK_MASKBIT_DTYPE);
      /* Append subsample ID. */
#if BRICKMASK_WFITS_SUBID == 1
      memcpy(chunk + idx, data->subid + i - 1, sizeof(unsigned char));
      idx++;
#endif
    }

    /* Write the columns at once. */
    idx = 0;
    for (long i = nread; i < nread + nrow; i++) {
      if (fits_write_tblbytes(fp, i, 1, owidth, chunk + idx, &status))
        FITS_WRITE_ABORT;
      idx += owidth;
    }

    nread += nrow;
    nrest -= nrow;
  }

  free(chunk);
  chunk = NULL;
  if (fits_close_file(ifp, &status)) FITS_WRITE_ABORT;
  ifp = NULL;

  /* Save the file in memory to an actual disk file. */
  if (memptr) {
    if (fits_create_file(&ifp, fname, &status)) FITS_WRITE_ABORT;
    if (fits_copy_file(fp, ifp, 1, 1, 1, &status)) FITS_WRITE_ABORT;
    if (fits_close_file(ifp, &status)) FITS_WRITE_ABORT;
    if (fits_close_file(fp, &status)) {
      P_ERR("cfitsio error: "); fits_report_error(stderr, status);
      free(memptr);
      return BRICKMASK_ERR_FILE;
    }
    free(memptr);
  }
  /* Close the file directly. */
  else if (fits_close_file(fp, &status)) {
    P_ERR("cfitsio error: "); fits_report_error(stderr, status);
    return BRICKMASK_ERR_FILE;
  }

  return 0;
}


#undef BRICKMASK_WFITS_MTYPE
#undef BRICKMASK_WFITS_SUBID
#undef BRICKMASK_WFITS_ALLCOL

#undef BRICKMASK_MASKBIT_DTYPE
#undef BRICKMASK_MASKBIT_TFORM
#undef FITS_WRITE_SUBID_NAME
#undef FITS_WRITE_ALLCOL_NAME

#undef FITS_WRITE_ABORT

#endif
