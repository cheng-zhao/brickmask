/*******************************************************************************
* save_res.c: this file is part of the brickmask program.

* brickmask: assign bit codes defined on Legacy Survey brick pixels
             to a catalogue with sky coordinates.

* Github repository:
        https://github.com/cheng-zhao/brickmask

* Copyright (c) 2020 -- 2021 Cheng Zhao <zhaocheng03@gmail.com>
  [MIT license]
 
*******************************************************************************/

#include "define.h"
#include "save_file.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <stdbool.h>
#include <limits.h>
#include <inttypes.h>

/* Shortcut for writing a line to the file. */
#define WRITE_LINE(...)                                         \
  if (output_writeline(__VA_ARGS__)) {                          \
    output_destroy(ofile); return BRICKMASK_ERR_FILE;           \
  }

/* Data structure for writing ASCII file by chunk. */
typedef struct {
  const char *fname;    /* name of the output file                    */
  FILE *fp;             /* pointer to the file stream                 */
  char *chunk;          /* buffer for writing file by chunks          */
  int size;             /* number of used characters of the buffer    */
  int max;              /* maximum number of characters of the buffer */
} OFILE;

/*============================================================================*\
                       Functions for writing ASCII files
\*============================================================================*/

/******************************************************************************
Function `output_init`:
  Initialise the interface for writing ASCII files.
Return:
  Address of the interface.
******************************************************************************/
static OFILE *output_init(void) {
  OFILE *ofile = malloc(sizeof *ofile);
  if (!ofile) {
    P_ERR("failed to initialise the interface for file writing\n");
    return NULL;
  }

  ofile->fname = NULL;
  ofile->fp = NULL;
  ofile->size = 0;
  ofile->max = BRICKMASK_FILE_CHUNK;
  ofile->chunk = calloc(ofile->max, sizeof(char));
  if (!ofile->chunk) {
    P_ERR("failed to allocate memory for writing file by chunk\n");
    free(ofile);
    return NULL;
  }

  return ofile;
}

/******************************************************************************
Function `output_flush`:
  Write the buffer string to file.
Arguments:
  * `ofile`:    structure for writing ASCII files.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int output_flush(OFILE *ofile) {
  if (!ofile) {
    P_ERR("the interface for file writing is not initialised\n");
    return BRICKMASK_ERR_INIT;
  }

  if (!ofile->size) return 0;

  if (fwrite(ofile->chunk, ofile->size * sizeof(char), 1, ofile->fp) != 1) {
    P_ERR("failed to write to the output file: `%s'\n", ofile->fname);
    return BRICKMASK_ERR_FILE;
  }

  ofile->size = 0;
  return 0;
}

/******************************************************************************
Function `output_destroy`:
  Deconstruct the interface for writing ASCII files.
Arguments:
  * `ofile`:    structure for writing ASCII files.
******************************************************************************/
static void output_destroy(OFILE *ofile) {
  if (!ofile) return;
  if (output_flush(ofile))
    P_WRN("closing the file with unsaved buffer\n");
  free(ofile->chunk);
  if (ofile->fp) {
    if (fclose(ofile->fp))
      P_WRN("failed to close file: `%s'\n", ofile->fname);
  }
  free(ofile);
}

/******************************************************************************
Function `output_newfile`:
  Flush the buffer to the existing file and open a new file.
Arguments:
  * `ofile`:    structure for writing ASCII files;
  * `fname`:    name of the file to be written to.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int output_newfile(OFILE *ofile, const char *fname) {
  if (!ofile) {
    P_ERR("the interface for file writing is not initialised\n");
    return BRICKMASK_ERR_INIT;
  }
  if (!fname || !(*fname)) {
    P_ERR("invalid output file name\n");
    return BRICKMASK_ERR_INIT;
  }

  if (output_flush(ofile)) {
    P_ERR("failed to flush the buffer to the opened file: `%s'\n",
        ofile->fname);
    return BRICKMASK_ERR_FILE;
  }
  if (ofile->fp && fclose(ofile->fp))
    P_WRN("failed to close file: `%s'\n", ofile->fname);

  ofile->fname = fname;
  if (!(ofile->fp = fopen(fname, "w"))) {
    P_ERR("failed to open the file for writing: `%s'\n", fname);
    return BRICKMASK_ERR_FILE;
  }
  return 0;
}

/******************************************************************************
Function `output_writeline`:
  Write a line to the buffer and save it to the file if necessary.
Arguments:
  * `ofile`:    structure for writing ASCII files;
  * `format`:   a string specifying how the data is interpreted;
  * `...`:      arguments specifying data to be saved.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int output_writeline(OFILE *ofile, const char *restrict format, ...) {
  if (!ofile) {
    P_ERR("the interface for file writing is not initialised\n");
    return BRICKMASK_ERR_INIT;
  }
  if (!ofile->fp) {
    P_ERR("no opened file for writing the line\n");
    return BRICKMASK_ERR_FILE;
  }
  if (!format || !(*format)) {
    P_ERR("the formatter string is not given for writing files\n");
    return BRICKMASK_ERR_INIT;
  }

  if (ofile->size == ofile->max) {
    if (output_flush(ofile)) return BRICKMASK_ERR_FILE;
  }

  va_list args;
  va_start(args, format);
  int size = vsnprintf(ofile->chunk + ofile->size, ofile->max - ofile->size,
      format, args);
  va_end(args);

  if (size < 0) {
    P_ERR("failed to save the result in format: `%s'\n", format);
    return BRICKMASK_ERR_FILE;
  }
  if (size >= BRICKMASK_MAX_CHUNK) {
    P_ERR("the line to be saved is too long\n");
    return BRICKMASK_ERR_FILE;
  }

  /* Check if the buffer is full. */
  if (size >= ofile->max - ofile->size) {
    bool enlarge = false;
    /* Enlarge the size of the buffer if necessary. */
    while (size >= ofile->max) {
      if (BRICKMASK_MAX_CHUNK / 2 < ofile->max) ofile->max = size + 1;
      else ofile->max <<= 1;
      enlarge = true;
    }

    if (enlarge) {
      char *tmp = realloc(ofile->chunk, ofile->max * sizeof(char));
      if (!tmp) {
        P_ERR("failed to allocate memory for saving the line\n");
        return BRICKMASK_ERR_MEMORY;
      }
      ofile->chunk = tmp;
    }

    if (output_flush(ofile)) return BRICKMASK_ERR_FILE;

    va_start(args, format);
    size = vsnprintf(ofile->chunk, ofile->max, format, args);
    va_end(args);

    if (size >= ofile->max) {
      P_ERR("failed to detect the size of the line\n");
      return BRICKMASK_ERR_UNKNOWN;
    }

    ofile->size = size;
  }
  else ofile->size += size;

  return 0;
}


/*============================================================================*\
                Interface for saving the ASCII-format catalogue
\*============================================================================*/

/******************************************************************************
Function `save_ascii`:
  Write the data catalogue to an ASCII file.
Arguments:
  * `conf`:     structure for storing configurations;
  * `data`:     structure for the the data catalogue;
  * `idx`:      index of the output catalogue.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int save_ascii(const CONF *conf, const DATA *data, const int idx) {
  if (!conf) {
    P_ERR("configuration parameters are not loaded\n");
    return BRICKMASK_ERR_INIT;
  }
  if (!data) {
    P_ERR("the data catalog is not read\n");
    return BRICKMASK_ERR_INIT;
  }

  /* Initialise the interface for writing files. */
  OFILE *ofile = output_init();
  if (!ofile) return BRICKMASK_ERR_FILE;
  if (output_newfile(ofile, conf->output[idx])) {
    output_destroy(ofile);
    return BRICKMASK_ERR_FILE;
  }

  /* Write the catalog. */
  const char *content = data->content;
  for (size_t i = data->iidx[idx]; i < data->iidx[idx + 1]; i++) {
    /* Write columns of the original file. */
    WRITE_LINE(ofile, "%s", content + data->cidx[i]);

    /* Write maskbits. */
    WRITE_LINE(ofile, "%" PRIu64, data->mask[i]);

    /* Write subsample IDs if applicable. */
    if (data->subid) WRITE_LINE(ofile, " %" PRIu8, data->subid[i]);

    WRITE_LINE(ofile, "\n");
  }

  output_destroy(ofile);
  return 0;
}
