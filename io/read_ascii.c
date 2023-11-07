/*******************************************************************************
* read_ascii.c: this file is part of the brickmask program.

* brickmask: assign bit codes defined on Legacy Survey brick pixels
             to a catalogue with sky coordinates.

* Github repository:
        https://github.com/cheng-zhao/brickmask

* Copyright (c) 2020 -- 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "read_file.h"
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* Data structure for information of the ASCII columns. */
typedef struct {
  int max;              /* the maximum number of columns to be read     */
  size_t *idx;          /* starting indices of different columns        */
  int c[2];             /* column indices for RA and Dec                */
  int ncol;             /* number of columns to be read                 */
  int *cid;             /* indices of columns to be read                */
} ASCII_COL_t;

/*============================================================================*\
                      Functions for reading ASCII columns
\*============================================================================*/

/******************************************************************************
Function `ascii_col_destroy`:
  Deconstruct the structure for ASCII columns.
Arguments:
  * `col`:      structure for ASCII columns.
******************************************************************************/
static void ascii_col_destroy(ASCII_COL_t *col) {
  if (!col) return;
  if (col->idx) free(col->idx);
  if (col->cid) free(col->cid);
  free(col);
}

/******************************************************************************
Function `ascii_col_init`:
  Initialise the structure for ASCII columns.
Arguments:
  * `conf`:     structure for storing configurations.
Return:
  Address of the structure for ASCII columns on success; NULL on error.
******************************************************************************/
static ASCII_COL_t *ascii_col_init(const CONF *conf) {
  ASCII_COL_t *col = calloc(1, sizeof(ASCII_COL_t));
  if (!col) {
    P_ERR("failed to allocate memory for ASCII columns\n");
    return NULL;
  }
  col->idx = NULL;
  col->cid = NULL;

  col->ncol = conf->ncol;
  col->c[0] = conf->cnum[0] - 1;        /* column index starting from 0 */
  col->c[1] = conf->cnum[1] - 1;
  col->max = (conf->cnum[0] > conf->cnum[1]) ? conf->cnum[0] : conf->cnum[1];

  if (conf->ncol) {     /* read columns required by the output */
    if (!(col->cid = malloc(col->ncol * sizeof(int)))) {
      P_ERR("failed to allocate memory for reading ASCII columns\n");
      ascii_col_destroy(col);
      return NULL;
    }
    for (int i = 0; i < conf->ncol; i++) {
      if (col->max < conf->onum[i]) col->max = conf->onum[i];
      col->cid[i] = conf->onum[i] - 1;
    }
  }

  if (!(col->idx = calloc(col->max + 1, sizeof(size_t)))) {
    P_ERR("failed to allocate memory for column indices\n");
    ascii_col_destroy(col);
    return NULL;
  }

  return col;
}

/******************************************************************************
Function `chunk_resize`:
  Enlarge the size of the chunk for reading files.
Arguments:
  * `chunk`:    address of the chunk;
  * `size`:     current size of the chunk.
Return:
  New address of the chunk on success; NULL on error.
******************************************************************************/
static inline char *chunk_resize(void *chunk, size_t *size) {
  if (BRICKMASK_MAX_CHUNK / 2 < *size) return NULL;
  *size *= 2;
  char *tmp = realloc(chunk, *size * sizeof(char));
  return tmp;
}

/******************************************************************************
Function `column_index`:
  Find column indices of a line.
Arguments:
  * `line`:     starting address of the line;
  * `num`:      number of characters in the line;
  * `col`:      structure for storing information of columns.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static inline int column_index(const char *line, const size_t num,
    ASCII_COL_t *col) {
  col->idx[0] = 0;
  bool incol = true;            /* now we are at the first column */
  int k = 1;
  for (size_t i = 1; i < num; i++) {
    if (incol) {
      if (!isspace(line[i])) continue;
      incol = false;
    }
    else {
      if (isspace(line[i])) continue;
      incol = true;
      col->idx[k++] = i;
      if (k == col->max) break;
    }
  }
  if (k != col->max || !incol) {
    P_ERR("too few columns of line:\n%s\n", line);
    return BRICKMASK_ERR_FILE;
  }

  /* Find the end of the last column. */
  col->idx[col->max] = num;
  if (col->ncol) {
    for (size_t i = col->idx[col->max - 1]; i < num; i++) {
      if (isspace(line[i])) {
        col->idx[col->max] = i;
        break;
      }
    }
  }
  return 0;
}

/******************************************************************************
Function `copy_column`:
  Copy a column from the input catalogue into the memory.
Arguments:
  * `content`:  address in memory for receiving columns;
  * `size`:     size of the memory that is already consumed;
  * `max`:      mamxium size of the memory allocated for receiving columns;
  * `p`:        starting address of the source to be copied;
  * `num`:      number of characters to be copied.
Return:
  Address for receiving columns on success; NULL on error.
******************************************************************************/
static inline char *copy_column(char *content, size_t *size, size_t *max,
    const char *p, const size_t num) {
  /* Enlarge the memory if necessary. */
  while (*size + num >= *max) {         /* >=: reserve one byte for '\0'. */
    if (*max < BRICKMASK_CONTENT_MAX_DOUBLE_SIZE) {
      if (BRICKMASK_CONTENT_MAX_SIZE / 2 < *max) return NULL;
      *max *= 2;
    }
    else {
      if (BRICKMASK_CONTENT_MAX_SIZE - BRICKMASK_CONTENT_MAX_DOUBLE_SIZE < *max)
        return NULL;
      *max += BRICKMASK_CONTENT_MAX_DOUBLE_SIZE;
    }
  }
  content = realloc(content, *max * sizeof(char));
  if (!content) return NULL;
  memcpy(content + *size, p, num);
  *size += num;
  return content;
}

/******************************************************************************
Function `read_ascii_col`:
  Read columns of an ASCII text file.
Arguments:
  * `fname`:    filename of the input catalogue;
  * `comment`:  symbol indicating lines to be skipped;
  * `col`:      structure for columns to be read;
  * `data`:     structure for storing the input data catalogue.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int read_ascii_col(const char *fname, const char comment,
    ASCII_COL_t *col, DATA *data) {
  /* Allocate memory for the chunk. */
  size_t cmax = BRICKMASK_FILE_CHUNK;
  char *chunk = malloc(cmax * sizeof(char));
  if (!chunk) {
    P_ERR("failed to allocate memory for reading file by chunks\n");
    return BRICKMASK_ERR_MEMORY;
  }

  /* Open the file for reading. */
  FILE *fp;
  if (!(fp = fopen(fname, "r"))) {
    P_ERR("cannot open file for reading: `%s'\n", fname);
    free(chunk);
    return BRICKMASK_ERR_FILE;
  }
  size_t nread, nrest;
  nrest = 0;

  /* Start reading the file by chunk. */
  while ((nread = fread(chunk + nrest, sizeof(char), cmax - nrest, fp))) {
    char *p = chunk;
    char *end = p + nrest + nread;
    char *endl;
    if (nread < cmax - nrest) *end = '\n';     /* append '\n' to last line */

    /* Process lines in the chunk. */
    while ((endl = memchr(p, '\n', end - p))) {
      *endl = '\0';             /* replace '\n' by string terminator '\0' */
      while (isspace(*p)) ++p;          /* omit leading whitespaces */
      if (*p == comment || *p == '\0') {        /* comment or empty */
        p = endl + 1;
        continue;
      }

      /* Find indices of columns. */
      if (column_index(p, endl - p, col)) {
        free(chunk); fclose(fp);
        return BRICKMASK_ERR_FILE;
      }

      /* Copy output columns into memory. */
      data->cidx[data->n] = data->csize;
      if (col->ncol) {          /* copy given columns */
        for (int i = 0; i < col->ncol; i++) {
          const int c = col->cid[i];            /* current column number */
          char *tmp = copy_column(data->content, &data->csize, &data->cmax,
              p + col->idx[c], col->idx[c + 1] - col->idx[c]);
          if (!tmp) {
            P_ERR("failed to save columns of the input file: `%s'\n", fname);
            free(chunk); fclose(fp);
            return BRICKMASK_ERR_MEMORY;
          }
          data->content = tmp;
          /* Append whitespace to the last column of the original file. */
          if (c == col->max - 1) {
            if (!(tmp = copy_column(data->content, &data->csize, &data->cmax,
                " ", 1))) {
              P_ERR("failed to save columns of the input file: `%s'\n", fname);
              free(chunk); fclose(fp);
              return BRICKMASK_ERR_MEMORY;
            }
            data->content = tmp;
          }
        }
      }
      else {                    /* copy all columns */
        char *tmp = copy_column(data->content, &data->csize, &data->cmax,
            p, endl - p);
        if (!tmp) {
          P_ERR("failed to save columns of the input file: `%s'\n", fname);
          free(chunk); fclose(fp);
          return BRICKMASK_ERR_MEMORY;
        }
        data->content = tmp;
        /* Append white space to the end of the line. */
        if (!(tmp = copy_column(data->content, &data->csize, &data->cmax,
            " ", 1))) {
          P_ERR("failed to save columns of the input file: `%s'\n", fname);
          free(chunk); fclose(fp);
          return BRICKMASK_ERR_MEMORY;
        }
        data->content = tmp;
      }
      /* Always append a '\0' at the end. */
      *((char *) data->content + data->csize++) = '\0';

      /* Parse RA and Dec. */
      if (sscanf(p + col->idx[col->c[0]], "%lf", data->ra + data->n) != 1 ||
          sscanf(p + col->idx[col->c[1]], "%lf", data->dec + data->n) != 1) {
        P_ERR("failed to read coordinates from file: `%s':\n%s\n", fname, p);
        free(chunk); fclose(fp);
        return BRICKMASK_ERR_FILE;
      }

      /* Enlarge memory for the data if necessary. */
      if (++data->n >= data->nmax) {
        if (SIZE_MAX / 2 < data->nmax) {
          P_ERR("too many objects in the file: `%s'\n", fname);
          free(chunk); fclose(fp);
          return BRICKMASK_ERR_FILE;
        }
        data->nmax *= 2;
        double *tmp = realloc(data->ra, data->nmax * sizeof(double));
        if (!tmp) {
          P_ERR("failed to enlarge memory for the input catalog\n");
          free(chunk); fclose(fp);
          return BRICKMASK_ERR_MEMORY;
        }
        data->ra = tmp;
        if (!(tmp = realloc(data->dec, data->nmax * sizeof(double)))) {
          P_ERR("failed to enlarge memory for the input catalog\n");
          free(chunk); fclose(fp);
          return BRICKMASK_ERR_MEMORY;
        }
        data->dec = tmp;
        size_t *stmp = realloc(data->cidx, data->nmax * sizeof(size_t));
        if (!stmp) {
          P_ERR("failed to enlarge memory for the input catalog\n");
          free(chunk); fclose(fp);
          return BRICKMASK_ERR_MEMORY;
        }
        data->cidx = stmp;
      }

      /* Continue with the next line. */
      p = endl + 1;
    }

    /* The chunk cannot hold a full line. */
    if (p == chunk) {
      char *tmp = chunk_resize(chunk, &cmax);
      if (!tmp) {
        P_ERR("failed to enlarge memory for reading file by chunk\n");
        free(chunk); fclose(fp);
        return BRICKMASK_ERR_MEMORY;
      }
      chunk = tmp;
      nrest += nread;
      continue;
    }

    /* Copy the remaining characters to the beginning of the chunk. */
    nrest = end - p;
    memmove(chunk, p, nrest);
  }

  free(chunk);
  if (!feof(fp)) {
    P_ERR("unexpected end of file: `%s'\n", fname);
    fclose(fp);
    return BRICKMASK_ERR_FILE;
  }
  if (fclose(fp)) P_WRN("failed to close file: `%s'\n", fname);

  return 0;
}


/*============================================================================*\
                       Interfaces for reading ASCII files
\*============================================================================*/

/******************************************************************************
Function `read_ascii`:
  Read data from the input ASCII catalogue.
Arguments:
  * `fname`:    filename of the input catalogue;
  * `conf`:     structure for storing configurations;
  * `data`:     structure for the input data catalogue.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
int read_ascii(const char *fname, const CONF *conf, DATA *data) {
  if (!conf) {
    P_ERR("configuration parameters are not loaded\n");
    return BRICKMASK_ERR_INIT;
  }
  if (!data) {
    P_ERR("structure for the input data is not initialised\n");
    return BRICKMASK_ERR_INIT;
  }

  /* Check columns to be read. */
  ASCII_COL_t *col = ascii_col_init(conf);
  if (!col) {
    P_ERR("failed to allocate memory for ASCII columns\n");
    return BRICKMASK_ERR_MEMORY;
  }

  /* Read the file. */
  if (read_ascii_col(fname, conf->comment, col, data)) {
    ascii_col_destroy(col);
    return BRICKMASK_ERR_FILE;
  }

  ascii_col_destroy(col);
  return 0;
}

/******************************************************************************
Function `read_fname`:
  Read filenames from a text file.
Arguments:
  * `fname`:    name of the file storing filenames;
  * `output`:   string array for storing filenames;
  * `num`:      number of filenames read from file.
Return:
  Number of bytes read on success; zero on error.
******************************************************************************/
size_t read_fname(const char *fname, char ***output, size_t *num) {
  if (!fname || !(*fname)) {
    P_ERR("the list for maskbit filenames is not available\n");
    return 0;
  }
  if (!output || !num) {
    P_ERR("the bricks are not initialised\n");
    return 0;
  }

  /* Allocate memory for filenames. */
  size_t nmax = BRICKMASK_DATA_INIT_NUM;
  size_t imax = BRICKMASK_DATA_INIT_NUM;
  size_t nsize = 0;
  char *names = NULL;           /* a single array for all filenames */
  size_t *idx = NULL;           /* indices of filenames in the array */
  if (!(names = malloc(nmax * sizeof(char))) ||
      !(idx = malloc(imax * sizeof(size_t)))) {
    P_ERR("failed to allocate memory for maskbit filenames\n");
    if (names) free(names);
    return 0;
  }

  /* Allocate memory for the chunk. */
  size_t cmax = BRICKMASK_FILE_CHUNK;
  char *chunk = malloc(cmax * sizeof(char));
  if (!chunk) {
    P_ERR("failed to allocate memory for reading file by chunks\n");
    free(names); free(idx);
    return 0;
  }

  /* Open the file for reading. */
  FILE *fp;
  if (!(fp = fopen(fname, "r"))) {
    P_ERR("cannot open file for reading: `%s'\n", fname);
    free(names); free(idx); free(chunk);
    return 0;
  }
  size_t n, nread, nrest;
  n = nrest = 0;

  /* Start reading the file by chunk. */
  while ((nread = fread(chunk + nrest, sizeof(char), cmax - nrest, fp))) {
    char *p = chunk;
    char *end = p + nrest + nread;
    char *endl;
    if (nread < cmax - nrest) *end = '\n';     /* append '\n' to last line */

    /* Process lines in the chunk. */
    while ((endl = memchr(p, '\n', end - p))) {
      *endl = '\0';             /* replace '\n' by string terminator '\0' */
      while (isspace(*p)) ++p;          /* omit leading whitespaces */
      if (*p == BRICKMASK_READ_COMMENT || *p == '\0') { /* comment or empty */
        p = endl + 1;
        continue;
      }

      /* Parse the line. */
      idx[n] = nsize;
      do {
        if (isspace(*p)) {
          /* Overwrite the escape character if applicable. */
          if (*(p - 1) == BRICKMASK_SPACE_ESCAPE) nsize--;
          else break;
        }
        names[nsize] = *p;
        /* Enlarge memory for the filenames if necessary. */
        if (++nsize >= nmax - 1) {      /* reserve one byte for '\0' */
          char *tmp = chunk_resize(names, &nmax);
          if (!tmp) {
            P_ERR("failed to enlarge memory for the filenames\n");
            free(names); free(idx); free(chunk); fclose(fp);
            return 0;
          }
          names = tmp;
        }
      }
      while (*(++p) != '\0');
      names[nsize++] = '\0';            /* null termination */

      /* Enlarge memory for the indices if necessary. */
      if (++n >= imax) {
        if (SIZE_MAX / 2 < imax) {
          P_ERR("too many entries in the file: `%s'\n", fname);
          free(names); free(idx); free(chunk); fclose(fp);
          return 0;
        }
        imax *= 2;
        size_t *tmp = realloc(idx, imax * sizeof(size_t));
        if (!tmp) {
          P_ERR("failed to enlarge memory for the filename indices\n");
          free(names); free(idx); free(chunk); fclose(fp);
          return 0;
        }
        idx = tmp;
      }

      /* Continue with the next line. */
      p = endl + 1;
    }

    /* The chunk cannot hold a full line. */
    if (p == chunk) {
      char *tmp = chunk_resize(chunk, &cmax);
      if (!tmp) {
        P_ERR("failed to enlarge memory for reading the file by chunk\n");
        free(names); free(idx); free(chunk); fclose(fp);
        return 0;
      }
      chunk = tmp;
      nrest += nread;
      continue;
    }

    /* Copy the remaining characters to the beginning of the chunk. */
    nrest = end - p;
    memmove(chunk, p, nrest);
  }

  free(chunk);
  if (!feof(fp)) {
    P_ERR("unexpected end of file: `%s'\n", fname);
    free(names); free(idx); fclose(fp);
    return 0;
  }
  if (fclose(fp)) P_WRN("failed to close file: `%s'\n", fname);

  if (!n) {
    P_ERR("no valid filename found in file: `%s'\n", fname);
    free(names); free(idx);
    return 0;
  }

  /* Reduce the memory cost if applicable. */
  char *tmp = realloc(names, nsize * sizeof(char));
  if (tmp) names = tmp;

  /* Allocate memory for the filename array. */
  if (!(*output = malloc(n * sizeof(char *)))) {
    P_ERR("failed to allocate memory for filenames\n");
    free(names); free(idx);
    return 0;
  }
  **output = names;
  for (size_t i = 1; i < n; i++) (*output)[i] = names + idx[i];
  *num = n;

  free(idx);
  return nsize;
}
