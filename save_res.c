/*******************************************************************************
* 
* brickmask: assign bit codes defined on Legacy Survey brick pixels
*            to a catalogue with sky coordinates.
*
* Github repository:
*       https://github.com/cheng-zhao/brickmask
*
* Copyright (c) 2020 Cheng Zhao <zhaocheng03@gmail.com>
* 
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
* 
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*
*******************************************************************************/

#include "brickmask.h"

/******************************************************************************
Function `write_ascii`:
  Write the input galaxy data with the information of masks to the output.

Arguments:
  * `fname`:    the filename for the output;
  * `data`:     a pointer to the data to be saved;
  * `num`:      the number of galaxies;
  * `verbose`:  0 for concise outputs, 1 for detailed outputs.
Return:
  A non-zero integer if there is problem.
******************************************************************************/
int write_ascii(const char *fname, const DATA *data, const size_t num,
    const int verbose) {
  FILE *fp;
  int n, k;
  size_t i;
  char *buf, *end;
  char stmp[LEN_OUT_LINE];

  if (verbose) printf("\n  Filename : %s.\n", fname);

  if (!(fp = fopen(fname, "w"))) {
    P_ERR("cannot write to file `%s'.\n", fname);
    return ERR_FILE;
  }

  MY_ALLOC(buf, char, CHUNK, writing outputs);
  end = buf;

  /* Write file by chunks. */
  for (i = 0; i < num; i++) {
#if NUMSUB == 1
    n = snprintf(stmp, LEN_OUT_LINE, OFMT_DBL " " OFMT_DBL " %s %" PRIu64 "\n",
        data[i].ra, data[i].dec, data[i].rest, data[i].mask);
    CHECKSTR(n, LEN_OUT_LINE, "output line too long: " OFMT_DBL " " OFMT_DBL
        " %s %" PRIu64 "\nPlease enlarge LEN_OUT_LINE in `define.h`.\n",
        data[i].ra, data[i].dec, data[i].rest, data[i].mask);
#else
    n = snprintf(stmp, LEN_OUT_LINE, OFMT_DBL " " OFMT_DBL " %s %" PRIu64
        " %" SUB_COL_FMT "\n",
        data[i].ra, data[i].dec, data[i].rest, data[i].mask, data[i].flag);
    CHECKSTR(n, LEN_OUT_LINE, "output line too long: " OFMT_DBL " " OFMT_DBL
        " %s %" PRIu64 " %" SUB_COL_FMT "\n"
        "Please enlarge LEN_OUT_LINE in `define.h`.\n",
        data[i].ra, data[i].dec, data[i].rest, data[i].mask, data[i].flag);
#endif

    if (end - buf + n < CHUNK) {        /* there is still space in buf */
      k = safe_strcpy(end, stmp, n + 1);
      CHECKSTR(k, n + 1, "unexpected error for writing: %s\n", stmp);
      end += k;
    }
    else {                              /* write buf to file */
      if (fwrite(buf, sizeof(char) * (end - buf), 1, fp) != 1) {
        P_ERR("failed to write to output: %s\n", stmp);
        return ERR_FILE;
      }
      k = safe_strcpy(buf, stmp, n + 1);
      CHECKSTR(k, n + 1, "unexpected error for writing: %s\n", stmp);
      end = buf + k;
    }
  }

  if ((n = end - buf) > 0) {
    if (fwrite(buf, sizeof(char) * n, 1, fp) != 1) {
      P_ERR("failed to write to output: %s\n", stmp);
      return ERR_FILE;
    }
  }

  fclose(fp);
  free(buf);
  return 0;
}


/******************************************************************************
Function `write_fits`:
  Append the mask and chunk information to the input galaxies, and write to a
  fits file.

Argument:
  * `input`:    filename of the input galaxy catalog;
  * `output`:   filename of the output file;
  * `data`:     a pointer to the data to be saved;
  * `num`:      the number of galaxies;
  * `verbose`:  0 for concise outputs, 1 for detailed outputs.
Return:
  A non-zero integer if there is problem.
******************************************************************************/
int write_fits(const char *input, const char *output, const DATA *data,
    const long num, const int masktype, const int verbose) {
  fitsfile *fin, *fptr;
  char fname[FLEN_FILENAME];
  int n, exist, col, status, nbyte;
  long i;
  unsigned char *mskcol;

#if NUMSUB == 1
  const int nnew = 1;
  char *ttype[1] = {MSK_COL};
  char *tform[1];
#else
  subid_t *subcol;
  const int nnew = 2;
  char *ttype[2] = {MSK_COL, SUB_COL};
  char *tform[2];

  tform[1] = SUB_COL_TYPE;
#endif

  switch (masktype) {
    case TBYTE:
      tform[0] = "B";
      nbyte = 1;
      break;
    case TSHORT:
      tform[0] = "UI";
      nbyte = 2;
      break;
    case TINT:
      tform[0] = "UK";
      nbyte = 4;
      break;
    case TLONG:
      tform[0] = "UJJ";
      nbyte = 8;
      break;
    default:
      P_EXT("wrong data type for maskbit codes.\n");
      return ERR_INPUT;
  }

  if (verbose) printf("\n  Filename: %s\n", output);

  status = 0;
  if (fits_file_exists(output, &exist, &status)) FITS_EXIT;
  if (exist > 0) n = snprintf(fname, FLEN_FILENAME, "!%s", output);
  else n = safe_strcpy(fname, output, FLEN_FILENAME);
  CHECKSTR(n, FLEN_FILENAME, "length of the output filename is too long.\n");

  status = 0;
  /* Copy input to output. */
  if (fits_open_data(&fin, input, READONLY, &status)) FITS_EXIT;
  if (fits_create_file(&fptr, fname, &status)) FITS_EXIT;
  if (fits_copy_file(fin, fptr, 1, 1, 1, &status)) FITS_EXIT;
  if (fits_close_file(fin, &status)) FITS_EXIT;

  /* Append columns to the output. */
  if (fits_get_num_cols(fptr, &col, &status)) FITS_EXIT;
  if (fits_insert_cols(fptr, col + 1, nnew, ttype, tform, &status)) FITS_EXIT;

  MY_ALLOC(mskcol, unsigned char, (num * nbyte), the output column);
#if NUMSUB != 1
  if (sizeof(subid_t) == nbyte) {
    MY_ALLOC(subcol, subid_t, num, the output column);
  }
  else subcol = (subid_t *) mskcol;
#endif

  /* Generate maskbit array, taking into account the reordering of data */
  switch (masktype) {
    case TBYTE:
      for (i = 0; i < num; i++)
        ((uint8_t *) mskcol)[data[i].oridx] = data[i].mask;
      break;
    case TSHORT:
      for (i = 0; i < num; i++)
        ((uint16_t *) mskcol)[data[i].oridx] = data[i].mask;
      break;
    case TINT:
      for (i = 0; i < num; i++)
        ((uint32_t *) mskcol)[data[i].oridx] = data[i].mask;
      break;
    case TLONG:
      for (i = 0; i < num; i++)
        ((uint64_t *) mskcol)[data[i].oridx] = data[i].mask;
      break;
  }
  if (fits_write_col(fptr, TBYTE, col + 1, 1, 1, num, mskcol, &status))
    FITS_EXIT;

#if NUMSUB != 1
  for (i = 0; i < num; i++) subcol[data[i].oridx] = data[i].flag;
  if (fits_write_col(fptr, TBYTE, col + 2, 1, 1, num, subcol, &status))
    FITS_EXIT;
#endif

  if (fits_close_file(fptr, &status)) FITS_EXIT;
  free(mskcol);
#if NUMSUB != 1
  if (sizeof(subid_t) == nbyte) {
    free(subcol);
  }
#endif
  return 0;
}

