/*******************************************************************************
* define.h: this file is part of the brickmask program.

* brickmask: assign bit codes defined on Legacy Survey brick pixels
             to a catalogue with sky coordinates.

* Github repository:
        https://github.com/cheng-zhao/brickmask

* Copyright (c) 2021 Cheng Zhao <zhaocheng03@gmail.com> and Andrei Variu

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

*******************************************************************************/

#ifndef __DEFINE_H__
#define __DEFINE_H__

/*============================================================================*\
                 Definitions of mathematical/physical constants
\*============================================================================*/
#define DEGREE_2_RAD    0x1.1df46a2529d39p-6    /* M_PI / 180 */
#define RAD_2_DEGREE    0x1.ca5dc1a63c1f8p+5    /* 180 / M_PI */

/*============================================================================*\
                         Definitions for configurations
\*============================================================================*/
/* Default value for unset parameters. */
#define DEFAULT_CONF_FILE               "brickmask.conf"
#define DEFAULT_FILE_TYPE               BRICKMASK_FFMT_ASCII
#define DEFAULT_ASCII_COMMENT           '\0'
#define DEFAULT_OVERWRITE               0
#define DEFAULT_VERBOSE                 true

#ifdef EBOSS
#define DEFAULT_MASK_NULL               0
#else
#define DEFAULT_MASK_NULL               1
#endif

#define BRICKMASK_MAX_SUBID             UCHAR_MAX
#define BRICKMASK_MAX_COLUMN            65536

/* Priority of parameters from different sources. */
#define BRICKMASK_PRIOR_CMD             5
#define BRICKMASK_PRIOR_FILE            1

/*============================================================================*\
                            Definitions for file IO
\*============================================================================*/
#define BRICKMASK_PATH_SEP      '/'     /* separator for file paths     */
#define BRICKMASK_FILE_CHUNK    1048576 /* chunk size for ASCII file IO */
#define BRICKMASK_MAX_CHUNK     INT_MAX /* maximum allowed chunk size   */
#define BRICKMASK_READ_COMMENT  '#'     /* comment symbol for reading   */
/* Initial number of objects allocated for the input catalogue.           */
#define BRICKMASK_DATA_INIT_NUM                 128
/* Initial size for columns other than coordinates to be read from file   */
#define BRICKMASK_CONTENT_INIT_SIZE             1024
/* Maximum content doubling size                                          */
#define BRICKMASK_CONTENT_MAX_DOUBLE_SIZE       INT_MAX
/* Maximum content size                                                   */
#define BRICKMASK_CONTENT_MAX_SIZE              SIZE_MAX

/*============================================================================*\
                            Other runtime constants
\*============================================================================*/
#define BRICKMASK_CODE_NAME     "BRICKMASK"     /* name of the program */
#define BRICKMASK_SPACE_ESCAPE  '\\'    /* escape character for spaces */
#define BRICKMASK_TOL   1e-9    /* tolerance for coordinate comparison */
/* Names of FITS columns. */
#define BRICKMASK_FITS_BRICKNAME        "BRICKNAME"
#define BRICKMASK_FITS_RAMIN            "RA1"
#define BRICKMASK_FITS_RAMAX            "RA2"
#define BRICKMASK_FITS_DECMIN           "DEC1"
#define BRICKMASK_FITS_DECMAX           "DEC2"
#define BRICKMASK_FITS_SUBID            "SUBID"
/* Maximum length of FITS columns. */
#define BRICKMASK_FITS_MAX_COLNAME      32
/* Case sensitivity of FITS columns. */
#define BRICKMASK_FITS_CASESEN          CASEINSEN
/* Minimum FITS memory file size. */
#define BRICKMASK_FITS_MIN_MEMSIZE      1048576
/* Number of revisions for showing progress. */
#define BRICKMASK_PROGRESS_NUM          100

#ifdef EBOSS
#define EBOSS_MASK_VALID(bit)           ((bit) & 1)
#define EBOSS_XYBUG_BIT                 4
#define EBOSS_XYBUG_VALID(bit)          ((bit) & (EBOSS_XYBUG_BIT))
#endif

#ifdef MPI
#define BRICKMASK_MPI_ROOT              0       /* root rank of MPI */
#define BRICKMASK_MAX_DATA              INT_MAX /* maximum number of objects */
#endif

/*============================================================================*\
                     Definitions for the format of outputs
\*============================================================================*/
#define FMT_WARN "\n\x1B[35;1mWarning:\x1B[0m"          /* Magenta "Warning" */
#define FMT_ERR  "\n\x1B[31;1mError:\x1B[0m"            /* Red "Error"       */
#define FMT_EXIT "\x1B[31;1mExit:\x1B[0m"               /* Red "Exit"        */
#define FMT_DONE "\r\x1B[70C[\x1B[32;1mDONE\x1B[0m]\n"  /* Green "DONE"      */
#define FMT_FAIL "\r\x1B[70C[\x1B[31;1mFAIL\x1B[0m]\n"  /* Red "FAIL"        */
#define FMT_KEY(key)    "\x1B[36;1m" #key "\x1B[0m"     /* Cyan keyword      */
#define OFMT_DBL "%.10lg"             /* Output format for double parameters */

/*============================================================================*\
                          Definitions for error codes
\*============================================================================*/
#define BRICKMASK_ERR_MEMORY            (-1)
#define BRICKMASK_ERR_CFG               (-4)
#define BRICKMASK_ERR_BRICK             (-2)
#define BRICKMASK_ERR_FILE              (-3)
#define BRICKMASK_ERR_INIT              (-5)
#define BRICKMASK_ERR_MASK              (-6)
#define BRICKMASK_ERR_MPI               (-7)
#define BRICKMASK_ERR_SAVE              (-12)
#define BRICKMASK_ERR_UNKNOWN           (-99)

/*============================================================================*\
                           Definitions for shortcuts
\*============================================================================*/
#define P_ERR(...) fprintf(stderr, FMT_ERR " " __VA_ARGS__)
#define P_WRN(...) fprintf(stderr, FMT_WARN " " __VA_ARGS__)
#define P_EXT(...) fprintf(stderr, FMT_EXIT " " __VA_ARGS__)

#endif

