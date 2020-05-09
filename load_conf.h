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

#ifndef _LOAD_CONF_H_
#define _LOAD_CONF_H_

#include "define.h"

typedef enum {
  parse_conf_start = 1,
  parse_conf_keyword = 2,
  parse_conf_equal = 3,
  parse_conf_value_start = 4,
  parse_conf_value_quote = 5,
  parse_conf_value = 6,
  parse_conf_done = 0,
  parse_conf_abort = -1
} PARSE_CONF_STATE;

typedef struct {
  int format;                           /* FILE_TYPE */
  int force;                            /* FORCE */
  int verbose;                          /* VERBOSE */
  char comment;                         /* COMMENT */
  char list[FLEN_FILENAME];             /* BRICK_LIST */
  char bdir[FLEN_FILENAME];             /* BRICK_DIR */
  char input[FLEN_FILENAME];            /* INPUT */
  char output[FLEN_FILENAME];           /* OUTPUT */
} CONF;


void init_conf(CONF *);

void read_opt(const int, char * const [], char *, CONF *);

int read_conf(const char *, CONF *);

int parse_conf(const char *, char *, char *, const size_t, const size_t,
    const size_t);

int check_conf(CONF *);

int check_input(const char *, const char *);

int check_output(const char *, const char *, const int);

void print_conf(const char *, const CONF *);

void temp_conf(void);

void usage(char *);

#endif

