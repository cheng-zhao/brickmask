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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <getopt.h>
#include "load_conf.h"
#include "brickmask.h"

/******************************************************************************
Function `init_conf`:
  Initialize the configuration by invalid values.

Arguments:
  * `conf`:     the structure for configuration.
******************************************************************************/
void init_conf(CONF *conf) {
  memset(conf->list, 0, FLEN_FILENAME);
  memset(conf->bdir, 0, FLEN_FILENAME);
  memset(conf->input, 0, FLEN_FILENAME);
  memset(conf->output, 0, FLEN_FILENAME);
  conf->comment = '\0';
  conf->format = conf->verbose = INIT_INT;
  conf->force = INIT_INT;
}


/******************************************************************************
Function `read_opt`:
  Read command line options.

Arguments:
  * `argc`:     number of command line options;
  * `argv`:     array of command line options;
  * `cfile`:    filename of the configuration file;
  * `conf`:     the structure for configuration.
******************************************************************************/
void read_opt(const int argc, char * const argv[], char *cfile, CONF *conf) {
  int opts, idx, n;
  const char *optstr = "htyc:l:d:i:f:s:o:";
  const struct option long_opts[] = {
    { "help",           no_argument,            NULL,   'h'     },
    { "template",       no_argument,            NULL,   't'     },
    { "force",          no_argument,            NULL,   'y'     },
    { "conf",           required_argument,      NULL,   'c'     },
    { "brick-list",     required_argument,      NULL,   'l'     },
    { "brick-dir",      required_argument,      NULL,   'd'     },
    { "input",          required_argument,      NULL,   'i'     },
    { "format",         required_argument,      NULL,   'f'     },
    { "comment",        required_argument,      NULL,   's'     },
    { "output",         required_argument,      NULL,   'o'     },
    { "verbose",        no_argument,  &conf->verbose,    1      },
    { "brief",          no_argument,  &conf->verbose,    0      },
    { 0, 0, 0, 0}
  };

  opts = idx = 0;
  while ((opts = getopt_long(argc, argv, optstr, long_opts, &idx)) != -1) {
    switch (opts) {
      case 0:
        break;
      case '?':
        P_EXT("please check the command line options.\n");
        exit(ERR_ARG);
      case 'h':
        usage(argv[0]);
        exit(0);
      case 't':
        temp_conf();
        exit(0);
      case 'y':
        conf->force = 1;
        break;
      case 'c':
        n = safe_strcpy(cfile, optarg, FLEN_FILENAME);
        CHECKSTR(n, FLEN_FILENAME, "argument too long: %s\n", optarg);
        break;
      case 'l':
        n = safe_strcpy(conf->list, optarg, FLEN_FILENAME);
        CHECKSTR(n, FLEN_FILENAME, "argument too long: %s\n", optarg);
        break;
      case 'd':
        n = safe_strcpy(conf->bdir, optarg, FLEN_FILENAME);
        CHECKSTR(n, FLEN_FILENAME, "argument too long: %s\n", optarg);
        break;
      case 'i':
        n = safe_strcpy(conf->input, optarg, FLEN_FILENAME);
        CHECKSTR(n, FLEN_FILENAME, "argument too long: %s\n", optarg);
        break;
      case 'f':
        sscanf(optarg, "%d", &conf->format);
        break;
      case 's':
        conf->comment = optarg[0];
        break;
      case 'o':
        n = safe_strcpy(conf->output, optarg, FLEN_FILENAME);
        CHECKSTR(n, FLEN_FILENAME, "argument too long: %s\n", optarg);
        break;
      default:
        break;
    }
  }

  if (optind < argc) {
    P_WRN("unknown command line options:\n ");
    while (optind < argc) printf(" %s", argv[optind++]);
    printf("\n");
  }
}


/******************************************************************************
Function `parse_conf`:
  Getting `KEYWORD` and `VALUE` from a line with the format
  `KEYWORD = VALUE # COMMENTS`.

Arguments:
  * `line`:     the input line;
  * `key`:      a string storing `KEYWORD`;
  * `value`:    a string storing `VALUE`;
  * `nline`:    maximum length of `line`;
  * `nkey`:     maximum length of `KEYWORD`;
  * `nval`:     maximum length of `VALUE`.
Return:
  Status of the parser (see enum `PARSE_CONF_STATE`).
******************************************************************************/
int parse_conf(const char *line, char *key, char *value,
    const size_t nline, const size_t nkey, const size_t nval) {
  int i, ik, iv;
  char c, quote;
  PARSE_CONF_STATE state;

  if (nline <= 1 || nkey <= 1 || nval <= 1) return parse_conf_abort;

  ik = iv = 0;
  quote = '\0';
  state = parse_conf_start;

  for (i = 0; i < nline; i++) {
    c = line[i];
    if (c == '\0' || c == '\n' || c == EOF) {
      key[ik] = value[iv] = '\0';
      if (state == parse_conf_value) state = parse_conf_done;
      break;
    }

    switch (state) {
      case parse_conf_start:
        if (isalpha(c) || c == '_') {
          key[ik++] = c;
          state = parse_conf_keyword;
        }
        else if (!isspace(c)) return state;
        break;
      case parse_conf_keyword:
        if (isspace(c)) {
          key[ik] = '\0';
          state = parse_conf_equal;
        }
        else if (c == '=') {
          state = parse_conf_value_start;
        }
        else if (isalnum(c) || c == '_') {
          if (ik >= nkey - 1) {
            key[ik] = '\0';
            return parse_conf_abort;
          }
          key[ik++] = c;
        }
        else return state;
        break;
      case parse_conf_equal:
        if (!isspace(c)) {
          if (c == '=') state = parse_conf_value_start;
          else return state;
        }
        break;
      case parse_conf_value_start:
        if (!isspace(c)) {
          if (c == '"' || c == '\'') {
            quote = c;
            state = parse_conf_value_quote;
          }
          else if (c == '#') return state;
          else if (isprint(c)) {
            value[iv++] = c;
            state = parse_conf_value;
          }
          else return state;
        }
        break;
      case parse_conf_value_quote:
        if (c == quote) {
          value[iv] = '\0';
          return parse_conf_done;
        }
        else {
          if (iv >= nval - 1) {
            value[iv] = '\0';
            return parse_conf_abort;
          }
          value[iv++] = c;
        }
        break;
      case parse_conf_value:
        if (isprint(c) && c != '#') {
          if (iv >= nval - 1) {
            value[iv] = '\0';
            return parse_conf_abort;
          }
          value[iv++] = c;
        }
        else {
          value[iv] = '\0';
          return parse_conf_done;
        }
        break;
      default:
        return state;
    }
  }

  return state;
}


/******************************************************************************
Function `read_conf`:
  Read the configuration file with a format of `KEYWORD = VALUE # COMMENTS`.

Arguments:
  * `fname`:    the filename of the configuration file;
  * `conf`:     the structure for configuration.
Return:
  A non-zero integer if there is problem.
******************************************************************************/
int read_conf(const char *fname, CONF *conf) {
  FILE *fp;
  char line[LEN_CONF_LINE];
  char keyword[LEN_CONF_KEY];
  char value[LEN_CONF_VALUE];
  int n;

  if (!(fp = fopen(fname, "r"))) {
    P_WRN("cannot open configuration file `%s'.\n", fname);
    return ERR_FILE;
  }
  while (fgets(line, LEN_CONF_LINE, fp) != NULL) {
    if (parse_conf(line, keyword, value, LEN_CONF_LINE, LEN_CONF_KEY,
          LEN_CONF_VALUE) != parse_conf_done) continue;

    if (!strncmp(keyword, "BRICK_LIST", LEN_CONF_KEY)
        && conf->list[0] == '\0') {
      n = safe_strcpy(conf->list, value, FLEN_FILENAME);
      CHECKSTR(n, FLEN_FILENAME, FMT_KEY(BRICK_LIST) " too long: %s", value);
    }
    else if (!strncmp(keyword, "BRICK_DIR", LEN_CONF_KEY)
        && conf->bdir[0] == '\0') {
      n = safe_strcpy(conf->bdir, value, FLEN_FILENAME);
      CHECKSTR(n, FLEN_FILENAME, FMT_KEY(BRICK_DIR) " too long: %s", value);
    }
    else if (!strncmp(keyword, "INPUT", LEN_CONF_KEY)
        && conf->input[0] == '\0') {
      n = safe_strcpy(conf->input, value, FLEN_FILENAME);
      CHECKSTR(n, FLEN_FILENAME, FMT_KEY(INPUT) " too long: %s", value);
    }
    else if (!strncmp(keyword, "OUTPUT", LEN_CONF_KEY)
        && conf->output[0] == '\0') {
      n = safe_strcpy(conf->output, value, FLEN_FILENAME);
      CHECKSTR(n, FLEN_FILENAME, FMT_KEY(OUTPUT) " too long: %s", value);
    }
    else if (!strncmp(keyword, "COMMENT", LEN_CONF_KEY)
        && conf->comment == '\0')
      conf->comment = value[0];
    else if (!strncmp(keyword, "FILE_TYPE", LEN_CONF_KEY)
        && conf->format == INIT_INT)
      sscanf(value, "%d", &conf->format);
    else if (!strncmp(keyword, "FORCE", LEN_CONF_KEY)
        && conf->force == INIT_INT)
      sscanf(value, "%d", &conf->force);
    else if (!strncmp(keyword, "VERBOSE", LEN_CONF_KEY)
        && conf->verbose == INIT_INT)
      sscanf(value, "%d", &conf->verbose);
  }

  fclose(fp);
  return 0;
}


/******************************************************************************
Function `check_conf`:
  Check the loaded configuration to see if the values are set correctly, and
  the input/output files are valid.

Arguments:
  * `conf`:     the structure for configuration.
Return:
  A non-zero integer if there is problem.
******************************************************************************/
int check_conf(CONF *conf) {
  int err;

  if ((err = check_input(conf->list, "BRICK_LIST"))) return err;
  if ((err = check_input(conf->bdir, "BRICK_DIR"))) return err;
  if ((err = check_input(conf->input, "INPUT"))) return err;

  if (conf->format < 0 || conf->format > 1) {
    P_WRN(FMT_KEY(FORMAT) " is not set correctly.\n"
        "Use the default value (%d) instead.\n", DEFAULT_FILE_TYPE);
    conf->format = DEFAULT_FILE_TYPE;
  }

  if (conf->format == 0) {
    if (conf->comment == '\0') {
      P_WRN(FMT_KEY(COMMENT) " is not set correctly.\n"
          "Use the default value (%c) instead.\n", DEFAULT_COMMENT);
      conf->comment = DEFAULT_COMMENT;
    }
  }

  if (conf->force == INIT_INT) conf->force = 0;
  if ((err = check_output(conf->output, "OUTPUT", conf->force))) return err;

  if (conf->verbose != 0 && conf->verbose != 1) {
    P_WRN(FMT_KEY(VERBOSE) " is not set correctly.\n"
        "Use the default value (%d) instead.\n", DEFAULT_VERBOSE);
    conf->verbose = DEFAULT_VERBOSE;
  }

  return 0;
}


/******************************************************************************
Function `check_input`:
  Check whether the input file can be read.

Arguments:
  * `fname`:    the filename of the input file;
  * `dscp`:     description of the input file.
Return:
  A non-zero integer if there is problem.
******************************************************************************/
int check_input(const char *fname, const char *dscp) {
  if (fname[0] == '\0' || fname[0] == ' ') {
    P_ERR("the input " FMT_KEY(%s) " is not set.\n", dscp);
    return ERR_FILE;
  }
  if (access(fname, R_OK)) {
    P_ERR("cannot open " FMT_KEY(%s) ": `%s'.\n", dscp, fname);
    return ERR_FILE;
  }
  return 0;
}


/******************************************************************************
Function `check_output`:
  Check whether the output file can be written. If the file already exists,
  prompt a choice for overwriting it.

Arguments:
  * `fname`:    the filename of the output file to be checked;
  * `dscp`:     description of the output file.
  * `force`:    non-zero for overwriting existing files without notifications.
Return:
  A non-zero integer if there is problem.
******************************************************************************/
int check_output(const char *fname, const char *dscp, const int force) {
  char confirm, *end;
  char path[FLEN_FILENAME];
  int cnt = 0;

  if (fname[0] == '\0' || fname[0] == ' ') {
    P_ERR("the output " FMT_KEY(%s) " is not set.\n", dscp);
    return ERR_FILE;
  }
  if (!access(fname, F_OK) && force == 0) {     /* if the output file exists */
    P_WRN("the output " FMT_KEY(%s) " `%s' exists.\n", dscp, fname);
    do {
      if ((++cnt) == TIMEOUT) {
        P_ERR("too many failed inputs.\n");
        return ERR_INPUT;
      }
      fprintf(stderr, "Are you going to overwrite it? (y/n): ");
      if (scanf("%c", &confirm) != 1) continue;
      while(getchar() != '\n');         /* Ignore invalid inputs. */
    }
    while (confirm != 'y' && confirm != 'n');
    if(confirm == 'n') {
      P_ERR("cannot write to the file.\n");
      return ERR_FILE;
    }
  }

  if (!access(fname, F_OK) && access(fname, W_OK)) {
    P_ERR("cannot write to file `%s'.\n", fname);
    return ERR_FILE;
  }

  /* Check the path of the output file. */
  safe_strcpy(path, fname, FLEN_FILENAME);
  if ((end = strrchr(path, '/')) != NULL) {
    *(end + 1) = '\0';
    if(access(path, X_OK)) {
      P_ERR("cannot access the path `%s'.\n", path);
      return ERR_FILE;
    }
  }

  return 0;
}


/******************************************************************************
Function `print_conf`:
  Print the relevant items for a given configuration.

Arguments:
  * `conf`:     the structure for configuration.
******************************************************************************/
void print_conf(const char *conf_file, const CONF *conf) {
  printf("\n  Filename: %s\n", conf_file);
  printf("  " FMT_KEY(BRICK_LIST) " = %s\n", conf->list);
  printf("  " FMT_KEY(BRICK_DIR) "  = %s\n", conf->bdir);
  printf("  " FMT_KEY(INPUT) "      = %s\n", conf->input);
  printf("  " FMT_KEY(FILE_TYPE) "  = %d\n", conf->format);
  if (conf->format == 0)
    printf("  " FMT_KEY(COMMENT) "    = %c\n", conf->comment);
  printf("  " FMT_KEY(OUTPUT) "     = %s\n", conf->output);
  printf("  " FMT_KEY(FORCE) "      = %d\n", conf->force);
  printf("  " FMT_KEY(VERBOSE) "    = %d\n", conf->verbose);
}


/******************************************************************************
Function `temp_conf`:
  Print a template configuration file.
******************************************************************************/
void temp_conf(void) {
  printf("# Configuration file (default: `%s').\n\
# NOTE that command line options have priority over this file.\n\
# FORMAT: KEYWORD = VALUE # COMMENT\n\
# Please surround VALUE by quotation marks if there are special characters.\n\
\n\
BRICK_LIST      = \n\
        # A fits file for the list of bricks, e.g.\n\
        # http://portal.nersc.gov/project/cosmo/data/legacysurvey/dr7/survey-bricks.fits.gz\n\
BRICK_DIR       = \n\
        # Directory for brick mask files.\n\
INPUT           = \n\
        # Input data catalog.\n\
FILE_TYPE       = %d\n\
        # Format of the input catalog (default: %d).\n\
        # The allowed values are:\n\
        # - 0: `ASCII` format, with the first two columns being RA and Dec;\n\
        # - 1: `fits` format;\n\
COMMENT         = '%c'\n\
        # Comment symbol for the input ASCII catalog (default: '%c').\n\
OUTPUT          = \n\
        # The output file with masks.\n\
FORCE           = %d\n\
        # Non-zero integer for force overwriting the output (default: %d).\n\
VERBOSE         = %d\n\
        # 0 for concise standard outputs; 1 for detailed outputs (default: %d).\
\n",
      DEFAULT_CONF_FILE, DEFAULT_FILE_TYPE, DEFAULT_FILE_TYPE,
      DEFAULT_COMMENT, DEFAULT_COMMENT, DEFAULT_OVERWRITE, DEFAULT_OVERWRITE,
      DEFAULT_VERBOSE, DEFAULT_VERBOSE);
}


/******************************************************************************
Function `usage`:
  Print the usage of command line options.

Arguments:
  * `pname`:    name of this program.
******************************************************************************/
void usage(char *pname) {
  printf("Usage: %s [OPTION [VALUE]]\n\
Apply veto masks to a input catalog according to mask bricks.\n\n\
  -h, --help\n\
        Display this message and exit.\n\
  -t, --template\n\
        Display a template configuration and exit.\n\
  -c, --conf=CONF_FILE\n\
        Set the configuration file to CONF_FILE.\n\
        The default configuration file is `%s'.\n\
  -l, --brick-list=BRICK_LIST\n\
        Set the file for brick lists to BRICK_LIST.\n\
  -d, --brick-dir=BRICK_DIR\n\
        Set the directory for brick masks to BRICK_DIR.\n\
  -i, --input=INPUT\n\
        Set the input data catalog to INPUT.\n\
  -f, --format=FILE_TYPE\n\
        Set the format of the input catalog to FILE_TYPE.\n\
  -s, --comment=COMMENT\n\
        Set the comment symbol for the input ASCII catalog to COMMENT.\n\
  -o, --output=OUTPUT\n\
        Set the output file to OUTPUT.\n\
  -y, --force\n\
        Force overwriting existing output files without notifications.\n\
      --verbose\n\
        Display detailed standard outputs.\n\
      --brief\n\
        Display concise standard outputs.\n\n\
Consult the -t option for more information on the configuraion.\n\
Report bugs to <zhaocheng03@gmail.com>.\n",
    pname, DEFAULT_CONF_FILE);
}

