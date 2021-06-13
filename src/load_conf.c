/*******************************************************************************
* load_conf.c: this file is part of the brickmask program.

* brickmask: assign bit codes defined on Legacy Survey brick pixels
             to a catalogue with sky coordinates.

* Github repository:
        https://github.com/cheng-zhao/brickmask

* Copyright (c) 2020 -- 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "load_conf.h"
#include "data_io.h"
#include "libcfg.h"
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>

/*============================================================================*\
                           Macros for error handling
\*============================================================================*/
/* Check existence of configuration parameters. */
#define CHECK_EXIST_PARAM(name, cfg, var)                       \
  if (!cfg_is_set((cfg), (var))) {                              \
    P_ERR(FMT_KEY(name) " is not set\n");                       \
    return BRICKMASK_ERR_CFG;                                   \
  }
#define CHECK_EXIST_ARRAY(name, cfg, var, num)                  \
  if (!((num) = cfg_get_size((cfg), (var)))) {                  \
    P_ERR(FMT_KEY(name) " is not set\n");                       \
    return BRICKMASK_ERR_CFG;                                   \
  }

/* Check length of array. */
#define CHECK_ARRAY_LENGTH(name, cfg, var, fmt, num, nexp)      \
  if ((num) < (nexp)) {                                         \
    P_ERR("too few elements of " FMT_KEY(name) "\n");           \
    return BRICKMASK_ERR_CFG;                                   \
  }                                                             \
  if ((num) > (nexp)) {                                         \
    P_WRN("omitting the following " FMT_KEY(name) ":");         \
    for (int i = (nexp); i < (num); i++)                        \
      fprintf(stderr, " " fmt, (var)[i]);                       \
    fprintf(stderr, "\n");                                      \
  }

#define CHECK_STR_ARRAY_LENGTH(name, cfg, var, num, nexp)       \
  if ((num) < (nexp)) {                                         \
    P_ERR("too few elements of " FMT_KEY(name) "\n");           \
    return BRICKMASK_ERR_CFG;                                   \
  }                                                             \
  if ((num) > (nexp)) {                                         \
    P_WRN("omitting the following " FMT_KEY(name) ":\n");       \
    for (int i = (nexp); i < (num); i++)                        \
      fprintf(stderr, "  %s\n", (var)[i]);                      \
  }

/* Release memory for configuration parameters. */
#define FREE_ARRAY(x)           {if(x) free(x);}
#define FREE_STR_ARRAY(x)       {if(x) {if (*(x)) free(*(x)); free(x);}}

/* Print the warning and error messages. */
#define P_CFG_WRN(cfg)  cfg_pwarn(cfg, stderr, FMT_WARN);
#define P_CFG_ERR(cfg)  {                                       \
  cfg_perror(cfg, stderr, FMT_ERR);                             \
  cfg_destroy(cfg);                                             \
  return NULL;                                                  \
}


/*============================================================================*\
                    Functions called via command line flags
\*============================================================================*/

/******************************************************************************
Function `usage`:
  Print the usage of command line options.
******************************************************************************/
static void usage(void *args) {
  (void) args;
  printf("Usage: " BRICKMASK_CODE_NAME " [OPTION]\n\
Assign maskbit codes to a catalog with sky coordinates.\n\
  -h, --help\n\
        Display this message and exit\n\
  -t, --template\n\
        Print a template configuration file to the standard output and exit\n\
  -c, --conf            " FMT_KEY(CONFIG_FILE) "     String\n\
        Specify the configuration file (default: `%s')\n\
  -l, --brick-list      " FMT_KEY(BRICK_LIST) "      String\n\
        Specify the FITS table listing all survey bricks\n\
  -m, --mask-file       " FMT_KEY(MASKBIT_FILES) "   String array\n\
        Specify text files with paths of maskbit files\n\
  -n, --mask-null       " FMT_KEY(MASKBIT_NULL) "    Integer\n\
        Set the maskbit code for objects outside all maskbit bricks\n\
  -s, --sample-id       " FMT_KEY(SUBSAMPLE_ID) "    Integer array\n\
        Set IDs of subsamples\n\
  -i, --input           " FMT_KEY(INPUT) "           String\n\
        Specify the input catalog\n\
  -f, --file-type       " FMT_KEY(FILE_TYPE) "       Integer\n\
        Specify the file type of the input catalog\n\
      --comment         " FMT_KEY(ASCII_COMMENT) "   Character\n\
        Specify the comment symbol for ASCII-format input catalog\n\
  -C, --coord-col       " FMT_KEY(COORD_COLUMN) "    String array\n\
        Specify columns for RA and Dec in the input catalog\n\
  -o, --output          " FMT_KEY(OUTPUT) "          String\n\
        Set the output catalog\n\
  -e, --output-col      " FMT_KEY(OUTPUT_COLUMN) "   String array\n\
        Set columns to be written to the output catalog\n\
  -O, --overwrite       " FMT_KEY(OVERWRITE) "       Integer\n\
        Indicate whether to overwrite existing output files\n\
  -v, --verbose         " FMT_KEY(VERBOSE) "         Boolean\n\
        Indicate whether to display detailed standard outputs\n\
Consult the -t option for more information on the parameters\n\
Github repository: https://github.com/cheng-zhao/brickmask\n\
Licence: MIT\n",
    DEFAULT_CONF_FILE);
  exit(0);
}

/******************************************************************************
Function `conf_template`:
  Print a template configuration file.
******************************************************************************/
static void conf_template(void *args) {
  (void) args;
  printf("# Configuration file for " BRICKMASK_CODE_NAME " (default: `"
DEFAULT_CONF_FILE "').\n\
# Format: keyword = value # comment\n\
#     or: keyword = [element1, element2]\n\
#    see: https://github.com/cheng-zhao/libcfg for details.\n\
# Some of the entries allow expressions, see\n\
#         https://github.com/cheng-zhao/libast for details.\n\
# NOTE that command line options have priority over this file.\n\
# Unnecessary entries can be left unset.\n\
\n\
BRICK_LIST      = \n\
    # Filename for the FITS table with the list of all bricks, see e.g.\n\
    # https://www.legacysurvey.org/dr9/files/#survey-bricks-fits-gz\n\
MASKBIT_FILES   = \n\
    # String or string array, ASCII files with the paths of maskbit files.\n\
    # Each element specifies maskbit files for a subsample, such as NGC or SGC.\n\
    # Each row of the ASCII files specifies the path of a maskbit file.\n\
    # Each space in the paths must be escaped by a leading '\\' character.\n\
    # Name of the bricks must present in the filenames.\n\
    # Lines starting with '#' are omitted.\n\
MASKBIT_NULL    = \n\
    # Integer, bit code for objects outisde all maskbit bricks (unset: %d).\n\
SUBSAMPLE_ID    = \n\
    # If set, the IDs of subsamples are saved to the output as an extra column.\n\
    # Integer or integer array, same dimension as `MASKBIT_FILES`.\n\
INPUT           = \n\
    # Filename of the input data catalog.\n\
FILE_TYPE       = \n\
    # Integer, format of the input catalog (default: %d).\n\
    # The allowed values are:\n\
    # * %d: ASCII text file;\n\
    # * %d: FITS table.\n\
ASCII_COMMENT   = \n\
    # Character indicating comment lines for ASCII-format catalog (unset: '%c%s.\n\
COORD_COLUMN    = \n\
    # 2-element integer or string array, columns of (RA,Dec) for `INPUT`.\n\
    # They must be integers indicating the column numbers (starting from 1) for\n\
    # an ASCII file, or strings indicating the column names for a FITS table.\n\
OUTPUT          = \n\
    # Filename for the output catalog, with the same format as `INPUT`.\n\
OUTPUT_COLUMN   = \n\
    # Integer or String arrays, columns to be saved to `OUTPUT`.\n\
    # If not set, all columns of `INPUT` are saved in the original order.\n\
    # Note that maskbits (and optionally subsample IDs) are always saved\n\
    # as the last column (or last two columns).\n\
OVERWRITE       = \n\
    # Flag indicating whether to overwrite existing files, integer (unset: %d).\n\
    # Allowed values are:\n\
    # * 0: quit the program when an output file exist;\n\
    # * positive: force overwriting output files whenever possible;\n\
    # * negative: notify at most this number of times for existing files.\n\
VERBOSE         = \n\
    # Boolean option, indicate whether to show detailed outputs (unset: %c).\n",
      DEFAULT_MASK_NULL, DEFAULT_FILE_TYPE,
      BRICKMASK_FFMT_ASCII, BRICKMASK_FFMT_FITS,
      DEFAULT_ASCII_COMMENT ? DEFAULT_ASCII_COMMENT : '\'',
      DEFAULT_ASCII_COMMENT ? "')" : ")",
      DEFAULT_OVERWRITE, DEFAULT_VERBOSE ? 'T' : 'F');
  exit(0);
}


/*============================================================================*\
                      Function for reading configurations
\*============================================================================*/

/******************************************************************************
Function `conf_init`:
  Initialise the structure for storing configurations.
Return:
  Address of the structure.
******************************************************************************/
static CONF *conf_init(void) {
  CONF *conf = calloc(1, sizeof *conf);
  if (!conf) return NULL;
  conf->fconf = conf->flist = conf->input = conf->output = NULL;
  conf->fmask = conf->cname = conf->ocol = NULL;
  conf->subid = conf->onum = NULL;
  return conf;
}

/******************************************************************************
Function `conf_read`:
  Read configurations.
Arguments:
  * `conf`:     structure for storing configurations;
  * `argc`:     number of arguments passed via command line;
  * `argv`:     array of command line arguments.
Return:
  Interface of libcfg.
******************************************************************************/
static cfg_t *conf_read(CONF *conf, const int argc, char *const *argv) {
  if (!conf) {
    P_ERR("the structure for configurations is not initialised\n");
    return NULL;
  }
  cfg_t *cfg = cfg_init();
  if (!cfg) P_CFG_ERR(cfg);

  /* Functions to be called via command line flags. */
  const cfg_func_t funcs[] = {
    {   'h',        "help",             usage,          NULL},
    {   't',    "template",     conf_template,          NULL}
  };

  /* Configuration parameters. */
  const cfg_param_t params[] = {
    {'c', "conf"        , "CONFIG_FILE"    , CFG_DTYPE_STR , &conf->fconf   },
    {'l', "brick-list"  , "BRICK_LIST"     , CFG_DTYPE_STR , &conf->flist   },
    {'m', "mask-file"   , "MASKBIT_FILES"  , CFG_ARRAY_STR , &conf->fmask   },
    {'n', "mask-null"   , "MASKBIT_NULL"   , CFG_DTYPE_INT , &conf->mnull   },
    {'s', "sample-id"   , "SUBSAMPLE_ID"   , CFG_ARRAY_INT , &conf->subid   },
    {'i', "input"       , "INPUT"          , CFG_DTYPE_STR , &conf->input   },
    {'f', "file-type"   , "FILE_TYPE"      , CFG_DTYPE_INT , &conf->ftype   },
    { 0 , "comment"     , "ASCII_COMMENT"  , CFG_DTYPE_CHAR, &conf->comment },
    {'C', "coord-col"   , "COORD_COLUMN"   , CFG_ARRAY_STR , &conf->cname   },
    {'o', "output"      , "OUTPUT"         , CFG_DTYPE_STR , &conf->output  },
    {'e', "output-col"  , "OUTPUT_COLUMN"  , CFG_ARRAY_STR , &conf->ocol    },
    {'O', "overwrite"   , "OVERWRITE"      , CFG_DTYPE_INT , &conf->ovwrite },
    {'v', "verbose"     , "VERBOSE"        , CFG_DTYPE_BOOL, &conf->verbose }
  };

  /* Register functions and parameters. */
  if (cfg_set_funcs(cfg, funcs, sizeof(funcs) / sizeof(funcs[0])))
    P_CFG_ERR(cfg);
  P_CFG_WRN(cfg);
  if (cfg_set_params(cfg, params, sizeof(params) / sizeof(params[0])))
    P_CFG_ERR(cfg);
  P_CFG_WRN(cfg);

  /* Read configurations from command line options. */
  int optidx;
  if (cfg_read_opts(cfg, argc, argv, BRICKMASK_PRIOR_CMD, &optidx))
    P_CFG_ERR(cfg);
  P_CFG_WRN(cfg);

  /* Read parameters from configuration file. */
  if (!cfg_is_set(cfg, &conf->fconf)) conf->fconf = DEFAULT_CONF_FILE;
  if (access(conf->fconf, R_OK))
    P_WRN("cannot access the configuration file: `%s'\n", conf->fconf);
  else if (cfg_read_file(cfg, conf->fconf, BRICKMASK_PRIOR_FILE))
    P_CFG_ERR(cfg);
  P_CFG_WRN(cfg);

  return cfg;
}


/*============================================================================*\
                      Functions for parameter verification
\*============================================================================*/

/******************************************************************************
Function `check_input`:
  Check whether an input file can be read.
Arguments:
  * `fname`:    filename of the input file;
  * `key`:      keyword of the input file.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static inline int check_input(const char *fname, const char *key) {
  if (!fname || *fname == '\0') {
    P_ERR(FMT_KEY(%s) " is not set\n", key);
    return BRICKMASK_ERR_CFG;
  }
  if (access(fname, R_OK)) {
    P_ERR("cannot access " FMT_KEY(%s) ": `%s'\n", key, fname);
    return BRICKMASK_ERR_FILE;
  }
  return 0;
}

/******************************************************************************
Function `check_output`:
  Check whether an output file can be written.
Arguments:
  * `fname`:    filename of the input file;
  * `key`:      keyword of the input file;
  * `ovwrite`:  option for overwriting exisiting files.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int check_output(char *fname, const char *key, const int ovwrite) {
  if (!fname || *fname == '\0') {
    P_ERR(FMT_KEY(%s) " is not set\n", key);
    return BRICKMASK_ERR_CFG;
  }

  /* Check if the file exists. */
  if (!access(fname, F_OK)) {
    if (ovwrite == 0) {         /* not overwriting */
      P_ERR(FMT_KEY(%s) " exists: `%s'\n", key, fname);
      return BRICKMASK_ERR_FILE;
    }
    else if (ovwrite > 0) {     /* force overwriting */
      P_WRN(FMT_KEY(%s) " will be overwritten: `%s'\n",
          key, fname);
    }
    else {                      /* ask for decision */
      P_WRN(FMT_KEY(%s) " exists: `%s'\n", key, fname);
      char confirm = 0;
      for (int i = 0; i != ovwrite; i--) {
        fprintf(stderr, "Are you going to overwrite it? (y/n): ");
        if (scanf("%c", &confirm) != 1) continue;
        int c;
        while((c = getchar()) != '\n' && c != EOF) continue;
        if (confirm == 'n') {
          P_ERR("cannot write to the file\n");
          return BRICKMASK_ERR_FILE;
        }
        else if (confirm == 'y') break;
      }
      if (confirm != 'y') {
        P_ERR("too many failed inputs\n");
        return BRICKMASK_ERR_FILE;
      }
    }

    /* Check file permission for overwriting. */
    if (access(fname, W_OK)) {
      P_ERR("cannot write to file `%s'\n", fname);
      return BRICKMASK_ERR_FILE;
    }
  }
  /* Check the path permission. */
  else {
    char *end;
    if ((end = strrchr(fname, BRICKMASK_PATH_SEP)) != NULL) {
      *end = '\0';
      if (access(fname, X_OK)) {
        P_ERR("cannot access the directory `%s'\n", fname);
        return BRICKMASK_ERR_FILE;
      }
      *end = BRICKMASK_PATH_SEP;
    }
  }
  return 0;
}

/******************************************************************************
Function `conf_verify`:
  Verify configuration parameters.
Arguments:
  * `cfg`:      interface of libcfg;
  * `conf`:     structure for storing configurations.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
static int conf_verify(const cfg_t *cfg, CONF *conf) {
  int e, num;

  /* BRICK_LIST */
  CHECK_EXIST_PARAM(BRICK_LIST, cfg, &conf->flist);
  if ((e = check_input(conf->flist, "BRICK_LIST"))) return e;

  /* MASKBIT_FILES */
  CHECK_EXIST_ARRAY(MASKBIT_FILES, cfg, &conf->fmask, conf->nsub);
  for (int i = 0; i < conf->nsub; i++) {
    if ((e = check_input(conf->fmask[i], "MASKBIT_FILES"))) return e;
  }

  /* MASKBIT_NULL */
  if (!cfg_is_set(cfg, &conf->mnull)) conf->mnull = DEFAULT_MASK_NULL;
  if (conf->mnull < 0) {
    P_ERR(FMT_KEY(MASKBIT_NULL) " must be non-negative\n");
    return BRICKMASK_ERR_CFG;
  }

  /* SUBSAMPLE_ID */
  if ((num = cfg_get_size(cfg, &conf->subid))) {
    CHECK_ARRAY_LENGTH(SUBSAMPLE_ID, cfg, conf->subid, "%d", num, conf->nsub);
    for (int i = 0; i < conf->nsub; i++) {
      if (conf->subid[i] < 0 || conf->subid[i] > BRICKMASK_MAX_SUBID) {
        P_ERR(FMT_KEY(SUBSAMPLE_ID) " must be between 0 and %d\n",
            BRICKMASK_MAX_SUBID);
        return BRICKMASK_ERR_CFG;
      }
    }
  }

  /* INPUT */
  CHECK_EXIST_PARAM(INPUT, cfg, &conf->input);
  if ((e = check_input(conf->input, "INPUT"))) return e;

  /* FILE_TYPE */
  if (!cfg_is_set(cfg, &conf->ftype)) conf->ftype = DEFAULT_FILE_TYPE;
  switch (conf->ftype) {
    case BRICKMASK_FFMT_ASCII:
      /* ASCII_COMMENT */
      if (!cfg_is_set(cfg, &conf->comment))
        conf->comment = DEFAULT_ASCII_COMMENT;
      if (conf->comment && !isgraph(conf->comment)) {
        P_ERR("invalid " FMT_KEY(ASCII_COMMENT) ": '%c' (ASCII code: %d)\n",
            conf->comment, conf->comment);
        return BRICKMASK_ERR_CFG;
      }
      break;
    case BRICKMASK_FFMT_FITS:
      break;
    default:
      P_ERR("invalid " FMT_KEY(FILE_TYPE) ": %d\n", conf->ftype);
      return BRICKMASK_ERR_CFG;
  }

  /* COORD_COLUMN */
  CHECK_EXIST_ARRAY(COORD_COLUMN, cfg, &conf->cname, num);
  CHECK_ARRAY_LENGTH(COORD_COLUMN, cfg, conf->cname, "%s", num, 2);
  if (conf->ftype == BRICKMASK_FFMT_ASCII) {
    /* Convert strings to integers. */
    for (int i = 0; i < 2; i++) {
      conf->cnum[i] = 0;
      if (sscanf(conf->cname[i], "%d", conf->cnum + i) != 1) {
        P_ERR(FMT_KEY(COORD_COLUMN) " must be integers\n");
        return BRICKMASK_ERR_CFG;
      }
      if (conf->cnum[i] <= 0 || conf->cnum[i] > BRICKMASK_MAX_COLUMN) {
        P_ERR(FMT_KEY(COORD_COLUMN) " must be postive and not larger than %d\n",
            BRICKMASK_MAX_COLUMN);
        return BRICKMASK_ERR_CFG;
      }
    }
    if (conf->cnum[0] == conf->cnum[1]) {
      P_ERR("identical RA and Dec columns: %d\n", conf->cnum[0]);
      return BRICKMASK_ERR_CFG;
    }
  }
  else if (!strcmp(conf->cname[0], conf->cname[1])) {
    P_ERR("identical RA and Dec columns: %s\n", conf->cname[0]);
    return BRICKMASK_ERR_CFG;
  }

  /* OVERWRITE */
  if (!cfg_is_set(cfg, &conf->ovwrite)) conf->ovwrite = DEFAULT_OVERWRITE;

  /* OUTPUT */
  CHECK_EXIST_PARAM(OUTPUT, cfg, &conf->output);
  if ((e = check_output(conf->output, "OUTPUT", conf->ovwrite))) return e;

  /* OUTPUT_COLUMN */
  if ((conf->ncol = cfg_get_size(cfg, &conf->ocol))) {
    if (conf->ftype == BRICKMASK_FFMT_ASCII) {
      if (!(conf->onum = calloc(conf->ncol, sizeof(int)))) {
        P_ERR("failed to allocate memory for " FMT_KEY(OUTPUT_COLUMN) "\n");
        return BRICKMASK_ERR_MEMORY;
      }
      bool find_ra, find_dec;
      find_ra = find_dec = false;
      for (int i = 0; i < conf->ncol; i++) {
        if (sscanf(conf->ocol[i], "%d", conf->onum + i) != 1) {
          P_ERR(FMT_KEY(OUTPUT_COLUMN) " must be integers\n");
          return BRICKMASK_ERR_CFG;
        }
        if (conf->onum[i] <= 0 || conf->onum[i] > BRICKMASK_MAX_COLUMN) {
          P_ERR(FMT_KEY(OUTPUT_COLUMN) " must be positive and not larger than "
              "%d\n", BRICKMASK_MAX_COLUMN);
          return BRICKMASK_ERR_CFG;
        }
        /* Check if RA and Dec are available. */
        if (conf->onum[i] == conf->cnum[0]) find_ra = true;
        else if (conf->onum[i] == conf->cnum[1]) find_dec = true;
      }
      if (!find_ra)
        P_WRN("Right ascension not in " FMT_KEY(OUTPUT_COLUMN) "\n");
      if (!find_dec) P_WRN("Declination not in " FMT_KEY(OUTPUT_COLUMN) "\n");

      /* Check duplicates. */
      for (int i = 1; i < conf->ncol; i++) {
        for (int j = 0; j < i; j++) {
          if (conf->onum[i] == conf->onum[j]) {
            P_ERR("duplicate " FMT_KEY(OUTPUT_COLUMN) ": %d\n", conf->onum[i]);
            return BRICKMASK_ERR_CFG;
          }
        }
      }
    }
    else {
      /* Check if RA and Dec are available. */
      bool find_ra, find_dec;
      find_ra = find_dec = false;
      for (int i = 0; i < conf->ncol; i++) {
        if (!strcmp(conf->ocol[i], conf->cname[0])) find_ra = true;
        else if (!strcmp(conf->ocol[i], conf->cname[1])) find_dec = true;
      }
      if (!find_ra)
        P_WRN("Right ascension not in " FMT_KEY(OUTPUT_COLUMN) "\n");
      if (!find_dec) P_WRN("Declination not in " FMT_KEY(OUTPUT_COLUMN) "\n");

      /* Check duplicates. */
      for (int i = 1; i < conf->ncol; i++) {
        for (int j = 0; j < i; j++) {
          if (!strcmp(conf->ocol[i], conf->ocol[j])) {
            P_ERR("duplicate " FMT_KEY(OUTPUT_COLUMN) ": %s\n", conf->ocol[i]);
            return BRICKMASK_ERR_CFG;
          }
        }
      }
    }
  }

  /* VERBOSE */
  if (!cfg_is_set(cfg, &conf->verbose)) conf->verbose = DEFAULT_VERBOSE;

  return 0;
}


/*============================================================================*\
                      Function for printing configurations
\*============================================================================*/

/******************************************************************************
Function `conf_print`:
  Print configuration parameters.
Arguments:
  * `conf`:     structure for storing configurations.
******************************************************************************/
static void conf_print(const CONF *conf) {
  /* Configuration file */
  printf("\n  CONFIG_FILE     = %s", conf->fconf);

  /* Input catalogs. */
  printf("\n  BRICK_LIST      = %s", conf->flist);
  printf("\n  MASKBIT_FILES   = %s", conf->fmask[0]);
  for (int i = 1; i < conf->nsub; i++)
    printf("\n                    %s", conf->fmask[i]);
  if (conf->subid) {
    printf("\n  SUBSAMPLE_ID    = %d", conf->subid[0]);
    for (int i = 1; i < conf->nsub; i++) printf(" , %d", conf->subid[i]);
  }
  printf("\n  INPUT           = %s", conf->input);

  const char *ftype[2] = {"ASCII", "FITS"};
  printf("\n  FILE_TYPE       = %d (%s)", conf->ftype, ftype[conf->ftype]);
  if (conf->ftype == BRICKMASK_FFMT_ASCII) {
    if (conf->comment == 0) printf("\n  ASCII_COMMENT   = ''");
    else printf("\n  ASCII_COMMENT   = '%c'", conf->comment);
    printf("\n  COORD_COLUMN    = %d , %d", conf->cnum[0], conf->cnum[1]);
  }
  else printf("\n  COORD_COLUMN    = %s , %s", conf->cname[0], conf->cname[1]);

  printf("\n  OUTPUT          = %s", conf->output);
  if (conf->ncol) {
    if (conf->ftype == BRICKMASK_FFMT_ASCII) {
      printf("\n  OUTPUT_COLUMN   = %d", conf->onum[0]);
      for (int i = 1; i < conf->ncol; i++) printf(" , %d", conf->onum[i]);
    }
    else {
      printf("\n  OUTPUT_COLUMN   = %s", conf->ocol[0]);
      for (int i = 1; i < conf->ncol; i++) printf(" , %s", conf->ocol[i]);
    }
  }

  printf("\n  OVERWRITE       = %d\n", conf->ovwrite);
}


/*============================================================================*\
                      Interface for loading configurations
\*============================================================================*/

/******************************************************************************
Function `load_conf`:
  Read, check, and print configurations.
Arguments:
  * `argc`:     number of arguments passed via command line;
  * `argv`:     array of command line arguments.
Return:
  The structure for storing configurations.
******************************************************************************/
CONF *load_conf(const int argc, char *const *argv) {
  CONF *conf = conf_init();
  if (!conf) return NULL;

  cfg_t *cfg = conf_read(conf, argc, argv);
  if (!cfg) {
    conf_destroy(conf);
    return NULL;
  }

  printf("Loading configurations ...");
  fflush(stdout);

  if (conf_verify(cfg, conf)) {
    if (cfg_is_set(cfg, &conf->fconf)) free(conf->fconf);
    conf_destroy(conf);
    cfg_destroy(cfg);
    return NULL;
  }

  if (conf->verbose) conf_print(conf);

  if (cfg_is_set(cfg, &conf->fconf)) free(conf->fconf);
  cfg_destroy(cfg);

  printf(FMT_DONE);
  return conf;
}

/******************************************************************************
Function `conf_destroy`:
  Release memory allocated for the configurations.
Arguments:
  * `conf`:     the structure for storing configurations.
******************************************************************************/
void conf_destroy(CONF *conf) {
  if (!conf) return;
  FREE_ARRAY(conf->flist);
  FREE_STR_ARRAY(conf->fmask);
  FREE_ARRAY(conf->subid);
  FREE_ARRAY(conf->input);
  FREE_STR_ARRAY(conf->cname);
  FREE_ARRAY(conf->output);
  FREE_STR_ARRAY(conf->ocol);
  FREE_ARRAY(conf->onum);
  free(conf);
}
