# Configuration parameters for `BRICKMASK`

A template configuration file can be generated using the `-t`/`--template` command line options, e.g. [`brickmask.conf`](brickmask.conf).

### `BRICK_LIST` (`-l` / `--brick-list`)

Filename of a FITS table with the information of all bricks, including their coordinate ranges and names. See e.g. [https://www.legacysurvey.org/dr9/files/#survey-bricks-fits-gz](https://www.legacysurvey.org/dr9/files/#survey-bricks-fits-gz).

In particular, the brick list file for the eBOSS ELG sample can be downloaded here: [http://portal.nersc.gov/project/cosmo/data/legacysurvey/dr7/survey-bricks.fits.gz](http://portal.nersc.gov/project/cosmo/data/legacysurvey/dr7/survey-bricks.fits.gz).

### `MASKBIT_FILES` (`-m` / `--mask-file`)

Files containing paths of all maskbits files. Each file indicates a collection of maskbits files for a certain subsample, such as northern (NGC) or southern (SGC) galactic caps.

Each row of the file sets the path of a maskbits file, which must contain the name of the corresponding brick. Note that each white space in the path should be escaped by a leading '`\`' character. For instance, the path `/path with space/mask-0001m002.fits` should be provided as `/path\ with\ space/mask-0001m002.fits`.

A simple way to generate the file list is through the `ls` command, e.g. (the eBOSS DR16 ELG masks are available at [https://data.sdss.org/sas/dr16/eboss/lss/catalogs/DR16/ELGmasks](https://data.sdss.org/sas/dr16/eboss/lss/catalogs/DR16/ELGmasks))
```bash
$ ls DR16/ELGmasks/mask-eboss21-*.fits.gz > masks-eboss21.txt
$ ls DR16/ELGmasks/mask-eboss22-*.fits.gz > masks-eboss22.txt
$ ls DR16/ELGmasks/mask-eboss23-*.fits.gz > masks-eboss23.txt
$ ls DR16/ELGmasks/mask-eboss25-*.fits.gz > masks-eboss25.txt
```

The `MASKBIT_FILES` parameter can then be set as
```nginx
MASKBIT_FILES = [masks-eboss21.txt, masks-eboss22.txt,\
                 masks-eboss23.txt, masks-eboss25.txt]
```

Moreover, the script [`legacy_mask_list.py`](scripts/legacy_mask_list.py) is for generating lists of Legacy Survey maskbits files for both NGC and SGC on [NERSC](https://www.nersc.gov/).

### `MASKBIT_NULL` (`-n` / `--mask-null`)

Bit code for objects that are not found in any of the given maskbits files. It is recommended to set it to `1` for eBOSS ELG masks, and `0` for Legacy Survey masks (e.g. [https://www.legacysurvey.org/dr9/bitmasks](https://www.legacysurvey.org/dr9/bitmasks/)).

### `SUBSAMPLE_ID` (`-s` / `--sample-id`)

Optional parameter for specifying IDs of subsamples. If it is set, the corresponding values will be appended to the input catalogue, to indicate the subsample that each object belongs to. For instance, to identify eBOSS chunks of ELGs, one may set

```nginx
MASKBIT_FILES = [masks-eboss21.txt, masks-eboss22.txt,\
                 masks-eboss23.txt, masks-eboss25.txt]
SUBSAMPLE_ID  = [0,1,2,3]
    # 0 -> eboss21; 1 -> eboss22; 2 -> eboss23; 3 -> eboss25.
```

### `INPUT` (`-i` / `--input`)

Filename of the input catalogue, which must contain right ascensions (RA) and declinations (Dec) of objects. Note that an ASCII-format catalogue can be passed via pipe, but in general a FITS-format catalogue cannot.

### `FILE_TYPE` (`-f` / `--file-type`)

An integer indicating the format of the input catalogue. Allowed values are:
-   `0`: for ASCII text file;
-   `1`: for FITS table.

### `ASCII_COMMENT` (`--comment`)

Comment symbol for the input catalogue (lines starting with this symbol are omitted), e.g.

```nginx
ASCII_COMMENT = '#'
```

### `COORD_COLUMN` (`-C` / `--coord-col`)

Column numbers (ASCII format, starting from 1) or names (FITS format, case insensitive) for the RA and Dec of each object in the input catalogue, e.g.

```nginx
COORD_COLUMN = [1,2]
    # ASCII: first column for RA, second column for Dec
COORD_COLUMN = [RA,DEC]
    # FITS: "RA" for RA, "Dec" for Dec
```

### `OUTPUT` (`-o` / `--output`)

Filename of the output catalogue. It is safe to specify the same filename as `INPUT`, though this is not recommended. It can be a pipe for the ASCII-format catalogue, but cannot for FITS format.

### `OUTPUT_COLUMN` (`-e` / `--output-col`)

Optional parameter for column numbers (ASCII format, starting from 1) or names (FITS format, case insensitive) to be saved to the output catalogue. The columns must be available in the input catalogue. Columns of the output catalogue will be in the same order as is specified here, and the evaluated maskbit and subsample ID (if applicable) columns are always at the end.

If `OUTPUT_COLUMN` is not set, all columns of the input catalogue will be saved in the original order to the output, in addition with the columns for maskbits and subsample IDs.

### `OVERWRITE` (`-O` / `--overwrite`)

An integer value indicating whether to overwrite existing files. Allowed values are

-   `0`: quit the program when an output file exist;
-   positive: force overwriting output files whenever possible;
-   negative: notify at most this number (absolute value) of times, for asking whether overwriting existing files.

### `VERBOSE` (`-v` / `--verbose`)

A boolean value indicating whether to show detailed standard outputs.

