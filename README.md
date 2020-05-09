# brickmask

![GitHub](https://img.shields.io/github/license/cheng-zhao/brickmask.svg)

## Table of Contents

- [Introduction](#introduction)
- [Dependencies](#dependencies)
- [Compilation](#compilation)
- [Running](#running)
- [Configurations](#configurations)
- [Todo list](#todo-list)
- [Acknowledgement](#acknowledgement)

## Introduction

brickmask is a tool for assigning bit codes defined on [Legacy Survey](http://legacysurvey.org) brick pixels<sup id="quote0">[1](#footnote1),[2](#footnote2)</sup> to a catalogue with sky coordinates: J2000 right ascension (RA) and declination (Dec) in degrees.

A common usage of brickmask is to mark objects to be removed based on the veto masks defined on brick pixels. The input catalogue can be either a plain ASCII file, or in the [fits format](https://fits.gsfc.nasa.gov/fits_home.html). And a collection of fits-format brick mask files has to be provided. The input objects are then located in the bricks and assigned the corresponding bit codes as an additional attribute. Moreover, if the bricks are separated as multiple subsamples, another column indicating the subsample ID is also appended to the input catalogue.

brickmask is compliant with the ISO C99 and IEEE POSIX.1-2008 standards. It is written by Cheng Zhao (&#36213;&#25104;), and is distributed under the MIT license (see [LICENSE.txt](LICENSE.txt) for the details).

If you use this tool in research work that results in publications, please cite the following paper:

> Zhao et al. in preparation.

<small>[\[TOC\]](#table-of-contents)</small>

## Dependencies

- [cfitsio](https://heasarc.gsfc.nasa.gov/fitsio): for fits file I/O.
- GLIBC `getopt.h`: for parsing command line options

<small>[\[TOC\]](#table-of-contents)</small>

## Compilation

If the cfitsio library is already installed in the standard library paths, such as `/usr/lib`, one should be able to compile brickmask simply with the command
```bash
make
```

Otherwise, the paths to the cfitsio header (`fitsio.h`) and library files (`libcfitsio.so`/`libcfitsio.a`) have to be supplemented to [Makefile](Makefile), e.g.
```makefile
LIBS = -lm -lcfitsio -L/custom_cfitsio_path/lib
INCL = -I/custom_cfitsio_path/include
```

Optionally, brickmask can be compiled with OpenMP, for reading the brick mask files in parallel with multiple threads. This feature can be enabled by uncommenting the `-DOMP -fopenmp` options in [Makefile](Makefile#L4). However, the bottleneck of the code is typically disk I/O. Thus, the overhead cannot be reduced with multiple OpenMP threads. That is why the OpenMP options are disabled by default.

<small>[\[TOC\]](#table-of-contents)</small>

## Running

If cfitsio is installed in a custom directory, the path to the library file has to be included in the environment variable `LD_LIBRARY_PATH` before running brickmask, e.g.
```bash
export LD_LIBRARY_PATH=/custom_cfitsio_path/lib:$LD_LIBRARY_PATH
```

Once brickmask is run successfully, it looks for the configuration file set via the `-c` command line option, or `brickmask.conf` by default, which sets mostly the input/output options (see [Configurations](#configurations) for details). If the options are unset or invalid, the code prints some error messages and aborts.

The configuration options can also be set via command line options, which override the entries in the configuration file. A list of all command line options can be found with the `-h`/`--help` option.

<small>[\[TOC\]](#table-of-contents)</small>

## Configurations

The entries of the configuration file is summarised as follows

| Parameter    | Description                                                                                                                                                                 | Default |
|--------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|---------|
| `BRICK_LIST` | FITS table<sup id="quote1">[1](#footnote1)</sup> with the RA, Dec bounds of all bricks                                                                                                |         |
| `BRICK_DIR`  | Directory for FITS table<sup id="quote2">[2](#footnote2)</sup> with bit codes for each brick                                                                                          |         |
| `INPUT`      | Input catalogue                                                                                                                                                             |         |
| `FILE_TYPE`  | Format of the input catalogue: <br />`0`: ASCII format, with RA and Dec being the first two columns<br />`1`: fits format, with RA and Dec being the `RA` and `DEC` columns | 0       |
| `COMMENT`    | Comment symbol for the input catalogue<br />(lines starting with this symbol are omitted)                                                                                   | '#'     |
| `OUTPUT`     | Output file<br />(in the same format as `INPUT`)                                                                                                                            |         |
| `FORCE`      | Non-zero: overwrite existing files without asking<br />`0`: ask before overwriting files                                                                                    | 0       |
| `VERBOSE`    | Non-zero: detailed standard outputs<br />`0`: concise standard outputs                                                                                                      | 1       |


Apart from the input/output options, there are also some runtime settings of the brickmask tool, defined as macros in [define.h](define.h). For instance, [`MSK_COL`](define.h#L53) and [`CHK_COL`](define.h#L54) indicate the column names for the bit codes and subsample IDs appended to the fits-format output file. And for large input files with not too many columns, one may want to adjust the [`LEN_IN_LINE`](define.h#L65) value, which specifies the maximum number of characters reserved for each line of the input catalogue, to reduce the memory consumption of the brickmask tool.

<small>[\[TOC\]](#table-of-contents)</small>

## Acknowledgement

I thank Dr. Anand Raichoor for helpful discussions on the development and debugging of this tool.

<small>[\[TOC\]](#table-of-contents)</small>

## Todo list

- Replace the GLIBC `getopt.h` reliance with the standalone [libcfg](https://github.com/cheng-zhao/libcfg) library.
- Smarter ASCII file reader.

<small>[\[TOC\]](#table-of-contents)</small>

<span id="footnote1">1</span>. see e.g. [http://legacysurvey.org/dr8/files/#survey-bricks-fits-gz](http://legacysurvey.org/dr8/files/#survey-bricks-fits-gz) [&#8617;](#quote0) [&#8617;](#quote1)

<span id="footnote2">2</span>. e.g. [http://legacysurvey.org/dr8/files/#region-survey-bricks-dr8-region-fits-gz](http://legacysurvey.org/dr8/files/#region-survey-bricks-dr8-region-fits-gz) [&#8617;](#quote0) [&#8617;](#quote2)
