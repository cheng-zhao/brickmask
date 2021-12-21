# brickmask

![GitHub](https://img.shields.io/github/license/cheng-zhao/brickmask.svg)
![Codacy grade](https://img.shields.io/codacy/grade/b780618f6c2144649f71de9814a36430.svg)

## Table of Contents

-   [Introduction](#introduction)
-   [Compilation](#compilation)
-   [Running](#running)
-   [Configurations](#configurations)
-   [Acknowledgements](#acknowledgements)
-   [References](#references)

## Introduction

brickmask is a tool for assigning bit codes defined on [Legacy Survey](http://legacysurvey.org) brick pixels<sup id="quote0">[1](#footnote1),[2](#footnote2)</sup> to a catalogue with sky coordinates: J2000 right ascension (RA) and declination (Dec) in degrees.

A common usage of brickmask is to mark objects to be removed based on the veto masks defined on brick pixels. The input catalogue can be either a plain ASCII file, or in the [FITS format](https://fits.gsfc.nasa.gov/fits_home.html). Besides, a collection of FITS-format maskbits files has to be provided. The input objects are then located in the bricks and assigned the corresponding bit codes. Optionally, an extra attribute specifying subsample IDs of the maskbits can be appended to the input catalogue as well.

This program is compliant with the ISO C99 and IEEE POSIX.1-2008 standards, and supports Message Passing Interface (MPI) parallelisation. It is written by Cheng Zhao (&#36213;&#25104;), and is distributed under the [MIT license](LICENSE.txt).

This tool was originally developed for the [*Extended Baryon Oscillation Spectroscopic Survey*](https://www.sdss.org/surveys/eboss) (eBOSS) Emission-line Galaxy (ELG) masks (see [\[1\]](#ref1),[\[2\]](#ref2)). It is also compatible with the veto masks used for the [*Dark Energy Spectroscopic Instrument*](https://www.desi.lbl.gov/) (DESI). If you use this tool in research work that results in publications, please cite the following paper:

> Zhao et al., 2021, [The completed SDSS-IV extended Baryon Oscillation Spectroscopic Survey: 1000 multi-tracer mock catalogues with redshift evolution and systematics for galaxies and quasars of the final data release](https://doi.org/10.1093/mnras/stab510), *MNRAS*, 503(1):1149&ndash;1173 ([arXiv:2007.08997](https://arxiv.org/abs/2007.08997))


<sub><span id="footnote1">1.</span> see e.g. [http://legacysurvey.org/dr9/files/#survey-bricks-fits-gz](http://legacysurvey.org/dr9/files/#survey-bricks-fits-gz) [&#8617;](#quote0)</sub><br />
<sub><span id="footnote2">2.</span> e.g. [https://www.legacysurvey.org/dr9/bitmasks](https://www.legacysurvey.org/dr9/bitmasks) [&#8617;](#quote0)</sub>

<sub>[\[TOC\]](#table-of-contents)</sub>

## Compilation

The only mandatory external library required by this program is [CFITSIO](https://heasarc.gsfc.nasa.gov/fitsio), for reading/writing FITS files. If CFITSIO is already installed in the standard library paths, such as `/usr/lib`, one should be able to compile brickmask simply with the command
```bash
make
```

Otherwise, the CFITSIO installation path has to be set via `CFITSIO_DIR` in [Makefile](Makefile#L12). In this case, the header file `fitsio.h` has to be available in `CFITSIO_DIR/include`, and the library file (`libcfitsio.so`) must be present in `CFITSIO_DIR/lib`.

To enable MPI support, a compiler wrapper for MPI programs (such as `mpicc`) must be available, and the option `USE_MPI` in [Makefile](Makefile#L7) should be set to `T`.

The other optional compilation flags (can be set via `CFLAGS`) are summarised below:

| Compilation Flag  | Usage                                                                        |
|:-----------------:|------------------------------------------------------------------------------|
| `-DEBOSS`         | for eBOSS ELG masks<sup id="quote1">[3](#footnote3)</sup>                    |
| `-DFAST_FITS_IMG` | enable low-level maskbits file reading<sup id="quote2">[4](#footnote4)</sup> |

<sub><span id="footnote3">3.</span> See [https://data.sdss.org/datamodel/files/EBOSS_LSS/catalogs/DR16/ELGmask/mask.html](https://data.sdss.org/datamodel/files/EBOSS_LSS/catalogs/DR16/ELGmask/mask.html). Note also that there are additional eBOSS ELG masks that should be set using the script [eBOSS_ELG_extra.py](scripts/eBOSS_ELG_extra.py). [&#8617;](#quote1)</sub><br />
<sub><span id="footnote4">4.</span> The low-level FITS image reader is &sim; 4 times faster than the default reader for plain images, but only marginally faster for gzipped images. Note that it should never be enabled for maskbits compressed with algorithms other than gzip (such as `.fits.fz` files). [&#8617;](#quote2)</sub>

<sub>[\[TOC\]](#table-of-contents)</sub>

## Running

If the CFITSIO library is installed in a custom directory, the path to the library file (`libcfitsio.so`) has to be included in the environment variable `LD_LIBRARY_PATH` before running this program, e.g.
```bash
export LD_LIBRARY_PATH=/Custom_CFITSIO_DIR/lib:$LD_LIBRARY_PATH
```

Moreover, if MPI is enabled, the code should be run with an MPI executable, such as `mpirun` or `srun`.

Once the program is run successfully, it looks for the configuration file set via the `-c` command line option, or `brickmask.conf` by default (see [CONFIG.md](CONFIG.md) for details). Configuration parameters can also be set via command line options, which override the entries in the configuration file. A list of all command line options can be found with the `-h`/`--help` option.

<sub>[\[TOC\]](#table-of-contents)</sub>

## Configurations

Detailed descriptions of all configuration parameters can be found in [CONFIG.md](CONFIG.md).

Apart from parameters that can be supplied via command line options and the configuration file, there are also some runtime settings defined as macros in [`define.h`](src/define.h). For instance, [`BRICKMASK_FITS_SUBID`](src/define.h#L93) indicates the column name of subsample IDs for the FITS-format output catalogue.

<sub>[\[TOC\]](#table-of-contents)</sub>

## Acknowledgements

I thank Dr. Anand Raichoor for helpful discussions on the development and debugging of this tool, and Dr. Arnaud de Mattia for sharing the original script for assigning extra maskbit codes for eBOSS ELGs.

<sub>[\[TOC\]](#table-of-contents)</sub>

## References

<span id="ref1">\[1\]</span> Raichoor et al., 2021, [The completed SDSS-IV extended Baryon Oscillation Spectroscopic Survey: large-scale structure catalogues and measurement of the isotropic BAO between redshift 0.6 and 1.1 for the Emission Line Galaxy Sample](https://doi.org/10.1093/mnras/staa3336), *MNRAS*, 500(3):3254&ndash;3274 ([arXiv:2007.09007](https://arxiv.org/abs/2007.09007))

<span id="ref2">\[2\]</span> Zhao et al., 2021, [The completed SDSS-IV extended Baryon Oscillation Spectroscopic Survey: 1000 multi-tracer mock catalogues with redshift evolution and systematics for galaxies and quasars of the final data release](https://doi.org/10.1093/mnras/stab510), *MNRAS*, 503(1):1149&ndash;1173 ([arXiv:2007.08997](https://arxiv.org/abs/2007.08997))

<sub>[\[TOC\]](#table-of-contents)</sub>

