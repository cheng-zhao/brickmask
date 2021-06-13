/*******************************************************************************
* brickmask.c: this file is part of the brickmask program.

* brickmask: assign bit codes defined on Legacy Survey brick pixels
             to a catalogue with sky coordinates.

* Github repository:
        https://github.com/cheng-zhao/brickmask

* Copyright (c) 2020 -- 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#include "define.h"
#include "load_conf.h"
#include "get_brick.h"
#include "data_io.h"
#include "assign_mask.h"
#include "save_file.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef MPI
#include "mpi_schedule.h"
#include <mpi.h>
#define BRICKMASK_QUIT(status)  MPI_Abort(MPI_COMM_WORLD, status); exit(status);
#else
#define BRICKMASK_QUIT(status)  return status;
#endif

int main(int argc, char *argv[]) {
#ifdef MPI
  if (MPI_Init(NULL, NULL)) {
    P_EXT("failed to initialise MPI\n");
    return BRICKMASK_ERR_MPI;
  }
  int rank;
  if (MPI_Comm_rank(MPI_COMM_WORLD, &rank)) {
    P_ERR("failed to obtain MPI ranks\n");
    BRICKMASK_QUIT(BRICKMASK_ERR_MPI);
  }
#endif

  bool verbose = false;
  CONF *conf = NULL;
  BRICK *brick = NULL;
  DATA *data = NULL;

#ifdef MPI
  if (rank == BRICKMASK_MPI_ROOT) {
#endif
    if (!(conf = load_conf(argc, argv))) {
      printf(FMT_FAIL);
      P_EXT("failed to load configuration parameters\n");
      BRICKMASK_QUIT(BRICKMASK_ERR_CFG);
    }
    else verbose = conf->verbose;

    if (!(brick = get_brick(conf))) {
      printf(FMT_FAIL);
      P_EXT("failed to get information of the bricks\n");
      conf_destroy(conf);
      BRICKMASK_QUIT(BRICKMASK_ERR_BRICK);
    }

    if (!(data = read_data(conf))) {
      printf(FMT_FAIL);
      P_EXT("failed to read the input data catalog\n");
      conf_destroy(conf); brick_destroy(brick);
      BRICKMASK_QUIT(BRICKMASK_ERR_FILE);
    }

    if (sort_data(brick, data, verbose)) {
      printf(FMT_FAIL);
      P_EXT("failed to sort the input data\n");
      conf_destroy(conf); brick_destroy(brick); data_destroy(data);
      BRICKMASK_QUIT(BRICKMASK_ERR_BRICK);
    }
#ifdef MPI
  }

  /* Broadcast verbose. */
  if (MPI_Bcast(&verbose, 1, MPI_C_BOOL, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD)) {
    P_ERR("failed to communicate between MPI tasks\n");
    BRICKMASK_QUIT(BRICKMASK_ERR_MPI);
  }
  mpi_init_worker(&brick, &data, verbose);
#endif

  if (assign_mask(brick, data, verbose)) {
    printf(FMT_FAIL);
    P_EXT("failed to assign maskbits to the data\n");
    conf_destroy(conf); brick_destroy(brick); data_destroy(data);
    BRICKMASK_QUIT(BRICKMASK_ERR_MASK);
  }

#ifdef MPI
  if (MPI_Barrier(MPI_COMM_WORLD)) {
    P_ERR("failed to set MPI barrier\n");
    BRICKMASK_QUIT(BRICKMASK_ERR_MPI);
  }
  mpi_gather_data(brick, data);

  if (rank == BRICKMASK_MPI_ROOT) {
#endif
    brick_destroy(brick);

    if (save_data(conf, data)) {
      printf(FMT_FAIL);
      P_EXT("failed to save the output data catalog\n");
      conf_destroy(conf); data_destroy(data);
      BRICKMASK_QUIT(BRICKMASK_ERR_SAVE);
    }

    conf_destroy(conf);
#ifdef MPI
  }

  if (MPI_Finalize()) {
    P_ERR("failed to finalize MPI tasks\n");
    BRICKMASK_QUIT(BRICKMASK_ERR_MPI);
  }
#endif
  return 0;
}
