/*******************************************************************************
* mpi_schedule.c: this file is part of the brickmask program.

* brickmask: assign bit codes defined on Legacy Survey brick pixels
             to a catalogue with sky coordinates.

* Github repository:
        https://github.com/cheng-zhao/brickmask

* Copyright (c) 2020 -- 2021 Cheng Zhao <zhaocheng03@gmail.com>  [MIT license]

*******************************************************************************/

#ifdef MPI
#include "define.h"
#include "mpi_schedule.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <mpi.h>

/* Define the MPI data type for size_t */
#if SIZE_MAX == UINT_MAX
  #define MY_MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
  #define MY_MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
  #define MY_MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
  #error "unsupported size of type \"size_t\""
#endif

#if BRICKMASK_MPI_ROOT != 0
  #error "MPI root rank must be 0"
#endif

/*============================================================================*\
                 Functions for sharing information with workers
\*============================================================================*/

/******************************************************************************
Function `mpi_bcast_brick`:
  Broadcast information of bricks.
Arguments:
  * `brick`:    structure for storing information of bricks.
******************************************************************************/
static void mpi_bcast_brick(BRICK **brick) {
  int rank;
  if (MPI_Comm_rank(MPI_COMM_WORLD, &rank))
    MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);

  if (rank != BRICKMASK_MPI_ROOT && !(*brick = calloc(1, sizeof(BRICK)))) {
    P_ERR("failed to allocate memory for task-private brick information\n");
    MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MEMORY);
  }
  BRICK *b = *brick;
  size_t len;           /* length of brick path */
  if (rank != BRICKMASK_MPI_ROOT) {
    b->ra1 = b->ra2 = b->dec1 = b->dec2 = NULL;
    b->name = NULL;
    b->photsys = NULL;
  }
  else {
    len = strlen(b->bpath) + 1;
  }

  /* Setup custom MPI datatype for broadcasting. */
  MPI_Datatype dtypes[5] =
      {MY_MPI_SIZE_T, MY_MPI_SIZE_T, MPI_UINT64_T, MPI_UINT64_T, MPI_INT};
  int lengths[5] = {1, 1, 1, 1, 1};
  MPI_Aint addrs[5];
  MPI_Datatype send_type;

  /* Broadcast number of bricks and length of brick names. */
  if (MPI_Get_address(&len, addrs) ||
      MPI_Get_address(&(b->n), addrs + 1) ||
      MPI_Get_address(&(b->mnull), addrs + 2) ||
      MPI_Get_address(&(b->enull), addrs + 3) ||
      MPI_Get_address(&(b->nlen), addrs + 4) ||
      MPI_Type_create_struct(5, lengths, addrs, dtypes, &send_type) ||
      MPI_Type_commit(&send_type) ||
      MPI_Bcast(MPI_BOTTOM, 1, send_type, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD) ||
      MPI_Type_free(&send_type)) {
    P_ERR("failed to broadcast brick information\n");
    MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);
  }

  /* Allocate memory for photsys, bpath, and brick names. */
  if (rank != BRICKMASK_MPI_ROOT) {
    if (!(b->photsys = malloc(b->n * sizeof(unsigned char))) ||
        !(b->bpath = calloc(len, sizeof(char))) ||
        !(b->name = malloc(b->n * sizeof(char *)))) {
      P_ERR("failed to allocate memory for task-private brick information\n");
    }
    for (size_t i = 0; i < b->n; i++) b->name[i] = NULL;
    if (!(b->name[0] = malloc(b->n * (b->nlen + 1) * sizeof(char)))) {
      P_ERR("failed to allocate memory for task-private brick information\n");
      MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MEMORY);
    }
  }

  /* Broadcast photsys and brick names. */
  dtypes[0] = MPI_UNSIGNED_CHAR;
  dtypes[1] = dtypes[2] = MPI_CHAR;
  lengths[0] = b->n;
  lengths[1] = len;
  lengths[2] = b->n * (b->nlen + 1);
  if (MPI_Get_address(b->photsys, addrs) ||
      MPI_Get_address(b->bpath, addrs + 1) ||
      MPI_Get_address(b->name[0], addrs + 2) ||
      MPI_Type_create_struct(3, lengths, addrs, dtypes, &send_type) ||
      MPI_Type_commit(&send_type) ||
      MPI_Bcast(MPI_BOTTOM, 1, send_type, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD) ||
      MPI_Type_free(&send_type)) {
    P_ERR("failed to broadcast brick information\n");
    MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);
  }

  /* Define brick names, and allocate memory for the workers. */
  if (rank != BRICKMASK_MPI_ROOT) {
    for (size_t i = 1; i < b->n; i++)
      b->name[i] = b->name[0] + i * (b->nlen + 1);
  }
}

/******************************************************************************
Function `mpi_scatter_data`:
  Scatter parts of the data to the workers.
Arguments:
  * `data`:     structure for storing the input data.
******************************************************************************/
static void mpi_scatter_data(DATA **data) {
  int size, rank;
  size = rank = 0;
  if (MPI_Comm_size(MPI_COMM_WORLD, &size) ||
      MPI_Comm_rank(MPI_COMM_WORLD, &rank))
    MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);

  int *nsend, *disp;    /* Scatterv only accepts ints. */
  if (!(nsend = malloc(size * sizeof(int))) ||
      !(disp = malloc((size + 1) * sizeof(int)))) {
    P_ERR("failed to allocate memory for sharing data information\n");
    MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MEMORY);
  }

  DATA *d;
  /* Split the data by brick IDs. */
  if (rank == BRICKMASK_MPI_ROOT) {
    d = *data;
    const size_t pnum = d->nbrick / size;   /* number of bricks per worker */
    const int rem = d->nbrick % size;
    /* Go to the next worker if the desired number of bricks are assigned. */
    size_t goal = (0 < rem) ? pnum + 1 : pnum;
    int wid = 0;
    size_t cnt = 0;
    disp[0] = 0;
    long prev = d->id[0];

    for (size_t i = 1; i < d->n; i++) {
      if (d->id[i] != prev) {
        prev = d->id[i];
        if (++cnt == goal) {
          if (wid++ == size) {
            P_ERR("unexpected number of MPI tasks\n");
            MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_UNKNOWN);
          }
          disp[wid] = i;
          goal = (wid < rem) ? (pnum + 1) * (wid + 1) : pnum * (wid + 1) + rem;
        }
      }
    }
    for (int i = wid + 1; i <= size; i++) disp[i] = d->n;
    for (int i = 0; i < size; i++) nsend[i] = disp[i + 1] - disp[i];
  }
  /* Allocate memory for the task-private data. */
  else {
    if (!(*data = calloc(1, sizeof(DATA)))) {
      P_ERR("failed to allocate memory for task-private brick information\n");
      MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MEMORY);
    }
    d = *data;
    d->ra = d->dec = NULL;
    d->idx = d->cidx = d->iidx = NULL;
    d->id = NULL;
    d->mask = d->nexp[0] = d->nexp[1] = d->nexp[2] = NULL;
    d->content = NULL;
  }

  /* Broadcast the length of data and number of bricks for each task. */
  MPI_Datatype dtypes[3], send_type;
  int lengths[3];
  MPI_Aint addrs[3];

  dtypes[0] = dtypes[1] = MPI_INT;
  dtypes[2] = MY_MPI_SIZE_T;
  lengths[0] = lengths[1] = size;
  lengths[2] = 1;
  if (MPI_Get_address(nsend, addrs) ||
      MPI_Get_address(disp, addrs + 1) ||
      MPI_Get_address(&(d->nbrick), addrs + 2) ||
      MPI_Type_create_struct(3, lengths, addrs, dtypes, &send_type) ||
      MPI_Type_commit(&send_type) ||
      MPI_Bcast(MPI_BOTTOM, 1, send_type, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD) ||
      MPI_Type_free(&send_type)) {
    P_ERR("failed to share data information\n");
    MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);
  }

  /* Set length of data and number of bricks for each task. */
  d->n = nsend[rank];
  const size_t pnum = d->nbrick / size;
  const int rem = d->nbrick % size;
  d->nbrick = (rank < rem) ? pnum + 1 : pnum;

  /* Allocate memory for the data. */
  if (rank != BRICKMASK_MPI_ROOT && d->n) {
    if (!(d->ra = malloc(d->n * sizeof(double))) ||
        !(d->dec = malloc(d->n * sizeof(double))) ||
        !(d->id = malloc(d->n * sizeof(long))) ||
        !(d->mask = calloc(d->n, sizeof(uint64_t))) ||
        !(d->nexp[0] = calloc(d->n, sizeof(uint64_t))) ||
        !(d->nexp[1] = calloc(d->n, sizeof(uint64_t))) ||
        !(d->nexp[2] = calloc(d->n, sizeof(uint64_t)))) {
      P_ERR("failed to allocate memory for the task-private data\n");
      MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MEMORY);
    }
  }

  /* Scatter the data. */
  MPI_Request req[3];
  if (rank == BRICKMASK_MPI_ROOT) {     /* in-place scatter from the root */
    if (MPI_Iscatterv(d->ra, nsend, disp, MPI_DOUBLE, MPI_IN_PLACE, d->n,
        MPI_DOUBLE, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD, req) ||
        MPI_Iscatterv(d->dec, nsend, disp, MPI_DOUBLE, MPI_IN_PLACE, d->n,
        MPI_DOUBLE, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD, req + 1) ||
        MPI_Iscatterv(d->id, nsend, disp, MPI_LONG, MPI_IN_PLACE, d->n,
        MPI_LONG, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD, req + 2)) {
      P_ERR("failed to share data information\n");
      MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);
    }
  }
  else {        /* receive only for the workers */
    if (MPI_Iscatterv(NULL, nsend, disp, MPI_DOUBLE, d->ra, d->n,
        MPI_DOUBLE, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD, req) ||
        MPI_Iscatterv(NULL, nsend, disp, MPI_DOUBLE, d->dec, d->n,
        MPI_DOUBLE, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD, req + 1) ||
        MPI_Iscatterv(NULL, nsend, disp, MPI_LONG, d->id, d->n,
        MPI_LONG, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD, req + 2)) {
      P_ERR("failed to share data information\n");
      MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);
    }
  }

  if (MPI_Waitall(3, req, MPI_STATUSES_IGNORE)) {
    P_ERR("failed to share data information\n");
    MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);
  }

  free(nsend);
  free(disp);
}


/*============================================================================*\
                         Interfaces for MPI schedulers
\*============================================================================*/

/******************************************************************************
Function `mpi_init_worker`:
  Initialise MPI workers with brick lists and the input data.
Arguments:
  * `brick`:    structure for bricks;
  * `data`:     structure for the data catalogue;
  * `verbose`:  indicate whether to show detailed standard outputs.
Return:
  Zero on success; non-zero on error.
******************************************************************************/
void mpi_init_worker(BRICK **brick, DATA **data, const bool verbose) {
  int rank;
  if (MPI_Comm_rank(MPI_COMM_WORLD, &rank))
    MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);

  if (rank == BRICKMASK_MPI_ROOT) {
    printf("Distributing data to MPI tasks ...");
    if (verbose) printf("\n");
    fflush(stdout);
  }

  if (!brick || !data) {
    P_ERR("the brick information or input data is not initialised\n");
    MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_INIT);
  }

  /* Send brick information to workers. */
  mpi_bcast_brick(brick);

  /* Send data information to workers. */
  mpi_scatter_data(data);

  if (verbose) {
    printf("  Task %d: %zu objects in %zu bricks\n", rank,
        (*data)->n, (*data)->nbrick);
    fflush(stdout);
  }

  if (MPI_Barrier(MPI_COMM_WORLD)) {
    P_ERR("failed to set MPI barrier\n");
    MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);
  }

  if (rank == BRICKMASK_MPI_ROOT) {
    printf(FMT_DONE);
    fflush(stdout);
  }
}

/******************************************************************************
Function `mpi_gather_data`:
  Gather data from different MPI tasks.
Arguments:
  * `brick`:    structure for bricks;
  * `data`:     structure for the data catalogue.
******************************************************************************/
void mpi_gather_data(BRICK *brick, DATA *data) {
  int size, rank;
  size = rank = 0;
  if (MPI_Comm_size(MPI_COMM_WORLD, &size) ||
      MPI_Comm_rank(MPI_COMM_WORLD, &rank))
    MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);

  if (rank == BRICKMASK_MPI_ROOT) {
    printf("Gathering data from MPI tasks ...");
    fflush(stdout);
  }

  if (!brick || !data) {
    P_ERR("the brick information or input data is not initialised\n");
    MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_INIT);
  }

  /* Gather the number of objects assigned to each task. */
  int n = data->n;
  int *nrecv, *disp;
  nrecv = disp = NULL;
  if (rank == BRICKMASK_MPI_ROOT) {
    if (!(nrecv = calloc(size, sizeof(int))) ||
        !(disp = calloc(size, sizeof(int)))) {
      P_ERR("failed to allocate memory for gathering data from tasks\n");
      MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MEMORY);
    }

    /* Gather the length of data segments. */
    if (MPI_Gather(&n, 1, MPI_INT, nrecv, 1, MPI_INT, BRICKMASK_MPI_ROOT,
        MPI_COMM_WORLD)) {
      P_ERR("failed to gather data information\n");
      MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);
    }
  }
  else {
    /* Send the length of data segments. */
    if (MPI_Gather(&n, 1, MPI_INT, NULL, 1, MPI_INT, BRICKMASK_MPI_ROOT,
        MPI_COMM_WORLD)) {
      P_ERR("failed to gather data information\n");
      MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);
    }
  }

  /* Check if there is no effective worker. */
  bool skip = true;
  if (rank == BRICKMASK_MPI_ROOT) {
    for (int i = 1; i < size; i++) {
      if (nrecv[i]) {
        skip = false;
        break;
      }
    }
    if (MPI_Bcast(&skip, 1, MPI_C_BOOL, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD)) {
      P_ERR("failed to gather data information\n");
      MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);
    }
    if (skip) {
      free(nrecv); free(disp);
      printf(FMT_DONE);
      fflush(stdout);
      return;
    }
  }
  else {
    if (MPI_Bcast(&skip, 1, MPI_C_BOOL, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD)) {
      P_ERR("failed to gather data information\n");
      MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);
    }
    if (skip) {
      brick_destroy(brick);
      data_destroy(data);
      return;
    }
  }

  if (rank == BRICKMASK_MPI_ROOT) {
    /* Compute displacements. */
    for (int i = 1; i < size; i++) disp[i] = disp[i - 1] + nrecv[i - 1];

    /* Gather data. */
    if (MPI_Gatherv(MPI_IN_PLACE, n, MPI_DOUBLE, data->ra, nrecv, disp,
        MPI_DOUBLE, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD) ||
        MPI_Gatherv(MPI_IN_PLACE, n, MPI_DOUBLE, data->dec, nrecv, disp,
        MPI_DOUBLE, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD) ||
        MPI_Gatherv(MPI_IN_PLACE, n, MPI_UINT64_T, data->mask, nrecv, disp,
        MPI_UINT64_T, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD) ||
        MPI_Gatherv(MPI_IN_PLACE, n, MPI_UINT64_T, data->nexp[0], nrecv, disp,
        MPI_UINT64_T, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD) ||
        MPI_Gatherv(MPI_IN_PLACE, n, MPI_UINT64_T, data->nexp[1], nrecv, disp,
        MPI_UINT64_T, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD) ||
        MPI_Gatherv(MPI_IN_PLACE, n, MPI_UINT64_T, data->nexp[2], nrecv, disp,
        MPI_UINT64_T, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD)) {
      P_ERR("failed to gather data information\n");
      MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);
    }

    /* Get the total length of data. */
    data->n = 0;
    for (int i = 0; i < size; i++) data->n += nrecv[i];
    free(nrecv);
    free(disp);

    /* Get the largest mask widths. */
    if (MPI_Reduce(MPI_IN_PLACE, &data->mtype, 1, MPI_INT, MPI_MAX,
        BRICKMASK_MPI_ROOT, MPI_COMM_WORLD) ||
        MPI_Reduce(MPI_IN_PLACE, data->etype, 1, MPI_INT, MPI_MAX,
        BRICKMASK_MPI_ROOT, MPI_COMM_WORLD) ||
        MPI_Reduce(MPI_IN_PLACE, data->etype + 1, 1, MPI_INT, MPI_MAX,
        BRICKMASK_MPI_ROOT, MPI_COMM_WORLD) ||
        MPI_Reduce(MPI_IN_PLACE, data->etype + 2, 1, MPI_INT, MPI_MAX,
        BRICKMASK_MPI_ROOT, MPI_COMM_WORLD)) {
      P_ERR("failed to gather data information\n");
      MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);
    }
  }
  else {
    /* Gather data. */
    if (MPI_Gatherv(data->ra, n, MPI_DOUBLE, NULL, NULL, NULL,
        MPI_DOUBLE, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD) ||
        MPI_Gatherv(data->dec, n, MPI_DOUBLE, NULL, NULL, NULL,
        MPI_DOUBLE, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD) ||
        MPI_Gatherv(data->mask, n, MPI_UINT64_T, NULL, NULL, NULL,
        MPI_UINT64_T, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD) ||
        MPI_Gatherv(data->nexp[0], n, MPI_UINT64_T, NULL, NULL, NULL,
        MPI_UINT64_T, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD) ||
        MPI_Gatherv(data->nexp[1], n, MPI_UINT64_T, NULL, NULL, NULL,
        MPI_UINT64_T, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD) ||
        MPI_Gatherv(data->nexp[2], n, MPI_UINT64_T, NULL, NULL, NULL,
        MPI_UINT64_T, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD)) {
      P_ERR("failed to gather data information\n");
      MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);
    }

    /* Send the largest mask width. */
    if (MPI_Reduce(&data->mtype, NULL, 1, MPI_INT, MPI_MAX,
        BRICKMASK_MPI_ROOT, MPI_COMM_WORLD) ||
        MPI_Reduce(data->etype, NULL, 1, MPI_INT, MPI_MAX,
        BRICKMASK_MPI_ROOT, MPI_COMM_WORLD) ||
        MPI_Reduce(data->etype + 1, NULL, 1, MPI_INT, MPI_MAX,
        BRICKMASK_MPI_ROOT, MPI_COMM_WORLD) ||
        MPI_Reduce(data->etype + 2, NULL, 1, MPI_INT, MPI_MAX,
        BRICKMASK_MPI_ROOT, MPI_COMM_WORLD)) {
      P_ERR("failed to gather data information\n");
      MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);
    }
  }

  if (rank == BRICKMASK_MPI_ROOT) {
    printf(FMT_DONE);
    fflush(stdout);
  }
  else {
    brick_destroy(brick);
    data_destroy(data);
  }

  if (MPI_Barrier(MPI_COMM_WORLD)) {
    P_ERR("failed to set MPI barrier\n");
    MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);
  }
}

#endif
