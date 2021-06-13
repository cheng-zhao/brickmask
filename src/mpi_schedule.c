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
  if (rank != BRICKMASK_MPI_ROOT) {
    b->ra1 = b->ra2 = b->dec1 = b->dec2 = NULL;
    b->name = NULL;
    b->subid = NULL;
    b->nmask = NULL;
    b->fmask = NULL;
    b->mlen = NULL;
  }

  /* Broadcast number of bricks and length of brick names. */
  MPI_Request req[3];
  if (MPI_Ibcast(&b->n, 1, MY_MPI_SIZE_T, BRICKMASK_MPI_ROOT,
      MPI_COMM_WORLD, req) || MPI_Ibcast(&b->nlen, 1, MPI_INT,
      BRICKMASK_MPI_ROOT, MPI_COMM_WORLD, req + 1) ||
      MPI_Waitall(2, req, MPI_STATUSES_IGNORE)) {
    P_ERR("failed to broadcast brick information\n");
    MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);
  }

  /* Allocate memory for brick names. */
  if (rank != BRICKMASK_MPI_ROOT) {
    if (!(b->name = malloc(b->n * sizeof(char *)))) {
      P_ERR("failed to allocate memory for task-private brick information\n");
    }
    for (size_t i = 0; i < b->n; i++) b->name[i] = NULL;
    if (!(b->name[0] = malloc(b->n * (b->nlen + 1) * sizeof(char)))) {
      P_ERR("failed to allocate memory for task-private brick information\n");
      MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MEMORY);
    }
  }

  /* Broadcast brick names and number of subsamples. */
  if (MPI_Ibcast(b->name[0], b->n * (b->nlen + 1), MPI_CHAR,
      BRICKMASK_MPI_ROOT, MPI_COMM_WORLD, req) || MPI_Ibcast(&b->nsp, 1,
      MPI_INT, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD, req + 1) ||
      MPI_Waitall(2, req, MPI_STATUSES_IGNORE)) {
    P_ERR("failed to broadcast brick information\n");
    MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);
  }

  /* Define brick names, and allocate memory for the workers. */
  if (rank != BRICKMASK_MPI_ROOT) {
    for (size_t i = 1; i < b->n; i++)
      b->name[i] = b->name[0] + i * (b->nlen + 1);

    if (!(b->subid = malloc(b->nsp * sizeof(int))) ||
        !(b->nmask = malloc(b->nsp * sizeof(size_t))) ||
        !(b->mlen = malloc(b->nsp * sizeof(size_t))) ||
        !(b->fmask = malloc(b->nsp * sizeof(char **)))) {
      P_ERR("failed to allocate memory for task-private brick information\n");
      MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MEMORY);
    }
    for (int i = 0; i < b->nsp; i++) b->fmask[i] = NULL;
  }

  /* Broadcast subsample IDs, and number & length of maskbit files. */
  if (MPI_Ibcast(b->subid, b->nsp, MPI_INT, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD,
      req) || MPI_Ibcast(b->nmask, b->nsp, MY_MPI_SIZE_T, BRICKMASK_MPI_ROOT,
      MPI_COMM_WORLD, req + 1) || MPI_Ibcast(b->mlen, b->nsp, MY_MPI_SIZE_T,
      BRICKMASK_MPI_ROOT, MPI_COMM_WORLD, req + 2) ||
      MPI_Waitall(3, req, MPI_STATUSES_IGNORE)) {
    P_ERR("failed to broadcast brick information\n");
    MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);
  }

  /* Allocate memory for maskbit filenames. */
  if (rank != BRICKMASK_MPI_ROOT) {
    for (int i = 0; i < b->nsp; i++) {
      if (b->nmask[i] &&
          !(b->fmask[i] = malloc(b->nmask[i] * sizeof(char *)))) {
        P_ERR("failed to allocate memory for task-private brick information\n");
        MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MEMORY);
      }
      for (size_t j = 0; j < b->nmask[i]; j++) b->fmask[i][j] = NULL;
      if (b->nmask[i] &&
          !(b->fmask[i][0] = malloc(b->mlen[i] * sizeof(char)))) {
        P_ERR("failed to allocate memory for task-private brick information\n");
        MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MEMORY);
      }
    }
  }

  /* Broadcast maskfit filenames. */
  for (int i = 0; i < b->nsp; i++) {
    if (MPI_Bcast(b->fmask[i][0], b->mlen[i], MPI_CHAR, BRICKMASK_MPI_ROOT,
        MPI_COMM_WORLD)) {
      P_ERR("failed to broadcase brick information\n");
      MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);
    }
    /* Search for individual filenames. */
    char *end = b->fmask[i][0] + b->mlen[i];
    for (size_t j = 1; j < b->nmask[i]; j++) {
      char *p = b->fmask[i][j - 1];
      b->fmask[i][j] = memchr(p, '\0', end - p) + 1;
    }
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

  bool subid = false;
  if (rank == BRICKMASK_MPI_ROOT && (*data)->subid) subid = true;
  if (MPI_Bcast(&subid, 1, MPI_C_BOOL, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD)) {
    P_ERR("failed to share data information\n");
    MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);
  }

  int *nsend, *disp;
  if (!(nsend = calloc(size, sizeof(int))) ||
      !(disp = calloc(size, sizeof(int)))) {
    P_ERR("failed to allocate memory for sharing data information\n");
    MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MEMORY);
  }

  DATA *d;
  /* Split the data by brick IDs. */
  if (rank == BRICKMASK_MPI_ROOT) {
    d = *data;
    long num = d->nbrick / size;        /* number of bricks per worker */
    long num0 = d->nbrick - num * (size - 1);   /* number of bricks for root */

    /* Count the number of data associated with each brick. */
    long prev = -1;
    long cnt = -1;
    size_t i_prev = 0;
    for (size_t i = 0; i < d->n; i++) {
      if (d->id[i] != prev) {
        prev = d->id[i];
        if (++cnt == num0) {
          i_prev = i;
          break;
        }
      }
    }
    if (!i_prev) i_prev = d->n;
    int wid = 0;
    nsend[wid++] = i_prev;      /* number of data to be processed by root */
    for (size_t i = i_prev; i < d->n; i++) {
      if (d->id[i] != prev) {
        prev = d->id[i];
        if (++cnt == num0 + num * wid) {
          /* Compute the length of data for this worker. */
          nsend[wid] = i - i_prev;
          disp[wid] = disp[wid - 1] + nsend[wid - 1];
          if (wid++ == size) {
            P_ERR("unexpected number of MPI tasks\n");
            MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_UNKNOWN);
          }
          i_prev = i;
        }
      }
    }
    if (wid != size) {
      nsend[wid] = d->n - i_prev;
      disp[wid] = disp[wid - 1] + nsend[wid - 1];
    }
  }
  /* Allocate memory for the task-private data. */
  else {
    if (!(*data = calloc(1, sizeof(DATA)))) {
      P_ERR("failed to allocate memory for task-private brick information\n");
      MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MEMORY);
    }
    d = *data;
    d->ra = d->dec = NULL;
    d->idx = NULL;
    d->id = NULL;
    d->mask = NULL;
    d->subid = NULL;
    d->content = NULL;
  }

  /* Broadcast the length of data and number of bricks for each task. */
  MPI_Request req[3];
  if (MPI_Ibcast(nsend, size, MPI_INT, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD,
      req) || MPI_Ibcast(disp, size, MPI_INT, BRICKMASK_MPI_ROOT,
      MPI_COMM_WORLD, req + 1) || MPI_Ibcast(&d->nbrick, 1, MY_MPI_SIZE_T,
      BRICKMASK_MPI_ROOT, MPI_COMM_WORLD, req + 2) ||
      MPI_Waitall(3, req, MPI_STATUSES_IGNORE)) {
    P_ERR("failed to share data information\n");
    MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);
  }

  /* Set length of data and number of bricks for each task. */
  d->n = nsend[rank];
  size_t nbrick = d->nbrick / size;
  if (rank == BRICKMASK_MPI_ROOT) d->nbrick -= nbrick * (size - 1);
  else d->nbrick = nbrick;

  /* Allocate memory for the data. */
  if (rank != BRICKMASK_MPI_ROOT && d->n) {
    if (!(d->ra = malloc(d->n * sizeof(double))) ||
        !(d->dec = malloc(d->n * sizeof(double))) ||
        !(d->id = malloc(d->n * sizeof(long))) ||
        !(d->mask = calloc(d->n, sizeof(uint64_t))) ||
        (subid && !(d->subid = calloc(d->n, sizeof(unsigned char))))) {
      P_ERR("failed to allocate memory for the task-private data\n");
      MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MEMORY);
    }
  }

  /* Scatter the data. */
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
  MPI_Request req[4];
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
    if (MPI_Igatherv(MPI_IN_PLACE, n, MPI_DOUBLE, data->ra, nrecv, disp,
        MPI_DOUBLE, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD, req) ||
        MPI_Igatherv(MPI_IN_PLACE, n, MPI_DOUBLE, data->dec, nrecv, disp,
        MPI_DOUBLE, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD, req + 1) ||
        MPI_Igatherv(MPI_IN_PLACE, n, MPI_UINT64_T, data->mask, nrecv, disp,
        MPI_DOUBLE, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD, req + 2) ||
        (data->subid && MPI_Igatherv(MPI_IN_PLACE, n, MPI_UNSIGNED_CHAR,
        data->subid, nrecv, disp, MPI_UNSIGNED_CHAR, BRICKMASK_MPI_ROOT,
        MPI_COMM_WORLD, req + 3))) {
      P_ERR("failed to gather data information\n");
      MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);
    }

    int nreq = (data->subid) ? 4 : 3;
    if (MPI_Waitall(nreq, req, MPI_STATUSES_IGNORE)) {
      P_ERR("failed to gather data information\n");
      MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);
    }

    /* Get the total length of data. */
    data->n = 0;
    for (int i = 0; i < size; i++) data->n += nrecv[i];
    free(nrecv);
    free(disp);

    /* Get the largest mask width. */
    if (MPI_Reduce(MPI_IN_PLACE, &data->mtype, 1, MPI_INT, MPI_MAX,
        BRICKMASK_MPI_ROOT, MPI_COMM_WORLD)) {
      P_ERR("failed to gather data information\n");
      MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);
    }
  }
  else {
    /* Gather data. */
    if (MPI_Igatherv(data->ra, n, MPI_DOUBLE, NULL, NULL, NULL,
        MPI_DOUBLE, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD, req) ||
        MPI_Igatherv(data->dec, n, MPI_DOUBLE, NULL, NULL, NULL,
        MPI_DOUBLE, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD, req + 1) ||
        MPI_Igatherv(data->mask, n, MPI_UINT64_T, NULL, NULL, NULL,
        MPI_DOUBLE, BRICKMASK_MPI_ROOT, MPI_COMM_WORLD, req + 2) ||
        (data->subid && MPI_Igatherv(data->subid, n, MPI_UNSIGNED_CHAR,
        NULL, NULL, NULL, MPI_UNSIGNED_CHAR, BRICKMASK_MPI_ROOT,
        MPI_COMM_WORLD, req + 3))) {
      P_ERR("failed to gather data information\n");
      MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);
    }

    int nreq = (data->subid) ? 4 : 3;
    if (MPI_Waitall(nreq, req, MPI_STATUSES_IGNORE)) {
      P_ERR("failed to gather data information\n");
      MPI_Abort(MPI_COMM_WORLD, BRICKMASK_ERR_MPI);
    }

    /* Send the largest mask width. */
    if (MPI_Reduce(&data->mtype, NULL, 1, MPI_INT, MPI_MAX, BRICKMASK_MPI_ROOT,
        MPI_COMM_WORLD)) {
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
