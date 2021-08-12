
  // ******************************************************************************************
  //                                                                                          *
  //                    Copyright (C) 2021 C. H. Patterson and A.-M. Elena                    *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

void mpi_begin_end(int*, int*, int, int, JOB_PARAM*, FILES );

void mpi_begin_end_lower_triangle(int*, int*, int*, int*, int, int, JOB_PARAM*, FILES );

void mpi_receive_offset(int*, int*, int*, int*, int ,JOB_PARAM*, FILES );

void mpi_receive_offset_pairs(int*, int*, int*, int*, PAIR_TRAN*, ATOM*, JOB_PARAM*, FILES);

void mpi_receive_offset_momentum(int*, int*, int*, int*, PAIR_TRAN*, ATOM*, JOB_PARAM*, FILES);

void mpi_receive_offset_nblocks(int*, int*, int*, int*, PAIR_TRAN*, int*, ATOM*, JOB_PARAM*, FILES);

void mpi_local_begin_end(int, int[], int[], int[], int[], JOB_PARAM*, FILES);

void CreateCounter(int, int, int** ,int, MPI_Win*, MPI_Comm);

int GetCounter(int, int, MPI_Win*);

void DestroyCounter(int, int, MPI_Win*, int*);

void CreateLongCounter(int, int, long** , long, MPI_Win*, MPI_Comm);

long GetLongCounter(int, long, MPI_Win*);

void DestroyLongCounter(int, int, MPI_Win*, long*);

