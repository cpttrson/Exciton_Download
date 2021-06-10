#ifndef PARALLELH
#define PARALLELH

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

#endif
