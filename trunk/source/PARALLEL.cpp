

  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 C. H. Patterson                             *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

#include <mpi.h>
#include <cstdlib>
#include "USER_DATA.h"
#include "PARALLEL.h"

using namespace std;

void mpi_begin_end(int *begin, int *end, int numval, int numtasks, JOB_PARAM *job, FILES file)

{

  int i, remk;

 for (i = 0; i < numtasks; i++) {
   begin[i] = i * (numval / numtasks);
   end[i] = (i + 1) * (numval / numtasks);
   remk = numval - (numval / numtasks) * numtasks;
   if (i < remk) {
     begin[i] += i;
     end[i] += i + 1;
   }
   if (i >= remk) {
     begin[i] += remk;
     end[i] += remk;
   }
   if (job->taskid == 0 && job->verbosity > 1)
   fprintf(file.out, "begin[%d] %d end[%d] %d\n", i, begin[i], i, end[i]);
 }

}

void mpi_begin_end_begin_end(int *begin1, int *end1, int *begin2, int *end2, int numval1, int numval2, int numtasks, JOB_PARAM *job, FILES file)

{

  int i, remk;

 for (i = 0; i < numtasks; i++) {
   begin1[i] = i * (numval1 / numtasks);
   end1[i] = (i + 1) * (numval1 / numtasks);
   remk = numval1 - (numval1 / numtasks) * numtasks;
   if (i < remk) {
     begin1[i] += i;
     end1[i] += i + 1;
   }
   if (i >= remk) {
     begin1[i] += remk;
     end1[i] += remk;
   }
   if (job->taskid == 0 && job->verbosity > 1)
   fprintf(file.out, "%d begin1[%d] %d end1[%d]\n", i, begin1[i], i, end1[i]);
   if (job->taskid == 0) {
   int cores = end1[i] - begin1[i];
   begin2[i] = i * (numval2 / cores);
   end2[i] = (i + 1) * (numval2 / cores);
   remk = numval2 - (numval2 / cores) * cores;
   if (i >= remk) {
     begin2[i] += remk;
     end2[i] += remk;
   }
  }
   fprintf(file.out, "%d begin2[%d] %d end2[%d]\n", i, begin2[i], i, end2[i]);
 }

}

void mpi_begin_end_lower_triangle(int *begin_row, int *end_row, int *begin_array, int *end_array, int numrows, int numtasks, JOB_PARAM *job, FILES file)

{

  int i, remrows, numrows2;
  numrows2 = (numrows * (numrows + 1)) / 2;
  for (i = 0; i < numtasks; i++) {
    begin_row[i] = (int) sqrt((i * 2 * numrows2) / numtasks);
    end_row[i] = (int) sqrt(((i + 1) * 2 * numrows2) / numtasks);
    remrows = numrows - (int) sqrt((numtasks * 2 * numrows2) / numtasks);
    begin_array[i] = (begin_row[i] * (begin_row[i] + 1)) / 2;
    end_array[i] = (end_row[i] * (end_row[i] + 1)) / 2;
    if (job->taskid == 0 && job->verbosity > 1)
    fprintf(file.out, "begin_row[%d] %d end_row[%d] %d begin_array[%d] %d end_array[%d] %d\n", \
    i, begin_row[i], i, end_row[i], i, begin_array[i], i, end_array[i]);
   }

}

void mpi_receive_offset(int *begin, int *end, int *receive, int *offset, int val, JOB_PARAM *job, FILES file)

{

  int i, j;

  for (i = 0; i < job->numtasks; i++) {
    offset[i] = 0;
    receive[i] = 0;
    for (j = begin[i]; j < end[i]; j++) {
      receive[i] += val;
    }
  }
  for (i = 0; i < job->numtasks; i++) {
    for (j = 0; j < i; j++) {
      offset[i] += receive[j];
    }
  }

  if (job->taskid == 0 && job->verbosity > 1)
  for (i = 0; i < job->numtasks; i++) {
    fprintf(file.out, "offset %d receive %d %d %d\n", offset[i], receive[i], begin[i], end[i]);
  }

}

void mpi_receive_offset_pairs(int *begin, int *end, int *receive, int *offset, PAIR_TRAN *pair_p, ATOM *atoms, JOB_PARAM *job, FILES file)

{

  int i, j;

  for (i = 0; i < job->numtasks; i++) {
    offset[i] = 0;
    receive[i] = 0;
    for (j = begin[i]; j < end[i]; j++) {
      receive[i] += atoms->bfnnumb_sh[pair_p->cell1[pair_p->posn[j]]] * atoms->bfnnumb_sh[pair_p->cell2[pair_p->posn[j]]];
    }
  }
  for (i = 0; i < job->numtasks; i++) {
    for (j = 0; j < i; j++) {
      offset[i] += receive[j];
    }
  }

  if (job->taskid == 0 && job->verbosity > 1)
  for (i = 0; i < job->numtasks; i++) {
    fprintf(file.out, "offset %d receive %d %d %d\n", offset[i], receive[i], begin[i], end[i]);
  }

}

void mpi_receive_offset_momentum(int *begin, int *end, int *receive, int *offset, PAIR_TRAN *pair_p, ATOM *atoms, JOB_PARAM *job, FILES file)

{

  int i, j;

  for (i = 0; i < job->numtasks; i++) {
    offset[i] = 0;
    receive[i] = 0;
    for (j = begin[i]; j < end[i]; j++) {
      receive[i] += 3 * atoms->bfnnumb_sh[pair_p->cell1[pair_p->posn[j]]] * atoms->bfnnumb_sh[pair_p->cell2[pair_p->posn[j]]];
    }
  }
  for (i = 0; i < job->numtasks; i++) {
    for (j = 0; j < i; j++) {
      offset[i] += receive[j];
    }
  }

  if (job->taskid == 0 && job->verbosity > 1)
  for (i = 0; i < job->numtasks; i++) {
    fprintf(file.out, "offset %d receive %d %d %d\n", offset[i], receive[i], begin[i], end[i]);
  }

}

void mpi_receive_offset_nblocks(int *begin, int *end, int *receive, int *offset, PAIR_TRAN *pair_p, int *n, ATOM *atoms, JOB_PARAM *job, FILES file)

{

  int i, j;

  for (i = 0; i < job->numtasks; i++) {
    offset[i] = 0;
    receive[i] = 0;
    for (j = begin[i]; j < end[i]; j++) {
      receive[i] += *n * atoms->bfnnumb_sh[pair_p->cell1[pair_p->posn[j]]] * atoms->bfnnumb_sh[pair_p->cell2[pair_p->posn[j]]];
    }
  }
  for (i = 0; i < job->numtasks; i++) {
    for (j = 0; j < i; j++) {
      offset[i] += receive[j];
    }
  }

  if (job->taskid == 0 && job->verbosity > 1)
  for (i = 0; i < job->numtasks; i++) {
    fprintf(file.out, "nblocks %d offset %d receive %d %d %d\n", *n, offset[i], receive[i], begin[i], end[i]);
  }

}

void mpi_local_begin_end(int nloop, int begin_k[], int end_k[], int begin_k_local[], int end_k_local[], JOB_PARAM *job, FILES file)

{

  int i, remk, nklocal;

      i=job->taskid;
      nklocal = end_k[i] - begin_k[i];
      fprintf(file.out,"begink[%d] %d endk[%d] %d\n",i,begin_k[i],i,end_k[i]);

      for(i=0;i<nloop;i++) {
      begin_k_local[i] =     i*(nklocal/nloop) ;
      end_k_local[i]   =    (i+1)*(nklocal/nloop) ;
      remk             =     nklocal - (nklocal/nloop)*nloop ;
      if (i<remk)         { begin_k_local[i] +=    i ; end_k_local[i] += i+1  ; }
      if (i>=remk)        { begin_k_local[i] += remk ; end_k_local[i] += remk ; }
      if (job->verbosity > 1)
      fprintf(file.out,"begin_k_local[%d] %d end_k_local[%d] %d\n",i,begin_k_local[i],i,end_k_local[i]);
      }

      for(i=0;i<nloop;i++) {
      begin_k_local[i] += begin_k[job->taskid];
      end_k_local[i]   += begin_k[job->taskid];
      }

}

/*!   \fn void DestroyCounter(int myrank,int winRank, MPI_Win *win, int *counter)
      \brief Destroys the pointer and the window from the memory of a process
      \param myrank the rank of the process
      \param winRank the rank of the process that holds the window
      \param win Pointer to the window
      \param counter Pointer to the counter

  */

void DestroyCounter(int myrank,int winRank, MPI_Win *win, int *counter)

{

  MPI_Win_free(win);

  if (winRank==myrank) {
    MPI_Free_mem(counter);
  }

}

/*!   \fn int GetCounter(int winRank, int increment , MPI_Win *win)
      \brief Reads and returns the value of the counter from window and then increments it by increment
      \param winRank the rank of the process that holds the window
      \param increment the value by which the counter is incremented
      \param win Pointer to the window
    */

int GetCounter(int winRank, int increment , MPI_Win *win) 

{

  int myCounter;

  MPI_Win_lock(MPI_LOCK_EXCLUSIVE,winRank,0,*win);
  MPI_Get(&myCounter, 1, MPI_INT, winRank, 0, 1, MPI_INT, *win);
  MPI_Accumulate(&increment, 1, MPI_INT, winRank, 0, 1, MPI_INT, MPI_SUM, *win);
  MPI_Win_unlock(winRank,*win);
  //printf("Process %d has the following value for counter: %d \n", winRank,myCounter);

  return myCounter;

}

/*!   \fn void CreateCounter(int rank, int winRank, int **counter, int initialValue, MPI_Win *win, MPI_Comm communicator){
      \brief Creates the counter and the window
      \param rank rank of the current process
      \param winRank the rank of the process that holds the window and the counter
      \param counter double pointer to the counter
      \param initialValue the initial value for the counter
      \param win Pointer to the window
      \param communicator the communicator that would have access to the window.
    */

void CreateCounter(int rank, int winRank, int **counter, int initialValue, MPI_Win *win, MPI_Comm communicator)

{

 int winSize=0;

 if (rank==winRank) {
    winSize=1;
    MPI_Alloc_mem(sizeof(int)*winSize, MPI_INFO_NULL, counter);
    **counter=initialValue;
    //printf("Process %d has the following value for counter: %d \n", winRank,**counter);
  }

  MPI_Win_create(*counter, winSize, 1, MPI_INFO_NULL,communicator, win);

}

void DestroyLongCounter(int myrank,int winRank, MPI_Win *win, long *counter)

{

  MPI_Win_free(win);

  if (winRank==myrank) {
    MPI_Free_mem(counter);
  }

}

/*!   \fn long GetCounter(int winRank, int increment , MPI_Win *win)
      \brief Reads and returns the value of the counter from window and then increments it by increment
      \param winRank the rank of the process that holds the window
      \param increment the value by which the counter is incremented
      \param win Pointer to the window
    */

long GetLongCounter(int winRank, long increment , MPI_Win *win)

{

  long myCounter;

  MPI_Win_lock(MPI_LOCK_EXCLUSIVE,winRank,0,*win);
  MPI_Get(&myCounter, 1, MPI_LONG, winRank, 0, 1, MPI_LONG, *win);
  MPI_Accumulate(&increment, 1, MPI_LONG, winRank, 0, 1, MPI_LONG, MPI_SUM, *win);
  MPI_Win_unlock(winRank,*win);
  //printf("Process %d has the following value for counter: %d \n", winRank,myCounter);

  return myCounter;

}

/*!   \fn void CreateCounter(int rank, int winRank, int **counter, int initialValue, MPI_Win *win, MPI_Comm communicator){
      \brief Creates the counter and the window
      \param rank rank of the current process
      \param winRank the rank of the process that holds the window and the counter
      \param counter double pointer to the counter
      \param initialValue the initial value for the counter
      \param win Pointer to the window
      \param communicator the communicator that would have access to the window.
    */

void CreateLongCounter(int rank, int winRank, long **counter, long initialValue, MPI_Win *win, MPI_Comm communicator)

{

 int winSize=0;

 if (rank==winRank) {
    winSize=1;
    MPI_Alloc_mem(sizeof(long)*winSize, MPI_INFO_NULL, counter);
    **counter=initialValue;
    //printf("Process %d has the following value for counter: %d \n", winRank,**counter);
  }

  MPI_Win_create(*counter, winSize, 1, MPI_INFO_NULL, communicator, win);

}
