
  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 C. H. Patterson                             *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

#include <cstring>
#include <cstdlib>
#include <mkl_scalapack.h>
#include "myconstants.h"
#include "USER_DATA.h"
#include "MATRIX_UTIL.h"
#include "PARALLEL.h"
#include "ALLOCATE_MEMORY_MOLECULE.h"

extern "C" int   Cblacs_pnum(int, int, int); 
extern "C" void  Cblacs_get( int context, int request, int* value);
extern "C" int   Cblacs_gridinit( int* context, char * order, int np_row, int np_col);
extern "C" void  Cblacs_gridinfo( int context, int*  np_row, int* np_col, int*  my_row, int*  my_col);
//extern "C" void  Cblacs_pinfo( int*, int*);
//extern "C" int   Cblacs_gridmap( int* context, int *imap, int ldimap , int np_row, int np_col);
//extern "C" void  Cblacs_gridexit( int context);
//extern "C" void  Cblacs_exit( int error_code);
//extern "C" int   numroc_( int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
//extern "C" void  descinit_( int *desc, int *m, int *n, int *mb, int *nb, int *irsrc, int *icsrc, int *ictxt, int *lld, int *info);

using namespace std;

void initialise_spk_grid(int *gridsize, int *ictxt, int *nbsize, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Initialise scalapack processor grid                                                    *
  // ******************************************************************************************
  
int i, j;
int info = 0, izero = 0, ione = 1;
int cblacs_taskid, itemp, nprow, npcol, myrow, mycol, mpA, nqA;
int descA[9];
double ratio;
char row[2] = "R";

  nprow = 1; 
  npcol = job->numtasks;
  ratio = (double) job->numtasks;

  for (i = 1; i < job->numtasks; i++) {
    if (i * (job->numtasks / i) == job->numtasks) { 
      //printf("%3d %3d %10.4f %10.4f\n",i,job->numtasks / i,(double) i / (double) (job->numtasks / i), ratio);
      if (fabs((double) i / (double) (job->numtasks / i) - k_one) < ratio) {
        j = job->numtasks / i; // i and j are factors of job->numtasks
        ratio = fabs((double) i / (double) j - k_one);
       } 
      } 
     } 

  nprow = job->numtasks / j;
  npcol = j;

  if (job->numtasks == 1) { nprow = 1; npcol = 1; }

  //else if (job->numtasks > 1) { npcol = min(job->numtasks, 14); nprow = job->numtasks / npcol; } // set npcol up to 14 
  //*nbsize = min(1 + ((*gridsize - 1)/ npcol), 256);
  *nbsize = min(1 + ((*gridsize - 1)/ npcol), 64);

  if (job->taskid == 0) printf("ROWS COLS NBSIZE %3d %3d %3d\n",nprow,npcol,*nbsize);
  if (nprow * npcol > job->numtasks) { 
  fprintf(file.out,"nprow * npcol (%d) > numtasks (%d) \n", nprow * npcol, job->numtasks); 
  MPI_Finalize(); 
  exit(1); 
 }
  Cblacs_get(-1, 0, ictxt);
  Cblacs_gridinit(ictxt, row, nprow, npcol);
  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  cblacs_taskid = -1;
  if ((myrow < nprow && myrow >= 0) && (mycol < npcol && mycol >= 0)) {
    mpA = numroc_(gridsize, nbsize, &myrow, &izero, &nprow);
    nqA = numroc_(gridsize, nbsize, &mycol, &izero, &npcol);
    if (mpA == 0 || nqA == 0) {
    fprintf(file.out,"numroc generated a zero row %3d or column %3d dimension for process %3d\n",mpA,nqA,job->taskid); 
    printf("nprow * npcol (%d) > numtasks (%d) %d %d\n", nprow * npcol, job->numtasks,*gridsize,*nbsize); 
    MPI_Finalize(); 
    exit(1); 
   }
    itemp = max(1, mpA);
    descinit_(descA, gridsize, gridsize, nbsize, nbsize, &izero, &izero, ictxt, &itemp, &info);
    cblacs_taskid = Cblacs_pnum(*ictxt,myrow,mycol);
   }

  else { mpA = 1; nqA = 1; }

}

void initialise_spk_grid_crystal(int *gridsize, FERMI *fermi, int nocc_nvir, int *ictxt, int *nbsize_row, int *nbsize_col, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Initialise scalapack processor grid                                                    *
  // ******************************************************************************************
  
int kpts_per_proc_row, kpts_per_proc_col;
int info = 0, izero = 0, ione = 1;
int cblacs_taskid, itemp, nprow, npcol, myrow, mycol, mpA, nqA;
int descA[9];
char row[2] = "R";

  if (job->numtasks == 1) { nprow = 1; npcol = 1; }
  else if (job->numtasks > 1) { npcol = min(job->numtasks, 14); nprow = job->numtasks / npcol; }

  *nbsize_col = *gridsize / fermi->nktot;
  *nbsize_row = *nbsize_col;
  if (job->taskid == 0) printf("ROWS COLS %3d %3d FIXXX NBSIZE ROW %3d COL %3d\n",nprow,npcol,*nbsize_row,*nbsize_col);

  if (nprow * npcol > job->numtasks) { 
  fprintf(file.out,"nprow * npcol (%d) > numtasks (%d) \n", nprow * npcol, job->numtasks); 
  MPI_Finalize(); 
  exit(1); 
 }
  Cblacs_get(-1, 0, ictxt);
  Cblacs_gridinit(ictxt, row, nprow, npcol);
  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  cblacs_taskid = -1;
  if ((myrow < nprow && myrow >= 0) && (mycol < npcol && mycol >= 0)) {
    mpA = numroc_(gridsize, nbsize_row, &myrow, &izero, &nprow);
    nqA = numroc_(gridsize, nbsize_col, &mycol, &izero, &npcol);
    if (mpA == 0 || nqA == 0) {
    printf("numroc generated a zero row %3d or column %3d dimension for process %3d  %2d/%2d %2d/%2d  %3d %3d %3d\n",\
    mpA,nqA,job->taskid,myrow,nprow,mycol,npcol,*gridsize,*nbsize_row,*nbsize_col); 
    MPI_Finalize(); 
    exit(1); 
   }
    itemp = max(1, mpA);
    descinit_(descA, gridsize, gridsize, nbsize_row, nbsize_col, &izero, &izero, ictxt, &itemp, &info);
    cblacs_taskid = Cblacs_pnum(*ictxt,myrow,mycol);
   }

  else { mpA = 1; nqA = 1; }

}

void setup_procs(int *begin_j, int *end_j, int *begin_q, int *end_q, int *MPA, int *NQA, int *dim_ham, int *dim1, int *ictxt, int *ntransitions, int *nbsize_row, int *nbsize_col, FERMI *fermi, FILES file, JOB_PARAM *job)

{

int i, row, col, mpa, nqa;
int cblacs_taskid, itemp, nprow, npcol, myrow, mycol, mpA, nqA, izero = 0, info = 0, ione = 1;
int nprow_nbsize, nprow_myrow, npcol_nbsize, npcol_mycol;
int descA[9];
int num_proc, local_rank;

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(ntransitions, nbsize_row, &myrow, &izero, &nprow);
  nqA = numroc_(ntransitions, nbsize_col, &mycol, &izero, &npcol);
  itemp = max(1, mpA);
  descinit_(descA, ntransitions, ntransitions, nbsize_row, nbsize_col, &izero, &izero, ictxt, &itemp, &info);
  if (job->taskid == 0) printf("mpA %7d nqA %7d blocksizes %7d %7d ntransitions %7d\n",mpA,nqA,*nbsize_row,*nbsize_col,*ntransitions);
  nprow_nbsize = nprow * *nbsize_row;
  npcol_nbsize = npcol * *nbsize_col;
  nprow_myrow = ((nprow + myrow) % nprow) * *nbsize_row;
  npcol_mycol = ((npcol + mycol) % npcol) * *nbsize_col;
  num_proc = job->numtasks / fermi->nkunique; // if num_proc > 1, one or more q points are split over cores
  if (num_proc == 0)   num_proc = 1;
  if (num_proc > *dim1) num_proc = *dim1;
  local_rank = job->taskid % num_proc;
  mpi_begin_end(begin_j, end_j, fermi->nkunique * num_proc, job->numtasks, job, file);
  for (i = 0; i < job->numtasks; i++) { 
    begin_q[i]  = begin_j[i] / num_proc;
    end_q[i]    = end_j[i] / num_proc;
    begin_j[i] *= *dim1;  
    begin_j[i] /= num_proc;
    end_j[i]   *= *dim1;  
    end_j[i]   /= num_proc;
   }
  for (row = 0; row < nprow; row++) {
    mpa = numroc_(ntransitions, nbsize_row, &row, &izero, &nprow);
    for (col = 0; col < npcol; col++) {
      nqa = numroc_(ntransitions, nbsize_col, &col, &izero, &npcol);
      cblacs_taskid = Cblacs_pnum(*ictxt,row,col);
      dim_ham[cblacs_taskid] = mpa * nqa;
      MPA[cblacs_taskid] = mpa;
      NQA[cblacs_taskid] = nqa;
     }
    }

}

void block_cyclic_to_linear(int *nt, int *ictxt, int *nbsize, double *eigvec_buffer1, char *xx, JOB_PARAM *job, FILES file)

{

int i, j, k;
int I1, I2, il, jl, j1;
int local_pointer, myrow_pointer, mycol_pointer, global_pointer, blocksize;
int mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(nt, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(nt, nbsize, &mycol, &izero, &npcol);

  char buf4[120];
  strcpy(buf4,file.scf_eigvec);
  strcat(buf4,xx);
  MPI_File fh;
  MPI_File_open(MPI_COMM_WORLD,buf4,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&fh) ;

  if (job->bse_lim == 0 || job->bse_lim > *nt) job->bse_lim = *nt; // if number of vectors is not set, set to max value

  if ((myrow < nprow && myrow >= 0) && (mycol < npcol && mycol >= 0)) {
  //printf("rank mpA nqA ntrans %3d %3d %3d %3d myrow %3d mycol %3d nb %3d\n",job->taskid,mpA,nqA,*nt,myrow,mycol,*nbsize);
  for (j = 0; j < nqA / *nbsize ; j++) { 
    for (k = 0; k < *nbsize; k++) {
      local_pointer = (j * *nbsize + k) * mpA;
      for (i = 0; i < mpA / *nbsize; i++) { 
        myrow_pointer = i * nprow * *nbsize + ((nprow + myrow) % nprow) * *nbsize;
        mycol_pointer = j * npcol * *nbsize + ((npcol + mycol) % npcol) * *nbsize;
        global_pointer = (mycol_pointer - *nt + job->bse_lim + k) * *nt + myrow_pointer;
        //global_pointer = (mycol_pointer + k) * *nt + myrow_pointer;
        //if (mycol_pointer + k < job->bse_lim) {
        if (mycol_pointer + k >= *nt - job->bse_lim) {
        MPI_File_seek(fh, global_pointer * sizeof(double), MPI_SEEK_SET) ;
        MPI_File_write(fh, &eigvec_buffer1[local_pointer], *nbsize, MPI_DOUBLE, MPI_STATUS_IGNORE);
        //printf("task0 %3d %3d %3d %3d %3d %3d %3d %3d %10.4f\n",\
        job->taskid, j, k, myrow_pointer,mycol_pointer,global_pointer,local_pointer,blocksize,eigvec_buffer1[local_pointer]);
       }
        //printf("task %3d %3d %3d %3d %3d %3d\n",job->taskid, i, j, myrow_pointer,mycol_pointer,global_pointer);
        local_pointer += *nbsize;
       }
      }
     }

  if ((mpA / *nbsize) * *nbsize * nprow + ((nprow + myrow) % nprow) * *nbsize + mpA % *nbsize == *nt) {
    local_pointer = 0;
    blocksize = mpA % *nbsize;
    for (j = 0; j < nqA / *nbsize ; j++) {  // column
      for (k = 0; k < *nbsize; k++) {       // sweep bottom rows
        local_pointer = (j * *nbsize + k) * mpA + (mpA / *nbsize) * *nbsize;
        myrow_pointer = (mpA / *nbsize) * nprow * *nbsize + ((nprow + myrow) % nprow) * *nbsize;
        mycol_pointer = j * npcol * *nbsize + ((npcol + mycol) % npcol) * *nbsize;
        //global_pointer = (mycol_pointer + k) * *nt + myrow_pointer;
        global_pointer = (mycol_pointer - *nt + job->bse_lim + k) * *nt + myrow_pointer;
        if (mycol_pointer + k >= *nt - job->bse_lim) {
        //if (mycol_pointer + k < job->bse_lim) {
        MPI_File_seek(fh, global_pointer * sizeof(double), MPI_SEEK_SET) ;
        MPI_File_write(fh, &eigvec_buffer1[local_pointer], blocksize, MPI_DOUBLE, MPI_STATUS_IGNORE);
        //printf("task1 %3d %3d %3d %3d %3d %3d %3d %3d %10.4f\n",\
        job->taskid, j, k, myrow_pointer,mycol_pointer,global_pointer,local_pointer,blocksize,eigvec_buffer1[local_pointer]);
       }
      }
     }
    }

  if ((nqA / *nbsize) * *nbsize * npcol + ((npcol + mycol) % npcol) * *nbsize + nqA % *nbsize == *nt) {
    blocksize = *nbsize;
    for (j1 = 0; j1 < nqA % *nbsize; j1++) {
      local_pointer = ((nqA / *nbsize) * *nbsize + j1) * mpA;
      for (i = 0; i < mpA / *nbsize ; i++) { // row
        myrow_pointer = i * nprow * *nbsize + ((nprow + myrow) % nprow) * *nbsize;
        mycol_pointer = (nqA / *nbsize) * npcol * *nbsize + ((npcol + mycol) % npcol) * *nbsize + j1;
        //global_pointer = mycol_pointer * *nt + myrow_pointer;
        global_pointer = (mycol_pointer - *nt + job->bse_lim) * *nt + myrow_pointer;
        if (mycol_pointer >= *nt - job->bse_lim) {
        //if (mycol_pointer < job->bse_lim) {
        MPI_File_seek(fh, global_pointer * sizeof(double), MPI_SEEK_SET) ;
        MPI_File_write(fh, &eigvec_buffer1[local_pointer], blocksize, MPI_DOUBLE, MPI_STATUS_IGNORE);
        //printf("task2 %3d %3d %3d %3d %3d %3d %3d %10.4f\n",\
        job->taskid, j, myrow_pointer,mycol_pointer,global_pointer,local_pointer,blocksize,eigvec_buffer1[local_pointer]);
       }
        local_pointer += *nbsize;
       }
      }
     }

  if ((mpA / *nbsize) * *nbsize * nprow + ((nprow + myrow) % nprow) * *nbsize + mpA % *nbsize == *nt  && \
      (nqA / *nbsize) * *nbsize * npcol + ((npcol + mycol) % npcol) * *nbsize + nqA % *nbsize == *nt) {
    blocksize = mpA % *nbsize;
    for (j1 = 0; j1 < nqA % *nbsize; j1++) {
      local_pointer = ((nqA / *nbsize) * *nbsize + j1) * mpA + (mpA / *nbsize) * *nbsize;
      myrow_pointer = (mpA / *nbsize) * nprow * *nbsize + ((nprow + myrow) % nprow) * *nbsize;
      mycol_pointer = (nqA / *nbsize) * npcol * *nbsize + ((npcol + mycol) % npcol) * *nbsize + j1;
      //global_pointer = mycol_pointer * *nt + myrow_pointer;
      global_pointer = (mycol_pointer - *nt + job->bse_lim) * *nt + myrow_pointer;
      if (mycol_pointer >= *nt - job->bse_lim) {
      //if (mycol_pointer < job->bse_lim) {
      MPI_File_seek(fh, global_pointer * sizeof(double), MPI_SEEK_SET) ;
      MPI_File_write(fh, &eigvec_buffer1[local_pointer], blocksize, MPI_DOUBLE, MPI_STATUS_IGNORE);
      //printf("task3 %3d %3d %3d %3d %3d %3d %3d %10.4f\n",\
      job->taskid, j, myrow_pointer,mycol_pointer,global_pointer,local_pointer,blocksize,eigvec_buffer1[local_pointer]);
     }
     }
    }

    } // close if (

  MPI_File_close(&fh);

}

void block_cyclic_to_linear_limit(int *nt, int *ictxt, int *nbsize, int limit, double *eigvec_buffer1, char *xx, JOB_PARAM *job, FILES file)

{

int i, j, k;
int I1, I2, il, jl, j1;
int blocksize;
int mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;
long long local_pointer, global_pointer, myrow_pointer, mycol_pointer;

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(nt, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(nt, nbsize, &mycol, &izero, &npcol);

  char buf4[120];
  strcpy(buf4,file.scf_eigvec);
  strcat(buf4,xx);
  MPI_File fh;
  MPI_File_open(MPI_COMM_WORLD,buf4,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&fh) ;
  //if (job->bse_lim == 0 || job->bse_lim > *nt) job->bse_lim = *nt; // if number of vectors is not set, set to max value
  if (limit == 0 || limit > *nt) limit = *nt; // if number of vectors is not set, set to max value

  if ((myrow < nprow && myrow >= 0) && (mycol < npcol && mycol >= 0)) {
  for (j = 0; j < nqA / *nbsize ; j++) { 
    for (k = 0; k < *nbsize; k++) {
      local_pointer = (j * *nbsize + k) * mpA;
      for (i = 0; i < mpA / *nbsize; i++) { 
        myrow_pointer = i * nprow * *nbsize + ((nprow + myrow) % nprow) * *nbsize;
        mycol_pointer = j * npcol * *nbsize + ((npcol + mycol) % npcol) * *nbsize;
        global_pointer = ((mycol_pointer + k) * *nt + myrow_pointer) * sizeof(double);
        if (mycol_pointer + k < limit) {
        //if (mycol_pointer + k < job->bse_lim) {
        MPI_File_seek(fh, global_pointer, MPI_SEEK_SET) ;
        MPI_File_write(fh, &eigvec_buffer1[local_pointer], *nbsize, MPI_DOUBLE, MPI_STATUS_IGNORE);
       }
        //printf("task %3d %3d %3d %3d %3d %3d\n",job->taskid, i, j, myrow_pointer,mycol_pointer,global_pointer);
        local_pointer += *nbsize;
       }
      }
        //printf("task %3d %3d %3d %lld %lld %lld\n",job->taskid, i, j, myrow_pointer,mycol_pointer,global_pointer);
     }

  if ((mpA / *nbsize) * *nbsize * nprow + ((nprow + myrow) % nprow) * *nbsize + mpA % *nbsize == *nt) {
    local_pointer = 0;
    blocksize = mpA % *nbsize;
    for (j = 0; j < nqA / *nbsize ; j++) {  // column
      for (k = 0; k < *nbsize; k++) {       // sweep bottom rows
        local_pointer = (j * *nbsize + k) * mpA + (mpA / *nbsize) * *nbsize;
        myrow_pointer = (mpA / *nbsize) * nprow * *nbsize + ((nprow + myrow) % nprow) * *nbsize;
        mycol_pointer = j * npcol * *nbsize + ((npcol + mycol) % npcol) * *nbsize;
        global_pointer = ((mycol_pointer + k) * *nt + myrow_pointer) * sizeof(double);
        //global_pointer = (mycol_pointer + k) * *nt + myrow_pointer;
        //if (mycol_pointer + k < job->bse_lim) {
        if (mycol_pointer + k < limit && blocksize > 0) {
        MPI_File_seek(fh, global_pointer, MPI_SEEK_SET) ;
        MPI_File_write(fh, &eigvec_buffer1[local_pointer], blocksize, MPI_DOUBLE, MPI_STATUS_IGNORE);
       }
      }
     }
        //printf("task %3d %3d %lld %lld %lld\n",job->taskid, j, myrow_pointer,mycol_pointer,global_pointer);
    }

  if ((nqA / *nbsize) * *nbsize * npcol + ((npcol + mycol) % npcol) * *nbsize + nqA % *nbsize == *nt) {
    blocksize = *nbsize;
    for (j1 = 0; j1 < nqA % *nbsize; j1++) {
      local_pointer = ((nqA / *nbsize) * *nbsize + j1) * mpA;
      for (i = 0; i < mpA / *nbsize ; i++) { // row
        myrow_pointer = i * nprow * *nbsize + ((nprow + myrow) % nprow) * *nbsize;
        mycol_pointer = (nqA / *nbsize) * npcol * *nbsize + ((npcol + mycol) % npcol) * *nbsize + j1;
        global_pointer = (mycol_pointer * *nt + myrow_pointer) * sizeof(double);
        //global_pointer = mycol_pointer * *nt + myrow_pointer;
        //if (mycol_pointer < job->bse_lim) {
        if (mycol_pointer < limit) {
        MPI_File_seek(fh, global_pointer, MPI_SEEK_SET) ;
        MPI_File_write(fh, &eigvec_buffer1[local_pointer], blocksize, MPI_DOUBLE, MPI_STATUS_IGNORE);
       }
        local_pointer += *nbsize;
       }
      }
        //printf("j1task %3d %3d %lld %lld %lld\n",job->taskid, j1, myrow_pointer,mycol_pointer,global_pointer);
     }

  if ((mpA / *nbsize) * *nbsize * nprow + ((nprow + myrow) % nprow) * *nbsize + mpA % *nbsize == *nt  && \
      (nqA / *nbsize) * *nbsize * npcol + ((npcol + mycol) % npcol) * *nbsize + nqA % *nbsize == *nt) {
    blocksize = mpA % *nbsize;
    for (j1 = 0; j1 < nqA % *nbsize; j1++) {
      local_pointer = ((nqA / *nbsize) * *nbsize + j1) * mpA + (mpA / *nbsize) * *nbsize;
      myrow_pointer = (mpA / *nbsize) * nprow * *nbsize + ((nprow + myrow) % nprow) * *nbsize;
      mycol_pointer = (nqA / *nbsize) * npcol * *nbsize + ((npcol + mycol) % npcol) * *nbsize + j1;
      global_pointer = (mycol_pointer * *nt + myrow_pointer) * sizeof(double);
      //global_pointer = mycol_pointer * *nt + myrow_pointer;
      //if (mycol_pointer < job->bse_lim) {
      if (mycol_pointer < limit) {
      MPI_File_seek(fh, global_pointer, MPI_SEEK_SET) ;
      MPI_File_write(fh, &eigvec_buffer1[local_pointer], blocksize, MPI_DOUBLE, MPI_STATUS_IGNORE);
     }
    }
        //printf("J1task %3d %3d %lld %lld %lld\n",job->taskid, j1, myrow_pointer,mycol_pointer,global_pointer);
   }

  } // close if (

     MPI_File_close(&fh);

}

void block_cyclic_to_linear_limit_complex(int *nt, int *ictxt, int *nbsize, int limit, Complex *eigvec_buffer1, char *xx, JOB_PARAM *job, FILES file)

{

int i, j, k;
int I1, I2, il, jl, j1;
int blocksize;
int mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;
long long local_pointer, global_pointer, myrow_pointer, mycol_pointer;

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  mpA = numroc_(nt, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(nt, nbsize, &mycol, &izero, &npcol);

  char buf4[120];
  strcpy(buf4,file.scf_eigvec);
  strcat(buf4,xx);
  MPI_File fh;
  MPI_File_open(MPI_COMM_WORLD,buf4,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&fh) ;

  //if (job->bse_lim == 0 || job->bse_lim > *nt) job->bse_lim = *nt; // if number of vectors is not set, set to max value
  if (limit == 0 || limit > *nt) limit = *nt; // if number of vectors is not set, set to max value

  if (job->taskid == 0) printf("limit %3d nt %6d \n",limit,*nt);
  if ((myrow < nprow && myrow >= 0) && (mycol < npcol && mycol >= 0)) {
  for (j = 0; j < nqA / *nbsize ; j++) { 
    for (k = 0; k < *nbsize; k++) {
      local_pointer = (j * *nbsize + k) * mpA;
      for (i = 0; i < mpA / *nbsize; i++) { 
        myrow_pointer = i * nprow * *nbsize + ((nprow + myrow) % nprow) * *nbsize;
        mycol_pointer = j * npcol * *nbsize + ((npcol + mycol) % npcol) * *nbsize;
        global_pointer = ((mycol_pointer + k) * *nt + myrow_pointer) * sizeof(Complex);
        if (mycol_pointer + k < limit) {
        //if (mycol_pointer + k < job->bse_lim) {
        MPI_File_seek(fh, global_pointer, MPI_SEEK_SET) ;
        MPI_File_write(fh, &eigvec_buffer1[local_pointer], 2 * *nbsize, MPI_DOUBLE, MPI_STATUS_IGNORE);
       }
        //printf("task %3d %3d %3d %3d %3d %3d\n",job->taskid, i, j, myrow_pointer,mycol_pointer,global_pointer);
        local_pointer += *nbsize;
       }
      }
        //printf("task %3d %3d %3d %lld %lld %lld\n",job->taskid, i, j, myrow_pointer,mycol_pointer,global_pointer);
     }

  if ((mpA / *nbsize) * *nbsize * nprow + ((nprow + myrow) % nprow) * *nbsize + mpA % *nbsize == *nt) {
    local_pointer = 0;
    blocksize = mpA % *nbsize;
    for (j = 0; j < nqA / *nbsize ; j++) {  // column
      for (k = 0; k < *nbsize; k++) {       // sweep bottom rows
        local_pointer = (j * *nbsize + k) * mpA + (mpA / *nbsize) * *nbsize;
        myrow_pointer = (mpA / *nbsize) * nprow * *nbsize + ((nprow + myrow) % nprow) * *nbsize;
        mycol_pointer = j * npcol * *nbsize + ((npcol + mycol) % npcol) * *nbsize;
        global_pointer = ((mycol_pointer + k) * *nt + myrow_pointer) * sizeof(Complex);
        //global_pointer = (mycol_pointer + k) * *nt + myrow_pointer;
        //if (mycol_pointer + k < job->bse_lim) {
        if (mycol_pointer + k < limit && blocksize > 0) {
        MPI_File_seek(fh, global_pointer, MPI_SEEK_SET) ;
        MPI_File_write(fh, &eigvec_buffer1[local_pointer], 2 * blocksize, MPI_DOUBLE, MPI_STATUS_IGNORE);
       }
      }
     }
        //printf("j0task %3d %3d %lld %lld %lld\n",job->taskid, j, myrow_pointer,mycol_pointer,global_pointer);
    }

  if ((nqA / *nbsize) * *nbsize * npcol + ((npcol + mycol) % npcol) * *nbsize + nqA % *nbsize == *nt) {
    blocksize = *nbsize;
    for (j1 = 0; j1 < nqA % *nbsize; j1++) {
      local_pointer = ((nqA / *nbsize) * *nbsize + j1) * mpA;
      for (i = 0; i < mpA / *nbsize ; i++) { // row
        myrow_pointer = i * nprow * *nbsize + ((nprow + myrow) % nprow) * *nbsize;
        mycol_pointer = (nqA / *nbsize) * npcol * *nbsize + ((npcol + mycol) % npcol) * *nbsize + j1;
        global_pointer = (mycol_pointer * *nt + myrow_pointer) * sizeof(Complex);
        //global_pointer = mycol_pointer * *nt + myrow_pointer;
        //if (mycol_pointer < job->bse_lim) {
        if (mycol_pointer < limit) {
        MPI_File_seek(fh, global_pointer, MPI_SEEK_SET) ;
        MPI_File_write(fh, &eigvec_buffer1[local_pointer], 2 * blocksize, MPI_DOUBLE, MPI_STATUS_IGNORE);
       }
        local_pointer += *nbsize;
       }
      }
        //printf("j1task %3d %3d %lld %lld %lld\n",job->taskid, j1, myrow_pointer,mycol_pointer,global_pointer);
     }

  if ((mpA / *nbsize) * *nbsize * nprow + ((nprow + myrow) % nprow) * *nbsize + mpA % *nbsize == *nt  && \
      (nqA / *nbsize) * *nbsize * npcol + ((npcol + mycol) % npcol) * *nbsize + nqA % *nbsize == *nt) {
    blocksize = mpA % *nbsize;
    for (j1 = 0; j1 < nqA % *nbsize; j1++) {
      local_pointer = ((nqA / *nbsize) * *nbsize + j1) * mpA + (mpA / *nbsize) * *nbsize;
      myrow_pointer = (mpA / *nbsize) * nprow * *nbsize + ((nprow + myrow) % nprow) * *nbsize;
      mycol_pointer = (nqA / *nbsize) * npcol * *nbsize + ((npcol + mycol) % npcol) * *nbsize + j1;
      global_pointer = (mycol_pointer * *nt + myrow_pointer) * sizeof(Complex);
      //global_pointer = mycol_pointer * *nt + myrow_pointer;
      //if (mycol_pointer < job->bse_lim) {
      if (mycol_pointer < limit) {
      MPI_File_seek(fh, global_pointer, MPI_SEEK_SET) ;
      MPI_File_write(fh, &eigvec_buffer1[local_pointer], 2 * blocksize, MPI_DOUBLE, MPI_STATUS_IGNORE);
     }
    }
        //printf("J1task %3d %3d %lld %lld %lld\n",job->taskid, j1, myrow_pointer,mycol_pointer,global_pointer);
   }

 } // close if (

     MPI_File_close(&fh);

}

void block_cyclic_zero_triangle(char *uplo, int *nt, int *ictxt, int *nbsize, double *buffer)

{

int il, jl;
int mpA, nqA, nprow, npcol, myrow, mycol, izero = 0;

  Cblacs_gridinfo(*ictxt, &nprow, &npcol, &myrow, &mycol);
  if ((myrow < nprow && myrow >= 0) && (mycol < npcol && mycol >= 0) != 1) return; 
  mpA = numroc_(nt, nbsize, &myrow, &izero, &nprow);
  nqA = numroc_(nt, nbsize, &mycol, &izero, &npcol);

  if (*uplo == 'U') {
  for (il = 0; il < mpA; il++) {
    for (jl = 0; jl < nqA; jl++) {
      if (nprow * *nbsize * (il / *nbsize) + il % *nbsize + ((nprow + myrow) % nprow) * *nbsize < \
          npcol * *nbsize * (jl / *nbsize) + jl % *nbsize + ((npcol + mycol) % npcol) * *nbsize) 
          buffer[il + jl * mpA] = k_zero;
         }
        }
       }

  else if (*uplo == 'L') {
  for (il = 0; il < mpA; il++) {
    for (jl = 0; jl < nqA; jl++) {
      if (nprow * *nbsize * (il / *nbsize) + il % *nbsize + ((nprow + myrow) % nprow) * *nbsize > \
          npcol * *nbsize * (jl / *nbsize) + jl % *nbsize + ((npcol + mycol) % npcol) * *nbsize) 
          buffer[il + jl * mpA] = k_zero;
         }
        }
       }

}

