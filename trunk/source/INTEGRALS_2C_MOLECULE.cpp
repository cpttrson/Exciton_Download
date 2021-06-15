#include <mpi.h>
#include "myconstants.h"
#include "USER_DATA.h"
#include "PARALLEL.h"
#include "MATRIX_UTIL.h"
#include "RECURSION.h"
#include "CARTESIAN_TO_SH.h"
#include "E_COEFFICIENTS.h"
#include "INCOMPLETE_GAMMA.h"
#include "ALLOCATE_MEMORY_MOLECULE.h"
#include "MCMURCHIE_DAVIDSON.h"
#include "INTEGRALS_2C_MOLECULE.h"

using namespace std;

void integrals_molecule_ij(double *Coulomb, PAIR_TRAN *pair_p, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, REAL_LATTICE *R, JOB_PARAM *job, FILES file)

{

int ip, jp, gi, gj, p;
int i, j, i4, j4, count;
int index_i, index_j;
int dim, dim3;
int bfposi, bfposj, bfposi1, bfposj1;
int imax, jmax, im, jm;
int gausposi, gausposj, shelposi, shelposj;
int sheli, shelj, sheli1, shelj1;
int nd1, nd2, nd3, nd4, nd12, nd34;
int begin_p[job->numtasks], end_p[job->numtasks];
double R_AB_1esqrd;
double E1_max;
double time1, time2;
VECTOR_DOUBLE R_AB_1e;

  time1 = MPI_Wtime();

  mpi_begin_end(begin_p,end_p,pair_p->nump,job->numtasks,job,file);
  //printf("process int_1e %d begin %d end %d\n", job->taskid, begin_p[job->taskid],end_p[job->taskid]);

  im = 0; jm = 0;
  dim3 = 0;
  for (p = 0; p < begin_p[job->taskid]; p++) {
  ip = pair_p->cell1[pair_p->posn[p]];
  jp = pair_p->cell2[pair_p->posn[p]];
  nd3 = atoms->bfnnumb_sh[ip];
  nd4 = atoms->bfnnumb_sh[jp];
  dim3 += nd3 * nd4;
 }

  dim = 0;
  //dim3 = 0;
  //for (p = 0; p < pair_p->nump; p++) {
  for (p = begin_p[job->taskid]; p < end_p[job->taskid]; p++) {
    ip = pair_p->cell1[pair_p->posn[p]];
    jp = pair_p->cell2[pair_p->posn[p]];
    gi = 0;
    gj = pair_p->latt2[pair_p->posn[p]];
    nd1 = atoms->bfnnumb[ip];
    nd2 = atoms->bfnnumb[jp];
    nd3 = atoms->bfnnumb_sh[ip];
    nd4 = atoms->bfnnumb_sh[jp];
    nd12 = nd1 * nd2;
    nd34 = nd3 * nd4;
    double coulomb_cart[nd12]; 
    for (i = 0; i < nd12; i++) coulomb_cart[i] = k_zero;
    shelposi = atoms->shelposn_sh[ip];
    gausposi = atoms->gausposn_sh[ip];
    bfposi  = 0;
    bfposi1 = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
      sheli  = shells->type_sh[index_i];
      sheli1 = shells->type1_sh[index_i];
      imax   = shells->imax_sh[index_i];
      shelposj = atoms->shelposn_sh[jp];
      gausposj = atoms->gausposn_sh[jp];
      bfposj  = 0;
      bfposj1 = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
        R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1;
        R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2;
        R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3;
        R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);
        shelj  = shells->type_sh[index_j];
        shelj1 = shells->type1_sh[index_j];
        jmax   = shells->imax_sh[index_j];
        double sab[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double E1x[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+im+jm+1) * (imax+im+1) * (jmax+jm+1)];
        double E1y[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+im+jmax+jm+1) * (imax+im+1) * (jmax+jm+1)];
        double E1z[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+im+jmax+jm+1) * (imax+im+1) * (jmax+jm+1)];
        double pab_inv[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        VECTOR_DOUBLE R_AB[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        E_coefficients1(ip,jp,gi,gj,index_i,index_j,gausposi,gausposj,im,jm,E1x,E1y,E1z,&E1_max,sab,pab_inv,&R_AB_1esqrd,\
        R_AB,R,atoms,shells,gaussians,job,file);
        two_centre_coulomb(coulomb_cart,index_i,index_j,bfposi1,bfposj1,gausposi,gausposj,nd2,im,jm,E1x,E1y,E1z,&R_AB_1e,\
        atoms,shells,gaussians,job,file);
        cartesian_to_sh_ij(coulomb_cart,&Coulomb[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,bfposj1,\
        nd2,nd4,shells,job,file);
        //two_centre_cartesian_to_sh_ij(coulomb_cart,&Coulomb[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,bfposj1,\
        nd2,nd4,shells,job,file);
        bfposj   += shelj;
        bfposj1  += shelj1;
        gausposj += shells->ng_sh[index_j];
       } // close loop over index_j
      bfposi   += sheli;
      bfposi1  += sheli1;
      gausposi += shells->ng_sh[index_i];
     } // close loop over index_i

 if (job->verbosity > 1) {
  count = 0;
   fprintf(file.out,"Coulomb\n");
      fprintf(file.out,"p %3d ip %3d jp %3d gj %3d\n",p,ip,jp,gj);
      for(i=0;i<atoms->bfnnumb_sh[ip];i++) {
        for(j=0;j<atoms->bfnnumb_sh[jp];j++) {
          fprintf(file.out,"%6.2lf",Coulomb[dim3 + dim + count]);
           count++;
          }
         fprintf(file.out,"\n");
        }
        fprintf(file.out,"\n");
       }

    dim += nd34;

  } // close loop on p

  time2 = MPI_Wtime();

}

void nuclear_repulsion_energy(double *nuc_nuc, REAL_LATTICE *R, RECIPROCAL_LATTICE *G, ATOM *atoms, SYMMETRY *symmetry, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  /******************************************************************************************
   * Nuclear-nuclear contribution to total energy                                           *
   ******************************************************************************************/

  int count;
  int index_R, index_S, index_G;
  int atm_no, atm_n1;
  double result1, result2, result3;
  double gamma_0_inv = pow(crystal->primitive_cell_volume , two_thirds) / eight;
  double gamma_0 = k_one / G->gamma_0_inv;
  double root_gamma_0 = sqrt(gamma_0);
  double GdotR;
  double Rsqrd, root_Rsqrd;
  double shell_sum;
  double h, z;
  VECTOR_DOUBLE s_12;

  *nuc_nuc = k_zero;

  for (atm_no = 0; atm_no < atoms->number_of_atoms_in_unit_cell; atm_no++) {
    for (atm_n1 = 0; atm_n1 < atoms->number_of_atoms_in_unit_cell; atm_n1++) {

      result1 = k_zero;
      result2 = k_zero;

      s_12.comp1 = atoms->cell_vector[atm_no].comp1 - atoms->cell_vector[atm_n1].comp1;
      s_12.comp2 = atoms->cell_vector[atm_no].comp2 - atoms->cell_vector[atm_n1].comp2;
      s_12.comp3 = atoms->cell_vector[atm_no].comp3 - atoms->cell_vector[atm_n1].comp3;

         switch (crystal->type[0]) {

           case 'C':

        s_12.comp1 = atoms->cell_vector[atm_no].comp1 - atoms->cell_vector[atm_n1].comp1;
        s_12.comp2 = atoms->cell_vector[atm_no].comp2 - atoms->cell_vector[atm_n1].comp2;
        s_12.comp3 = atoms->cell_vector[atm_no].comp3 - atoms->cell_vector[atm_n1].comp3;

        count = 1;
        for (index_S = 1; index_S < G->number_of_shells; index_S++) {
        for (index_G = 0; index_G < G->num[index_S]; index_G++) {
          GdotR = double_vec_dot(&G->vec_b2[count], &s_12);
          count++;
          if (crystal->type[0] == 'C')
            result1 += G->EXPFAC[index_S] * cos(GdotR);
       }
      }
        count = 0;
      for (index_S = 0; index_S < R->number_of_ewald_shells; index_S++) {
        shell_sum = k_zero;
        for (index_R = 0; index_R < R->num[index_S]; index_R++) {
        s_12.comp1 = atoms->cell_vector[atm_no].comp1 - atoms->cell_vector[atm_n1].comp1 + R->vec_ai[count].comp1;
        s_12.comp2 = atoms->cell_vector[atm_no].comp2 - atoms->cell_vector[atm_n1].comp2 + R->vec_ai[count].comp2;
        s_12.comp3 = atoms->cell_vector[atm_no].comp3 - atoms->cell_vector[atm_n1].comp3 + R->vec_ai[count].comp3;
        Rsqrd = double_vec_dot(&s_12, &s_12);
        root_Rsqrd = sqrt(Rsqrd);
        if (atm_n1 == atm_no && Rsqrd < 0.00001) {
          result2 -= two * root_gamma_0 / rtpi;
        } else {
          result2 += erfc(root_gamma_0 * root_Rsqrd) / root_Rsqrd;
          shell_sum += erfc(root_gamma_0 * root_Rsqrd) / root_Rsqrd;
        } // close else { result2
            count++;
           } // end loop on index_R
              } // end loop over index_S

        result3 = -pi * G->gamma_0_inv / crystal->primitive_cell_volume;

       break;

       case 'S':

        result1 = k_zero;
        result2 = k_zero;

        s_12.comp1 = atoms->cell_vector[atm_no].comp1 - atoms->cell_vector[atm_n1].comp1;
        s_12.comp2 = atoms->cell_vector[atm_no].comp2 - atoms->cell_vector[atm_n1].comp2;
        s_12.comp3 = k_zero;
        z = atoms->cell_vector[atm_no].comp3 - atoms->cell_vector[atm_n1].comp3;

        count = 1;
        for (index_S = 1; index_S < G->number_of_shells; index_S++) {
          h = sqrt(G->sqr[index_S]);
          for (index_G = 0; index_G < G->num[index_S]; index_G++) {
            GdotR = double_vec_dot(&G->vec_b2[count], &s_12);
            result1 += two * pi / crystal->primitive_cell_volume * cos(GdotR) * (exp(h * z) * \
            erfc(z * root_gamma_0 + h / two / root_gamma_0) + \
            exp(-h * z) * erfc(-z * root_gamma_0 + h / two / root_gamma_0)) / h;
            count++;
          }
         }

        count = 0;
        for (index_S = 0; index_S < R->number_of_ewald_shells; index_S++) {
          shell_sum = k_zero;
          for (index_R = 0; index_R < R->num[index_S]; index_R++) {
            s_12.comp1 = atoms->cell_vector[atm_no].comp1 - atoms->cell_vector[atm_n1].comp1 + R->vec_ai[count].comp1;
            s_12.comp2 = atoms->cell_vector[atm_no].comp2 - atoms->cell_vector[atm_n1].comp2 + R->vec_ai[count].comp2;
            s_12.comp3 = atoms->cell_vector[atm_no].comp3 - atoms->cell_vector[atm_n1].comp3 + R->vec_ai[count].comp3;
            Rsqrd = double_vec_dot(&s_12, &s_12);
            root_Rsqrd = sqrt(Rsqrd);
            if (atm_n1 == atm_no && Rsqrd < 0.00001) {
            result2 -= two * root_gamma_0 / rtpi;
           } else {
            result2 += erfc(root_gamma_0 * root_Rsqrd) / root_Rsqrd;
            shell_sum += erfc(root_gamma_0 * root_Rsqrd) / root_Rsqrd;
           } // close else { result2
            count++;
           } // end loop on index_R
          } // end loop over index_S

           result3 = -two * pi / crystal->primitive_cell_volume * \
           (z * erf(root_gamma_0 * z) + exp(-gamma_0 * z * z) / rtpi / root_gamma_0);

       break;

       case 'P':

       result1 =  k_zero;
       result2 =  k_zero;
       result3 =  k_zero;

       if (job->taskid == 0)
       fprintf(file.out,"Nuclear-nuclear repulsion not coded for 1-D systems\n");

       break;

       case 'M':

         result1 =  k_zero;
         result2 =  k_zero;
         result3 =  k_zero;

         s_12.comp1 = atoms->cell_vector[atm_no].comp1 - atoms->cell_vector[atm_n1].comp1;
         s_12.comp2 = atoms->cell_vector[atm_no].comp2 - atoms->cell_vector[atm_n1].comp2;
         s_12.comp3 = atoms->cell_vector[atm_no].comp3 - atoms->cell_vector[atm_n1].comp3;

         if (atm_no != atm_n1) {
         Rsqrd = double_vec_dot(&s_12, &s_12);
         result2 =  k_one / sqrt(Rsqrd);
        }

       break;

      } // close switch

      *nuc_nuc += double(atoms->atomic_number[atm_no] * atoms->atomic_number[atm_n1]) * (result1 + result2 + result3) / two;

    }
  }

  job->nuc_nuc = *nuc_nuc;

}

void fock_element_1e1(INT_1E *one_ints, int dim, PAIR_TRAN *pair_p, int num_p, int Function[8], REAL_LATTICE *R, RECIPROCAL_LATTICE *G, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  int ip, jp, gi, gj, p;
  int i, j, i4, j4, k, l, mm, nn, atm0, count;
  int index_i, index_j;
  int dim1, dim2, dim3, dim4, dim5, dim6, dim7, dim8, dim9, atm1, atm2;
  int bfposi, bfposj, bfposi1, bfposj1;
  int imax, jmax, im, jm;
  int gausposi, gausposj, shelposi, shelposj;
  int sheli, shelj, sheli1, shelj1;
  int nd1, nd2, nd3, nd4, nd12, nd34;
  int nsheli, nshelj;
  int begin_p[job->numtasks], end_p[job->numtasks];
  double R_AB_1esqrd;
  double E1_max;
  double time1, time2;
  VECTOR_DOUBLE R_AB_1e;

  // Function[0] Fock Operator
  // Function[1] Kinetic Energy Operator
  // Function[2] Electron-Nuclear Potential Energy = 1 (finite systems) = 2 (extended systems) = 3 (atoms)
  // Function[3] Momentum Operator Matrix Elements
  // Function[4] Overlap
  // Function[5] 
  // Function[6] Dipole Operator Matrix Elements
  // Function[7] Coulomb Operator Matrix Elements

  time1 = MPI_Wtime();

  //mpi_begin_end(begin_p,end_p,num_p,job->numtasks,job,file);
  //printf("process int_1e %d begin %d end %d\n", job->taskid, begin_p[job->taskid],end_p[job->taskid]);

  //dim3 = 0;
  //for (p = 0; p < begin_p[job->taskid]; p++) {
  //ip = pair_p->cell1[pair_p->posn[p]];
  //jp = pair_p->cell2[pair_p->posn[p]];
  //nd3 = atoms->bfnnumb_sh[ip];
  //nd4 = atoms->bfnnumb_sh[jp];
  //dim3 += nd3 * nd4;
 //}
    //dim = 0;
    //dim4 = 0;
  
  //for (p = begin_p[job->taskid]; p < end_p[job->taskid]; p++) {


  int *counter, winRank, myCounter, increment;
  const int izero = 0;
  MPI_Win win;
  winRank = 0;
  CreateCounter(job->taskid, winRank, &counter, izero, &win, MPI_COMM_WORLD);
  myCounter = 0;

  if (job->taskid > 0 || job->numtasks == 1)
  while (myCounter < pair_p->nump) {

    increment = (pair_p->nump - myCounter) / job->numtasks / 2;
    if (increment <  1) increment =  1;
    if (increment > 64) increment = 64;
    myCounter = GetCounter(winRank, increment , &win);
    //printf("proc %d myCounter %d\n",job->taskid,myCounter);
    int begin_p1 = (myCounter             < pair_p->nump) ? myCounter : pair_p->nump;
    int end_p1   = (myCounter + increment < pair_p->nump) ? myCounter + increment : pair_p->nump;
    //printf("proc %3d Fock begin %6d end %6d out of %6d\n",job->taskid,begin_p1,end_p1,pair_p->nump);
    dim = 0;
    dim3 = 0;
    dim4 = 0;
    for (p = 0; p < begin_p1; p++) {
    ip = pair_p->cell1[pair_p->posn[p]];
    jp = pair_p->cell2[pair_p->posn[p]];
    nd3 = atoms->bfnnumb_sh[ip];
    nd4 = atoms->bfnnumb_sh[jp];
    dim3 += nd3 * nd4;
   }

  for (p = begin_p1; p < end_p1; p++) {

    ip = pair_p->cell1[pair_p->posn[p]];
    jp = pair_p->cell2[pair_p->posn[p]];
    gi = 0;
    gj = pair_p->latt2[pair_p->posn[p]];
    nd1 = atoms->bfnnumb[ip];
    nd2 = atoms->bfnnumb[jp];
    nd3 = atoms->bfnnumb_sh[ip];
    nd4 = atoms->bfnnumb_sh[jp];
    nd12 = nd1 * nd2;
    nd34 = nd3 * nd4;
    double kinetic_cart[nd12]; 
    double elecnuc_cart[nd12]; 
    double momentum_cart[3 * nd12]; 
    double overlap_cart[nd12]; 
    double dipole_cart[3 * nd12]; 
    double coulomb_cart[nd12]; 
    // order matters for these loops to set im and jm
    im = 0;
    jm = 0;
    if (Function[0] || Function[2]) {
    for (i = 0; i < nd12; i++) elecnuc_cart[i] = k_zero;
    im = 0; jm = 0;
   }
    if (Function[4]) {
    for (i = 0; i < nd12; i++) overlap_cart[i] = k_zero;
    im = 0; jm = 0;
   }
    if (Function[7]) {
    for (i = 0; i < nd12; i++) coulomb_cart[i] = k_zero;
    im = 0; jm = 0;
   }
    if (Function[3]) {
    for (i = 0; i < 3 * nd12; i++) momentum_cart[i] = k_zero;
    im = 1; jm = 1;
   }
    if (Function[6]) {
    for (i = 0; i < 3 * nd12; i++) dipole_cart[i] = k_zero;
    im = 1; jm = 1;
   }
    if (Function[0] || Function[1]) {
    for (i = 0; i < nd12; i++) kinetic_cart[i] = k_zero;
    im = 2; jm = 2;
   }
    shelposi = atoms->shelposn_sh[ip];
    gausposi = atoms->gausposn_sh[ip];
    bfposi  = 0;
    bfposi1 = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
      sheli  = shells->type_sh[index_i];
      sheli1 = shells->type1_sh[index_i];
      imax   = shells->imax_sh[index_i];
      shelposj = atoms->shelposn_sh[jp];
      gausposj = atoms->gausposn_sh[jp];
      bfposj  = 0;
      bfposj1 = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
        R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1 - R->vec_ai[gj].comp1;
        R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2 - R->vec_ai[gj].comp2;
        R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3 - R->vec_ai[gj].comp3;
        R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);
        shelj  = shells->type_sh[index_j];
        shelj1 = shells->type1_sh[index_j];
        jmax   = shells->imax_sh[index_j];
        double sab[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double E1x[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+im+jm+1) * (imax+im+1) * (jmax+jm+1)];
        double E1y[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+im+jmax+jm+1) * (imax+im+1) * (jmax+jm+1)];
        double E1z[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+im+jmax+jm+1) * (imax+im+1) * (jmax+jm+1)];
        double pab_inv[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        VECTOR_DOUBLE R_AB[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        E_coefficients1(ip,jp,gi,gj,index_i,index_j,gausposi,gausposj,im,jm,E1x,E1y,E1z,&E1_max,sab,pab_inv,&R_AB_1esqrd,\
        R_AB,R,atoms,shells,gaussians,job,file);
    if (Function[0] || Function[1]) {
        Kinetic(kinetic_cart, index_i, index_j, bfposi1, bfposj1, gausposi, gausposj, nd2, im, jm, E1x, E1y, E1z, sab, \
        shells, gaussians, job, file);
        //two_centre_cartesian_to_sh_shell_ij(kinetic_cart,&one_ints->Kinetic[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,\
        bfposj1,nd2,nd4,shells,job,file);
        //two_centre_cartesian_to_sh_ij(kinetic_cart,&one_ints->Kinetic[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,\
        bfposj1,nd2,nd4,shells,job,file);
        cartesian_to_sh_ij(kinetic_cart,&one_ints->Kinetic[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,\
        bfposj1,nd2,nd4,shells,job,file);
       }
    if (Function[0] || Function[2]) {
        ElecNuc(elecnuc_cart, index_i, index_j, bfposi1, bfposj1, nd2,  im, jm, E1x, E1y, E1z, R_AB, pab_inv, sab, R, G, \
        crystal, atoms,shells,job,file);
        //two_centre_cartesian_to_sh_shell_ij(elecnuc_cart,&one_ints->ElecNuc[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,\
        bfposj1,nd2,nd4,shells,job,file);
        //two_centre_cartesian_to_sh_ij(elecnuc_cart,&one_ints->ElecNuc[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,\
        bfposj1,nd2,nd4,shells,job,file);
        cartesian_to_sh_ij(elecnuc_cart,&one_ints->ElecNuc[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,\
        bfposj1,nd2,nd4,shells,job,file);
       }
    if (Function[3]) {
        Momentum(momentum_cart, index_i, index_j, bfposi1, bfposj1, gausposi, gausposj, nd1, nd2, im, jm, E1x, E1y, E1z, sab, \
        shells, gaussians, job,file);
        //two_centre_vector_cartesian_to_sh_shell(momentum_cart,&one_ints->Momentum[3 * dim3 + dim4],index_i,index_j,bfposi,bfposj,\
        bfposi1,bfposj1,nd1,nd2,nd3,nd4,shells,job,file);
        //two_centre_cartesian_to_sh_ij_vector(momentum_cart,&one_ints->Momentum[3 * (dim3 + dim)],index_i,index_j,bfposi,bfposj,\
        bfposi1,bfposj1,nd1,nd2,nd3,nd4,shells,job,file);
        cartesian_to_sh_ij_vector(momentum_cart,&one_ints->Momentum[3 * (dim3 + dim)],index_i,index_j,bfposi,bfposj,\
        bfposi1,bfposj1,nd1,nd2,nd3,nd4,shells,job,file);
       }
    if (Function[4]) {
        Overlap(overlap_cart, index_i, index_j, bfposi1, bfposj1, nd2, im, jm, E1x, E1y, E1z, sab, shells, job, file);
        //two_centre_cartesian_to_sh_shell_ij(overlap_cart,&one_ints->Overlap[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,\
        bfposj1,nd2,nd4,shells,job,file);
        //two_centre_cartesian_to_sh_ij(overlap_cart,&one_ints->Overlap[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,\
        bfposj1,nd2,nd4,shells,job,file);
        cartesian_to_sh_ij(overlap_cart,&one_ints->Overlap[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,\
        bfposj1,nd2,nd4,shells,job,file);
       }
    if (Function[6]) {
        Dipole(dipole_cart,index_i,index_j,bfposi1,bfposj1,gausposi,gausposj,nd1,nd2,ip,im,jm,E1x,E1y,E1z,sab,R_AB,atoms,shells,\
        gaussians,job,file);
        //fix dim
        //two_centre_vector_cartesian_to_sh_shell_ij(dipole_cart,&one_ints->Dipole[3 * (dim3 + dim)],index_i,index_j,bfposi,bfposj,\
        bfposi1,bfposj1,nd1,nd2,nd3,nd4,shells,job,file);
        //two_centre_cartesian_to_sh_ij_vector(dipole_cart,&one_ints->Dipole[3 * (dim3 + dim)],index_i,index_j,bfposi,bfposj,\
        bfposi1,bfposj1,nd1,nd2,nd3,nd4,shells,job,file);
        cartesian_to_sh_ij_vector(dipole_cart,&one_ints->Dipole[3 * (dim3 + dim)],index_i,index_j,bfposi,bfposj,\
        bfposi1,bfposj1,nd1,nd2,nd3,nd4,shells,job,file);
       }
    if (Function[7]) {
        two_centre_coulomb(coulomb_cart,index_i,index_j,bfposi1,bfposj1,gausposi,gausposj,nd2,im,jm,E1x,E1y,E1z,&R_AB_1e,\
        atoms,shells,gaussians,job,file);
        //two_centre_cartesian_to_sh_ij(coulomb_cart,&one_ints->Coulomb[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,bfposj1,\
        nd2,nd4,shells,job,file);
        cartesian_to_sh_ij(coulomb_cart,&one_ints->Coulomb[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,bfposj1,\
        nd2,nd4,shells,job,file);
       }
        bfposj   += shelj;
        bfposj1  += shelj1;
        gausposj += shells->ng_sh[index_j];
       } // close loop over index_j
      bfposi   += sheli;
      bfposi1  += sheli1;
      gausposi += shells->ng_sh[index_i];
     } // close loop over index_i

 if (job->verbosity > 1) {
  count = 0;
   fprintf(file.out,"ElecNuc\n");
      fprintf(file.out,"XYZ %3d ip %3d jp %3d gj %3d\n",p,ip,jp,gj);
      for(i=0;i<atoms->bfnnumb_sh[ip];i++) {
        for(j=0;j<atoms->bfnnumb_sh[jp];j++) {
          fprintf(file.out,"%6.2lf",one_ints->Overlap[dim3 + dim + count]);
           count++;
          }
         fprintf(file.out,"\n");
        }
        fprintf(file.out,"\n");
       }

    dim += nd34;
    dim4 += 3 * nd34;

  } // close loop on p

    if (Function[0]) {
      for (i = dim3; i < dim3 + dim; i++) {
        one_ints->Fock[i] = one_ints->Kinetic[i] + one_ints->ElecNuc[i];
       }
      }

  } // close while (
        
   DestroyCounter(job->taskid, winRank, &win, counter);

   time2 = MPI_Wtime();

}

void fock_element_1e2(INT_1E *one_ints, PAIR_TRAN *pair_p, int Function[8], REAL_LATTICE *R, RECIPROCAL_LATTICE *G, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

int ip, jp, gi, gj, p;
int i, j, i4, j4, count;
int index_i, index_j;
int dim, dim3, dim4, dimall, dimall3, dimg;
int bfposi, bfposj, bfposi1, bfposj1;
int imax, jmax, im, jm;
int gausposi, gausposj, shelposi, shelposj;
int sheli, shelj, sheli1, shelj1;
int nd1, nd2, nd3, nd4, nd12, nd34;
int nsheli, nshelj;
int begin_p[job->numtasks], end_p[job->numtasks];
double R_AB_1esqrd;
double E1_max;
double time1, time2;
VECTOR_DOUBLE R_AB_1e;
INT_1E one_ints_buffer;

  // Function[0] Fock Operator
  // Function[1] Kinetic Energy Operator
  // Function[2] Electron-Nuclear Potential Energy = 1 (finite systems) = 2 (extended systems) = 3 (atoms)
  // Function[3] Momentum Operator Matrix Elements
  // Function[4] Overlap
  // Function[5] 
  // Function[6] Dipole Operator Matrix Elements
  // Function[7] Coulomb Operator Matrix Elements

  time1 = MPI_Wtime();

  mpi_begin_end(begin_p,end_p,pair_p->nump,job->numtasks,job,file);
  //printf("process int_1e %d begin %d end %d\n", job->taskid, begin_p[job->taskid],end_p[job->taskid]);
  array_dimensions(&dimall, &dimg, pair_p, atoms, job, file); // don't change
  allocate_INT_1E(&one_ints_buffer, dimall, Function, job, file);
  dimall3 = 3 * dimall;

  dim3 = 0;
  for (p = 0; p < begin_p[job->taskid]; p++) {
  ip = pair_p->cell1[pair_p->posn[p]];
  jp = pair_p->cell2[pair_p->posn[p]];
  nd3 = atoms->bfnnumb_sh[ip];
  nd4 = atoms->bfnnumb_sh[jp];
  dim3 += nd3 * nd4;
 }

  dim = 0;
  dim4 = 0;
  for (p = begin_p[job->taskid]; p < end_p[job->taskid]; p++) {
    ip = pair_p->cell1[pair_p->posn[p]];
    jp = pair_p->cell2[pair_p->posn[p]];
    gi = 0;
    gj = pair_p->latt2[pair_p->posn[p]];
    nd1 = atoms->bfnnumb[ip];
    nd2 = atoms->bfnnumb[jp];
    nd3 = atoms->bfnnumb_sh[ip];
    nd4 = atoms->bfnnumb_sh[jp];
    nd12 = nd1 * nd2;
    nd34 = nd3 * nd4;
    double kinetic_cart[nd12]; 
    double elecnuc_cart[nd12]; 
    double momentum_cart[3 * nd12]; 
    double overlap_cart[nd12]; 
    double dipole_cart[3 * nd12]; 
    double coulomb_cart[nd12]; 
    // order matters for these loops to set im and jm
    im = 0;
    jm = 0;
    if (Function[0] || Function[2]) {
    for (i = 0; i < nd12; i++) elecnuc_cart[i] = k_zero;
    im = 0; jm = 0;
   }
    if (Function[4]) {
    for (i = 0; i < nd12; i++) overlap_cart[i] = k_zero;
    im = 0; jm = 0;
   }
    if (Function[7]) {
    for (i = 0; i < nd12; i++) coulomb_cart[i] = k_zero;
    im = 0; jm = 0;
   }
    if (Function[3]) {
    for (i = 0; i < 3 * nd12; i++) momentum_cart[i] = k_zero;
    im = 1; jm = 1;
   }
    if (Function[6]) {
    for (i = 0; i < 3 * nd12; i++) dipole_cart[i] = k_zero;
    im = 1; jm = 1;
   }
    if (Function[0] || Function[1]) {
    for (i = 0; i < nd12; i++) kinetic_cart[i] = k_zero;
    im = 2; jm = 2;
   }
    shelposi = atoms->shelposn_sh[ip];
    gausposi = atoms->gausposn_sh[ip];
    bfposi  = 0;
    bfposi1 = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
      sheli  = shells->type_sh[index_i];
      sheli1 = shells->type1_sh[index_i];
      imax   = shells->imax_sh[index_i];
      shelposj = atoms->shelposn_sh[jp];
      gausposj = atoms->gausposn_sh[jp];
      bfposj  = 0;
      bfposj1 = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
        R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1 - R->vec_ai[gj].comp1;
        R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2 - R->vec_ai[gj].comp2;
        R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3 - R->vec_ai[gj].comp3;
        R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);
        shelj  = shells->type_sh[index_j];
        shelj1 = shells->type1_sh[index_j];
        jmax   = shells->imax_sh[index_j];
        double sab[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double E1x[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+im+jm+1) * (imax+im+1) * (jmax+jm+1)];
        double E1y[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+im+jmax+jm+1) * (imax+im+1) * (jmax+jm+1)];
        double E1z[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+im+jmax+jm+1) * (imax+im+1) * (jmax+jm+1)];
        double pab_inv[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        VECTOR_DOUBLE R_AB[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        E_coefficients1(ip,jp,gi,gj,index_i,index_j,gausposi,gausposj,im,jm,E1x,E1y,E1z,&E1_max,sab,pab_inv,&R_AB_1esqrd,\
        R_AB,R,atoms,shells,gaussians,job,file);
    if (Function[0] || Function[1]) {
        Kinetic(kinetic_cart, index_i, index_j, bfposi1, bfposj1, gausposi, gausposj, nd2, im, jm, E1x, E1y, E1z, sab, \
        shells, gaussians, job, file);
        //two_centre_cartesian_to_sh_ij(kinetic_cart,&one_ints_buffer.Kinetic[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,\
        bfposj1,nd2,nd4,shells,job,file);
        cartesian_to_sh_ij(kinetic_cart,&one_ints_buffer.Kinetic[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,\
        bfposj1,nd2,nd4,shells,job,file);
       }
    if (Function[0] || Function[2]) {
        ElecNuc(elecnuc_cart, index_i, index_j, bfposi1, bfposj1, nd2,  im, jm, E1x, E1y, E1z, R_AB, pab_inv, sab, R, G, \
        crystal, atoms,shells,job,file);
        //two_centre_cartesian_to_sh_ij(elecnuc_cart,&one_ints_buffer.ElecNuc[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,\
        bfposj1,nd2,nd4,shells,job,file);
        cartesian_to_sh_ij(elecnuc_cart,&one_ints_buffer.ElecNuc[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,\
        bfposj1,nd2,nd4,shells,job,file);
       }
    if (Function[3]) {
        Momentum(momentum_cart, index_i, index_j, bfposi1, bfposj1, gausposi, gausposj, nd1, nd2, im, jm, E1x, E1y, E1z, sab, \
        shells, gaussians, job,file);
        //two_centre_cartesian_to_sh_ij_vector(momentum_cart,&one_ints_buffer.Momentum[3 * dim3 + dim4],index_i,index_j,bfposi,bfposj,\
        bfposi1,bfposj1,nd1,nd2,nd3,nd4,shells,job,file);
        cartesian_to_sh_ij_vector(momentum_cart,&one_ints_buffer.Momentum[3 * dim3 + dim4],index_i,index_j,bfposi,bfposj,\
        bfposi1,bfposj1,nd1,nd2,nd3,nd4,shells,job,file);
       }
    if (Function[4]) {
        Overlap(overlap_cart, index_i, index_j, bfposi1, bfposj1, nd2, im, jm, E1x, E1y, E1z, sab, shells, job, file);
        //two_centre_cartesian_to_sh_ij(overlap_cart,&one_ints_buffer.Overlap[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,\
        bfposj1,nd2,nd4,shells,job,file);
        cartesian_to_sh_ij(overlap_cart,&one_ints_buffer.Overlap[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,\
        bfposj1,nd2,nd4,shells,job,file);
       }
    if (Function[6]) {
        Dipole(dipole_cart,index_i,index_j,bfposi1,bfposj1,gausposi,gausposj,nd1,nd2,ip,im,jm,E1x,E1y,E1z,sab,R_AB,atoms,shells,\
        gaussians,job,file);
        //two_centre_cartesian_to_sh_ij_vector(dipole_cart,&one_ints_buffer.Dipole[3 * dim3 + dim4],index_i,index_j,bfposi,bfposj,\
        bfposi1,bfposj1,nd1,nd2,nd3,nd4,shells,job,file);
        cartesian_to_sh_ij_vector(dipole_cart,&one_ints_buffer.Dipole[3 * dim3 + dim4],index_i,index_j,bfposi,bfposj,\
        bfposi1,bfposj1,nd1,nd2,nd3,nd4,shells,job,file);
       }
    if (Function[7]) {
        two_centre_coulomb(coulomb_cart,index_i,index_j,bfposi1,bfposj1,gausposi,gausposj,nd2,im,jm,E1x,E1y,E1z,&R_AB_1e,\
        atoms,shells,gaussians,job,file);
        //two_centre_cartesian_to_sh_ij(coulomb_cart,&one_ints_buffer.Coulomb[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,bfposj1,\
        nd2,nd4,shells,job,file);
        cartesian_to_sh_ij(coulomb_cart,&one_ints_buffer.Coulomb[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1,bfposj1,\
        nd2,nd4,shells,job,file);
       }
        bfposj   += shelj;
        bfposj1  += shelj1;
        gausposj += shells->ng_sh[index_j];
       } // close loop over index_j
      bfposi   += sheli;
      bfposi1  += sheli1;
      gausposi += shells->ng_sh[index_i];
     } // close loop over index_i

 if (job->verbosity > 1) {
  count = 0;
   fprintf(file.out,"Coulomb\n");
      fprintf(file.out,"p %3d ip %3d jp %3d gj %3d\n",p,ip,jp,gj);
      for(i=0;i<atoms->bfnnumb_sh[ip];i++) {
        for(j=0;j<atoms->bfnnumb_sh[jp];j++) {
          fprintf(file.out,"%6.2lf",one_ints_buffer.Coulomb[dim3 + dim + count]);
           count++;
          }
         fprintf(file.out,"\n");
        }
        fprintf(file.out,"\n");
       }

    dim += nd34;
    dim4 += 3 * nd34;

  } // close loop on p

  // Function[0] Fock Operator
  // Function[1] Kinetic Energy Operator
  // Function[2] Electron-Nuclear Potential Energy = 1 (finite systems) = 2 (extended systems) = 3 (atoms)
  // Function[3] Momentum Operator Matrix Elements
  // Function[4] Overlap
  // Function[5] 
  // Function[6] Dipole Operator Matrix Elements
  // Function[7] Coulomb Operator Matrix Elements

  if (Function[1]) MPI_Allreduce(&one_ints_buffer.Kinetic[0],  &one_ints->Kinetic[0],  dimall,  MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  if (Function[2]) MPI_Allreduce(&one_ints_buffer.ElecNuc[0],  &one_ints->ElecNuc[0],  dimall,  MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  if (Function[3]) MPI_Allreduce(&one_ints_buffer.Momentum[0], &one_ints->Momentum[0], dimall3, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  if (Function[4]) MPI_Allreduce(&one_ints_buffer.Overlap[0],  &one_ints->Overlap[0],  dimall,  MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  if (Function[6]) MPI_Allreduce(&one_ints_buffer.Dipole[0],   &one_ints->Dipole[0],   dimall,  MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  if (Function[7]) MPI_Allreduce(&one_ints_buffer.Coulomb[0],  &one_ints->Coulomb[0],  dimall  ,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  free_INT_1E(&one_ints_buffer, Function, job, file);

  if (Function[0]) {
    for (i = 0; i < dimall; i++) {
      one_ints->Fock[i] = one_ints->Kinetic[i] + one_ints->ElecNuc[i];
     }
    }

   time2 = MPI_Wtime();

}

void Kinetic(double *Kinetic, int index_i, int index_j, int bfposi, int bfposj, int gausposi, int gausposj, int nd2,  int im, int jm, double *E1x, double *E1y, double *E1z, double *sab, SHELL *shells, GAUSSIAN *gaussians, JOB_PARAM *job, FILES file)

{

  int i, j, i1, i4, j4;
  int n1, n2, n3, n4, n5, n6;
  int sheli, shelj;
  int imax, jmax;
  int off1, off2, off3;
  int count;

  imax = shells->imax_sh[index_i];
  jmax = shells->imax_sh[index_j];
  off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2 = (imax + im + jmax + jm + 1) * off3;
  off1 = (jmax + jm + 1) * off2;
  sheli = shells->type1_sh[index_i];
  shelj = shells->type1_sh[index_j];
  for (i = 0; i < sheli; i++) {
    for (j = 0; j < shelj; j++) {
      n1 = shells->tuv[imax][i][0];
      n2 = shells->tuv[imax][i][1];
      n3 = shells->tuv[imax][i][2];
      n4 = shells->tuv[jmax][j][0];
      n5 = shells->tuv[jmax][j][1];
      n6 = shells->tuv[jmax][j][2];
      i1 = 0;
      count = (bfposi + i) * nd2 + bfposj + j;
      for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
      for (j4 = 0; j4 < shells->ng_sh[index_j]; j4++) {
      if (n4 > 1)
      Kinetic[count] -= (double)(n4 * (n4 - 1)) * E1x[n1 * off1 + (n4 - 2) * off2 + i1] * E1y[n2 * off1 + n5 * off2 + i1] * \
                         E1z[n3 * off1 + n6 * off2 + i1] * sab[i1] / two;
      if (n5 > 1)
      Kinetic[count] -= (double)(n5 * (n5 - 1)) * E1x[n1 * off1 + n4 * off2 + i1] * E1y[n2 * off1 + (n5 - 2) * off2 + i1] * \
                         E1z[n3 * off1 + n6 * off2 + i1] * sab[i1] / two;
      if (n6 > 1)
      Kinetic[count] -= (double)(n6 * (n6 - 1)) * E1x[n1 * off1 + n4 * off2 + i1] * E1y[n2 * off1 + n5 * off2 + i1] * \
                         E1z[n3 * off1 + (n6 - 2) * off2 + i1] * sab[i1] / two;
      Kinetic[count] -= (-two * gaussians->expo_sh[gausposj + j4] * (double)(2 * (n4 + n5 + n6) + 3) * \
                    E1x[n1 * off1 + n4 * off2 + i1] * E1y[n2 * off1 + n5 * off2 + i1] * E1z[n3 * off1 + n6 * off2 + i1] + \
                    four * gaussians->expo_sh[gausposj + j4] * gaussians->expo_sh[gausposj + j4] * \
                   (E1x[n1 * off1 + (n4 + 2) * off2 + i1] * E1y[n2 * off1 + n5 * off2 + i1] * E1z[n3 * off1 + n6 * off2 + i1] + \
                    E1x[n1 * off1 + n4 * off2 + i1] * E1y[n2 * off1 + (n5 + 2) * off2 + i1] * E1z[n3 * off1 + n6 * off2 + i1] + \
                    E1x[n1 * off1 + n4 * off2 + i1] * E1y[n2 * off1 + n5 * off2 + i1] * E1z[n3 * off1 + (n6 + 2) * off2 + i1])) * sab[i1] / two;
       i1++;
      } // close loop over j4
     } // close loop over i4
    } // close loop over i
   } // close loop over i

}

void ElecNuc(double *ElecNuc, int index_i, int index_j, int bfposi, int bfposj, int nd2,  int im, int jm, double *E1x, double *E1y, double *E1z, VECTOR_DOUBLE *R_AB, double *pab_inv, double *sab, REAL_LATTICE *R, RECIPROCAL_LATTICE *G, CRYSTAL *crystal, ATOM *atoms, SHELL *shells, JOB_PARAM *job, FILES file)

{

  int i, j, k, i4, nn, t, u, v;
  int n1, n2, n3, n4, n5, n6;
  int tmax, umax, vmax;
  int sheli, shelj;
  int off1, off2, off3;
  int imax, jmax;
  int count;
  int atm0;
  int count1, index_G, index_R, index_S;
  double Rsqrd, fn[55], gn[55];
  double gamma_0, root_gamma_0;
  double r1, r2, r3, r4, GdotR;
  VECTOR_DOUBLE r_12, s_12, t_12, Rvec_tmp;

  double em[1][55];
  double f[13][13][13][13];
  double temp, fgtuv_max, fgtuv_temp;
  double *fgtuv;
  int mm, mm0, mm1, mm2, mm3, mm4;

  double Az;

  imax = shells->imax_sh[index_i];
  jmax = shells->imax_sh[index_j];
  off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2 = (imax + im + jmax + jm + 1) * off3;
  off1 = (jmax + jm + 1) * off2;
  sheli = shells->type1_sh[index_i];
  shelj = shells->type1_sh[index_j];

  mm  = imax + jmax;
  mm3 = off3;
  mm2 = (mm + 1) * mm3;
  mm1 = (mm + 1) * mm2;
  mm0 = (mm + 1) * mm1;

  fgtuv = (double *) malloc(mm0 * sizeof(double));
  if (fgtuv == NULL) {
  if (job->taskid == 0)
  fprintf(stderr, "ERROR: not enough memory for double fgtuv! \n");
  MPI_Finalize();
  exit(1);
 }

int n;
double r5, Rzsqrd, sum5[mm + 1], sum2[mm + 1], x;
double fac1, fac2, expfac, x1, x2, sign;
double D_erf[16 + 1], D_exp[16 + 1], derivative_1[16 + 1], derivative_2[16 + 1];

  switch (crystal->type[0]) {

  case 'M':

  for (atm0 = 0; atm0 < atoms->number_of_atoms_in_unit_cell; atm0++) {
    for (k = 0; k < mm0; k++) fgtuv[k] = k_zero;
    for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
      s_12.comp1 = R_AB[i4].comp1 - atoms->cell_vector[atm0].comp1;
      s_12.comp2 = R_AB[i4].comp2 - atoms->cell_vector[atm0].comp2;
      s_12.comp3 = R_AB[i4].comp3 - atoms->cell_vector[atm0].comp3;
      Rsqrd = double_vec_dot(&s_12, &s_12);
      f000m(&em[0][0], Rsqrd / pab_inv[i4], k_one / pab_inv[i4], mm);
      non_recursive_ftuvn(mm, 0, f, em, &s_12);
      temp = sab[i4] * atoms->atomic_number[atm0];
      for (t = 0; t <= mm; t++) {
        for (u = 0; u <= mm; u++) {
          for (v = 0; v <= mm; v++) {
            if (t + u + v > mm) break;
            fgtuv_temp = f[t][u][v][0] * temp;
            fgtuv[t * mm1 + u * mm2 + v * mm3  + i4] += fgtuv_temp;
            fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
            fgtuv_max = k_one;
           }
          }
         }
        }

  for (i = 0; i < sheli; i++) {
    for (j = 0; j < shelj; j++) {
      count = (bfposi + i) * nd2 + bfposj + j;
      tmax = shells->tuv[imax][i][0] + shells->tuv[jmax][j][0];
      umax = shells->tuv[imax][i][1] + shells->tuv[jmax][j][1];
      vmax = shells->tuv[imax][i][2] + shells->tuv[jmax][j][2];
      n1 = shells->tuv[imax][i][0];
      n2 = shells->tuv[imax][i][1];
      n3 = shells->tuv[imax][i][2];
      n4 = shells->tuv[jmax][j][0];
      n5 = shells->tuv[jmax][j][1];
      n6 = shells->tuv[jmax][j][2];
      for (t = 0; t <= tmax; t++) {
        for (u = 0; u <= umax; u++) {
          for (v = 0; v <= vmax; v++) {
            for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
               ElecNuc[count] -= E1x[n1 * off1 + n4 * off2 + t * off3 + i4] * E1y[n2 * off1 + n5 * off2 + u * off3 + i4] * \
                                 E1z[n3 * off1 + n6 * off2 + v * off3 + i4] * fgtuv[t * mm1 + u * mm2 + v * mm3  + i4];
              }
             }
            }
           } // end t u v loop
          } // close loop over j
         } // close loop over i
        } // end loop over atm0

        free(fgtuv);

  break;

  case 'C':

  gamma_0 = k_one / G->gamma_0_inv;
  root_gamma_0 = sqrt(gamma_0);

  imax = shells->imax_sh[index_i];
  jmax = shells->imax_sh[index_j];
  off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2 = (imax + im + jmax + jm + 1) * off3;
  off1 = (jmax + jm + 1) * off2;
  sheli = shells->type1_sh[index_i];
  shelj = shells->type1_sh[index_j];
  for (i = 0; i < sheli; i++) {
    for (j = 0; j < shelj; j++) {
      count = (bfposi + i) * nd2 + bfposj + j;
      tmax = shells->tuv[imax][i][0] + shells->tuv[jmax][j][0];
      umax = shells->tuv[imax][i][1] + shells->tuv[jmax][j][1];
      vmax = shells->tuv[imax][i][2] + shells->tuv[jmax][j][2];
      n1 = shells->tuv[imax][i][0];
      n2 = shells->tuv[imax][i][1];
      n3 = shells->tuv[imax][i][2];
      n4 = shells->tuv[jmax][j][0];
      n5 = shells->tuv[jmax][j][1];
      n6 = shells->tuv[jmax][j][2];
      r3 = k_zero;
      for (atm0 = 0; atm0 < atoms->number_of_atoms_in_unit_cell; atm0++) {
        for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
          s_12.comp1 = R_AB[i4].comp1 - atoms->cell_vector[atm0].comp1;
          s_12.comp2 = R_AB[i4].comp2 - atoms->cell_vector[atm0].comp2;
          s_12.comp3 = R_AB[i4].comp3 - atoms->cell_vector[atm0].comp3;
          Rsqrd = double_vec_dot(&s_12, &s_12);
          r4 = atoms->atomic_number[atm0] * sab[i4];
          for (t = 0; t <= tmax; t++) {
            for (u = 0; u <= umax; u++) {
              for (v = 0; v <= vmax; v++) {
                nn = t + u + v;
                count1 = 1;
                r1 = r4 * E1x[n1 * off1 + n4 * off2 + t * off3 + i4] * E1y[n2 * off1 + n5 * off2 + u * off3 + i4] * \
                E1z[n3 * off1 + n6 * off2 + v * off3+i4];
                for (index_S = 1; index_S < G->number_of_shells; index_S++) {
                  r2 = r1 * G->EXPFAC[index_S]; 
                  for (index_G = 0; index_G < G->num[index_S]; index_G++) {
                    GdotR = double_vec_dot(&G->vec_b2[count1], &s_12);
                    r3 += r2 * G->x[count1 * 9 + t] * G->y[count1 * 9 + u] * G->z[count1 * 9 + v] * cosfactor(nn, GdotR);
                    count1++;
                   } // end loop over index_G
                  } // end loop over index_S
                 }
                }
               } // end t u v loop
              } // end loop over i4
             } // end loop over atm0
            ElecNuc[count] -= r3;
           } // close loop over j
          } // close loop over i


  off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2 = (imax + im + jmax + jm + 1) * off3;
  off1 = (jmax + jm + 1) * off2;
  sheli = shells->type1_sh[index_i];
  shelj = shells->type1_sh[index_j];
  for (i = 0; i < sheli; i++) {
    for (j = 0; j < shelj; j++) {
      count = (bfposi + i) * nd2 + bfposj + j;
      tmax = shells->tuv[imax][i][0] + shells->tuv[jmax][j][0];
      umax = shells->tuv[imax][i][1] + shells->tuv[jmax][j][1];
      vmax = shells->tuv[imax][i][2] + shells->tuv[jmax][j][2];
      n1 = shells->tuv[imax][i][0];
      n2 = shells->tuv[imax][i][1];
      n3 = shells->tuv[imax][i][2];
      n4 = shells->tuv[jmax][j][0];
      n5 = shells->tuv[jmax][j][1];
      n6 = shells->tuv[jmax][j][2];
      for (atm0 = 0; atm0 < atoms->number_of_atoms_in_unit_cell; atm0++) {
        for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
          r4 = sab[i4] * (double)-atoms->atomic_number[atm0];
          r1 = r4 * -pi * (G->gamma_0_inv - pab_inv[i4]) / crystal->primitive_cell_volume;
          s_12.comp1 = R_AB[i4].comp1 - atoms->cell_vector[atm0].comp1;
          s_12.comp2 = R_AB[i4].comp2 - atoms->cell_vector[atm0].comp2;
          s_12.comp3 = R_AB[i4].comp3 - atoms->cell_vector[atm0].comp3;
          map_to_wigner(crystal,&s_12, &t_12, &Rvec_tmp); // map s_12 back to primitive cell to accelerate convergence
          count1 = 0;
          for (index_S = 0; index_S < R->number_of_ewald_shells; index_S++) {
            double shell_sum = k_zero;
            for (index_R = 0; index_R < R->num[index_S]; index_R++) {
              r_12.comp1 = t_12.comp1 + R->vec_ai[count1].comp1;
              r_12.comp2 = t_12.comp2 + R->vec_ai[count1].comp2;
              r_12.comp3 = t_12.comp3 + R->vec_ai[count1].comp3;
              Rsqrd = double_vec_dot(&r_12, &r_12);
              for (t = 0; t <= tmax; t++) {
                for (u = 0; u <= umax; u++) {
                  for (v = 0; v <= vmax; v++) {
                    nn = t + u + v;
                    f000m(fn, Rsqrd / pab_inv[i4], k_one / pab_inv[i4], nn);
                    f000m(gn, gamma_0 * Rsqrd, gamma_0, nn);
		    shell_sum += E1x[n1 * off1 + n4 * off2 + t * off3 + i4] * E1y[n2 * off1 + n5 * off2 + u * off3 + i4] * \
                    E1z[n3 * off1 + n6 * off2 + v * off3 + i4] * (ftuvn(t, u, v, 0, fn, r_12) - ftuvn(t, u, v, 0, gn, r_12)) * r4;
		    //ElecNuc[count] += E1x[n1 * off1 + n4 * off2 + t * off3 + i4] * E1y[n2 * off1 + n5 * off2 + u * off3 + i4] * \
                    E1z[n3 * off1 + n6 * off2 + v * off3 + i4] * (ftuvn(t, u, v, 0, fn, r_12) - ftuvn(t, u, v, 0, gn, r_12)) * r4;
                   }
                  }
                 } // end t u v loop
                count1++;
               } // end index_R loop
               ElecNuc[count] += shell_sum;
               if (fabs(shell_sum) < 1e-13 && sab[i4] < 1e-13) break;
              } // end index_S loop
          ElecNuc[count] += E1x[n1 * off1 + n4 * off2 + i4] * E1y[n2 * off1 + n5 * off2 + i4] * E1z[n3 * off1 + n6 * off2 + i4] * r1;
            } // end loop over i4
           } // end loop over atm0
          } // close loop over j
         } // close loop over i

  break;

  case 'S':

  gamma_0 = k_one / G->gamma_0_inv;
  root_gamma_0 = sqrt(gamma_0);
  imax = shells->imax_sh[index_i];
  jmax = shells->imax_sh[index_j];
  off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2 = (imax + im + jmax + jm + 1) * off3;
  off1 = (jmax + jm + 1) * off2;
  sheli = shells->type1_sh[index_i];
  shelj = shells->type1_sh[index_j];
  for (i = 0; i < sheli; i++) {
    for (j = 0; j < shelj; j++) {
      count = (bfposi + i) * nd2 + bfposj + j;
      tmax = shells->tuv[imax][i][0] + shells->tuv[jmax][j][0];
      umax = shells->tuv[imax][i][1] + shells->tuv[jmax][j][1];
      vmax = shells->tuv[imax][i][2] + shells->tuv[jmax][j][2];
      n1 = shells->tuv[imax][i][0];
      n2 = shells->tuv[imax][i][1];
      n3 = shells->tuv[imax][i][2];
      n4 = shells->tuv[jmax][j][0];
      n5 = shells->tuv[jmax][j][1];
      n6 = shells->tuv[jmax][j][2];
      r3 = k_zero;
      for (atm0 = 0; atm0 < atoms->number_of_atoms_in_unit_cell; atm0++) {
        for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
          s_12.comp1 = R_AB[i4].comp1 - atoms->cell_vector[atm0].comp1;
          s_12.comp2 = R_AB[i4].comp2 - atoms->cell_vector[atm0].comp2;
          s_12.comp3 = R_AB[i4].comp3 - atoms->cell_vector[atm0].comp3;
          Rzsqrd = s_12.comp3 * s_12.comp3;
          Rsqrd = double_vec_dot(&s_12, &s_12);
          r4 = -atoms->atomic_number[atm0] * sab[i4];
          //r5 = r4 * -two * pi / crystal->primitive_cell_volume * \
         (exp(-gamma_0 * Rzsqrd) / root_gamma_0 / rtpi + s_12.comp3 * erf(root_gamma_0 * s_12.comp3));
          x = root_gamma_0 * s_12.comp3;
          erf_exp_derivative(D_erf, D_exp, mm, root_gamma_0, x, G);
          sum5[0] = D_exp[0] / rtpi / root_gamma_0 + s_12.comp3 * D_erf[0];
          for (n = 1; n <= mm; n++) {
          sum5[n] = D_exp[n] / rtpi / root_gamma_0 + s_12.comp3 * D_erf[n] + (double) n * D_erf[n - 1];
         }
          r5 = k_zero;
          for (v = 0; v <= vmax; v++) {
            r5 += -two * pi / crystal->primitive_cell_volume * sum5[v] * r4 * \
            E1x[n1 * off1 + n4 * off2 + i4] * E1y[n2 * off1 + n5 * off2 + i4] * E1z[n3 * off1 + n6 * off2 + v * off3 + i4];
           }
	    ElecNuc[count] += r5;
          for (t = 0; t <= tmax; t++) {
            for (u = 0; u <= umax; u++) {
              for (v = 0; v <= vmax; v++) {
                nn = t + u;
                count1 = 1;
              r1 = r4 * E1x[n1 * off1 + n4 * off2 + t * off3 + i4] * E1y[n2 * off1 + n5 * off2 + u * off3 + i4] * \
              E1z[n3 * off1 + n6 * off2 + v * off3 + i4];
                for (index_S = 1; index_S < G->number_of_shells; index_S++) {
                  x1 = G->shell_mag[9 * index_S + 1] / two / root_gamma_0 + root_gamma_0 * s_12.comp3;
                  erfc_derivative(derivative_1, mm, root_gamma_0, x1, G);
                  x2 = G->shell_mag[9 * index_S + 1] / two / root_gamma_0 - root_gamma_0 * s_12.comp3;
                  erfc_derivative(derivative_2, mm, -root_gamma_0, x2, G);
                  expfac = exp(G->shell_mag[9 * index_S + 1] * s_12.comp3);
                  for (n = 0; n <= mm; n++) {
                    sign = k_one;
                    sum2[n] = k_zero;
                    for (k = 0; k <= n; k++) {
                      sum2[n] += G->A->a[n][k] * G->shell_mag[9 * index_S + k] * (expfac * derivative_1[n - k] + sign / expfac * \
                      derivative_2[n - k]) / G->shell_mag[9 * index_S + 1];
                      sign *= -k_one;
                     }
                    }
                  //fac2 = fac1 * G->EXPFAC[index_S];
                  //double shell_sum = k_zero;
                  r2 = r1 *  two * pi / crystal->primitive_cell_volume * sum2[v];
                  for (index_G = 0; index_G < G->num[index_S]; index_G++) {
                    GdotR = double_vec_dot(&G->vec_b2[count1], &s_12);
                    r3 += r2 * G->x[count1 * 9 + t] * G->y[count1 * 9 + u] * cosfactor(nn, GdotR);
                    count1++;
                   } // end loop over index_G
                  //r3 += shell_sum;
                  //if (fabs(shell_sum) < 1e-24 && sab[i4] < 1e-24) continue;
                  } // end loop over index_S
                 }
                }
               } // end t u v loop
              } // end loop over i4
             } // end loop over atm0
            ElecNuc[count] += r3;
           } // close loop over j
          } // close loop over i


  off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2 = (imax + im + jmax + jm + 1) * off3;
  off1 = (jmax + jm + 1) * off2;
  sheli = shells->type1_sh[index_i];
  shelj = shells->type1_sh[index_j];
  for (i = 0; i < sheli; i++) {
    for (j = 0; j < shelj; j++) {
      count = (bfposi + i) * nd2 + bfposj + j;
      tmax = shells->tuv[imax][i][0] + shells->tuv[jmax][j][0];
      umax = shells->tuv[imax][i][1] + shells->tuv[jmax][j][1];
      vmax = shells->tuv[imax][i][2] + shells->tuv[jmax][j][2];
      n1 = shells->tuv[imax][i][0];
      n2 = shells->tuv[imax][i][1];
      n3 = shells->tuv[imax][i][2];
      n4 = shells->tuv[jmax][j][0];
      n5 = shells->tuv[jmax][j][1];
      n6 = shells->tuv[jmax][j][2];
      for (atm0 = 0; atm0 < atoms->number_of_atoms_in_unit_cell; atm0++) {
        for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
          r4 = sab[i4] * (double)-atoms->atomic_number[atm0];
          s_12.comp1 = R_AB[i4].comp1 - atoms->cell_vector[atm0].comp1;
          s_12.comp2 = R_AB[i4].comp2 - atoms->cell_vector[atm0].comp2;
          s_12.comp3 = R_AB[i4].comp3 - atoms->cell_vector[atm0].comp3;
          map_to_wigner(crystal,&s_12, &t_12, &Rvec_tmp); // map s_12 back to primitive cell to accelerate convergence
          count1 = 0;
          for (index_S = 0; index_S < R->number_of_ewald_shells; index_S++) {
            double shell_sum = k_zero;
            for (index_R = 0; index_R < R->num[index_S]; index_R++) {
              r_12.comp1 = t_12.comp1 + R->vec_ai[count1].comp1;
              r_12.comp2 = t_12.comp2 + R->vec_ai[count1].comp2;
              r_12.comp3 = t_12.comp3 + R->vec_ai[count1].comp3;
              Rsqrd = double_vec_dot(&r_12, &r_12);
              for (t = 0; t <= tmax; t++) {
                for (u = 0; u <= umax; u++) {
                  for (v = 0; v <= vmax; v++) {
                    nn = t + u + v;
                    f000m(fn, Rsqrd / pab_inv[i4], k_one / pab_inv[i4], nn);
                    f000m(gn, gamma_0 * Rsqrd, gamma_0, nn);
		    shell_sum += E1x[n1 * off1 + n4 * off2 + t * off3 + i4] * E1y[n2 * off1 + n5 * off2 + u * off3 + i4] * \
                    E1z[n3 * off1 + n6 * off2 + v * off3 + i4] * (ftuvn(t, u, v, 0, fn, r_12) - ftuvn(t, u, v, 0, gn, r_12)) * r4;
                    //fprintf(file.out,"shellsum %3d %10.3e\n",index_S,r4*(ftuvn(t, u, v, 0, fn, r_12) - ftuvn(t, u, v, 0, gn, r_12)));
		    //ElecNuc[count] += E1x[n1 * off1 + n4 * off2 + t * off3 + i4] * E1y[n2 * off1 + n5 * off2 + u * off3 + i4] * \
                    E1z[n3 * off1 + n6 * off2 + v * off3 + i4] * (ftuvn(t, u, v, 0, fn, r_12) - ftuvn(t, u, v, 0, gn, r_12)) * r4;
                   }
                  }
                 } // end t u v loop
                count1++;
               } // end index_R loop
               ElecNuc[count] += shell_sum;
               //fprintf(file.out," %20.14lf\n", shell_sum);
               //if (fabs(shell_sum) < 1e-13 && sab[i4] < 1e-13) break;
              } // end index_S loop
	     //ElecNuc[count] += E1x[n1 * off1 + n4 * off2 + i4] * E1y[n2 * off1 + n5 * off2 + i4] * E1z[n3 * off1 + n6 * off2 + i4] * r1;
             //fprintf(file.out,"Recip real other all %10.3e  %10.3e %10.3e %16.8e\n",r3,ElecNuc[count]-r3,r1,ElecNuc[count]);
            } // end loop over i4
           } // end loop over atm0
          } // close loop over j
         } // close loop over i

  break;

  case 'P':

  if (job->taskid == 0)

  fprintf(file.out,"Electron-nuclear attraction for 1-D not coded\n");

  break;

      } // close switch

}

void Momentum(double *Momentum, int index_i, int index_j, int bfposi, int bfposj, int gausposi, int gausposj, int nd1, int nd2, int im, int jm, double *E1x, double *E1y, double *E1z, double *sab, SHELL *shells, GAUSSIAN *gaussians, JOB_PARAM *job, FILES file)

{

  int i, j, i1, i4, j4;
  int n1, n2, n3, n4, n5, n6;
  int sheli, shelj;
  int off1, off2, off3;
  int imax, jmax;
  int countx, county, countz, nd12;

  imax = shells->imax_sh[index_i];
  jmax = shells->imax_sh[index_j];
  off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2 = (imax + im + jmax + jm + 1) * off3;
  off1 = (jmax + jm + 1) * off2;
  sheli = shells->type1_sh[index_i];
  shelj = shells->type1_sh[index_j];
  nd12 = nd1 * nd2;
  for (i = 0; i < sheli; i++) {
    for (j = 0; j < shelj; j++) {
      countx = (bfposi + i) * nd2 + bfposj + j;
      county = nd12 + countx;
      countz = nd12 + county;
      n1 = shells->tuv[imax][i][0];
      n2 = shells->tuv[imax][i][1];
      n3 = shells->tuv[imax][i][2];
      n4 = shells->tuv[jmax][j][0];
      n5 = shells->tuv[jmax][j][1];
      n6 = shells->tuv[jmax][j][2];
      i1 = 0;
      for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
      for (j4 = 0; j4 < shells->ng_sh[index_j]; j4++) {
        if (n4 > 0)
        Momentum[countx] += double(n4) * E1x[n1 * off1 + (n4 - 1) * off2 + i1] * E1y[n2 * off1 + n5 * off2 + i1] * E1z[n3 * off1 + n6 * off2 + i1] * sab[i1];
        Momentum[countx] -= two * gaussians->expo_sh[gausposj + j4] * \
                            E1x[n1 * off1 + (n4 + 1) * off2 + i1] * E1y[n2 * off1 + n5 * off2 + i1] * E1z[n3 * off1 + n6 * off2 + i1] * sab[i1];
        if (n5 > 0)
        Momentum[county] += double(n5) * E1x[n1 * off1 + n4 * off2 + i1] * E1y[n2 * off1 + (n5 - 1) * off2 + i1] * E1z[n3 * off1 + n6 * off2 + i1] * sab[i1];
        Momentum[county] -= two * gaussians->expo_sh[gausposj + j4] * \
                            E1x[n1 * off1 + n4 * off2 + i1] * E1y[n2 * off1 + (n5 + 1) * off2 + i1] * E1z[n3 * off1 + n6 * off2 + i1] * sab[i1];
        if (n6 > 0)
        Momentum[countz] += double(n6) * E1x[n1 * off1 + n4 * off2 + i1] * E1y[n2 * off1 + n5 * off2 + i1] * E1z[n3 * off1 + (n6 - 1) * off2 + i1] * sab[i1];
        Momentum[countz] -= two * gaussians->expo_sh[gausposj + j4] * \
                            E1x[n1 * off1 + n4 * off2 + i1] * E1y[n2 * off1 + n5 * off2 + i1] * E1z[n3 * off1 + (n6 + 1) * off2 + i1] * sab[i1];
        i1++;
       } // close loop over j4
      } // close loop over i4
     } // close loop over i
    } // close loop over i

}

void Overlap(double *Overlap, int index_i, int index_j, int bfposi, int bfposj, int nd2,  int im, int jm, double *E1x, double *E1y, double *E1z, double *sab, SHELL *shells, JOB_PARAM *job, FILES file)

{

  int i, j, i1, i4, j4;
  int n1, n2, n3, n4, n5, n6;
  int tmax, umax, vmax;
  int sheli, shelj;
  int off1, off2, off3;
  int imax, jmax;
  int count;

  imax = shells->imax_sh[index_i];
  jmax = shells->imax_sh[index_j];
  off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2 = (imax + im + jmax + jm + 1) * off3;
  off1 = (jmax + jm + 1) * off2;
  sheli = shells->type1_sh[index_i];
  shelj = shells->type1_sh[index_j];
  for (i = 0; i < sheli; i++) {
    for (j = 0; j < shelj; j++) {
      count = (bfposi + i) * nd2 + bfposj + j;
      tmax = shells->tuv[imax][i][0] + shells->tuv[jmax][j][0];
      umax = shells->tuv[imax][i][1] + shells->tuv[jmax][j][1];
      vmax = shells->tuv[imax][i][2] + shells->tuv[jmax][j][2];
      n1 = shells->tuv[imax][i][0];
      n2 = shells->tuv[imax][i][1];
      n3 = shells->tuv[imax][i][2];
      n4 = shells->tuv[jmax][j][0];
      n5 = shells->tuv[jmax][j][1];
      n6 = shells->tuv[jmax][j][2];
      for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
      Overlap[count] += E1x[n1 * off1 + n4 * off2 + i4] * E1y[n2 * off1 + n5 * off2 + i4] * E1z[n3 * off1 + n6 * off2 + i4] * sab[i4];
     } // close loop over i4
    } // close loop over j
   } // close loop over i

}

void Dipole(double *Dipole, int index_i, int index_j, int bfposi, int bfposj, int gausposi, int gausposj, int nd1, int nd2, int ip, int im, int jm, double *E1x, double *E1y, double *E1z, double *sab, VECTOR_DOUBLE *R_AB, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, JOB_PARAM *job, FILES file)

{

  int i, j, i1, i4, j4;
  int n1, n2, n3, n4, n5, n6;
  int sheli, shelj;
  int off1, off2, off3;
  int imax, jmax;
  int countx, county, countz, nd12;
  double PAx, PAy, PAz;

  imax = shells->imax_sh[index_i];
  jmax = shells->imax_sh[index_j];
  off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2 = (imax + im + jmax + jm + 1) * off3;
  off1 = (jmax + jm + 1) * off2;
  sheli = shells->type1_sh[index_i];
  shelj = shells->type1_sh[index_j];
  nd12 = nd1 * nd2;
  for (i = 0; i < sheli; i++) {
    for (j = 0; j < shelj; j++) {
      countx = (bfposi + i) * nd2 + bfposj + j;
      county = nd12 + countx;
      countz = nd12 + county;
      n1 = shells->tuv[imax][i][0];
      n2 = shells->tuv[imax][i][1];
      n3 = shells->tuv[imax][i][2];
      n4 = shells->tuv[jmax][j][0];
      n5 = shells->tuv[jmax][j][1];
      n6 = shells->tuv[jmax][j][2];
      i1 = 0;
      for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
      for (j4 = 0; j4 < shells->ng_sh[index_j]; j4++) {
        Dipole[countx] += (E1x[n1 * off1 + n4 * off2 + off3 + i1] + R_AB[i1].comp1 * E1x[n1 * off1 + n4 * off2 + i1]) * \
                           E1y[n2 * off1 + n5 * off2 + i1] * E1z[n3 * off1 + n6 * off2 + i1] * sab[i1];
        Dipole[county] += (E1y[n2 * off1 + n5 * off2 + off3 + i1] + R_AB[i1].comp2 * E1y[n2 * off1 + n5 * off2 + i1]) * \
                           E1x[n1 * off1 + n4 * off2 + i1] * E1z[n3 * off1 + n6 * off2 + i1] * sab[i1];
        Dipole[countz] += (E1z[n3 * off1 + n6 * off2 + off3 + i1] + R_AB[i1].comp3 * E1z[n3 * off1 + n6 * off2 + i1]) * \
                           E1x[n1 * off1 + n4 * off2 + i1] * E1y[n2 * off1 + n5 * off2 + i1] * sab[i1];
        i1++;
       } // close loop over j4
      } // close loop over i4
     } // close loop over i
    } // close loop over i

}

void two_centre_coulomb(double *Coulomb, int index_i, int index_j, int bfposi1, int bfposj1, int gausposi, int gausposj, int nd2,  int im, int jm, double *E1x, double *E1y, double *E1z, VECTOR_DOUBLE *R_AB_1e, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, JOB_PARAM *job, FILES file)

{

int t, u, v, i4, j4;
int imax, jmax;
int count;
int mm, mm0, mm1, mm2, mm3;
double Rsqrd;
double em[1][55];
double fgtuv_max, fgtuv_temp;
double *fgtuv;

  double prefac[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
  double ab[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
  count = 0;
  for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
  for (j4 = 0; j4 < shells->ng_sh[index_j]; j4++) {
  prefac[count] = gaussians->c_sh[gausposi + i4] * gaussians->c_sh[gausposj + j4] * \
  pi32 * pi32 / gaussians->expo_sh[gausposi + i4] / 
  gaussians->expo_sh[gausposj + j4] / sqrt(gaussians->expo_sh[gausposi + i4] * gaussians->expo_sh[gausposj + j4]); 
  ab[count] = k_one / gaussians->expo_sh[gausposi + i4] + k_one /  gaussians->expo_sh[gausposj + j4];
  count++; 
 }
}

  imax = shells->imax_sh[index_i];
  jmax = shells->imax_sh[index_j];
  mm  = imax + jmax;
  mm3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  mm2 = (mm + 1) * mm3;
  mm1 = (mm + 1) * mm2;
  mm0 = (mm + 1) * mm1;

  fgtuv = (double *) calloc(mm0, sizeof(double));
  if (fgtuv == NULL) {
  if (job->taskid == 0)
  fprintf(stderr, "ERROR: not enough memory for double fgtuv! \n");
  MPI_Finalize();
  exit(1);
 }

  Rsqrd = double_vec_dot(R_AB_1e, R_AB_1e);
  for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
  f000m(&em[0][0], Rsqrd / ab[i4], k_one / ab[i4], mm);
  for (t = 0; t <= mm; t++) {
    for (u = 0; u <= mm; u++) {
      for (v = 0; v <= mm; v++) {
        if (t + u + v > mm) break;
          fgtuv_temp = ftuvn(t,u,v,0,&em[0][0],R_AB_1e[0]) * prefac[i4];
          fgtuv[t * mm1 + u * mm2 + v * mm3  + i4] += fgtuv_temp;
          fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
          fgtuv_max = k_one;
         }
        }
       }
      }

  double C1_max, C2_max;
  double C1x[shells->ng_sh[index_i] * (imax+1) * (imax+1)];
  double C1y[shells->ng_sh[index_i] * (imax+1) * (imax+1)];
  double C1z[shells->ng_sh[index_i] * (imax+1) * (imax+1)];
  double C2x[shells->ng_sh[index_j] * (jmax+1) * (jmax+1)];
  double C2y[shells->ng_sh[index_j] * (jmax+1) * (jmax+1)];
  double C2z[shells->ng_sh[index_j] * (jmax+1) * (jmax+1)];
  E_coefficients_1c(index_i,gausposi,C1x,C1y,C1z,&C1_max,shells,gaussians,job,file);
  E_coefficients_1c(index_j,gausposj,C2x,C2y,C2z,&C2_max,shells,gaussians,job,file);
  mcmurchie_davidson_ij(Coulomb,index_i,index_j,bfposi1,bfposj1,nd2,C1x,C1y,C1z,C2x,C2y,C2z,fgtuv,shells,job,file);
  free(fgtuv);

}

/*
void fock_element_elecnuc(double *elecnuc_shar, PAIR_TRAN *pair_p, REAL_LATTICE *R, RECIPROCAL_LATTICE *G, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  int ip, jp, gi, gj, p;
  int i, j, i4, j4, k, l, mm, nn, atm0, count;
  int index_i, index_j;
  int dim, dim1, dim2, dim3, dim4, dim5, dim6, dim7, dim8, dim9, atm1, atm2;
  int bfposi, bfposj, bfposi1, bfposj1;
  int imax, jmax, im, jm;
  int gausposi, gausposj, shelposi, shelposj;
  int sheli, shelj, sheli1, shelj1;
  int nd1, nd2, nd3, nd4, nd12, nd34;
  int nsheli, nshelj;
  double R_AB_1esqrd;
  double E1_max;
  double time1, time2;
  VECTOR_DOUBLE R_AB_1e;

  // Function[0] Fock Operator
  // Function[1] Kinetic Energy Operator
  // Function[2] Electron-Nuclear Potential Energy = 1 (finite systems) = 2 (extended systems) = 3 (atoms)
  // Function[3] Momentum Operator Matrix Elements
  // Function[4] Overlap
  // Function[5] 
  // Function[6] Dipole Operator Matrix Elements
  // Function[7] Coulomb Operator Matrix Elements

  time1 = MPI_Wtime();

    dim = 0;
    for (p = 0; p < pair_p->nump; p++) {
    ip = pair_p->cell1[pair_p->posn[p]];
    jp = pair_p->cell2[pair_p->posn[p]];
    gi = 0;
    gj = pair_p->latt2[pair_p->posn[p]];
    nd1 = atoms->bfnnumb[ip];
    nd2 = atoms->bfnnumb[jp];
    nd3 = atoms->bfnnumb_sh[ip];
    nd4 = atoms->bfnnumb_sh[jp];
    nd12 = nd1 * nd2;
    nd34 = nd3 * nd4;

VECTOR_DOUBLE Rvec_tmp;
double Rsqrd;
Rvec_tmp.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1 - R->vec_ai[gj].comp1;
Rvec_tmp.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2 - R->vec_ai[gj].comp2;
Rvec_tmp.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3 - R->vec_ai[gj].comp3;
Rsqrd = double_vec_dot(&Rvec_tmp,&Rvec_tmp);
if (Rsqrd > job->itol1 * job->itol1) { dim += nd34;  continue; }

    //double kinetic_cart[nd12]; 
    double elecnuc_cart[nd12]; 
    //double momentum_cart[3 * nd12]; 
    //double overlap_cart[nd12]; 
    //double dipole_cart[3 * nd12]; 
    //double coulomb_cart[nd12]; 
    // order matters for these loops to set im and jm
    //im = 0;
    //jm = 0;
    //if (Function[0] || Function[2]) {
    for (i = 0; i < nd12; i++) elecnuc_cart[i] = k_zero;

    //fprintf(file.out,"ip jp gj numb %3d %3d %3d %3d %3d\n",p,ip,jp,gj,pair_p->numb[p]) ;
    shelposi = atoms->shelposn_sh[ip];
    gausposi = atoms->gausposn_sh[ip];
    bfposi  = 0;
    bfposi1 = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
      sheli  = shells->type_sh[index_i];
      sheli1 = shells->type1_sh[index_i];
      imax   = shells->imax_sh[index_i];
      shelposj = atoms->shelposn_sh[jp];
      gausposj = atoms->gausposn_sh[jp];
      bfposj  = 0;
      bfposj1 = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
        R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1 - R->vec_ai[gj].comp1;
        R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2 - R->vec_ai[gj].comp2;
        R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3 - R->vec_ai[gj].comp3;
        R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);
        shelj  = shells->type_sh[index_j];
        shelj1 = shells->type1_sh[index_j];
        jmax   = shells->imax_sh[index_j];
        double sab[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double E1x[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+im+jm+1) * (imax+im+1) * (jmax+jm+1)];
        double E1y[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+im+jmax+jm+1) * (imax+im+1) * (jmax+jm+1)];
        double E1z[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+im+jmax+jm+1) * (imax+im+1) * (jmax+jm+1)];
        double pab_inv[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        VECTOR_DOUBLE R_AB[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        E_coefficients1(ip,jp,gi,gj,index_i,index_j,gausposi,gausposj,im,jm,E1x,E1y,E1z,&E1_max,sab,pab_inv,&R_AB_1esqrd,\
        R_AB,R,atoms,shells,gaussians,job,file);
        ElecNuc(elecnuc_cart, index_i, index_j, bfposi1, bfposj1, nd2,  im, jm, E1x, E1y, E1z, R_AB, pab_inv, sab, R, G, \
        crystal, atoms,shells,job,file);
        two_center_cartesian_to_sh_shell(elecnuc_cart,&elecnuc_shar[dim],index_i,index_j,bfposi,bfposj,bfposi1,\
        bfposj1,nd2,nd4,shells,job,file);
        bfposj   += shelj;
        bfposj1  += shelj1;
        gausposj += shells->ng_sh[index_j];
       } // close loop over index_j
      bfposi   += sheli;
      bfposi1  += sheli1;
      gausposi += shells->ng_sh[index_i];
     } // close loop over index_i

 if (job->verbosity > 1) {
  count = 0;
   fprintf(file.out,"elecnuc\n");
      fprintf(file.out,"pair %3d ip %3d jp %3d gj %3d\n",p,ip,jp,gj);
      for(i=0;i<atoms->bfnnumb_sh[ip];i++) {
        for(j=0;j<atoms->bfnnumb_sh[jp];j++) {
          fprintf(file.out,"%6.2lf",elecnuc_shar[dim + count]);
           count++;
          }
         fprintf(file.out,"\n");
        }
        fprintf(file.out,"\n");
       }

    dim += nd34;

  } // close loop on p


   time2 = MPI_Wtime();

}
*/

/*
void two_centre_coulomb1(INT_1E *one_ints, PAIR_TRAN *pair_p, REAL_LATTICE *R, RECIPROCAL_LATTICE *G, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  int i, j, m, n, p, t, u, v;
  int ip, jp, gi, gj;
  int i4, j4;
  int n1, n2, n3, n4, n5, n6;
  int index_i, index_j;
  int bfposi, bfposj, bfposi1, bfposj1;
  int imax, jmax, im, jm;
  int gausposi, gausposj, shelposi, shelposj;
  int sheli, shelj, sheli1, shelj1;
  int nd1, nd2, nd3, nd4;
  int nsheli, nshelj;
  int tmax, umax, vmax;
  int off1, off2, off3, off4;
  int count;
  int mm, mm0, mm1, mm2, mm3, mm4;
  double Rsqrd, p_inv_ex, fn[55];
  double gamma_0, root_gamma_0, gamma_1;
  double en[R->last_ewald_vector][55], em[R->last_ewald_vector][55];
  double f[13][13][13][13];
  double fgtuv_max, fgtuv_temp;
  double *fgtuv;
  VECTOR_DOUBLE R_AB_1e;
  double R_AB_1esqrd, C1_max, C2_max;

int begin_p1 = 0, end_p1 = 1;
int dim, dimg, dim3=0;

int i1, index_R, index_S, index_G, tuv;
double fac2, dot_product;
double shell_sum, pi_vol = pi / crystal->primitive_cell_volume, gamma_1_inv = G->gamma_0_inv;
gamma_1 = k_one / G->gamma_0_inv;
VECTOR_DOUBLE r_12, t_12, Rvec_tmp;

  array_dimensions(&dim, &dimg, pair_p, atoms, job, file); // don't change

  for (i = 0; i < dim; i++) one_ints->Coulomb[i] = k_zero;

  switch (crystal->type[0]) {

  case 'M':

    dim = 0;

  for (p = 0; p < pair_p->nump; p++) {

    ip = pair_p->cell1[pair_p->posn[p]];
    jp = pair_p->cell2[pair_p->posn[p]];
    gi = 0;
    gj = pair_p->latt2[pair_p->posn[p]];
    nd1 = atoms->bfnnumb[ip];
    nd2 = atoms->bfnnumb[jp];
    nd3 = atoms->bfnnumb_sh[ip];
    nd4 = atoms->bfnnumb_sh[jp];
    double coulomb_cart[nd1 * nd2];
    for (i = 0; i < nd1 * nd2; i++) coulomb_cart[i] = k_zero;
    shelposi = atoms->shelposn_sh[ip];
    gausposi = atoms->gausposn_sh[ip];
    bfposi  = 0;
    bfposi1 = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
      sheli  = shells->shar[index_i];
      sheli1 = shells->cart[index_i];
      imax   = shells->imax_sh[index_i];
      shelposj = atoms->shelposn_sh[jp];
      gausposj = atoms->gausposn_sh[jp];
      bfposj  = 0;
      bfposj1 = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
        R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1 - R->vec_ai[gj].comp1;
        R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2 - R->vec_ai[gj].comp2;
        R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3 - R->vec_ai[gj].comp3;
        R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);
        shelj  = shells->shar[index_j];
        shelj1 = shells->cart[index_j];
        jmax   = shells->imax_sh[index_j];
        double sab_fac[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double ab_inv[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double C1x[shells->ng_sh[index_i] * (imax+1) * (imax+1)];
        double C1y[shells->ng_sh[index_i] * (imax+1) * (imax+1)];
        double C1z[shells->ng_sh[index_i] * (imax+1) * (imax+1)];
        double C2x[shells->ng_sh[index_j] * (jmax+1) * (jmax+1)];
        double C2y[shells->ng_sh[index_j] * (jmax+1) * (jmax+1)];
        double C2z[shells->ng_sh[index_j] * (jmax+1) * (jmax+1)];
        count = 0;
        for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
          for (j4 = 0; j4 < shells->ng_sh[index_j]; j4++) {
            ab_inv[count] = k_one / gaussians->expo_sh[gausposi + i4] + k_one / gaussians->expo_sh[gausposj + j4];
            sab_fac[count] = pi32 * pi32 / sqrt(gaussians->expo_sh[gausposi + i4] * gaussians->expo_sh[gausposj + j4]) / \
            gaussians->expo_sh[gausposi + i4] / gaussians->expo_sh[gausposj + j4] * gaussians->c_sh[gausposi + i4] * \
            gaussians->c_sh[gausposj + j4];
            count++;
           }
          }
        E_coefficients_1c(index_i,gausposi,C1x,C1y,C1z,&C1_max,shells,gaussians,job,file);
        E_coefficients_1c(index_j,gausposj,C2x,C2y,C2z,&C2_max,shells,gaussians,job,file);
        mm  = imax + jmax;
        mm4 = shells->ng_sh[index_j];
        mm3 = shells->ng_sh[index_i] * mm4;
        mm2 = (mm + 1) * mm3;
        mm1 = (mm + 1) * mm2;
        mm0 = (mm + 1) * mm1;
        fgtuv = (double *) calloc(mm0, sizeof(double));
        if (fgtuv == NULL) {
        if (job->taskid == 0)
        fprintf(stderr, "ERROR: not enough memory for double fgtuv! \n");
        MPI_Finalize();
        exit(1);
       }
        Rsqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);
        for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
          f000m(&em[0][0], Rsqrd / ab_inv[i4], k_one / ab_inv[i4], mm);
          //CHANGE BACK non_recursive_ftuvn(mm, 0, f, em, &R_AB_1e);
          for (t = 0; t <= mm; t++) {
            for (u = 0; u <= mm; u++) {
              for (v = 0; v <= mm; v++) {
                if (t + u + v > mm) break;
                  fgtuv_temp = ftuvn(t,u,v,0,&em[0][0],R_AB_1e) * sab_fac[i4];
                  //fgtuv_temp = f[t][u][v][0]* sab_fac[i4];
                  //printf("fgtuv2 %3d %3d %3d %10.4lf %10.4lf %10.4lf %10.4lf\n",t,u,v,fgtuv_temp,f[t][u][v][0],Rsqrd,ab_inv[i4]);
                  fgtuv[t * mm1 + u * mm2 + v * mm3  + i4] += fgtuv_temp;
                  fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
                  fgtuv_max = k_one;
                 }
                }
               }
              }
              mcmurchie_davidson_2c(coulomb_cart,index_i,index_j,bfposi1,bfposj1,nd2,C1x,C1y,C1z,C2x,C2y,C2z,fgtuv,shells,job,file);
              free(fgtuv);
              two_center_cartesian_to_sh_shell(coulomb_cart,&one_ints->Coulomb[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1, \
              bfposj1,nd2,nd4,shells,job,file);
              bfposj   += shelj;
              bfposj1  += shelj1;
              gausposj += shells->ng_sh[index_j];
             } // close loop over index_j
            bfposi   += sheli;
            bfposi1  += sheli1;
            gausposi += shells->ng_sh[index_i];
           } // close loop over index_i
            dim += nd3 * nd4;
              //for (int jj = 0; jj < nd1*nd2; jj++) fprintf(file.out,"2C CCART finite jj %3d %10.4lf\n",jj,coulomb_cart[jj]);
          } // close loop on p

   break;

  case 'C':

    dim = 0;

  for (p = 0; p < pair_p->nump; p++) {

    ip = pair_p->cell1[pair_p->posn[p]];
    jp = pair_p->cell2[pair_p->posn[p]];
    gi = 0;
    gj = pair_p->latt2[pair_p->posn[p]];
    nd1 = atoms->bfnnumb[ip];
    nd2 = atoms->bfnnumb[jp];
    nd3 = atoms->bfnnumb_sh[ip];
    nd4 = atoms->bfnnumb_sh[jp];
    double coulomb_cart[nd1 * nd2];
    for (i = 0; i < nd1 * nd2; i++) coulomb_cart[i] = k_zero;
    shelposi = atoms->shelposn_sh[ip];
    gausposi = atoms->gausposn_sh[ip];
    bfposi  = 0;
    bfposi1 = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
      sheli  = shells->shar[index_i];
      sheli1 = shells->cart[index_i];
      imax   = shells->imax_sh[index_i];
      //E_coefficients_1c(index_i,gausposi,C1x,C1y,C1z,&C1_max,shells,gaussians,job,file);
      shelposj = atoms->shelposn_sh[jp];
      gausposj = atoms->gausposn_sh[jp];
      bfposj  = 0;
      bfposj1 = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
        R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1 - R->vec_ai[gj].comp1;
        R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2 - R->vec_ai[gj].comp2;
        R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3 - R->vec_ai[gj].comp3;
        R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);
        shelj  = shells->shar[index_j];
        shelj1 = shells->cart[index_j];
        jmax   = shells->imax_sh[index_j];
        double sab_fac[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double ab_inv[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double pab[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double C1x[shells->ng_sh[index_i] * (imax+1) * (imax+1)];
        double C1y[shells->ng_sh[index_i] * (imax+1) * (imax+1)];
        double C1z[shells->ng_sh[index_i] * (imax+1) * (imax+1)];
        double C2x[shells->ng_sh[index_j] * (jmax+1) * (jmax+1)];
        double C2y[shells->ng_sh[index_j] * (jmax+1) * (jmax+1)];
        double C2z[shells->ng_sh[index_j] * (jmax+1) * (jmax+1)];
        count = 0;
        for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
          for (j4 = 0; j4 < shells->ng_sh[index_j]; j4++) {
            //pab[count] = gaussians->expo_sh[gausposi + i4] + gaussians->expo_sh[gausposj + j4];
            ab_inv[count] = k_one / gaussians->expo_sh[gausposi + i4] + k_one / gaussians->expo_sh[gausposj + j4];
            sab_fac[count] = pi32 * pi32 / sqrt(gaussians->expo_sh[gausposi + i4] * gaussians->expo_sh[gausposj + j4]) / \
            gaussians->expo_sh[gausposi + i4] / gaussians->expo_sh[gausposj + j4] * gaussians->c_sh[gausposi + i4] * \
            gaussians->c_sh[gausposj + j4];
            count++;
           }
          }
        E_coefficients_1c(index_i,gausposi,C1x,C1y,C1z,&C1_max,shells,gaussians,job,file);
        E_coefficients_1c(index_j,gausposj,C2x,C2y,C2z,&C2_max,shells,gaussians,job,file);
        mm  = imax + jmax;
        mm4 = shells->ng_sh[index_j];
        mm3 = shells->ng_sh[index_i] * mm4;
        mm2 = (mm + 1) * mm3;
        mm1 = (mm + 1) * mm2;
        mm0 = (mm + 1) * mm1;
        fgtuv = (double *) calloc(mm0, sizeof(double));
        if (fgtuv == NULL) {
        if (job->taskid == 0)
        fprintf(stderr, "ERROR: not enough memory for double fgtuv! \n");
        MPI_Finalize();
        exit(1);
       }

      fgtuv_max = k_zero;
      for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
        fgtuv[i4] -= pi_vol * sab_fac[i4] * (gamma_1_inv - ab_inv[i4]);
        //fprintf(file.out,"pi vol fac real %10.4f %10.4f %10.4f %10.4f %10.4f\n",pi_vol, sab_fac[i4], gamma_1_inv, - ab_inv[i4],\
             - pi_vol * sab_fac[i4] * (gamma_1_inv - ab_inv[i4]));
        count = 1;
        for (index_S = 1; index_S < G->number_of_shells; index_S++) {
          fac2 = sab_fac[i4] * G->EXPFAC[index_S];
          for (index_G = 0; index_G < G->num[index_S]; index_G++) {
            dot_product = G->vec_b2[count].comp1 * R_AB_1e.comp1 + G->vec_b2[count].comp2 * R_AB_1e.comp2 + \
            G->vec_b2[count].comp3 * R_AB_1e.comp3;
            for (t = 0; t <= mm; t++) {
              for (u = 0; u <= mm; u++) {
                for (v = 0; v <= mm; v++) {
                  tuv = t+ u + v;
                  if (tuv > mm) break;
                  fgtuv_temp = fac2 * G->x[count * 9 + t] * G->y[count * 9 + u] * G->z[count * 9 + v] * cosfactor(tuv, dot_product);
                  fgtuv[t * mm1 + u * mm2 + v * mm3  + i4] += fgtuv_temp;
                  fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
                  fgtuv_max = k_one;
                 }
                }
               }
              count++;
             }
            }
           }
                 //for (i=0;i<mm0;i++) fprintf(file.out,"2C FGTUV0 real %3d %10.4f\n",i,fgtuv[i]);

    for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
      map_to_wigner(crystal,&R_AB_1e, &t_12, &Rvec_tmp); // map s_12 back to primitive cell to accelerate convergence
      count = 0;
      for (index_S = 0; index_S < R->number_of_ewald_shells; index_S++) {
        shell_sum = k_zero;
        for (index_R = 0; index_R < R->num[index_S]; index_R++) {
          r_12.comp1 = t_12.comp1 + R->vec_ai[count].comp1;
          r_12.comp2 = t_12.comp2 + R->vec_ai[count].comp2;
          r_12.comp3 = t_12.comp3 + R->vec_ai[count].comp3;
          Rsqrd = r_12.comp1 * r_12.comp1 + r_12.comp2 * r_12.comp2 + r_12.comp3 * r_12.comp3;
          f000m(&em[count][0], gamma_1 * Rsqrd, gamma_1, mm);
          f000m(&en[count][0], Rsqrd / ab_inv[i4], k_one / ab_inv[i4], mm);
          for (i1 = 0; i1 <= mm; i1++) {
            en[count][i1] -= em[count][i1];
           }
            shell_sum += fabs(en[count][0]);
            count++;
           } // end loop on index_R
             } // end loop on index_S
              for (index_R = 0; index_R < count; index_R++) {
                r_12.comp1 = t_12.comp1 + R->vec_ai[index_R].comp1;
                r_12.comp2 = t_12.comp2 + R->vec_ai[index_R].comp2;
                r_12.comp3 = t_12.comp3 + R->vec_ai[index_R].comp3;
                non_recursive_ftuvn(mm, index_R, f, en, &r_12);
                //replace for G functions fgtuv_temp = sab_fac[i4] * ftuvn(m,n,pm,0,&em[0][0],s_12);
                for (t = 0; t <= mm; t++) {
                  for (u = 0; u <= mm; u++) {
                    for (v = 0; v <= mm; v++) {
                      tuv = t+ u + v;
                      if (tuv > mm) break;
                      fgtuv_temp = sab_fac[i4] * f[t][u][v][0];
                  //2019fgtuv_temp = sab_fac[i4] * ftuvn(t,u,v,0,&en[index_R][0],r_12);
                      fgtuv[t * mm1 + u * mm2 + v * mm3  + i4] += fgtuv_temp;
                  //fprintf(file.out,"fgtuv1 real %3d %20.10lf\n",index_R,fgtuv_temp);
                      fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
                      fgtuv_max = k_one;
                     }
                    }
                   }
                  } 
                 } 
              mcmurchie_davidson_2c(coulomb_cart,index_i,index_j,bfposi1,bfposj1,nd2,C1x,C1y,C1z,C2x,C2y,C2z,fgtuv,shells,job,file);
              free(fgtuv);
              two_center_cartesian_to_sh_shell(coulomb_cart,&one_ints->Coulomb[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1, \
              bfposj1,nd2,nd4,shells,job,file);
              bfposj   += shelj;
              bfposj1  += shelj1;
              gausposj += shells->ng_sh[index_j];
             } // close loop over index_j
            bfposi   += sheli;
            bfposi1  += sheli1;
            gausposi += shells->ng_sh[index_i];
           } // close loop over index_i
            dim += nd3 * nd4;
          } // close loop on p

  break;

  case 'S':

  gamma_0 = k_one / G->gamma_0_inv;
  root_gamma_0 = sqrt(gamma_0);
  gamma_1 = k_one / G->gamma_0_inv;

    dim = 0;

  for (p = 0; p < pair_p->nump; p++) {

    ip = pair_p->cell1[pair_p->posn[p]];
    jp = pair_p->cell2[pair_p->posn[p]];
    gi = 0;
    gj = pair_p->latt2[pair_p->posn[p]];
    nd1 = atoms->bfnnumb[ip];
    nd2 = atoms->bfnnumb[jp];
    nd3 = atoms->bfnnumb_sh[ip];
    nd4 = atoms->bfnnumb_sh[jp];
    double coulomb_cart[nd1 * nd2];
    for(i = 0; i < nd1 * nd2; i++) coulomb_cart[i] = k_zero;
    shelposi = atoms->shelposn_sh[ip];
    gausposi = atoms->gausposn_sh[ip];
    bfposi  = 0;
    bfposi1 = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
      sheli  = shells->shar[index_i];
      sheli1 = shells->cart[index_i];
      imax   = shells->imax_sh[index_i];
      //E_coefficients_1c(index_i,gausposi,C1x,C1y,C1z,&C1_max,shells,gaussians,job,file);
      shelposj = atoms->shelposn_sh[jp];
      gausposj = atoms->gausposn_sh[jp];
      bfposj  = 0;
      bfposj1 = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
        R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1 - R->vec_ai[gj].comp1;
        R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2 - R->vec_ai[gj].comp2;
        R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3 - R->vec_ai[gj].comp3;
        R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);
        shelj  = shells->shar[index_j];
        shelj1 = shells->cart[index_j];
        jmax   = shells->imax_sh[index_j];
        double sab_fac[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double ab_inv[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double C1x[shells->ng_sh[index_i] * (imax+1) * (imax+1)];
        double C1y[shells->ng_sh[index_i] * (imax+1) * (imax+1)];
        double C1z[shells->ng_sh[index_i] * (imax+1) * (imax+1)];
        double C2x[shells->ng_sh[index_j] * (jmax+1) * (jmax+1)];
        double C2y[shells->ng_sh[index_j] * (jmax+1) * (jmax+1)];
        double C2z[shells->ng_sh[index_j] * (jmax+1) * (jmax+1)];
        count = 0;
        for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
          for (j4 = 0; j4 < shells->ng_sh[index_j]; j4++) {
            ab_inv[count] = k_one / gaussians->expo_sh[gausposi + i4] + k_one / gaussians->expo_sh[gausposj + j4];
            sab_fac[count] = pi32 * pi32 / sqrt(gaussians->expo_sh[gausposi + i4] * gaussians->expo_sh[gausposj + j4]) / \
            gaussians->expo_sh[gausposi + i4] / gaussians->expo_sh[gausposj + j4] * gaussians->c_sh[gausposi + i4] * \
	    gaussians->c_sh[gausposj + j4];
            count++;
           }
          }
        E_coefficients_1c(index_i,gausposi,C1x,C1y,C1z,&C1_max,shells,gaussians,job,file);
        E_coefficients_1c(index_j,gausposj,C2x,C2y,C2z,&C2_max,shells,gaussians,job,file);
        mm  = imax + jmax;
        mm4 = shells->ng_sh[index_j];
        mm3 = shells->ng_sh[index_i] * mm4;
        mm2 = (mm + 1) * mm3;
        mm1 = (mm + 1) * mm2;
        mm0 = (mm + 1) * mm1;
        fgtuv = (double *) calloc(mm0, sizeof(double));
        if (fgtuv == NULL) {
        if (job->taskid == 0)
        fprintf(stderr, "ERROR: not enough memory for double fgtuv! \n");
        MPI_Finalize();
        exit(1);
       }

int ijk, k, i1;
int index_R, index_S, index_G;
VECTOR_DOUBLE r_12, s_12, t_12, Rvec_tmp;
double en[R->last_ewald_vector][55], em[R->last_ewald_vector][55], in[55], hn[55];
double sum5[16 + 1], sum2[16 + 1], x;
double fac1, fac2, Rzsqrd, expfac, x1, x2, sign, dot_product;
double D_erf[16 + 1], D_exp[16 + 1], derivative_1[16 + 1], derivative_2[16 + 1];
double shell_sum;
double *p_fgtuv;

            fgtuv_max = k_zero;
            for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
              //C_max = C1x_max[i4] * C1y_max[i4] * C1z_max[i4];
              //for (j4 = 0; j4 < shells->ng_sh[index_k] * shells->ng_sh[index_l]; j4++) {
                //C_max *= C2x_max[j4] * C2y_max[j4] * C2z_max[j4];
                s_12.comp1 = R_AB_1e.comp1;
                s_12.comp2 = R_AB_1e.comp2;
                s_12.comp3 = R_AB_1e.comp3;

                Rsqrd = double_vec_dot(&s_12, &s_12);
                Rzsqrd = s_12.comp3 * s_12.comp3;
                fac1 = sab_fac[i4];
                x = root_gamma_0 * s_12.comp3;
                erf_exp_derivative(D_erf, D_exp, mm, root_gamma_0, x, G);
                sum5[0] = D_exp[0] / rtpi / root_gamma_0 + s_12.comp3 * D_erf[0];
                for (n = 1; n <= mm; n++) {
                  sum5[n] = D_exp[n] / rtpi / root_gamma_0 + s_12.comp3 * D_erf[n] + (double) n * D_erf[n - 1];
                 }
                for (k = 0; k <= mm; k++) {
                  p_fgtuv = fgtuv + k * mm3  + i4;
                 *p_fgtuv -= fac1 * two * pi / crystal->primitive_cell_volume * sum5[k]; 
                }
                count = 1;
                for (index_S = 1; index_S < G->number_of_shells; index_S++) {
                  x1 = G->shell_mag[9 * index_S + 1] / two / root_gamma_0 + root_gamma_0 * s_12.comp3;
                  erfc_derivative(derivative_1, mm, root_gamma_0, x1, G);
                  x2 = G->shell_mag[9 * index_S + 1] / two / root_gamma_0 - root_gamma_0 * s_12.comp3;
                  erfc_derivative(derivative_2, mm, -root_gamma_0, x2, G);
                  expfac = exp(G->shell_mag[9 * index_S + 1] * s_12.comp3);
                  for (n = 0; n <= mm; n++) {
                    sign = k_one;
                    sum2[n] = k_zero;
                    for (k = 0; k <= n; k++) {
                      sum2[n] += G->A->a[n][k] * G->shell_mag[9 * index_S + k] * (expfac * derivative_1[n - k] + sign / expfac \
                      * derivative_2[n - k]) / G->shell_mag[9 * index_S + 1];
                      sign *= -k_one;
                     }
                    }
                  //double test = G->EXPFAC[index_S] * fac1 * \
                  (k_one > G->shell_mag[index_S * 9 + mm] ? k_one : G->shell_mag[index_S * 9 + mm]);
                  //if (fabs(test) * C_max < 1.0e-18) break;
                  for (index_G = 0; index_G < G->num[index_S]; index_G++) {
                    dot_product = G->vec_b2[count].comp1 * s_12.comp1 + G->vec_b2[count].comp2 * s_12.comp2;
                    for (i = 0; i <= mm; i++) {
                      for (j = 0; j <= mm; j++) {
                        for (k = 0; k <= mm; k++) {
                          ijk = i + j + k;
                          if (ijk > mm) break;
                          fac2 = fac1 * two * pi / crystal->primitive_cell_volume * sum2[k];
                          fgtuv_temp = fac2 * G->x[count * 9 + i] * G->y[count * 9 + j] * cosfactor(i + j, dot_product);
                          //fprintf(file.out,"Reci sum %3d %3d %14.6e %3d\n",i4,j4,fgtuv_temp,index_S);
                          p_fgtuv = fgtuv + i * mm1 + j * mm2 + k * mm3  + i4;
                         *p_fgtuv += fgtuv_temp;
                          fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
                          fgtuv_max = k_one;
                         }
                        }
                       }
                      count++;
                     }
                    }
                   }
                  //}

  for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
      p_inv_ex = ab_inv[i4];
      s_12.comp1 = R_AB_1e.comp1;
      s_12.comp2 = R_AB_1e.comp2;
      s_12.comp3 = R_AB_1e.comp3;
      map_to_wigner(crystal,&s_12, &t_12, &Rvec_tmp); // map s_12 back to primitive cell to accelerate convergence
      fac1 = sab_fac[i4];
      count = 0;
      for (index_S = 0; index_S < R->number_of_ewald_shells; index_S++) {
        shell_sum = k_zero;
        for (index_R = 0; index_R < R->num[index_S]; index_R++) {
          r_12.comp1 = t_12.comp1 + R->vec_ai[count].comp1;
          r_12.comp2 = t_12.comp2 + R->vec_ai[count].comp2;
          r_12.comp3 = t_12.comp3 + R->vec_ai[count].comp3;
          Rsqrd = r_12.comp1 * r_12.comp1 + r_12.comp2 * r_12.comp2 + r_12.comp3 * r_12.comp3;
          f000m(&em[count][0], gamma_1 * Rsqrd, gamma_1, mm);
          f000m(&en[count][0], Rsqrd / p_inv_ex, k_one / p_inv_ex, mm);
          for (i1 = 0; i1 <= mm; i1++) {
            en[count][i1] -= em[count][i1];
           }
            shell_sum += fabs(en[count][0]);
            count++;
           } // end loop on index_R
            //if (fabs(fac1 * shell_sum * C_max) < 1.0e-14)
              //break;
      //fprintf(file.out,"Real sum %3d %3d %3d %3d %14.6e %14.6e %14.6e %14.6e %3d %3d\n",i4,j4,k4,l4,1.0/p_inv_ex,Rsqrd/p_inv_ex,en[count][0],shell_sum,index_S,index_R);
              //fprintf(file.out,"cc shell sum %3d %12.5e %12.5e\n",index_S,shell_sum,fac1);
             } // end loop on index_S
              for (index_R = 0; index_R < count; index_R++) {
                r_12.comp1 = t_12.comp1 + R->vec_ai[index_R].comp1;
                r_12.comp2 = t_12.comp2 + R->vec_ai[index_R].comp2;
                r_12.comp3 = t_12.comp3 + R->vec_ai[index_R].comp3;
                non_recursive_ftuvn(mm, index_R, f, en, &r_12);
                for (i = 0; i <= mm; i++) {
                  for (j = 0; j <= mm; j++) {
                    for (k = 0; k <= mm; k++) {
                      ijk = i + j + k;
                      if (ijk > mm) break;
                      fgtuv_temp = fac1 * f[i][j][k][0];
                      p_fgtuv  = fgtuv + i * mm1 + j * mm2 + k * mm3  + i4;
                     *p_fgtuv += fgtuv_temp;
                     }
                    }
                   }
                  } 
                 } 

        //Rsqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);
        //for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
          //f000m(&em[0][0], Rsqrd / ab_inv[i4], k_one / ab_inv[i4], mm);
          //non_recursive_ftuvn(mm, 0, f, em, &R_AB_1e);
          //for (t = 0; t <= mm; t++) {
            //for (u = 0; u <= mm; u++) {
              //for (v = 0; v <= mm; v++) {
                //if (t + u + v > mm) break;
                  //fgtuv_temp = f[t][u][v][0]* sab_fac[i4];
                  //printf("fgtuv2 %3d %3d %3d %10.4lf %10.4lf %10.4lf %10.4lf\n",t,u,v,fgtuv_temp,f[t][u][v][0],Rsqrd,ab_inv[i4]);
                  //fgtuv[t * mm1 + u * mm2 + v * mm3  + i4] += fgtuv_temp;
                  //fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
                  //fgtuv_max = k_one;
                 //}
                //}
               //}
              //}

              mcmurchie_davidson_2c(coulomb_cart,index_i,index_j,bfposi1,bfposj1,nd2,C1x,C1y,C1z,C2x,C2y,C2z,fgtuv,shells,job,file);
              free(fgtuv);
              two_center_cartesian_to_sh_shell(coulomb_cart,&one_ints->Coulomb[dim3+dim],index_i,index_j,bfposi,bfposj,bfposi1, \
              bfposj1,nd2,nd4,shells,job,file);
              bfposj   += shelj;
              bfposj1  += shelj1;
              gausposj += shells->ng_sh[index_j];
             } // close loop over index_j
            bfposi   += sheli;
            bfposi1  += sheli1;
            gausposi += shells->ng_sh[index_i];
           } // close loop over index_i
            dim += nd3 * nd4;
          } // close loop on p

  case 'P':

  imax = shells->imax_sh[index_i];
  jmax = shells->imax_sh[index_j];
  gamma_0 = k_one / G->gamma_0_inv;
  root_gamma_0 = sqrt(gamma_0);

  break;

      } // close switch

}

void two_centre_coulomb1_crystal(ComplexMatrix *V_q, REAL_LATTICE *R, Q_LATTICE *q_G, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  int i, j, k, m, n, p, t, u, v;
  int ip, jp, gi, gj;
  int i1, i4, j4;
  int index_i, index_j, index_G, index_R, index_S;
  int bfposi, bfposj, bfposi1, bfposj1;
  int bfposip, bfposjp;
  int imax, jmax, im, jm;
  int gausposi, gausposj, shelposi, shelposj;
  int sheli, shelj, sheli1, shelj1;
  int nd1, nd2, nd3, nd4;
  int count, size;
  int mm, mm0, mm1, mm2, mm3, mm4, ijk;
  double Rsqrd, fn[55];
  double expfac, fac2, gamma_1, dot_product;
  //double em[1][55];
  //double en[1][55];
double en[R->last_ewald_vector][55], em[R->last_ewald_vector][55];
  double f[13][13][13][13];
  double R_AB_1esqrd, C1_max, C2_max;
  double shell_sum, pi_vol = pi / crystal->primitive_cell_volume;
  Complex fgtuv_max, fgtuv_temp, temp;
  Complex *fgtuv, *p_fgtuv;
  VECTOR_DOUBLE R_AB_1e, r_12, s_12;

  ResetComplexMatrix(V_q);
  
  size = 4 * job->l_max + 1;
  gamma_1 = k_one / q_G->gamma_0_inv;
  //fprintf(file.out,"gamma_1 %10.4f pi_vol %10.4f\n",gamma_1,pi_vol);

  switch (crystal->type[0]) {

  case 'C':

  for (ip = 0; ip < atoms->number_of_atoms_in_unit_cell; ip++) {
  //for (jp = 0; jp < atoms->number_of_atoms_in_unit_cell; jp++) {
    for (jp = 0; jp <= ip; jp++) {

    nd1 = atoms->bfnnumb[ip];
    nd2 = atoms->bfnnumb[jp];
    nd3 = atoms->bfnnumb_sh[ip];
    nd4 = atoms->bfnnumb_sh[jp];
    Complex coulomb_cart[nd1 * nd2];
    Complex coulomb_sh[nd3 * nd4];
    for (i = 0; i < nd1 * nd2; i++) coulomb_cart[i] = Complex(k_zero, k_zero);
    for (i = 0; i < nd3 * nd4; i++) coulomb_sh[i]   = Complex(k_zero, k_zero);
    shelposi = atoms->shelposn_sh[ip];
    gausposi = atoms->gausposn_sh[ip];
    bfposi  = 0;
    bfposi1 = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
      sheli  = shells->shar[index_i];
      sheli1 = shells->cart[index_i];
      imax   = shells->imax_sh[index_i];
      //E_coefficients_1c(index_i,gausposi,C1x,C1y,C1z,&C1_max,shells,gaussians,job,file);
      shelposj = atoms->shelposn_sh[jp];
      gausposj = atoms->gausposn_sh[jp];
      bfposj  = 0;
      bfposj1 = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
        R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1;
        R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2;
        R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3;
        R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);
        shelj  = shells->shar[index_j];
        shelj1 = shells->cart[index_j];
        jmax   = shells->imax_sh[index_j];
        double sab_fac[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double ab_inv[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double C1x[shells->ng_sh[index_i] * (imax+1) * (imax+1)];
        double C1y[shells->ng_sh[index_i] * (imax+1) * (imax+1)];
        double C1z[shells->ng_sh[index_i] * (imax+1) * (imax+1)];
        double C2x[shells->ng_sh[index_j] * (jmax+1) * (jmax+1)];
        double C2y[shells->ng_sh[index_j] * (jmax+1) * (jmax+1)];
        double C2z[shells->ng_sh[index_j] * (jmax+1) * (jmax+1)];
        count = 0;
        for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
          for (j4 = 0; j4 < shells->ng_sh[index_j]; j4++) {
            ab_inv[count] = k_one / gaussians->expo_sh[gausposi + i4] + k_one / gaussians->expo_sh[gausposj + j4];
            sab_fac[count] = pi32 * pi32 / sqrt(gaussians->expo_sh[gausposi + i4] * gaussians->expo_sh[gausposj + j4]) / \
            gaussians->expo_sh[gausposi + i4] / gaussians->expo_sh[gausposj + j4] * gaussians->c_sh[gausposi + i4] * \
            gaussians->c_sh[gausposj + j4];
            count++;
           }
          }
        E_coefficients_1c(index_i,gausposi,C1x,C1y,C1z,&C1_max,shells,gaussians,job,file);
        E_coefficients_1c(index_j,gausposj,C2x,C2y,C2z,&C2_max,shells,gaussians,job,file);
        mm  = imax + jmax;
        mm4 = shells->ng_sh[index_j];
        mm3 = shells->ng_sh[index_i] * mm4;
        mm2 = (mm + 1) * mm3;
        mm1 = (mm + 1) * mm2;
        mm0 = (mm + 1) * mm1;
        fgtuv = (Complex *) calloc(mm0, sizeof(Complex));
        //fgtuv = (double *) calloc(mm0, sizeof(double));
        if (fgtuv == NULL) {
        if (job->taskid == 0)
        fprintf(stderr, "ERROR: not enough memory for Complex fgtuv! \n");
        MPI_Finalize();
        exit(1);
       }

          fgtuv_max = k_zero;
          expfac = four * pi / crystal->primitive_cell_volume;
          for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
             s_12.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1;
             s_12.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2;
             s_12.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3;
             //p_inv_ex = ab_inv[i4] + pc_inv[j4];
             //fac1 = sab[i4] * sc[j4];
             // *p_fgtuv -= pi_vol * fac1 * (gamma_1_inv - pab_inv[i4] - pc_inv[j4]);
             //p_fgtuv = fgtuv + i4 * mm4 + j4;
             if (fabs(q_G->sqr[0]) < 0.00001)
             fgtuv[i4] -= pi_vol * sab_fac[i4] * (q_G->gamma_0_inv - ab_inv[i4]);
   //fprintf(file.out,"pi vol fac %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",pi_vol, sab_fac[i4], q_G->gamma_0_inv, - ab_inv[i4],\
             - pi_vol * sab_fac[i4] * (q_G->gamma_0_inv - ab_inv[i4]),fgtuv[i4]);
int start = 0;
if (fabs(q_G->sqr[0]) < 0.00001) start = 1;

             //for (index_G = 1; index_G < q_G->last_vector; index_G++) {
             for (index_G = start; index_G < q_G->last_vector; index_G++) {
               fac2 = sab_fac[i4] * expfac * exp(-q_G->sqr[index_G] / four * q_G->gamma_0_inv) / q_G->sqr[index_G];
               // uncomment for G only fac2 = sab_fac[i4] * expfac * exp(-q_G->sqr[index_G] / four * ab_inv[i4]) / q_G->sqr[index_G];
               dot_product = q_G->vec[index_G].comp1 * s_12.comp1 + q_G->vec[index_G].comp2 * s_12.comp2 + \
               q_G->vec[index_G].comp3 * s_12.comp3;
               //fprintf(file.out,"fac2 dot %f %f\n", fac2, dot_product);
               for (i = 0; i <= mm; i++) {
                 for (j = 0; j <= mm; j++) {
                   for (k = 0; k <= mm; k++) {
                     ijk = i + j + k;
                     if (ijk > mm) break;
                     fgtuv_temp = fac2 * q_G->x[index_G * size + i] * q_G->y[index_G * size + j] * q_G->z[index_G * size + k] * \
                     conj(cosfactor_complex(ijk, dot_product));
                     //check conj(cosfactor_complex(ijk, dot_product));
                     ////cosfactor_complex(ijk, dot_product);
//fgtuv_temp = Complex(k_zero,k_zero);
                     //fprintf(file.out,"fgtuv %10.4f %10.4lf\n", (fgtuv_temp).real(),(fgtuv_temp).imag());
                     p_fgtuv = fgtuv + i * mm1 + j * mm2 + k * mm3  + i4;
                    *p_fgtuv += fgtuv_temp;
                     //fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
                     //fgtuv_max = k_one;
                    }
                   }
                  }
                 }
                }
                 //for (i=0;i<mm0;i++) fprintf(file.out,"2C FGTUV0 %3d %10.4f %10.4f\n",i,fgtuv[i].real(),fgtuv[i].imag());
                 //for (i=0;i<mm0;i++) fprintf(file.out,"2C FGTUV0 %3d %14.8f %14.8f\n",i,(fgtuv[i]).real(),(fgtuv[i]).imag());

  Complex fac;
                //printf("q_G %10.4lf %10.4lf %10.4lf\n",q_G->vec[0].comp1,q_G->vec[0].comp2,q_G->vec[0].comp3);

  for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
      s_12.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1;
      s_12.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2;
      s_12.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3;
      count = 0;
      for (index_S = 0; index_S < R->number_of_ewald_shells; index_S++) {
        shell_sum = k_zero;
        for (index_R = 0; index_R < R->num[index_S]; index_R++) {
//fprintf(file.out,"%3d %3d %3d %3d\n",count,index_S,index_R,R->num[index_S]); fflush(file.out);
          r_12.comp1 = s_12.comp1 + R->vec_ai[count].comp1;
          r_12.comp2 = s_12.comp2 + R->vec_ai[count].comp2;
          r_12.comp3 = s_12.comp3 + R->vec_ai[count].comp3;
          Rsqrd = r_12.comp1 * r_12.comp1 + r_12.comp2 * r_12.comp2 + r_12.comp3 * r_12.comp3;
          f000m(&em[count][0], gamma_1 * Rsqrd, gamma_1, mm);
          f000m(&en[count][0], Rsqrd / ab_inv[i4], k_one / ab_inv[i4], mm);
          for (i1 = 0; i1 <= mm; i1++) {
            en[count][i1] -= em[count][i1]; 
           }
            shell_sum += fabs(en[count][0]);
            count++;
           } // end loop on index_R
             } // end loop on index_S
              for (index_R = 0; index_R < count; index_R++) {
                r_12.comp1 = s_12.comp1 + R->vec_ai[index_R].comp1;
                r_12.comp2 = s_12.comp2 + R->vec_ai[index_R].comp2;
                r_12.comp3 = s_12.comp3 + R->vec_ai[index_R].comp3;
                dot_product = q_G->vec[0].comp1 * R->vec_ai[index_R].comp1 + q_G->vec[0].comp2 * R->vec_ai[index_R].comp2 + \
                q_G->vec[0].comp3 * R->vec_ai[index_R].comp3;
                //dot_product = q_G->vec[0].comp1 * r_12.comp1 + q_G->vec[0].comp2 * r_12.comp2 + \
                q_G->vec[0].comp3 * r_12.comp3;
                fac = sab_fac[i4] * Complex(cos(dot_product), sin(dot_product));
                //non_recursive_ftuvn(mm, index_R, f, en, &r_12);
                //replace for G functions fgtuv_temp = sab[i4] * sc[j4] * ftuvn(m,n,pm,0,&em[0][0],s_12);
                for (i = 0; i <= mm; i++) {
                  for (j = 0; j <= mm; j++) {
                    for (k = 0; k <= mm; k++) {
                      ijk = i + j + k;
                      if (ijk > mm) break;
                  fgtuv_temp = fac * ftuvn(i,j,k,0,&en[index_R][0],r_12);
                      //fgtuv_temp = fac * f[i][j][k][0];
                     //fprintf(file.out,"FGTUV %14.8lf %14.8lf %14.8lf  %14.8f %14.8lf\n", \
                     dot_product,(fac).real(),(fac).imag(),(fgtuv_temp).real(),(fgtuv_temp).imag());
                      //fgtuv_temp = sab_fac[i4] * f[i][j][k][0];
//fgtuv_temp = Complex(k_zero,k_zero);
                      p_fgtuv  = fgtuv + i * mm1 + j * mm2 + k * mm3  + i4;
                     *p_fgtuv += fgtuv_temp;
                      //fprintf(file.out,"fgtuv1 %3d %3d %20.10lf %20.10lf\n",k,index_R,fgtuv_temp.real(),fgtuv_temp.imag());
                     }
                    }
                   }
                  } 
                 } 

                 //for (i=0;i<mm0;i++) fprintf(file.out,"%3d  %18.12f %18.12f\n",\
                 i,fgtuv[i].real(),fgtuv[i].imag());

              //for (i=0;i<mm0;i++) fprintf(file.out,"2C FGTUV1 %3d %14.8f %14.8f\n",i,(fgtuv[i]).real(),(fgtuv[i]).imag());
            //mcmurchie_davidson_2c(coulomb_cart,index_i,index_j,bfposi1,bfposj1,nd2,C1x,C1y,C1z,C2x,C2y,C2z,fgtuv,shells,job,file);
              mcmurchie_davidson_2c_complex(coulomb_cart,fgtuv,index_i,index_j,bfposi1,bfposj1,nd2,C1x,C1y,C1z,C2x,C2y,C2z,\
              shells,job,file);
              free(fgtuv);
              //two_center_cartesian_to_sh_shell_complex(coulomb_cart,&V_q->a[0][0],index_i,index_j,bfposi,bfposj,bfposi1, \
              bfposj1,nd2,nd4,shells,job,file);
              two_center_cartesian_to_sh_shell_complex(coulomb_cart,coulomb_sh,index_i,index_j,bfposi,bfposj,bfposi1, \
              bfposj1,nd2,nd4,shells,job,file);
              //two_center_cartesian_to_sh_shell_complex(coulomb_cart,&V_q->a[0][0],index_i,index_j,bfposi,bfposj,bfposi1, \
              bfposj1,nd2,nd4,shells,job,file);
              bfposj   += shelj;
              bfposj1  += shelj1;
              gausposj += shells->ng_sh[index_j];
             } // close loop over index_j
            bfposi   += sheli;
            bfposi1  += sheli1;
            gausposi += shells->ng_sh[index_i];
           } // close loop over index_i

              bfposip = atoms->bfnposn_sh[ip];
              bfposjp = atoms->bfnposn_sh[jp];
              count = 0;
              for (i = 0; i < nd3;i++) {
                for (j = 0; j < nd4; j++) {
                  V_q->a[bfposip + i][bfposjp + j] = coulomb_sh[count];
                  count++;
                  //fprintf(file.out,"CCART %3d %3d  %3d %3d %10.4f %10.4lf\n",\
                  bfposip,bfposjp,bfposip+i,bfposjp+j,(V_q->a[bfposip+i][bfposjp+j]).real(),(V_q->a[bfposip+i][bfposjp+j]).imag());
                 }
                }
              count = 0;
              if (ip > jp) {
              //for (i = 0; i < nd4;i++) {
                //for (j = 0; j < nd3; j++) {
              for (i = 0; i < nd3;i++) {
                for (j = 0; j < nd4; j++) {
                  V_q->a[bfposjp + j][bfposip + i] = conj(coulomb_sh[count]);
                  count++;
                 }
                }
               }

          } // close loop on jp
         } // close loop on ip

  break;

  case 'S':

  break;

  case 'P':

  break;

  case 'M':

  break;

      } // close switch

}

void two_centre_exchange_crystal(ComplexMatrix *V_q, REAL_LATTICE *R, Q_LATTICE *q_G, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  int i, j, k, m, n, p, t, u, v;

  int n1, n2, n3, n4, n5, n6;
  int tmax, umax, vmax;
  int off1, off2, off3;

  int ip, jp, gi, gj;
  int i1, i4, j4;
  int index_i, index_j, index_G, index_R, index_S;
  int bfposi, bfposj, bfposi1, bfposj1;
  int imax, jmax, im, jm;
  int gausposi, gausposj, shelposi, shelposj;
  int sheli, shelj, sheli1, shelj1;
  int nd1, nd2, nd3, nd4;
  int count, size;
  int mm, mm0, mm1, mm2, mm3, mm4, ijk;
  double Rsqrd, fn[55];
  double expfac, sqrt_expfac,fac2, gamma_1, dot_product;
  double en[R->last_ewald_vector][55], em[R->last_ewald_vector][55];
  double f[13][13][13][13];
  double R_AB_1esqrd, C1_max, C2_max;
  double shell_sum, pi_vol = pi / crystal->primitive_cell_volume;
  double E1_max;
  Complex fgtuv_max, fgtuv_temp, temp;
  Complex *fgtuv, *p_fgtuv;
  VECTOR_DOUBLE R_AB_1e, r_12, s_12;

  ResetComplexMatrix(V_q);
  
  size = 4 * job->l_max + 1;
  gamma_1 = k_one / q_G->gamma_0_inv;
  //fprintf(file.out,"gamma_1 %10.4f pi_vol %10.4f\n",gamma_1,pi_vol);

gj = 1;

  switch (crystal->type[0]) {

  case 'C':

  for (ip = 0; ip < atoms->number_of_atoms_in_unit_cell; ip++) {
    for (jp = 0; jp <= ip; jp++) {

    nd1 = atoms->bfnnumb[ip];
    nd2 = atoms->bfnnumb[jp];
    nd3 = atoms->bfnnumb_sh[ip];
    nd4 = atoms->bfnnumb_sh[jp];
    Complex coulomb_cart[nd1 * nd2];
    for (i = 0; i < nd1 * nd2; i++) coulomb_cart[i] = Complex(k_zero, k_zero);
    shelposi = atoms->shelposn_sh[ip];
    gausposi = atoms->gausposn_sh[ip];
    bfposi  = 0;
    bfposi1 = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
      sheli  = shells->shar[index_i];
      sheli1 = shells->cart[index_i];
      imax   = shells->imax_sh[index_i];
      //E_coefficients_1c(index_i,gausposi,C1x,C1y,C1z,&C1_max,shells,gaussians,job,file);
      shelposj = atoms->shelposn_sh[jp];
      gausposj = atoms->gausposn_sh[jp];
      bfposj  = 0;
      bfposj1 = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
        R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1 - R->vec_ai[gj].comp1;
        R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2 - R->vec_ai[gj].comp2;
        R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3 - R->vec_ai[gj].comp3;
        R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);
        //shelj  = shells->type_sh[index_j];
        //shelj1 = shells->type1_sh[index_j];
        shelj  = shells->shar[index_j];
        shelj1 = shells->cart[index_j];
        jmax   = shells->imax_sh[index_j];
        int im = 0, jm = 0;
        double sab[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double E1x[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+im+jm+1) * (imax+im+1) * (jmax+jm+1)];
        double E1y[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+im+jmax+jm+1) * (imax+im+1) * (jmax+jm+1)];
        double E1z[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+im+jmax+jm+1) * (imax+im+1) * (jmax+jm+1)];
        double pab_inv[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        VECTOR_DOUBLE R_AB[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        E_coefficients1(ip,jp,gi,gj,index_i,index_j,gausposi,gausposj,im,jm,E1x,E1y,E1z,&E1_max,sab,pab_inv,&R_AB_1esqrd,\
        R_AB,R,atoms,shells,gaussians,job,file);

        mm  = imax + jmax;
        mm4 = shells->ng_sh[index_j];
        mm3 = shells->ng_sh[index_i] * mm4;
        mm2 = (mm + 1) * mm3;
        mm1 = (mm + 1) * mm2;
        mm0 = (mm + 1) * mm1;
        fgtuv = (Complex *) calloc(mm0, sizeof(Complex));
        if (fgtuv == NULL) {
        if (job->taskid == 0)
        fprintf(stderr, "ERROR: not enough memory for Complex fgtuv! \n");
        MPI_Finalize();
        exit(1);
       }

          Complex sum = Complex(k_zero, k_zero);
          fgtuv_max = k_zero;
          sqrt_expfac = sqrt(four * pi / crystal->primitive_cell_volume);
          for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
             //s_12.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1;
             //s_12.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2;
             //s_12.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3;
             s_12.comp1 = R_AB[i4].comp1;
             s_12.comp2 = R_AB[i4].comp2;
             s_12.comp3 = R_AB[i4].comp3;
             fprintf(file.out,"sqrt_expfac %14.8f sab %10.4f %10.4f %10.4f %10.4f\n",\
             sqrt_expfac,sab[i4],s_12.comp1,s_12.comp2,s_12.comp3);
             //p_inv_ex = ab_inv[i4] + pc_inv[j4];
             //fac1 = sab[i4] * sc[j4];
             // *p_fgtuv -= pi_vol * fac1 * (gamma_1_inv - pab_inv[i4] - pc_inv[j4]);
             //p_fgtuv = fgtuv + i4 * mm4 + j4;
             //if (fabs(q_G->sqr[0]) < 0.00001)
             //fgtuv[i4] -= pi_vol * sab_fac[i4] * (q_G->gamma_0_inv - ab_inv[i4]);
             //fprintf(file.out,"pi vol fac %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",pi_vol, sab_fac[i4], q_G->gamma_0_inv, \
             -ab_inv[i4], - pi_vol * sab_fac[i4] * (q_G->gamma_0_inv - ab_inv[i4]),fgtuv[i4]);
             int start = 0;
             if (fabs(q_G->sqr[0]) < 0.00001) start = 1;
             for (index_G = start; index_G < q_G->last_vector; index_G++) {
               //fac2 = sab_fac[i4] * expfac / q_G->sqr[index_G];
               fac2 = sqrt_expfac * sab[i4] *  exp(-q_G->sqr[index_G] * pab_inv[i4] / four) / sqrt(q_G->sqr[index_G]);
sum += fac2 * fac2;
               dot_product = q_G->vec[index_G].comp1 * s_12.comp1 + q_G->vec[index_G].comp2 * s_12.comp2 + \
               q_G->vec[index_G].comp3 * s_12.comp3;
               fprintf(file.out,"fac2 dot %f %10.4f %10.4f %14.8lf\n",fac2, dot_product, k_one / sqrt(q_G->sqr[index_G]),sum.real());
               for (i = 0; i <= mm; i++) {
                 for (j = 0; j <= mm; j++) {
                   for (k = 0; k <= mm; k++) {
                     ijk = i + j + k;
                     if (ijk > mm) break;
                     fgtuv_temp = fac2 * q_G->x[index_G * size + i] * q_G->y[index_G * size + j] * q_G->z[index_G * size + k] * \
                     cosfactor_complex(ijk, dot_product);
                     //fprintf(file.out,"fgtuv %f\n", fgtuv_temp);
                     p_fgtuv = fgtuv + i * mm1 + j * mm2 + k * mm3  + i4;
                    *p_fgtuv += fgtuv_temp;
                     //fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
                     //fgtuv_max = k_one;
                    }
                   }
                  }
                 }
                }
                 //for (i=0;i<mm0;i++) fprintf(file.out,"2C FGTUV0 %3d %10.4f %10.4f\n",i,fgtuv[i].real(),fgtuv[i].imag());

              //for (i=0;i<mm0;i++) fprintf(file.out,"2C FGTUV1 %3d %10.4f %10.4f\n",i,fgtuv[i].real(),fgtuv[i].imag());
  off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  off2 = (imax + im + jmax + jm + 1) * off3;
  off1 = (jmax + jm + 1) * off2;
  sheli = shells->type1_sh[index_i];
  shelj = shells->type1_sh[index_j];
  for (i = 0; i < sheli; i++) {
    for (j = 0; j < shelj; j++) {
      count = (bfposi + i) * nd2 + bfposj + j;
      tmax = shells->tuv[imax][i][0] + shells->tuv[jmax][j][0];
      umax = shells->tuv[imax][i][1] + shells->tuv[jmax][j][1];
      vmax = shells->tuv[imax][i][2] + shells->tuv[jmax][j][2];
      n1 = shells->tuv[imax][i][0];
      n2 = shells->tuv[imax][i][1];
      n3 = shells->tuv[imax][i][2];
      n4 = shells->tuv[jmax][j][0];
      n5 = shells->tuv[jmax][j][1];
      n6 = shells->tuv[jmax][j][2];
      for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
        for (t = 0; t <= tmax; t++) {
          for (u = 0; u <= umax; u++) {
            for (v = 0; v <= vmax; v++) {
              coulomb_cart[count] += E1x[n1 * off1 + n4 * off2 + t * off3 + i4] * E1y[n2 * off1 + n5 * off2 + u * off3 + i4] * \
                                     E1z[n3 * off1 + n6 * off2 + v * off3 + i4] * fgtuv[t * mm1 + u * mm2 + v * mm3 + i4];
             }
            }
           } // end t u v loop
          } // end loop over i4
         } // close loop over j
        } // close loop over i
          free(fgtuv);
          two_center_cartesian_to_sh_shell_complex(coulomb_cart,&V_q->a[0][0],index_i,index_j,bfposi,bfposj,bfposi1, \
          bfposj1,nd2,nd4,shells,job,file);
          bfposj   += shelj;
          bfposj1  += shelj1;
          gausposj += shells->ng_sh[index_j];
         } // close loop over index_j
        bfposi   += sheli;
        bfposi1  += sheli1;
        gausposi += shells->ng_sh[index_i];
       } // close loop over index_i
      } // close loop on jp
     } // close loop on ip

  break;

  case 'S':

  break;

  case 'P':

  break;

  case 'M':

  break;

      } // close switch

}

void two_centre_exchange_crystal1(int p, Complex *q_G_array, PAIR_TRAN *pair_p, REAL_LATTICE *R, Q_LATTICE *q_G, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  int i, j, k, m, n, t, u, v;
  int n1, n2, n3, n4, n5, n6;
  int tmax, umax, vmax;
  int off1, off2, off3;
  int ip, jp, gi, gj;
  int i1, i4, j4;
  int index_i, index_j, index_G, index_R, index_S;
  int bfposi, bfposj, bfposi1, bfposj1;
  int imax, jmax, im = 0, jm = 0;
  int gausposi, gausposj, shelposi, shelposj;
  int sheli, shelj, sheli1, shelj1;
  int nd1, nd2, nd3, nd4, nd12, nd34;
  int start, count, size;
  int mm, mm0, mm1, mm2, mm3, mm4, ijk; 
  int dim3;
  double Rsqrd, fn[55];
  double expfac, sqrt_expfac,fac2, gamma_1, dot_product;
  double en[R->last_ewald_vector][55], em[R->last_ewald_vector][55];
  double f[13][13][13][13];
  double R_AB_1esqrd, C1_max, C2_max;
  //double shell_sum, pi_vol = pi / crystal->primitive_cell_volume;
  double E1_max;
  Complex fgtuv_max, fgtuv_temp, temp;
  Complex *fgtuv, *p_fgtuv;
  VECTOR_DOUBLE R_AB_1e, r_12, s_12;

  size = 4 * job->l_max + 1;
  gamma_1 = k_one / q_G->gamma_0_inv;
  //fprintf(file.out,"gamma_1 %10.4f pi_vol %10.4f\n",gamma_1,pi_vol);

    ip = pair_p->cell1[p]; 
    jp = pair_p->cell2[p]; 
    gj = pair_p->latt2[p]; 

    nd1 = atoms->bfnnumb[ip];
    nd2 = atoms->bfnnumb[jp];
    nd3 = atoms->bfnnumb_sh[ip];
    nd4 = atoms->bfnnumb_sh[jp];
    nd12 = nd1 * nd2;
    nd34 = nd3 * nd4;
    Complex coulomb_cart[nd1 * nd2];

  switch (crystal->type[0]) {

  case 'C':

    sqrt_expfac = sqrt(four * pi / crystal->primitive_cell_volume);
    shelposi = atoms->shelposn_sh[ip];
    gausposi = atoms->gausposn_sh[ip];
    bfposi  = 0;
    bfposi1 = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
      sheli  = shells->shar[index_i];
      sheli1 = shells->cart[index_i];
      imax   = shells->imax_sh[index_i];
      //E_coefficients_1c(index_i,gausposi,C1x,C1y,C1z,&C1_max,shells,gaussians,job,file);
      shelposj = atoms->shelposn_sh[jp];
      gausposj = atoms->gausposn_sh[jp];
      bfposj  = 0;
      bfposj1 = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
        R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1 - R->vec_ai[gj].comp1;
        R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2 - R->vec_ai[gj].comp2;
        R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3 - R->vec_ai[gj].comp3;
        R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);
        //shelj  = shells->type_sh[index_j];
        //shelj1 = shells->type1_sh[index_j];
        shelj  = shells->shar[index_j];
        shelj1 = shells->cart[index_j];
        jmax   = shells->imax_sh[index_j];
        double sab[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double E1x[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+im+jm+1) * (imax+im+1) * (jmax+jm+1)];
        double E1y[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+im+jmax+jm+1) * (imax+im+1) * (jmax+jm+1)];
        double E1z[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+im+jmax+jm+1) * (imax+im+1) * (jmax+jm+1)];
        double pab_inv[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        VECTOR_DOUBLE R_AB[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        E_coefficients1(ip,jp,gi,gj,index_i,index_j,gausposi,gausposj,im,jm,E1x,E1y,E1z,&E1_max,sab,pab_inv,&R_AB_1esqrd,\
        R_AB,R,atoms,shells,gaussians,job,file);

        mm  = imax + jmax;
        mm4 = shells->ng_sh[index_j];
        mm3 = shells->ng_sh[index_i] * mm4;
        mm2 = (mm + 1) * mm3;
        mm1 = (mm + 1) * mm2;
        mm0 = (mm + 1) * mm1;
        fgtuv = (Complex *) malloc(mm0 * sizeof(Complex));
        if (fgtuv == NULL) {
        if (job->taskid == 0)
        fprintf(stderr, "ERROR: not enough memory for Complex fgtuv! \n");
        MPI_Finalize();
        exit(1);
       }

        dim3 = 0;
        for (index_G = 0; index_G < q_G->last_vector; index_G++) {
          ResetComplexArray(fgtuv,&mm0);
          ResetComplexArray(coulomb_cart,&nd12);
          if (index_G == 0 && fabs(q_G->sqr[0]) < 0.00001) {
          dim3 += nd34; 
          continue;
         }
          fgtuv_max = k_zero;
          for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
            s_12.comp1 = R_AB[i4].comp1;
            s_12.comp2 = R_AB[i4].comp2;
            s_12.comp3 = R_AB[i4].comp3;
            fac2 = sqrt_expfac * sab[i4] *  exp(-q_G->sqr[index_G] * pab_inv[i4] / four) / sqrt(q_G->sqr[index_G]);
            dot_product = q_G->vec[index_G].comp1 * s_12.comp1 + q_G->vec[index_G].comp2 * s_12.comp2 + \
            q_G->vec[index_G].comp3 * s_12.comp3;
            for (i = 0; i <= mm; i++) {
              for (j = 0; j <= mm; j++) {
                for (k = 0; k <= mm; k++) {
                  ijk = i + j + k;
                  if (ijk > mm) break;
                  fgtuv_temp = fac2 * q_G->x[index_G * size + i] * q_G->y[index_G * size + j] * q_G->z[index_G * size + k] * \
                  cosfactor_complex(ijk, dot_product);
                  p_fgtuv = fgtuv + i * mm1 + j * mm2 + k * mm3  + i4;
                 *p_fgtuv += fgtuv_temp;
                 }
                }
               }
              } // close loop on i4
            off3 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
            off2 = (imax + im + jmax + jm + 1) * off3;
            off1 = (jmax + jm + 1) * off2;
            for (i = 0; i < sheli; i++) {
              for (j = 0; j < shelj; j++) {
                count = (bfposi + i) * nd2 + bfposj + j;
                tmax = shells->tuv[imax][i][0] + shells->tuv[jmax][j][0];
                umax = shells->tuv[imax][i][1] + shells->tuv[jmax][j][1];
                vmax = shells->tuv[imax][i][2] + shells->tuv[jmax][j][2];
                n1 = shells->tuv[imax][i][0];
                n2 = shells->tuv[imax][i][1];
                n3 = shells->tuv[imax][i][2];
                n4 = shells->tuv[jmax][j][0];
                n5 = shells->tuv[jmax][j][1];
                n6 = shells->tuv[jmax][j][2];
                for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
                  for (t = 0; t <= tmax; t++) {
                    for (u = 0; u <= umax; u++) {
                      for (v = 0; v <= vmax; v++) {
                        coulomb_cart[count] += E1x[n1 * off1 + n4 * off2 + t * off3 + i4] * \
                                               E1y[n2 * off1 + n5 * off2 + u * off3 + i4] * \
                                               E1z[n3 * off1 + n6 * off2 + v * off3 + i4] * \
                                               fgtuv[t * mm1 + u * mm2 + v * mm3 + i4];
                       }
                      }
                     } // end t u v loop
                    } // end loop over i4
                   } // close loop over j
                  } // close loop over i
                   two_center_cartesian_to_sh_shell_complex(coulomb_cart,&q_G_array[dim3],index_i,index_j,bfposi,bfposj,bfposi1, \
                   bfposj1,nd2,nd4,shells,job,file);
                   dim3 += nd34;
                  } // close loop on index_G
                   free(fgtuv);
                   bfposj   += shelj;
                   bfposj1  += shelj1;
                   gausposj += shells->ng_sh[index_j];
                 } // close loop over index_j
                  bfposi   += sheli;
                  bfposi1  += sheli1;
                  gausposi += shells->ng_sh[index_i];
                } // close loop over index_i

  break;

  case 'S':

  break;

  case 'P':

  break;

  case 'M':

  break;

      } // close switch

}

void two_centre_overlap(INT_1E *one_ints, PAIR_TRAN *pair_p, REAL_LATTICE *R, RECIPROCAL_LATTICE *G, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  int ip, jp, gi, gj, p;
  int i, j, i4, j4, l, mm, nn, atm0, count;
  int index_i, index_j;
  int dim;
  int bfposi, bfposj, bfposi1, bfposj1;
  int imax, jmax, im, jm;
  int gausposi, gausposj, shelposi, shelposj;
  int sheli, shelj, sheli1, shelj1;
  int nd1, nd2, nd3, nd4, nd12, nd34;
  int nsheli, nshelj;
  double R_AB_1esqrd;
  double E1_max;
  double time1, time2;
  VECTOR_DOUBLE R_AB_1e;

  im = 0; 
  jm = 0;
  dim = 0;
  for (p = 0; p < pair_p->nump; p++) {
    ip = pair_p->cell1[pair_p->posn[p]];
    jp = pair_p->cell2[pair_p->posn[p]];
    gi = 0;
    gj = pair_p->latt2[pair_p->posn[p]];
    nd1 = atoms->bfnnumb[ip];
    nd2 = atoms->bfnnumb[jp];
    nd3 = atoms->bfnnumb_sh[ip];
    nd4 = atoms->bfnnumb_sh[jp];
    nd12 = nd1 * nd2;
    nd34 = nd3 * nd4;
    double overlap_cart[nd12];
    for (i = 0; i < nd12; i++) overlap_cart[i] = k_zero;
    shelposi = atoms->shelposn_sh[ip];
    gausposi = atoms->gausposn_sh[ip];
    bfposi  = 0;
    bfposi1 = 0;
    for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
      sheli  = shells->type_sh[index_i];
      sheli1 = shells->type1_sh[index_i];
      imax   = shells->imax_sh[index_i];
      shelposj = atoms->shelposn_sh[jp];
      gausposj = atoms->gausposn_sh[jp];
      bfposj  = 0;
      bfposj1 = 0;
      for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
        R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1 - R->vec_ai[gj].comp1;
        R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2 - R->vec_ai[gj].comp2;
        R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3 - R->vec_ai[gj].comp3;
        R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);
        shelj  = shells->type_sh[index_j];
        shelj1 = shells->type1_sh[index_j];
        jmax   = shells->imax_sh[index_j];
        double sab[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double E1x[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
        double E1y[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
        double E1z[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
        double pab_inv[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        VECTOR_DOUBLE R_AB[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        E_coefficients1(ip,jp,gi,gj,index_i,index_j,gausposi,gausposj,im,jm,E1x,E1y,E1z,&E1_max,sab,pab_inv,&R_AB_1esqrd,R_AB,\
        R,atoms,shells,gaussians,job,file);
        Overlap(overlap_cart, index_i, index_j, bfposi1, bfposj1, nd2, im, jm, E1x, E1y, E1z, sab, shells, job, file);
        two_center_cartesian_to_sh_shell(overlap_cart,&one_ints->Overlap[dim],index_i,index_j,bfposi,bfposj,bfposi1,bfposj1,\
        nd2,nd4,shells,job,file);
        bfposj   += shelj;
        bfposj1  += shelj1;
        gausposj += shells->ng_sh[index_j];
       } // close loop over index_j
      bfposi   += sheli;
      bfposi1  += sheli1;
      gausposi += shells->ng_sh[index_i];
     } // close loop over index_i
    dim += nd34;
    //dim4 += 3 * nd34;
 } // close loop on p

}

void E_coefficients1(int ip, int jp, int gi, int gj, int index_i, int index_j, int gausposi, int gausposj, int im, int jm, double *C1x, double *C1y, double *C1z, double *C1_max, double *sab, double *pab_inv, double *R_AB_1esqrd, VECTOR_DOUBLE *R_AB, REAL_LATTICE *R, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, JOB_PARAM *job, FILES file)
//CHANGES2015void E_coefficients1(int ip, int jp, int gi, int gj, int index_i, int index_j, int gausposi, int gausposj, int imax, int jmax, double *C1x, double *C1y, double *C1z, double *C1_max, double *sab, double *pab_inv, double *R_AB_1esqrd, VECTOR_DOUBLE *R_AB, REAL_LATTICE *R, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, JOB_PARAM *job, FILES file)

{

int i4, j4;
int t, m, n;
int count;
int off1, off2, off3;
int imax, jmax;
double ab, pab, p32, expnta, expntb;
double SAB, KAB;
double PAx, PAy, PAz, PBx, PBy, PBz;
double C1x_max, C1y_max, C1z_max;

  C1x_max = k_zero;
  C1y_max = k_zero;
  C1z_max = k_zero;
 *C1_max  = k_zero;
  imax = shells->imax_sh[index_i];
  jmax = shells->imax_sh[index_j];
  off3    = shells->ng_sh[index_i] * shells->ng_sh[index_j];
  //CHANGES2015off2    = (imax + jmax + 1) * off3;
  //off1    = (jmax + 1) * off2;
  off2    = (imax + im + jmax + jm + 1) * off3;
  off1    = (jmax + jm + 1) * off2;
  count = 0;
  for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
    for (j4 = 0; j4 < shells->ng_sh[index_j]; j4++) {
      pab = gaussians->expo_sh[gausposi + i4] + gaussians->expo_sh[gausposj + j4];
      ab  = gaussians->expo_sh[gausposi + i4] * gaussians->expo_sh[gausposj + j4];
      p32 = pab * sqrt(pab);
      KAB = ab * *R_AB_1esqrd / pab ;
      SAB = pi32 * exp(-KAB) / p32 ;
      pab_inv[count] = k_one / pab;
      sab[count]     = gaussians->c_sh[gausposi + i4] * gaussians->c_sh[gausposj + j4] * SAB;
      expnta = gaussians->expo_sh[gausposi + i4];
      expntb = gaussians->expo_sh[gausposj + j4];
      R_AB[count].comp1 = (expnta * (atoms->cell_vector[ip].comp1 + R->vec_ai[gi].comp1) + \
      expntb * (atoms->cell_vector[jp].comp1 + R->vec_ai[gj].comp1)) / pab;
      R_AB[count].comp2 = (expnta * (atoms->cell_vector[ip].comp2 + R->vec_ai[gi].comp2) + \
      expntb * (atoms->cell_vector[jp].comp2 + R->vec_ai[gj].comp2)) / pab;
      R_AB[count].comp3 = (expnta * (atoms->cell_vector[ip].comp3 + R->vec_ai[gi].comp3) + \
      expntb * (atoms->cell_vector[jp].comp3 + R->vec_ai[gj].comp3)) / pab;
      PAx = R_AB[count].comp1 - atoms->cell_vector[ip].comp1 - R->vec_ai[gi].comp1;
      PAy = R_AB[count].comp2 - atoms->cell_vector[ip].comp2 - R->vec_ai[gi].comp2;
      PAz = R_AB[count].comp3 - atoms->cell_vector[ip].comp3 - R->vec_ai[gi].comp3;
      PBx = R_AB[count].comp1 - atoms->cell_vector[jp].comp1 - R->vec_ai[gj].comp1;
      PBy = R_AB[count].comp2 - atoms->cell_vector[jp].comp2 - R->vec_ai[gj].comp2;
      PBz = R_AB[count].comp3 - atoms->cell_vector[jp].comp3 - R->vec_ai[gj].comp3;
      //CHANGES2015for (m = 0; m <= imax; m++) {
        //for (n = 0; n <= jmax; n++) {
          //for (t = 0; t <= imax + jmax; t++) {
      for (m = 0; m <= imax + im; m++) {
        for (n = 0; n <= jmax + jm; n++) {
          for (t = 0; t <= imax + im + jmax + jm; t++) {
            C1x[m * off1 + n * off2 + t * off3 + count] = e(m, n, t, pab, PAx, PBx);
            C1y[m * off1 + n * off2 + t * off3 + count] = e(m, n, t, pab, PAy, PBy);
            C1z[m * off1 + n * off2 + t * off3 + count] = e(m, n, t, pab, PAz, PBz);
            C1x_max = (C1x_max > fabs(C1x[m * off1 + n * off2 + t * off3 + count])) ? \
            C1x_max : fabs(C1x[m * off1 + n * off2 + t * off3 + count]);
            C1y_max = (C1y_max > fabs(C1y[m * off1 + n * off2 + t * off3 + count])) ? \
            C1y_max : fabs(C1y[m * off1 + n * off2 + t * off3 + count]);
            C1z_max = (C1z_max > fabs(C1z[m * off1 + n * off2 + t * off3 + count])) ? \
            C1z_max : fabs(C1z[m * off1 + n * off2 + t * off3 + count]);
          }
         }
        } // end loop to set up C1 factors
       *C1_max = (*C1_max > C1x_max * C1y_max * C1z_max) ? *C1_max : C1x_max * C1y_max * C1z_max;
        count++;
      }
     }

}

void E_coefficients_1c(int index_i, int gausposi, double *C1x, double *C1y, double *C1z, double *C1_max, SHELL *shells, GAUSSIAN *gaussians, JOB_PARAM *job, FILES file)

{

int i4, j4;
int m, t;
int count;
int off1, off2;
int imax, jmax;
double pa, C1x_max, C1y_max, C1z_max;

  C1x_max = k_zero;
  C1y_max = k_zero;
  C1z_max = k_zero;
 *C1_max  = k_zero;
  imax = shells->imax_sh[index_i];
  off2    = shells->ng_sh[index_i];
  off1    = (imax + 1) * off2;

  count = 0;
  for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
    pa = gaussians->expo_sh[gausposi + i4];
    for (m = 0; m <= imax; m++) {
      for (t = 0; t <= imax; t++) {
        C1x[m * off1 + t * off2 + count] = e(m, 0, t, pa, k_zero, k_zero);
        C1y[m * off1 + t * off2 + count] = e(m, 0, t, pa, k_zero, k_zero);
        C1z[m * off1 + t * off2 + count] = e(m, 0, t, pa, k_zero, k_zero);
       *C1_max = (*C1_max > C1x_max * C1y_max * C1z_max) ? *C1_max : C1x_max * C1y_max * C1z_max;
       }
      }
     count++;
     }

}

void E_coefficients_2c(int ip, int jp, int gi, int gj, int index_i, int index_j, int gausposi, int gausposj, double *C1x, double *C1y, double *C1z, double *C1_max, double *C2x, double *C2y, double *C2z, double *C2_max, double *sab_fac, double *ab_inv, double *R_AB_1esqrd, REAL_LATTICE *R, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, JOB_PARAM *job, FILES file)

{

int i4, j4;
int m, t;
int count;
int off1, off2, off3, off4;
int imax, jmax;
double ab, pab, p32, pa, pb;
double SAB, KAB;
double C1x_max, C1y_max, C1z_max;
double C2x_max, C2y_max, C2z_max;

  C1x_max = k_zero;
  C1y_max = k_zero;
  C1z_max = k_zero;
 *C1_max  = k_zero;
  C2x_max = k_zero;
  C2y_max = k_zero;
  C2z_max = k_zero;
 *C2_max  = k_zero;
  imax = shells->imax_sh[index_i];
  jmax = shells->imax_sh[index_j];
  off4    = shells->ng_sh[index_j];
  off3    = (jmax + 1) * off4;
  off2    = shells->ng_sh[index_i];
  off1    = (imax + 1) * off2;

  //count = 0;
  //for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
    //for (j4 = 0; j4 < shells->ng_sh[index_j]; j4++) {
      //pab = gaussians->expo_sh[gausposi + i4] + gaussians->expo_sh[gausposj + j4];
      //ab  = gaussians->expo_sh[gausposi + i4] * gaussians->expo_sh[gausposj + j4];
      //p32 = pab * sqrt(pab);
      //KAB = ab * *R_AB_1esqrd / pab ;
      //SAB = pi32 * exp(-KAB) / p32 ;
      //ab_inv[count] = k_one / gaussians->expo_sh[gausposi + i4] + k_one / gaussians->expo_sh[gausposj + j4];
      //sab_fac[count] =pi32*pi32/sqrt(gaussians->expo_sh[gausposi+i4]*gaussians->expo_sh[gausposj + j4]) / gaussians->expo_sh[gausposi + i4] / \
      gaussians->expo_sh[gausposj + j4] * gaussians->c_sh[gausposi + i4] * gaussians->c_sh[gausposj + j4];
      //count++;
     //}
    //}

  count = 0;
  for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
    pa = gaussians->expo_sh[gausposi + i4];
    for (m = 0; m <= imax; m++) {
      for (t = 0; t <= imax; t++) {
        C1x[m * off1 + t * off2 + count] = e(m, 0, t, pa, k_zero, k_zero);
        C1y[m * off1 + t * off2 + count] = e(m, 0, t, pa, k_zero, k_zero);
        C1z[m * off1 + t * off2 + count] = e(m, 0, t, pa, k_zero, k_zero);
       *C1_max = (*C1_max > C1x_max * C1y_max * C1z_max) ? *C1_max : C1x_max * C1y_max * C1z_max;
       }
      }
     count++;
     }

  count = 0;
  for (j4 = 0; j4 < shells->ng_sh[index_j]; j4++) {
    pb = gaussians->expo_sh[gausposj + j4];
    for (m = 0; m <= jmax; m++) {
      for (t = 0; t <= jmax; t++) {
        C2x[m * off3 + t * off4 + count] = e(m, 0, t, pb, k_zero, k_zero);
        C2y[m * off3 + t * off4 + count] = e(m, 0, t, pb, k_zero, k_zero);
        C2z[m * off3 + t * off4 + count] = e(m, 0, t, pb, k_zero, k_zero);
       *C2_max = (*C2_max > C2x_max * C2y_max * C2z_max) ? *C2_max : C2x_max * C2y_max * C2z_max;
       }
      }
     count++;
     }

}

*/
/*

void mcmurchie_davidson_2c_complex(Complex *F_cart, Complex *fgtuv, int index_i, int index_j, int bfposi, int bfposj, int nd2, double *C1x, double *C1y, double *C1z, double *C2x, double *C2y, double *C2z, SHELL *shells, JOB_PARAM *job, FILES file)

{

int i, j;
int i4, j4;
int n1, n2, n3, n4, n5, n6;
int off1, off2, off3, off4;
int mm, mm1, mm2, mm3, mm4;
int imax, jmax;
int t, u, v, tmax, umax, vmax;
int tp, up, vp, tpmax, upmax, vpmax;
int tpupvpsign, tpupvp2;
int sheli, shelj, sheli1, shelj1;
int nsheli, nshelj;
int count;
Complex c1fac[15][15][15];
//Complex *p_fgtuv;

//fprintf(file.out,"%3d %3d\n",index_i,index_j);

  imax   = shells->imax_sh[index_i];
  jmax   = shells->imax_sh[index_j];
  off2 = shells->ng_sh[index_i];
  off1 = (imax + 1) * off2;
  off4 = shells->ng_sh[index_j];
  off3 = (jmax + 1) * off4;
  sheli = shells->cart[index_i];
  shelj = shells->cart[index_j];
  mm  = imax + jmax;
  mm4 = shells->ng_sh[index_j];
  mm3 = shells->ng_sh[index_i] * mm4;
  mm2 = (mm + 1) * mm3;
  mm1 = (mm + 1) * mm2;
  for (i = 0; i < sheli; i++) {
    for (j = 0; j < shelj; j++) {
      tmax =  shells->tuv[imax][i][0];
      umax =  shells->tuv[imax][i][1];
      vmax =  shells->tuv[imax][i][2];
      tpmax = shells->tuv[jmax][j][0];
      upmax = shells->tuv[jmax][j][1];
      vpmax = shells->tuv[jmax][j][2];
      //printf("%3d %3d %3d %3d    %3d %3d %3d\n",index_i,index_j,i,j,tmax,umax,vmax);
      n1 = shells->tuv[imax][i][0];
      n2 = shells->tuv[imax][i][1];
      n3 = shells->tuv[imax][i][2];
      n4 = shells->tuv[jmax][j][0];
      n5 = shells->tuv[jmax][j][1];
      n6 = shells->tuv[jmax][j][2];
      count = (bfposi + i) * nd2 + bfposj + j;
      for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
        for (t = 0; t <= tmax; t++) {
          for (u = 0; u <= umax; u++) {
            for (v = 0; v <= vmax; v++) {
              c1fac[t][u][v] = Complex(k_zero, k_zero) ;
             }
            }
           }
            for (t = 0; t <= tmax; t++) {
              for (u = 0; u <= umax; u++) {
                for (v = 0; v <= vmax; v++) {
                  for (tp = 0; tp <= tpmax; tp++) {
                    for (up = 0; up <= upmax; up++) {
                      for (vp = 0; vp <= vpmax; vp++) {
                        tpupvpsign = -k_one;
                        tpupvp2 = tp + up + vp + 2;
                        if ((tpupvp2 / 2) * 2 == tpupvp2)
                          tpupvpsign = k_one;
                          for (j4 = 0; j4 < shells->ng_sh[index_j]; j4++) {
                            //p_fgtuv = fgtuv + (t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4;
                            c1fac[t][u][v] += C2x[n4 * off3 + tp * off4 + j4] * C2y[n5 * off3 + up * off4 + j4] * \
                            C2z[n6 * off3 + vp * off4 + j4] * tpupvpsign * \
                            fgtuv[(t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4];
                            //fprintf(file.out,"MM %14.8lf %14.8lf\n",\
                           (fgtuv[(t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4]).real(),\
                           (fgtuv[(t + tp) * mm1 + (u + up) * mm2 + (v + vp) * mm3 + i4 * mm4 + j4]).imag());
                            //c1fac[t][u][v] += C2x[n4 * off3 + tp * off4 + j4] * C2y[n5 * off3 + up * off4 + j4] * \
                            C2z[n6 * off3 + vp * off4 + j4] * *p_fgtuv * tpupvpsign;
                           } // end j4 loop
                          }
                         }
                        } // end tp up vp loop
                       }
                      }
                     } // end t u v loop
                      for (t = 0; t <= tmax; t++) {
                        for (u = 0; u <= umax; u++) {
                          for (v = 0; v <= vmax; v++) {
                            F_cart[count] += c1fac[t][u][v] * C1x[n1 * off1 + t * off2 + i4] * C1y[n2 * off1 + u * off2 + i4] * \
                            C1z[n3 * off1 + v * off2 + i4];
//printf("c1fac %3d %3d %3d %10.4lf %10.4lf\n",t,u,v,c1fac[t][u][v],F_cart[i * shelj + j]);
//j4 = 0;
//printf("%3d %3d %3d %10.4lf %10.4lf %10.4lf      %10.4lf\n",t,u,v,\
C1x[n1 * off1 + t * off2 + i4], C1y[n2 * off1 + u * off2 + i4],C1z[n3 * off1 + v * off2 + i4], c1fac[t][u][v]);
                           }
                          }
                         } // end t u v loop
                        } // end i4 loop
//fprintf(file.out,"%10.4lf ",F_cart[count]);
//printf("c2 %3d %3d %3d %10.4lf\n",bfposi,bfposj,count,F_cart[count]);
                       }
//fprintf(file.out,"\n");
                      }

}

void integrals_molecule_screen(double *F_sh, int ip, int jp, REAL_LATTICE *R, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, SYMMETRY *symmetry, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  // compute 4-centre integrals for a quad of atoms (ij|ij) using pair_p as indices
  // starting point is pair_p->posn[*p1]

  int index_i, index_j, i4, j4;
  int bfposi, bfposj;
  int bfposi1, bfposj1;
  int gausposi, gausposj;
  int count;
  int i, j, k, l, m, n, p, pm, q;
  int mm, mm0, mm1, mm2, mm3, mm4;
  int gi, gj;
  int imax, jmax;
  int shelposi, shelposj;
  int sheli, shelj;
  int sheli1, shelj1;
  int offset1, offset2;
  int dim1 = atoms->number_of_atoms_in_unit_cell;
  int dim2 = dim1 * dim1;
  int nd1 = atoms->bfnnumb[ip];
  int nd2 = atoms->bfnnumb[jp];
  int nd3 = atoms->bfnnumb_sh[ip];
  int nd4 = atoms->bfnnumb_sh[jp];
  double p_inv_ex, fac1;
  double R_AB_1esqrd, Rsqrd;
  double C1_max, C_max;
  double fgtuv_max, fgtuv_temp;
  double time1, time2;
  double em[1][55];
  double f[13][13][13][13];
  double *F_cart, *p_F_cart, *p_F_sh, *fgtuv, *p_fgtuv;
  VECTOR_DOUBLE R_AB_1e, s_12;

  time1 = MPI_Wtime();

    gi = 0;
    gj = 0;

    F_cart = (double *) calloc(nd1 * nd2 * nd1 * nd2, sizeof(double));
    if (F_cart == NULL) {
    fprintf(stderr, "ERROR: not enough memory for double F_cart! \n");
    exit(1);
   }
    
    R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1;
    R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2;
    R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3;
    R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);

    bfposi   = 0;
    bfposi1  = 0;
    shelposi = atoms->shelposn_sh[ip];
    gausposi = atoms->gausposn_sh[ip];
    for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
      sheli  = shells->type1_sh[index_i];
      sheli1 = shells->type_sh[index_i];
      imax   = shells->imax_sh[index_i];
      bfposj   = 0;
      bfposj1  = 0;
      shelposj = atoms->shelposn_sh[jp];
      gausposj = atoms->gausposn_sh[jp];
      for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
        shelj  = shells->type1_sh[index_j];
        shelj1 = shells->type_sh[index_j];
        jmax   = shells->imax_sh[index_j];
        double sab[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double C1x[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
        double C1y[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
        double C1z[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
        double pab_inv[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        VECTOR_DOUBLE R_AB[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        E_coefficients(ip,jp,gi,gj,index_i,index_j,gausposi,gausposj,C1x,C1y,C1z,&C1_max,sab,pab_inv,&R_AB_1esqrd,R_AB,R,atoms,\
        shells,gaussians,job,file);

        mm  = imax + jmax + imax + jmax;
        mm4 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
        mm3 = shells->ng_sh[index_i] * shells->ng_sh[index_j] * mm4;
        mm2 = (mm + 1) * mm3;
        mm1 = (mm + 1) * mm2;
        mm0 = (mm + 1) * mm1;

        fgtuv = (double *) calloc(mm0, sizeof(double));
        if (fgtuv == NULL) {
        fprintf(stderr, "ERROR: not enough memory for double fgtuv! \n");
        exit(1);
       }

        fgtuv_max = k_zero;
        for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
          for (j4 = 0; j4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; j4++) {
            s_12.comp1 = R_AB[i4].comp1 - R_AB[j4].comp1;
            s_12.comp2 = R_AB[i4].comp2 - R_AB[j4].comp2;
            s_12.comp3 = R_AB[i4].comp3 - R_AB[j4].comp3;
            fac1 = sab[i4] * sab[j4];
            p_inv_ex = pab_inv[i4] + pab_inv[j4];
            Rsqrd = double_vec_dot(&s_12, &s_12);
            f000m(&em[0][0], Rsqrd / p_inv_ex, k_one / p_inv_ex, mm);
            //non_recursive_ftuvn(mm, 0, f, em, &s_12);
            for (m = 0; m <= mm; m++) {
              for (n = 0; n <= mm; n++) {
                for (pm = 0; pm <= mm; pm++) {
                  if (m + n + pm > mm) break;
                  //CHANGES2015 CHANGEMEBACK
                  fgtuv_temp = fac1 * ftuvn(m,n,pm,0,&em[0][0],s_12);
                  //fgtuv_temp = fac1 * f[m][n][pm][0];
                  p_fgtuv = fgtuv + m * mm1 + n * mm2 + pm * mm3  + i4 * mm4 + j4;
                 *p_fgtuv += fgtuv_temp;
                  fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
                  fgtuv_max = k_one;
                 }
                }
               }
              }
             }
        mcmurchie_davidson_screen(F_cart, index_i, index_j, bfposi, bfposj, nd2, C1x, C1y, C1z, fgtuv, shells, job, file);
        free(fgtuv);
        two_center_cartesian_to_sh_shell1(F_cart,F_sh,index_i,index_j,bfposi,bfposj,bfposi1,bfposj1,nd2,nd4,shells,job,file);
        bfposj   += shelj;
        bfposj1  += shelj1;
        gausposj += shells->ng_sh[index_j];
       } // close loop over index_j
      bfposi   += sheli;
      bfposi1  += sheli1;
      gausposi += shells->ng_sh[index_i];
     } // close loop over index_i
//for(i=0;i<nd1;i++) { for(j=0;j<nd2;j++) { for(k=0;k<nd1;k++) { for(l=0;l<nd2;l++) { fprintf(file.out,"%3d %3d %3d %3d %10.4lf\n",\
i,j,k,l,F_cart[i*nd2*nd1*nd2+j*nd1*nd2+k*nd2+l]);}}}}
//for(i=0;i<nd3;i++) { for(j=0;j<nd4;j++) { for(k=0;k<nd3;k++) { for(l=0;l<nd4;l++) { fprintf(file.out,"%3d %3d %3d %3d %10.4lf\n",\
i,j,k,l,F_sh[i*nd4*nd3*nd4+j*nd3*nd4+k*nd4+l]);}}}}

      free(F_cart);

}

void integrals_crystal_exchange_screen(double *F_sh, int ip, int jp, int gj, REAL_LATTICE *R, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, SYMMETRY *symmetry, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  // compute 4-centre integrals for a quad of atoms (ij|ij) using pair_p as indices
  // starting point is pair_p->posn[*p1]

  int index_i, index_j, i4, j4;
  int bfposi, bfposj;
  int bfposi1, bfposj1;
  int gausposi, gausposj;
  int count;
  int i, j, k, l, m, n, p, pm, q;
  int mm, mm0, mm1, mm2, mm3, mm4;
  int gi;
  //int gi, gj;
  int imax, jmax;
  int shelposi, shelposj;
  int sheli, shelj;
  int sheli1, shelj1;
  int offset1, offset2;
  int dim1 = atoms->number_of_atoms_in_unit_cell;
  int dim2 = dim1 * dim1;
  int nd1 = atoms->bfnnumb[ip];
  int nd2 = atoms->bfnnumb[jp];
  int nd3 = atoms->bfnnumb_sh[ip];
  int nd4 = atoms->bfnnumb_sh[jp];
  double p_inv_ex, fac1;
  double R_AB_1esqrd, Rsqrd;
  double C1_max, C_max;
  double fgtuv_max, fgtuv_temp;
  double time1, time2;
  double em[1][55];
  double f[13][13][13][13];
  double *F_cart, *p_F_cart, *p_F_sh, *fgtuv, *p_fgtuv;
  VECTOR_DOUBLE R_AB_1e, s_12;

  time1 = MPI_Wtime();

    gi = 0;
    //gj = 0;

    F_cart = (double *) calloc(nd1 * nd2 * nd1 * nd2, sizeof(double));
    if (F_cart == NULL) {
    fprintf(stderr, "ERROR: not enough memory for double F_cart! \n");
    exit(1);
   }
    
    R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1 - R->vec_ai[gj].comp1;
    R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2 - R->vec_ai[gj].comp2;
    R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3 - R->vec_ai[gj].comp3;
    R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);

    bfposi   = 0;
    bfposi1  = 0;
    shelposi = atoms->shelposn_sh[ip];
    gausposi = atoms->gausposn_sh[ip];
    for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
      sheli  = shells->type1_sh[index_i];
      sheli1 = shells->type_sh[index_i];
      imax   = shells->imax_sh[index_i];
      bfposj   = 0;
      bfposj1  = 0;
      shelposj = atoms->shelposn_sh[jp];
      gausposj = atoms->gausposn_sh[jp];
      for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
        shelj  = shells->type1_sh[index_j];
        shelj1 = shells->type_sh[index_j];
        jmax   = shells->imax_sh[index_j];
        double sab[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double C1x[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
        double C1y[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
        double C1z[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
        double pab_inv[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        VECTOR_DOUBLE R_AB[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        E_coefficients(ip,jp,gi,gj,index_i,index_j,gausposi,gausposj,C1x,C1y,C1z,&C1_max,sab,pab_inv,&R_AB_1esqrd,R_AB,R,atoms,\
        shells,gaussians,job,file);

        mm  = imax + jmax + imax + jmax;
        mm4 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
        mm3 = shells->ng_sh[index_i] * shells->ng_sh[index_j] * mm4;
        mm2 = (mm + 1) * mm3;
        mm1 = (mm + 1) * mm2;
        mm0 = (mm + 1) * mm1;

        fgtuv = (double *) calloc(mm0, sizeof(double));
        if (fgtuv == NULL) {
        fprintf(stderr, "ERROR: not enough memory for double fgtuv! \n");
        exit(1);
       }
        fgtuv_max = k_zero;
        for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
          for (j4 = 0; j4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; j4++) {
            s_12.comp1 = R_AB[i4].comp1 - R_AB[j4].comp1;
            s_12.comp2 = R_AB[i4].comp2 - R_AB[j4].comp2;
            s_12.comp3 = R_AB[i4].comp3 - R_AB[j4].comp3;
            fac1 = sab[i4] * sab[j4];
            p_inv_ex = pab_inv[i4] + pab_inv[j4];
            Rsqrd = double_vec_dot(&s_12, &s_12);
            f000m(&em[0][0], Rsqrd / p_inv_ex, k_one / p_inv_ex, mm);
            //non_recursive_ftuvn(mm, 0, f, em, &s_12);
            for (m = 0; m <= mm; m++) {
              for (n = 0; n <= mm; n++) {
                for (pm = 0; pm <= mm; pm++) {
                  if (m + n + pm > mm) break;
                  //CHANGES2015 CHANGEMEBACK
                  fgtuv_temp = fac1 * ftuvn(m,n,pm,0,&em[0][0],s_12);
                  //fgtuv_temp = fac1 * f[m][n][pm][0];
                  p_fgtuv = fgtuv + m * mm1 + n * mm2 + pm * mm3  + i4 * mm4 + j4;
                 *p_fgtuv += fgtuv_temp;
                  fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
                  fgtuv_max = k_one;
                 }
                }
               }
              }
             }
        mcmurchie_davidson_screen(F_cart, index_i, index_j, bfposi, bfposj, nd2, C1x, C1y, C1z, fgtuv, shells, job, file);
        free(fgtuv);
        two_center_cartesian_to_sh_shell1(F_cart,F_sh,index_i,index_j,bfposi,bfposj,bfposi1,bfposj1,nd2,nd4,shells,job,file);
        bfposj   += shelj;
        bfposj1  += shelj1;
        gausposj += shells->ng_sh[index_j];
       } // close loop over index_j
      bfposi   += sheli;
      bfposi1  += sheli1;
      gausposi += shells->ng_sh[index_i];
     } // close loop over index_i
//for(i=0;i<nd1;i++) { for(j=0;j<nd2;j++) { for(k=0;k<nd1;k++) { for(l=0;l<nd2;l++) { fprintf(file.out,"%3d %3d %3d %3d %10.4lf\n",\
i,j,k,l,F_cart[i*nd2*nd1*nd2+j*nd1*nd2+k*nd2+l]);}}}}
//for(i=0;i<nd3;i++) { for(j=0;j<nd4;j++) { for(k=0;k<nd3;k++) { for(l=0;l<nd4;l++) { fprintf(file.out,"%3d %3d %3d %3d %10.4lf\n",\
i,j,k,l,F_sh[i*nd4*nd3*nd4+j*nd3*nd4+k*nd4+l]);}}}}

      free(F_cart);

}

void integrals_crystal_screen(double *F_sh, int ip, int jp, int gj, REAL_LATTICE *R, RECIPROCAL_LATTICE *G, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, SYMMETRY *symmetry, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  // compute 4-centre integrals for a quad of atoms (ij|ij) using pair_p as indices
  // starting point is pair_p->posn[*p1]

  int index_i, index_j, index_G, index_R, index_S, i4, j4;
  int bfposi, bfposj;
  int bfposi1, bfposj1;
  int gausposi, gausposj;
  int count;
  int i, j, k, l, m, n, p, pm, q;
  int mm, mm0, mm1, mm2, mm3, mm4;
  int i1, ijk;
  int gi;
  //int gi, gj;
  int imax, jmax;
  int shelposi, shelposj;
  int sheli, shelj;
  int sheli1, shelj1;
  int offset1, offset2;
  int dim1 = atoms->number_of_atoms_in_unit_cell;
  int dim2 = dim1 * dim1;
  int nd1 = atoms->bfnnumb[ip];
  int nd2 = atoms->bfnnumb[jp];
  int nd3 = atoms->bfnnumb_sh[ip];
  int nd4 = atoms->bfnnumb_sh[jp];
  double p_inv_ex, fac1, fac2, dot_product, shell_sum;
  double gamma_1;
  double R_AB_1esqrd, Rsqrd;
  double C1_max, C_max;
  double fgtuv_max, fgtuv_temp;
  double time1, time2;
  double em[R->last_ewald_vector][55];
  double en[R->last_ewald_vector][55];
  double pi_vol = pi / crystal->primitive_cell_volume;
  double gamma_1_inv = G->gamma_0_inv;
  //double f[13][13][13][13];
  double *F_cart, *p_F_cart, *p_F_sh, *fgtuv, *p_fgtuv;
  VECTOR_DOUBLE R_AB_1e, r_12, s_12, t_12, Rvec_tmp;

  gamma_1 = k_one / G->gamma_0_inv;

  time1 = MPI_Wtime();

    gi = 0;
    //gj = 0;

    F_cart = (double *) calloc(nd1 * nd2 * nd1 * nd2, sizeof(double));
    if (F_cart == NULL) {
    fprintf(stderr, "ERROR: not enough memory for double F_cart! \n");
    exit(1);
   }
    
    R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1 - R->vec_ai[gj].comp1;
    R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2 - R->vec_ai[gj].comp2;
    R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3 - R->vec_ai[gj].comp3;
    R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);

    //fprintf(file.out,"screen %3d %3d  %3d %10.4lf\n",ip,jp,gj,sqrt(R_AB_1esqrd));

    bfposi   = 0;
    bfposi1  = 0;
    shelposi = atoms->shelposn_sh[ip];
    gausposi = atoms->gausposn_sh[ip];
    for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
      sheli  = shells->type1_sh[index_i];
      sheli1 = shells->type_sh[index_i];
      imax   = shells->imax_sh[index_i];
      bfposj   = 0;
      bfposj1  = 0;
      shelposj = atoms->shelposn_sh[jp];
      gausposj = atoms->gausposn_sh[jp];
      for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
        shelj  = shells->type1_sh[index_j];
        shelj1 = shells->type_sh[index_j];
        jmax   = shells->imax_sh[index_j];
        double sab[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double C1x[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
        double C1y[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
        double C1z[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
        double pab_inv[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        VECTOR_DOUBLE R_AB[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        E_coefficients(ip,jp,gi,gj,index_i,index_j,gausposi,gausposj,C1x,C1y,C1z,&C1_max,sab,pab_inv,&R_AB_1esqrd,R_AB,R,atoms,\
        shells,gaussians,job,file);

        mm  = imax + jmax + imax + jmax;
        mm4 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
        mm3 = shells->ng_sh[index_i] * shells->ng_sh[index_j] * mm4;
        mm2 = (mm + 1) * mm3;
        mm1 = (mm + 1) * mm2;
        mm0 = (mm + 1) * mm1;

        fgtuv = (double *) calloc(mm0, sizeof(double));
        if (fgtuv == NULL) {
        fprintf(stderr, "ERROR: not enough memory for double fgtuv! \n");
        exit(1);
       }

        fgtuv_max = k_zero;

         switch (crystal->type[0]) {

           case 'C':

            fgtuv_max = k_zero;
        for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
          for (j4 = 0; j4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; j4++) {
            s_12.comp1 = R_AB[i4].comp1 - R_AB[j4].comp1;
            s_12.comp2 = R_AB[i4].comp2 - R_AB[j4].comp2;
            s_12.comp3 = R_AB[i4].comp3 - R_AB[j4].comp3;
            fac1 = sab[i4] * sab[j4];
            //fprintf(file.out,"i4 %3d j4 %3d  %10.2e %10.2e %10.2e\n",i4,j4,sab[i4],sab[j4],fac1);
            p_inv_ex = pab_inv[i4] + pab_inv[j4];
            p_fgtuv = fgtuv + i4 * mm4 + j4;
           *p_fgtuv -= pi_vol * fac1 * (gamma_1_inv - pab_inv[i4] - pab_inv[j4]);
            count = 1;
            for (index_S = 1; index_S < G->number_of_shells; index_S++) {
              fac2 = fac1 * G->EXPFAC[index_S];
              //double test = G->EXPFAC[index_S] * fac1 * \
              (k_one > G->shell_mag[index_S * 9 + mm] ? k_one : G->shell_mag[index_S * 9 + mm]);
              //if (fabs(test) * C_max < 1.0e-18) break;
              //fac2 = fac1 * expfac * exp(-G->sqr[index_S] / four * p_inv_ex) / G->sqr[index_S];
              for (index_G = 0; index_G < G->num[index_S]; index_G++) {
                dot_product = G->vec_b2[count].comp1 * s_12.comp1 + G->vec_b2[count].comp2 * s_12.comp2 + \
                G->vec_b2[count].comp3 * s_12.comp3;
                for (i = 0; i <= mm; i++) {
                  for (j = 0; j <= mm; j++) {
                    for (k = 0; k <= mm; k++) {
                      ijk = i + j + k;
                      if (ijk > mm) break;
                      fgtuv_temp = fac2 * G->x[count * 9 + i] * G->y[count * 9 + j] * G->z[count * 9 + k] * \
                      cosfactor(ijk, dot_product);
                      p_fgtuv = fgtuv + i * mm1 + j * mm2 + k * mm3  + i4 * mm4 + j4;
                     *p_fgtuv += fgtuv_temp;
                      fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
                      fgtuv_max = k_one;
                     }
                    }
                   }
                  count++;
                 }
                }
               }
              }

       break;

       case 'S':
       case 'P':

       break;

       case 'M':

       if (job->taskid == 0)
       fprintf(file.out,"integrals_crystal_screen routine is for periodic systems only\n");
       MPI_Finalize();
       exit(1);

       break;

      } // close switch

  for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
    for (j4 = 0; j4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; j4++) {
      s_12.comp1 = R_AB[i4].comp1 - R_AB[j4].comp1;
      s_12.comp2 = R_AB[i4].comp2 - R_AB[j4].comp2;
      s_12.comp3 = R_AB[i4].comp3 - R_AB[j4].comp3;
      map_to_wigner(crystal,&s_12, &t_12, &Rvec_tmp); // map s_12 back to primitive cell to accelerate convergence
      fac1 = sab[i4] * sab[j4];
      //fprintf(file.out,"i4 %3d j4 %3d  %10.2e %10.2e %10.2e\n",i4,j4,sab[i4],sab[j4],fac1);
      p_inv_ex = pab_inv[i4] + pab_inv[j4];
      count = 0;
      for (index_S = 0; index_S < R->number_of_ewald_shells; index_S++) {
        shell_sum = k_zero;
        for (index_R = 0; index_R < R->num[index_S]; index_R++) {
          r_12.comp1 = t_12.comp1 + R->vec_ai[count].comp1;
          r_12.comp2 = t_12.comp2 + R->vec_ai[count].comp2;
          r_12.comp3 = t_12.comp3 + R->vec_ai[count].comp3;
          Rsqrd = double_vec_dot(&r_12, &r_12);
          f000m(&em[count][0], gamma_1 * Rsqrd, gamma_1, mm);
          f000m(&en[count][0], Rsqrd / p_inv_ex, k_one / p_inv_ex, mm);
          for (i1 = 0; i1 <= mm; i1++) {
            en[count][i1] -= em[count][i1];
           }
            shell_sum += fabs(en[count][0]);
            count++;
           } // end loop on index_R
             } // end loop on index_S
              for (index_R = 0; index_R < count; index_R++) {
                r_12.comp1 = t_12.comp1 + R->vec_ai[index_R].comp1;
                r_12.comp2 = t_12.comp2 + R->vec_ai[index_R].comp2;
                r_12.comp3 = t_12.comp3 + R->vec_ai[index_R].comp3;
            for (m = 0; m <= mm; m++) {
              for (n = 0; n <= mm; n++) {
                for (pm = 0; pm <= mm; pm++) {
                  if (m + n + pm > mm) break;
                  fgtuv_temp = fac1 * ftuvn(m,n,pm,0,&en[index_R][0],r_12);
                  p_fgtuv = fgtuv + m * mm1 + n * mm2 + pm * mm3  + i4 * mm4 + j4;
                 *p_fgtuv += fgtuv_temp;
                  fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
                  fgtuv_max = k_one;
                 }
                }
               }
              }
             }
            }
        mcmurchie_davidson_screen(F_cart, index_i, index_j, bfposi, bfposj, nd2, C1x, C1y, C1z, fgtuv, shells, job, file);
        free(fgtuv);
        two_center_cartesian_to_sh_shell1(F_cart,F_sh,index_i,index_j,bfposi,bfposj,bfposi1,bfposj1,nd2,nd4,shells,job,file);
        bfposj   += shelj;
        bfposj1  += shelj1;
        gausposj += shells->ng_sh[index_j];
       } // close loop over index_j
      bfposi   += sheli;
      bfposi1  += sheli1;
      gausposi += shells->ng_sh[index_i];
     } // close loop over index_i
//for(i=0;i<nd1;i++) { for(j=0;j<nd2;j++) { for(k=0;k<nd1;k++) { for(l=0;l<nd2;l++) { fprintf(file.out,"%3d %3d %3d %3d %10.4lf\n",\
i,j,k,l,F_cart[i*nd2*nd1*nd2+j*nd1*nd2+k*nd2+l]);}}}}
//for(i=0;i<nd3;i++) { for(j=0;j<nd4;j++) { for(k=0;k<nd3;k++) { for(l=0;l<nd4;l++) { fprintf(file.out,"%3d %3d %3d %3d %10.4lf\n",\
i,j,k,l,F_sh[i*nd4*nd3*nd4+j*nd3*nd4+k*nd4+l]);}}}}

      free(F_cart);

}

void integrals_crystal_screen_complex(Complex *F_sh, int ip, int jp, int gj, REAL_LATTICE *R, RECIPROCAL_LATTICE *G, Q_LATTICE *q_G, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, SYMMETRY *symmetry, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  // compute 4-centre integrals for a quad of atoms (ij|ij) using pair_p as indices
  // starting point is pair_p->posn[*p1]

  int index_i, index_j, index_G, index_R, index_S, i4, j4;
  int bfposi, bfposj;
  int bfposi1, bfposj1;
  int gausposi, gausposj;
  int count;
  int i, j, k, l, m, n, p, pm, q;
  int mm, mm0, mm1, mm2, mm3, mm4;
  int i1, ijk;
  int gi;
  int imax, jmax;
  int shelposi, shelposj;
  int sheli, shelj;
  int sheli1, shelj1;
  int offset1, offset2;
  int dim1 = atoms->number_of_atoms_in_unit_cell;
  int dim2 = dim1 * dim1;
  int nd1 = atoms->bfnnumb[ip];
  int nd2 = atoms->bfnnumb[jp];
  int nd3 = atoms->bfnnumb_sh[ip];
  int nd4 = atoms->bfnnumb_sh[jp];
  double p_inv_ex, fac1, fac2, expfac, dot_product, shell_sum;
  double gamma_1;
  double R_AB_1esqrd, Rsqrd;
  double C1_max, C_max;
  //double fgtuv_max, fgtuv_temp;
  Complex fac, fgtuv_max, fgtuv_temp;
  double time1, time2;
  double em[R->last_ewald_vector][55];
  double en[R->last_ewald_vector][55];
  double pi_vol = pi / crystal->primitive_cell_volume;
  double gamma_1_inv = G->gamma_0_inv;
  Complex *F_cart, *p_F_cart, *p_F_sh, *fgtuv, *p_fgtuv;
  //double *F_cart, *p_F_cart, *p_F_sh, *fgtuv, *p_fgtuv;
  VECTOR_DOUBLE R_AB_1e, r_12, s_12, t_12, Rvec_tmp;

  int size;
  size = 4 * job->l_max + 1;
  gamma_1 = k_one / G->gamma_0_inv;

  time1 = MPI_Wtime();

    gi = 0;

    F_cart = (Complex *) calloc(nd1 * nd2 * nd1 * nd2, sizeof(Complex));
    if (F_cart == NULL) {
    fprintf(stderr, "ERROR: not enough memory for Complex F_cart! \n");
    exit(1);
   }
    
    R_AB_1e.comp1 = atoms->cell_vector[ip].comp1 - atoms->cell_vector[jp].comp1 - R->vec_ai[gj].comp1;
    R_AB_1e.comp2 = atoms->cell_vector[ip].comp2 - atoms->cell_vector[jp].comp2 - R->vec_ai[gj].comp2;
    R_AB_1e.comp3 = atoms->cell_vector[ip].comp3 - atoms->cell_vector[jp].comp3 - R->vec_ai[gj].comp3;
    R_AB_1esqrd = double_vec_dot(&R_AB_1e, &R_AB_1e);

    //fprintf(file.out,"screen %3d %3d  %3d %10.4lf\n",ip,jp,gj,sqrt(R_AB_1esqrd));

    bfposi   = 0;
    bfposi1  = 0;
    shelposi = atoms->shelposn_sh[ip];
    gausposi = atoms->gausposn_sh[ip];
    for (index_i = shelposi; index_i < shelposi + atoms->nshel_sh[ip]; index_i++) {
      sheli  = shells->type1_sh[index_i];
      sheli1 = shells->type_sh[index_i];
      imax   = shells->imax_sh[index_i];
      bfposj   = 0;
      bfposj1  = 0;
      shelposj = atoms->shelposn_sh[jp];
      gausposj = atoms->gausposn_sh[jp];
      for (index_j = shelposj; index_j < shelposj + atoms->nshel_sh[jp]; index_j++) {
        shelj  = shells->type1_sh[index_j];
        shelj1 = shells->type_sh[index_j];
        jmax   = shells->imax_sh[index_j];
        double sab[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        double C1x[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
        double C1y[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
        double C1z[shells->ng_sh[index_i] * shells->ng_sh[index_j] * (imax+jmax+1) * (imax+1) * (jmax+1)];
        double pab_inv[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        VECTOR_DOUBLE R_AB[shells->ng_sh[index_i] * shells->ng_sh[index_j]];
        E_coefficients(ip,jp,gi,gj,index_i,index_j,gausposi,gausposj,C1x,C1y,C1z,&C1_max,sab,pab_inv,&R_AB_1esqrd,R_AB,R,atoms,\
        shells,gaussians,job,file);

        mm  = imax + jmax + imax + jmax;
        mm4 = shells->ng_sh[index_i] * shells->ng_sh[index_j];
        mm3 = shells->ng_sh[index_i] * shells->ng_sh[index_j] * mm4;
        mm2 = (mm + 1) * mm3;
        mm1 = (mm + 1) * mm2;
        mm0 = (mm + 1) * mm1;

        fgtuv = (Complex *) calloc(mm0, sizeof(Complex));
        if (fgtuv == NULL) {
        fprintf(stderr, "ERROR: not enough memory for Complex fgtuv! \n");
        exit(1);
       }

        fgtuv_max = Complex(k_zero, k_zero);

         switch (crystal->type[0]) {

           case 'C':

          expfac = four * pi / crystal->primitive_cell_volume;
          for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
          for (j4 = 0; j4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; j4++) {
            s_12.comp1 = R_AB[i4].comp1 - R_AB[j4].comp1;
            s_12.comp2 = R_AB[i4].comp2 - R_AB[j4].comp2;
            s_12.comp3 = R_AB[i4].comp3 - R_AB[j4].comp3;
            fac1 = sab[i4] * sab[j4];
           //for (j4 = 0; j4 < shells_ax->ng_sh[index_k]; j4++) {
             //s_12.comp1 = R_AB[i4].comp1 - atoms_ax->cell_vector[kp].comp1;
             //s_12.comp2 = R_AB[i4].comp2 - atoms_ax->cell_vector[kp].comp2;
             //s_12.comp3 = R_AB[i4].comp3 - atoms_ax->cell_vector[kp].comp3;
             //p_inv_ex = pab_inv[i4] + pc_inv[j4];
             //fac1 = sab[i4] * sc[j4];
             p_fgtuv = fgtuv + i4 * mm4 + j4;
             if (fabs(q_G->sqr[0]) < 0.00001)
            *p_fgtuv -= pi_vol * fac1 * (gamma_1_inv - pab_inv[i4] - pab_inv[j4]);
             //fprintf(file.out,"pi_vol 3C %f\n",\
             -pi_vol * fac1 * (gamma_1_inv - pab_inv[i4] - pab_inv[j4]));
             int start = 0;
             if (fabs(q_G->sqr[0]) < 0.00001) start = 1;
             for (index_G = start; index_G < q_G->last_vector; index_G++) {
               fac2 = fac1 * expfac * exp(-q_G->sqr[index_G] / four * q_G->gamma_0_inv) / q_G->sqr[index_G];
               //uncomment for G only 
               //fac2 = fac1 * expfac * exp(-q_G->sqr[index_G] / four * (pab_inv[i4] + pc_inv[j4]) ) / q_G->sqr[index_G];
               dot_product = q_G->vec[index_G].comp1 * s_12.comp1 + q_G->vec[index_G].comp2 * s_12.comp2 + \
               q_G->vec[index_G].comp3 * s_12.comp3;
               for (i = 0; i <= mm; i++) {
                 for (j = 0; j <= mm; j++) {
                   for (k = 0; k <= mm; k++) {
                     ijk = i + j + k;
                     if (ijk > mm) break;
                     fgtuv_temp = fac2 * q_G->x[index_G * size + i] * q_G->y[index_G * size + j] * q_G->z[index_G * size + k] * \
                     conj(cosfactor_complex(ijk, dot_product)); // use complex conjugate for e^-iq.r
             //fgtuv_temp = Complex(k_zero,k_zero);
                     //fprintf(file.out,"fgtuv %f\n", fgtuv_temp.real());
                     p_fgtuv = fgtuv + i * mm1 + j * mm2 + k * mm3  + i4 * mm4 + j4;
                    *p_fgtuv += fgtuv_temp;
                    }
                   }
                  }
                 }
                }
               }
                 //for (i=0;i<mm0;i++) fprintf(file.out,"2C complex FGTUV0 %3d  %3d %3d %3d %10.4f %10.4f\n",\
                 i,index_i,index_j,index_k,fgtuv[i].real(),fgtuv[i].imag());
                 //fprintf(file.out,"\n");
                 //for (i=0;i<mm0;i++) fgtuv[i] = Complex(k_zero, k_zero);


 //           fgtuv_max = k_zero;
 //       for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
 //         for (j4 = 0; j4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; j4++) {
 //           s_12.comp1 = R_AB[i4].comp1 - R_AB[j4].comp1;
 //           s_12.comp2 = R_AB[i4].comp2 - R_AB[j4].comp2;
 //           s_12.comp3 = R_AB[i4].comp3 - R_AB[j4].comp3;
 //           fac1 = sab[i4] * sab[j4];
 //           //fprintf(file.out,"i4 %3d j4 %3d  %10.2e %10.2e %10.2e\n",i4,j4,sab[i4],sab[j4],fac1);
 //           p_inv_ex = pab_inv[i4] + pab_inv[j4];
 //           p_fgtuv = fgtuv + i4 * mm4 + j4;
 //          *p_fgtuv -= pi_vol * fac1 * (gamma_1_inv - pab_inv[i4] - pab_inv[j4]);
 //           count = 1;
 //           for (index_S = 1; index_S < G->number_of_shells; index_S++) {
 //             fac2 = fac1 * G->EXPFAC[index_S];
 //             //double test = G->EXPFAC[index_S] * fac1 * \
 //             (k_one > G->shell_mag[index_S * 9 + mm] ? k_one : G->shell_mag[index_S * 9 + mm]);
 //             //if (fabs(test) * C_max < 1.0e-18) break;
 //             //fac2 = fac1 * expfac * exp(-G->sqr[index_S] / four * p_inv_ex) / G->sqr[index_S];
 //             for (index_G = 0; index_G < G->num[index_S]; index_G++) {
 //               dot_product = G->vec_b2[count].comp1 * s_12.comp1 + G->vec_b2[count].comp2 * s_12.comp2 + \
 //               G->vec_b2[count].comp3 * s_12.comp3;
 //               for (i = 0; i <= mm; i++) {
 //                 for (j = 0; j <= mm; j++) {
 //                   for (k = 0; k <= mm; k++) {
 //                     ijk = i + j + k;
 //                     if (ijk > mm) break;
 //                     fgtuv_temp = fac2 * G->x[count * 9 + i] * G->y[count * 9 + j] * G->z[count * 9 + k] * \
 //                     cosfactor(ijk, dot_product);
 //                     p_fgtuv = fgtuv + i * mm1 + j * mm2 + k * mm3  + i4 * mm4 + j4;
 //                    *p_fgtuv += fgtuv_temp;
 //                     fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
 //                     fgtuv_max = k_one;
 //                    }
 //                   }
 //                  }
 //                 count++;
 //                }
 //               }
 //              }
 //             }
 //
 
       break;
 
       case 'S':
       case 'P':
 
       break;

       case 'M':

       if (job->taskid == 0)
       fprintf(file.out,"integrals_crystal_screen_complex routine is for periodic systems only\n");
       MPI_Finalize();
       exit(1);

       break;

      } // close switch

  for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
    for (j4 = 0; j4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; j4++) {
      s_12.comp1 = R_AB[i4].comp1 - R_AB[j4].comp1;
      s_12.comp2 = R_AB[i4].comp2 - R_AB[j4].comp2;
      s_12.comp3 = R_AB[i4].comp3 - R_AB[j4].comp3;
      map_to_wigner(crystal,&s_12, &t_12, &Rvec_tmp); // map s_12 back to primitive cell to accelerate convergence
      fac1 = sab[i4] * sab[j4];
      //fprintf(file.out,"i4 %3d j4 %3d  %10.2e %10.2e %10.2e\n",i4,j4,sab[i4],sab[j4],fac1);
      p_inv_ex = pab_inv[i4] + pab_inv[j4];
      count = 0;

    //for (j4 = 0; j4 < shells_ax->ng_sh[index_k]; j4++) {
      //s_12.comp1 = R_AB[i4].comp1 - atoms_ax->cell_vector[kp].comp1;
      //s_12.comp2 = R_AB[i4].comp2 - atoms_ax->cell_vector[kp].comp2;
      //s_12.comp3 = R_AB[i4].comp3 - atoms_ax->cell_vector[kp].comp3;
      //p_inv_ex = pab_inv[i4] + pc_inv[j4];
      //map_to_wigner(crystal,&s_12, &t_12, &Rvec_tmp); // map s_12 back to primitive cell to accelerate convergence
      //fac1 = sab[i4] * sc[j4];
      //count = 0;
      for (index_S = 0; index_S < R->number_of_ewald_shells; index_S++) {
        shell_sum = k_zero;
        for (index_R = 0; index_R < R->num[index_S]; index_R++) {
          //r_12.comp1 = t_12.comp1 + R->vec_ai[count].comp1;
          //r_12.comp2 = t_12.comp2 + R->vec_ai[count].comp2;
          //r_12.comp3 = t_12.comp3 + R->vec_ai[count].comp3;
          r_12.comp1 = s_12.comp1 + R->vec_ai[count].comp1;
          r_12.comp2 = s_12.comp2 + R->vec_ai[count].comp2;
          r_12.comp3 = s_12.comp3 + R->vec_ai[count].comp3;
          Rsqrd = r_12.comp1 * r_12.comp1 + r_12.comp2 * r_12.comp2 + r_12.comp3 * r_12.comp3;
          f000m(&em[count][0], gamma_1 * Rsqrd, gamma_1, mm);
          f000m(&en[count][0], Rsqrd / p_inv_ex, k_one / p_inv_ex, mm);
          for (i1 = 0; i1 <= mm; i1++) {
            en[count][i1] -= em[count][i1];
           }
            shell_sum += fabs(en[count][0]);
            count++;
           } // end loop on index_R
             } // end loop on index_S
              for (index_R = 0; index_R < count; index_R++) {
                //r_12.comp1 = t_12.comp1 + R->vec_ai[index_R].comp1;
                //r_12.comp2 = t_12.comp2 + R->vec_ai[index_R].comp2;
                //r_12.comp3 = t_12.comp3 + R->vec_ai[index_R].comp3;
                r_12.comp1 = s_12.comp1 + R->vec_ai[index_R].comp1;
                r_12.comp2 = s_12.comp2 + R->vec_ai[index_R].comp2;
                r_12.comp3 = s_12.comp3 + R->vec_ai[index_R].comp3;
                dot_product = q_G->vec[0].comp1 * R->vec_ai[index_R].comp1 + q_G->vec[0].comp2 * R->vec_ai[index_R].comp2 + \
                q_G->vec[0].comp3 * R->vec_ai[index_R].comp3;
                //fac = fac1 * Complex(cos(dot_product), -sin(dot_product));
                fac = fac1 * Complex(cos(dot_product), sin(dot_product));  // sign of sin(dot ??
                //non_recursive_ftuvn(mm, index_R, f, en, &r_12);
                //replace for G functions fgtuv_temp = sab[i4] * sc[j4] * ftuvn(m,n,pm,0,&em[0][0],s_12);
                for (i = 0; i <= mm; i++) {
                  for (j = 0; j <= mm; j++) {
                    for (k = 0; k <= mm; k++) {
                      ijk = i + j + k;
                      if (ijk > mm) break;
                      fgtuv_temp = fac * ftuvn(i,j,k,0,&en[index_R][0],r_12);
           //fgtuv_temp = Complex(k_zero,k_zero);
                      //fgtuv_temp = fac * f[i][j][k][0];
                      //fgtuv_temp = fac1 * f[i][j][k][0];
                      p_fgtuv  = fgtuv + i * mm1 + j * mm2 + k * mm3  + i4 * mm4 + j4;
                     *p_fgtuv += fgtuv_temp;
                     }
                    }
                   }
                  } 
                 } 
                } 
                 //for (i=0;i<mm0;i++) fprintf(file.out,"3C complex FGTUV1 %3d  %3d %3d %3d %10.4f %10.4f\n",\
                 i,index_i,index_j,index_k,fgtuv[i].real(),fgtuv[i].imag());
                 //for (i=0;i<mm0;i++) fprintf(file.out,"%3d  %3d %3d %3d %18.12f %18.12f\n",\
                 i,index_i,index_j,index_k,fgtuv[i].real(),fgtuv[i].imag());
                 //for (i=0;i<mm0;i++) fprintf(file.out,"%3d  %3d %3d %3d %10.4f %10.4f\n",\
                 i,index_i,index_j,index_k,fgtuv[i].real(),fgtuv[i].imag());
                 //fprintf(file.out,"\n");

 // for (i4 = 0; i4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; i4++) {
 //   for (j4 = 0; j4 < shells->ng_sh[index_i] * shells->ng_sh[index_j]; j4++) {
 //     s_12.comp1 = R_AB[i4].comp1 - R_AB[j4].comp1;
 //     s_12.comp2 = R_AB[i4].comp2 - R_AB[j4].comp2;
 //     s_12.comp3 = R_AB[i4].comp3 - R_AB[j4].comp3;
 //     map_to_wigner(crystal,&s_12, &t_12, &Rvec_tmp); // map s_12 back to primitive cell to accelerate convergence
 //     fac1 = sab[i4] * sab[j4];
 //     //fprintf(file.out,"i4 %3d j4 %3d  %10.2e %10.2e %10.2e\n",i4,j4,sab[i4],sab[j4],fac1);
 //     p_inv_ex = pab_inv[i4] + pab_inv[j4];
 //     count = 0;
 //     for (index_S = 0; index_S < R->number_of_ewald_shells; index_S++) {
 //       shell_sum = k_zero;
 //       for (index_R = 0; index_R < R->num[index_S]; index_R++) {
 //         r_12.comp1 = t_12.comp1 + R->vec_ai[count].comp1;
 //         r_12.comp2 = t_12.comp2 + R->vec_ai[count].comp2;
 //         r_12.comp3 = t_12.comp3 + R->vec_ai[count].comp3;
 //         Rsqrd = double_vec_dot(&r_12, &r_12);
 //         f000m(&em[count][0], gamma_1 * Rsqrd, gamma_1, mm);
 //         f000m(&en[count][0], Rsqrd / p_inv_ex, k_one / p_inv_ex, mm);
 //         for (i1 = 0; i1 <= mm; i1++) {
 //           en[count][i1] -= em[count][i1];
 //          }
 //           shell_sum += fabs(en[count][0]);
 //           count++;
 //          } // end loop on index_R
 //            } // end loop on index_S
 //             for (index_R = 0; index_R < count; index_R++) {
 //               r_12.comp1 = t_12.comp1 + R->vec_ai[index_R].comp1;
 //               r_12.comp2 = t_12.comp2 + R->vec_ai[index_R].comp2;
 //               r_12.comp3 = t_12.comp3 + R->vec_ai[index_R].comp3;
 //           for (m = 0; m <= mm; m++) {
 //             for (n = 0; n <= mm; n++) {
 //               for (pm = 0; pm <= mm; pm++) {
 //                 if (m + n + pm > mm) break;
 //                 fgtuv_temp = fac1 * ftuvn(m,n,pm,0,&en[index_R][0],r_12);
 //                 p_fgtuv = fgtuv + m * mm1 + n * mm2 + pm * mm3  + i4 * mm4 + j4;
 //                *p_fgtuv += fgtuv_temp;
 //                 fgtuv_max = (fgtuv_max > fabs(fgtuv_temp)) ? fgtuv_max : fabs(fgtuv_temp);
 //                 fgtuv_max = k_one;
 //                }
 //               }
 //              }
 //             }
 //            }
 //           }
 //      mcmurchie_davidson_3c_reversed_complex(Coulomb_cart,fgtuv,index_i,index_j,index_k,bfposi,bfposj,bfposk,nd2,nd3,\
 //      C1x,C1y,C1z,C2x,C2y,C2z,shells,shells_ax,job,file);
 //      //for (int ggg = 0; ggg < dimtp; ggg++) fprintf(file.out,"ggg %3d %10.4lf\n",ggg,Coulomb_cart[ggg]);
 //      three_center_cartesian_to_sh_shell_ax_reversed_complex(Coulomb_cart,Coulomb,index_i,index_j,index_k,bfposi1,bfposj1,bfposk1,\
 //      bfposi,bfposj,bfposk,nd2,nd3,nd5,nd6,shells,shells_ax,job,file);
 //      //for (int ggg = 0; ggg < dimtp; ggg++) fprintf(file.out,"ggg %3d %10.4lf\n",ggg,Coulomb[ggg]);
 //      free(fgtuv);
        mcmurchie_davidson_screen_complex(F_cart, index_i, index_j, bfposi, bfposj, nd2, C1x, C1y, C1z, fgtuv, shells, job, file);

//int num = sheli1 * shelj1 * sheli1 * shelj1;

   //for (int ggg = 0; ggg < num; ggg++) fprintf(file.out,"ggg %3d %16.10f %16.10f\n",ggg,(F_cart[ggg]).real(),(F_cart[ggg]).imag());
        free(fgtuv);
        two_center_cartesian_to_sh_shell1_complex(F_cart,F_sh,index_i,index_j,bfposi,bfposj,bfposi1,bfposj1,nd2,nd4,shells,job,file);
        bfposj   += shelj;
        bfposj1  += shelj1;
        gausposj += shells->ng_sh[index_j];
       } // close loop over index_j
      bfposi   += sheli;
      bfposi1  += sheli1;
      gausposi += shells->ng_sh[index_i];
     } // close loop over index_i
//for(i=0;i<nd1;i++){for(j=0;j<nd2;j++) { for(k=0;k<nd1;k++) { for(l=0;l<nd2;l++) { fprintf(file.out,"F %3d %3d %3d %3d %10.4lf\n",\
i,j,k,l,F_sh[i*nd2*nd1*nd2+j*nd1*nd2+k*nd2+l]);}}}}
//for(i=0;i<nd3;i++) { for(j=0;j<nd4;j++) { for(k=0;k<nd3;k++) { for(l=0;l<nd4;l++) { fprintf(file.out,"%3d %3d %3d %3d %10.4lf\n",\
i,j,k,l,F_sh[i*nd4*nd3*nd4+j*nd3*nd4+k*nd4+l]);}}}}

      free(F_cart);

}

*/
