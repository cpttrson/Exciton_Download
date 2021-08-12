
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
#include <iostream>
#include <cstring>
#include <cstdlib>
#include "mycomplex.h"
#include "conversion_factors.h"
#include "myconstants.h"
#include "mylogical.h"
#include "USER_DATA.h"
#include "TOOLS.h"
#include "LIMITS.h"
#include "MATRIX_UTIL.h"
#include "ALLOCATE_MEMORY_MOLECULE.h"
#include "SCF_MOLECULE.h"
#include "OPTICAL_SPECTRUM_MOLECULE.h"
#include "PLOTTING_MOLECULE.h"
#include "GW_BSE_MOLECULE.h"
#include "ERRORS.h"
#include "INPUT_MOLECULE.h"

using namespace std;

int startjob(ATOM *atoms, ATOM *atoms_ax, ATOM_TRAN *atom_p, SHELL *shells, GAUSSIAN *gaussians, SHELL *shells_ax, GAUSSIAN *gaussians_ax, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

  int i, j, k;

  // set default values for job parameters

  job->itol1         =  7.0; // real space cutoff for product of two Gaussians in Bohr (periodic systems only)
  job->itol2         =  7.0; // real space cutoff for product of two Gaussians in Bohr (periodic systems only)
  job->itol3         =  7.0; // real space cutoff for product of two Gaussians in Bohr (periodic systems only)
  job->itol4         = 16.0; // real space cutoff for product of two Gaussians in Bohr (periodic systems only)
  job->itol5         = 22.0; // real space cutoff for product of two Gaussians in Bohr (periodic systems only)

  job->spin_polarisation = 0;

  job->spin_index     = 0;      // not currently used
  job->mixing_type    = 0;      // not currently used

  job->type           = 0;      // default job type is scf
  job->xc_grd         = 1;      // default dft_grid is LARGE
  job->xc_rad         = 96;     // default number of radial points in DFT integration
  job->xc_lmx         = 12;     // default L max in Lebedev grid for DFT integration is L = 12
  job->xc_num         = 0;      // default is no DFT
  job->ham_type       = 0;      // default is DFT, not hybrid DFT nor Hartree-Fock
  job->ham_dft        = 0;      // default is DFT, not DFT-GGA
  job->ham_dft_num    = 1;      // number of functionals needed
  job->ham_dft_exc[0] = 1;      // default exchange functional is LDA
  job->ham_dft_exc[1] = 0;      // default exchange functional
  job->ham_dft_cor[0] = 1;      // default correlation functional is VWN
  job->ham_dft_cor[1] = 0;      // default correlation functional
  job->ham_hyb_wt     = k_zero; // default weight for Fock exchange in hybrid DFT Hamiltonian

  job->coul_int = 0;            // default is calculate Coulomb integrals with multipole approximation
  job->exch_int = 1;            // default is calculate Exchange integrals without multipole approximation

  job->kpoints = 0;             // periodic systems only
  job->values  = 0;

  job->fix_occ = 0;             // Occupancy for density matrix calculation determined by Fermi level

  job->linewidth = 0.1;         // default linewidth for plotting 0.1 eV

  switch (job->spin_polarisation) {
    case TRUE:
      job->spin_dim = 2;
      job->spin_fac = 1;
      job->spin_pol = 1;
      job->spin_orb = 0;
      job->spin_tot = 0;
      break;
    case FALSE:
      job->spin_dim = 1;
      job->spin_fac = 2;
      job->spin_pol = 0;
      job->spin_orb = 0;
      job->spin_tot = 0;
  }

      job->lmax = 1;            // default max l value for multipole expansion of 2e integrals
      job->lmax_fac = 3;        // default maximum number of monomials for given l value plus 1 for spheropole moment
      job->int_exist = 0;       // default is two electron integrals must be calculated rather than existing on disk
      job->int_exist_no_sym = 0;// default is two electron integrals must be calculated rather than existing on disk
      job->scf_coulomb = 1;     // default is to calculate SCF 2e coulomb integrals
      job->scf_exchange = 1;    // default is to calculate SCF 2e exchange integrals
      job->scf_direct = 0;      // default is to store SCF 2e integrals on disk
      job->scf_denfit = 0;      // default is to calculate four centre 2e integrals exactly without density fitting
      job->scf_dencou = 1;      // default is to use 3-centre coulomb coefficients for density fitting when it is on
      job->scf_trans  = 1;      // default is Cholesky transformation for basis orthogonalisation
      job->bse_denfit = 0;      // default is to calculate four centre 2e integrals for BSE exactly without density fitting
      job->bse_screxc = 0;      // default is to calculate unscreened 2e integrals for BSE 
      job->bse_spk = 0;         // default is not to use SCALAPACK for BSE and TDHF 
      job->bse_cou = 0;         // default is not to test calculation of 3-centre Coulomb integrals
      job->bse_exc = 0;         // default is not to test calculation of 3-centre exchange integrals
      job->gw = 0;              // default is not to calculate GW correction for BSE 

  int int_value[2], self[10];
  int nbands, nocc, nvir, ntransitions;
  int numfrag = 2, max_atoms = 40;
  int natoms[numfrag], nat1[max_atoms], nat[max_atoms][2];
  char sgs[5];
  char pms[5];
  char trs[5];
  char kss[5];
  char jobname[17];
  char jobname1[17];
  char taskname[15];
  char title[100];
  char title1[100];
  char title2[15];
  char title3[15];
  char fermi_homo[12], fermi_bands[12], fermi_bands1[12], nspec[14], spectrum[14], monkhorst_pack[9], ntrans[12];
  char hamiltonian[7], field[80], tamm_dancoff[12], spin_state[12], linear_algebra[12], scissor_shift[9], scale_factor[9];
  char int_range[12], double_range[13], plot_type[12];
  double double_value[2], energy_range[2];
  FERMI fermi;

  do {

    read_line(file.job, title, 99);
    sscanf(title, "%s", jobname);
    //printf("%s\n",jobname); fflush(stdout);

    // *****JOB: SPIN POLARISATION ***************************************************

    if (!strcmp(jobname, "SPIN_POLARISED") || !strcmp(jobname, "UHF")) {

      read_line(file.job, title, 99);
      sscanf(title, "%d %d", &job->spin_polarisation,&job->spin_orb);

       switch (job->spin_polarisation) {
         case TRUE:
           job->spin_dim = 2;
           job->spin_fac = 1;
           job->spin_pol = 1;
            if (job->taskid == 0) 
             fprintf(file.out, "SPIN POLARISED CALCULATION\n\n");
            if (job->taskid == 0 && job->spin_orb == 1) 
             fprintf(file.out, "SPIN ORBIT COUPLING ON\n\n");  // not programmed
           break;
         case FALSE:
           job->spin_dim = 1;
           job->spin_fac = 2;
           job->spin_pol = 0;
           job->spin_orb = 0;  // spin orbit always off for an unpolarised calculation
             //if (job->taskid == 0)
             //fprintf(file.out, "CLOSED SHELL CALCULATION\n\n");
       }
    
    }

    // *****JOB: SCF CALCULATION ***************************************************

    else if (!strcmp(jobname, "SCF")) {

      int ion, spin, number_of_magnetic_ions;
      int nexc, ncor;
      char exch_type[20];
      char ham_dft1[10], ham_dft2[10], ham_type[7];
      char hamiltonian[7], spin_pol[6], x_functional[15], c_functional[12], dft_grid[9], integrals[10], direct[4];
      char guess[13], intexist[13], monkhorst_pack[9];
      double w[2];

      hamiltonian[0] = 0; // initialise char arrays to be printed
      x_functional[0] = 0;
      c_functional[0] = 0;
      exch_type[0] = 0;
      ham_dft1[0] = 0;
      ham_dft2[0] = 0;
      ham_type[0] = 0;
      dft_grid[0] = 0;
      spin_pol[0] = 0;
      guess[0] = 0;
      direct[0] = 0;
      integrals[0] = 0;
      intexist[0] = 0;
      monkhorst_pack[0] = 0;
     
      // Set Job Defaults

      read_line(file.job, title, 99);
      sscanf(title, "%s", file.scf_eigvec); // Directory for MPI file IO 

      job->mpp = 0;             // Default is to write SCF eigenvectors to disk

      job->type = 0;            // type = 0 scf,  = 1 bse

      job->vectors = 2;
      job->values  = 2;
      job->density = 2;
      job->kpoints = 0;

      job->xc_hfx  = 1;          // Default is Hartree-Fock
      job->xc_num  = 0;          // Default is no DFT functionals
      job->xc_typ[0] = -1;          
      job->xc_typ[1] = -1;          

      job->spin_dim = 1;         // Default is no spin polarisation
      job->spin_fac = 2;
      job->spin_pol = 0;
      job->spin_orb = 0;

      for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++) {
         atoms->magnetic[i] = 1;
         atoms->spin[i]     = 0;
        }

      job->guess_type = 0;       // Default is new calculation from atomic densities

      job->max_cycle  = 50;      // Default is 50 cycles
      job->scf_tol    = 1.0e-06; // Default is 1e-06

      if (atoms->number_of_sh_bfns_in_unit_cell <= 3000) {
      job->scf_direct = 0;  
      sprintf(direct,"%s","OFF");
     }
      else if (atoms->number_of_sh_bfns_in_unit_cell >  3000) {
      job->scf_direct = 1;       // Default is On for more than 3000 basis functions per unit cell
      sprintf(direct,"%s","ON");
     }

      job->coul_int = 0;         // Default is no approximation of integrals

      job->fock_mixing = 0.7;    // Default is 0.7
      job->mixing_order = 1;     // Default is simple mixing with weight 0.7
      job->diis = 0;             // Default is off

      if (job->sgs == 1) sprintf(sgs,"%s","Y");
      if (job->pms == 1) sprintf(pms,"%s","Y");
      if (job->trs == 1) sprintf(trs,"%s","Y");
      if (job->kss == 1) sprintf(kss,"%s","Y");
      if (job->sgs == 0) sprintf(sgs,"%s","N");
      if (job->pms == 0) sprintf(pms,"%s","N");
      if (job->trs == 0) sprintf(trs,"%s","N");
      if (job->kss == 0) sprintf(kss,"%s","N");

      fermi.bands[0] = 1;
      fermi.bands[1] = atoms->number_of_sh_bfns_in_unit_cell;
      fermi.bands[2] = 1;
      fermi.bands[3] = atoms->number_of_sh_bfns_in_unit_cell;

        switch (crystal->type[0]) {

    case 'C':
      fermi.is[0] = 4;
      fermi.is[1] = 4;
      fermi.is[2] = 4;
      sprintf(monkhorst_pack,"%d %d %d",fermi.is[0], fermi.is[1], fermi.is[2]);
      break;

    case 'S':
      fermi.is[0] = 4;
      fermi.is[1] = 4;
      fermi.is[2] = 1;
      sprintf(monkhorst_pack,"%d %d %d",fermi.is[0], fermi.is[1], fermi.is[2]);
      break;

    case 'P':
      fermi.is[0] = 1;
      fermi.is[1] = 1;
      fermi.is[2] = 4;
      sprintf(monkhorst_pack,"%d %d %d",fermi.is[0], fermi.is[1], fermi.is[2]);
      break;

    case 'M':
      fermi.is[0] = 1;
      fermi.is[1] = 1;
      fermi.is[2] = 1;
      sprintf(monkhorst_pack,"%s","OFF");
      break;

    } // close switch

      // read input file

   do {

      read_line(file.job, title, 99);
      sscanf(title, "%s", jobname1);
      //printf("%s\n",jobname1);

        if (!strcmp(jobname1, "ITOL")) {  // applies to periodic systems only
          read_line(file.job, title, 99);
          sscanf(title, "%lf %lf %lf %lf %lf", &job->itol1, &job->itol2, &job->itol3, &job->itol4, &job->itol5);
         }

        if (!strcmp(jobname1, "MPP")) {
          job->mpp = 1;  // write eigenvectors to disk at end of SCF only in scf_evec MPI File
         }

        if (!strcmp(jobname1, "DIIS")) {
          job->diis = 1; 
         }

        if (!strcmp(jobname1, "DFT")) {

          job->xc_grd = 0;           // Default DFT grid is STANDARD
          job->xc_lmx = 13;
          job->xc_rad = 64;
          sprintf(dft_grid,"%s","STANDARD");
          sprintf(hamiltonian,"%s","DFT");

            read_line(file.job, title, 99);
            sscanf(title, "%d", &job->xc_num);
            if (job->xc_num < 0 || job->xc_num > 2 && job->taskid == 0) {
              fprintf(file.out,"ERROR: Number of DFT XC potentials %3d incorrect\n",job->xc_num);
              MPI_Finalize();
              exit(1);
             }

          for (i = 0; i < job->xc_num; i++) {
            read_line(file.job, title, 99);
            sscanf(title, "%s", exch_type);
            if (exch_type[0] != 'X' || exch_type[1] != 'C' && job->taskid == 0) {
              fprintf(file.out,"DFT Hamiltonian requires XC potential input\n");
              MPI_Finalize();
              exit(1);
             }

                 job->xc_hfx = 0;

              if (!strcmp(exch_type, "XC_LDA_X")) {
                 job->xc_typ[i]   = 1;      
                 sprintf(x_functional,"%s","LDA");
                }
              else if (!strcmp(exch_type, "XC_LDA_C_WIGNER")) {
                 job->xc_typ[i]   = 2;      
                 sprintf(c_functional,"%s","LDA_WIGNER");
                }
              else if (!strcmp(exch_type, "XC_LDA_C_RPA")) {
                 job->xc_typ[i]   = 3;      
                 sprintf(c_functional,"%s","LDA_RPA");
                }
              else if (!strcmp(exch_type, "XC_LDA_C_HL")) {
                 job->xc_typ[i]   = 4;      
                 sprintf(c_functional,"%s","LDA_HL");
                }
              else if (!strcmp(exch_type, "XC_LDA_C_GL")) {
                 job->xc_typ[i]   = 5;      
                }
              else if (!strcmp(exch_type, "XC_LDA_C_XALPHA")) {
                 job->xc_typ[i]   = 6;      
                }
              else if (!strcmp(exch_type, "XC_LDA_VWN")) {
                 job->xc_typ[i]   = 7;      
                }
              else if (!strcmp(exch_type, "XC_LDA_VWN_RPA")) {
                 job->xc_typ[i]   = 8;      
                }
              else if (!strcmp(exch_type, "XC_LDA_PZ")) {
                 job->xc_typ[i]   = 9;      
                }
              else if (!strcmp(exch_type, "XC_LDA_PZ_MOD")) {
                 job->xc_typ[i]   = 10;      
                }
              else if (!strcmp(exch_type, "XC_LDA_OB_PZ")) {
                 job->xc_typ[i]   = 11;      
                }
              else if (!strcmp(exch_type, "XC_LDA_PW")) {
                 job->xc_typ[i]   = 12;      
                }
              else if (!strcmp(exch_type, "XC_LDA_PW_MOD")) {
                 job->xc_typ[i]   = 13;      
                }
              else if (!strcmp(exch_type, "XC_LDA_OB_PW")) {
                 job->xc_typ[i]   = 14;      
                }
              else if (!strcmp(exch_type, "XC_LDA_C_2D_AMGB")) {
                 job->xc_typ[i]   = 15;      
                }
              else if (!strcmp(exch_type, "XC_LDA_C_2D_PRM")) {
                 job->xc_typ[i]   = 16;      
                }
              else if (!strcmp(exch_type, "XC_LDA_C_vBH")) {
                 job->xc_typ[i]   = 17;      
                }
              else if (!strcmp(exch_type, "XC_LDA_C_1D_CSC")) {
                 job->xc_typ[i]   = 18;      
                }
              else if (!strcmp(exch_type, "XC_LDA_X_2D")) {
                 job->xc_typ[i]   = 19;      
                }
              else if (!strcmp(exch_type, "XC_LDA_TETER93")) {
                 job->xc_typ[i]   = 20;      
                }
              else if (!strcmp(exch_type, "XC_LDA_X_1D")) {
                 job->xc_typ[i]   = 21;      
                }
              else if (!strcmp(exch_type, "XC_LDA_C_ML1")) {
                 job->xc_typ[i]   = 22;      
                }
              else if (!strcmp(exch_type, "XC_LDA_C_ML2")) {
                 job->xc_typ[i]   = 23;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_PBE")) {
                 job->xc_typ[i]   = 101;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_PBE_R")) {
                 job->xc_typ[i]   = 102;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_B86")) {
                 job->xc_typ[i]   = 103;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_B86_R")) {
                 job->xc_typ[i]   = 104;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_B86_MGC")) {
                 job->xc_typ[i]   = 105;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_B88")) {
                 job->xc_typ[i]   = 106;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_G96")) {
                 job->xc_typ[i]   = 107;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_PW86")) {
                 job->xc_typ[i]   = 108;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_PW91")) {
                 job->xc_typ[i]   = 109;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_OPTX")) {
                 job->xc_typ[i]   = 110;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_DK87_R1")) {
                 job->xc_typ[i]   = 111;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_DK87_R2")) {
                 job->xc_typ[i]   = 112;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_LG93")) {
                 job->xc_typ[i]   = 113;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_FT97_A")) {
                 job->xc_typ[i]   = 114;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_FT97_B")) {
                 job->xc_typ[i]   = 115;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_PBE_SOL")) {
                 job->xc_typ[i]   = 116;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_RPBE")) {
                 job->xc_typ[i]   = 117;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_WC")) {
                 job->xc_typ[i]   = 118;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_mPW91")) {
                 job->xc_typ[i]   = 119;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_AM05")) {
                 job->xc_typ[i]   = 120;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_PBEA")) {
                 job->xc_typ[i]   = 121;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_MPBE")) {
                 job->xc_typ[i]   = 122;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_XPBE")) {
                 job->xc_typ[i]   = 123;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_2D_B86_MGC")) {
                 job->xc_typ[i]   = 124;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_BAYESIAN")) {
                 job->xc_typ[i]   = 125;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_PBE_JSJR")) {
                 job->xc_typ[i]   = 126;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_2D_B88")) {
                 job->xc_typ[i]   = 127;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_2D_B86")) {
                 job->xc_typ[i]   = 128;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_2D_PBE")) {
                 job->xc_typ[i]   = 129;      
                }
              else if (!strcmp(exch_type, "XC_GGA_C_PBE")) {
                 job->xc_typ[i]   = 130;      
                }
              else if (!strcmp(exch_type, "XC_GGA_C_LYP")) {
                 job->xc_typ[i]   = 131;      
                }
              else if (!strcmp(exch_type, "XC_GGA_C_P86")) {
                 job->xc_typ[i]   = 132;      
                }
              else if (!strcmp(exch_type, "XC_GGA_C_PBE_SOL")) {
                 job->xc_typ[i]   = 133;      
                }
              else if (!strcmp(exch_type, "XC_GGA_C_PW91")) {
                 job->xc_typ[i]   = 134;      
                }
              else if (!strcmp(exch_type, "XC_GGA_C_AM05")) {
                 job->xc_typ[i]   = 135;      
                }
              else if (!strcmp(exch_type, "XC_GGA_C_XPBE")) {
                 job->xc_typ[i]   = 136;      
                }
              else if (!strcmp(exch_type, "XC_GGA_C_LM")) {
                 job->xc_typ[i]   = 137;      
                }
              else if (!strcmp(exch_type, "XC_GGA_C_PBE_JRGX")) {
                 job->xc_typ[i]   = 138;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_OPTB88_VDW")) {
                 job->xc_typ[i]   = 139;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_PBEK1_VDW")) {
                 job->xc_typ[i]   = 140;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_OPTPBE_VDW")) {
                 job->xc_typ[i]   = 141;      
                }
              else if (!strcmp(exch_type, "XC_GGA_X_RGE2")) {
                 job->xc_typ[i]   = 142;      
                }
              else if (!strcmp(exch_type, "XC_GGA_C_RGE2")) {
                 job->xc_typ[i]   = 143;      
                }
              else if (!strcmp(exch_type, "XC_GGA_XC_LB")) {
                 job->xc_typ[i]   = 160;      
                }
              else if (!strcmp(exch_type, "XC_GGA_XC_HCTH_93")) {
                 job->xc_typ[i]   = 161;      
                }
              else if (!strcmp(exch_type, "XC_GGA_XC_HCTH_120")) {
                 job->xc_typ[i]   = 162;      
                }
              else if (!strcmp(exch_type, "XC_GGA_XC_HCTH_147")) {
                 job->xc_typ[i]   = 163;      
                }
              else if (!strcmp(exch_type, "XC_GGA_XC_HCTH_407")) {
                 job->xc_typ[i]   = 164;      
                }
              else if (!strcmp(exch_type, "XC_GGA_XC_EDF1")) {
                 job->xc_typ[i]   = 165;      
                }
              else if (!strcmp(exch_type, "XC_GGA_XC_XLYP")) {
                 job->xc_typ[i]   = 166;      
                }
              else if (!strcmp(exch_type, "XC_GGA_XC_B97")) {
                 job->xc_typ[i]   = 167;      
                }
              else if (!strcmp(exch_type, "XC_GGA_XC_B97_1")) {
                 job->xc_typ[i]   = 168;      
                }
              else if (!strcmp(exch_type, "XC_GGA_XC_B97_2")) {
                 job->xc_typ[i]   = 169;      
                }
              else if (!strcmp(exch_type, "XC_GGA_XC_B97_D")) {
                 job->xc_typ[i]   = 170;      
                }
              else if (!strcmp(exch_type, "XC_GGA_XC_B97_K")) {
                 job->xc_typ[i]   = 171;      
                }
              else if (!strcmp(exch_type, "XC_GGA_XC_B97_3")) {
                 job->xc_typ[i]   = 172;      
                }
              else if (!strcmp(exch_type, "XC_GGA_XC_PBE1W")) {
                 job->xc_typ[i]   = 173;      
                }
              else if (!strcmp(exch_type, "XC_GGA_XC_MPWLYP1W")) {
                 job->xc_typ[i]   = 174;      
                }
              else if (!strcmp(exch_type, "XC_GGA_XC_PBELYP1W")) {
                 job->xc_typ[i]   = 175;      
                }
              else if (!strcmp(exch_type, "XC_GGA_XC_SB98_1a")) {
                 job->xc_typ[i]   = 176;      
                }
              else if (!strcmp(exch_type, "XC_GGA_XC_SB98_1b")) {
                 job->xc_typ[i]   = 177;      
                }
              else if (!strcmp(exch_type, "XC_GGA_XC_SB98_1c")) {
                 job->xc_typ[i]   = 178;      
                }
              else if (!strcmp(exch_type, "XC_GGA_XC_SB98_2a")) {
                 job->xc_typ[i]   = 179;      
                }
              else if (!strcmp(exch_type, "XC_GGA_XC_SB98_2b")) {
                 job->xc_typ[i]   = 180;      
                }
              else if (!strcmp(exch_type, "XC_GGA_XC_SB98_2c")) {
                 job->xc_typ[i]   = 181;      
                }
              else {
              //   xc_func_type func;
              //   if(xc_func_init(&func, job->xc_typ[i], 1) != 0){
                 if (job->taskid == 0) {
                   fprintf(file.out, "Functional %s not found\n",exch_type);
                   fprintf(file.out,"Choices are HARTREE_FOCK or at www.tddft.org/programs/octopus/wiki/index.php/Libxc:manual\n");
                  }
                   MPI_Finalize();
                   exit(0);
              //    }
              //   xc_func_end(&func);
                }
               } // close loop on i
   
     	     } // close DFT input

        if (!strcmp(jobname1, "HYBRID")) {
          job->xc_grd = 0;           // Default HYBRID DFT grid is STANDARD
	
          job->xc_lmx = 13;
          job->xc_rad = 64;
          job->xc_hfx = 1;
          sprintf(hamiltonian,"%s","HYBRID");

            read_line(file.job, title, 99);
            sscanf(title, "%d", &job->xc_num);
            if (job->xc_num < 0 || job->xc_num > 2 && job->taskid == 0) {
              fprintf(file.out,"ERROR: HYBRID Hamiltonian requires number of XC potentials\n");
              MPI_Finalize();
              exit(1);
             }

          for (i = 0; i < job->xc_num; i++) {
            read_line(file.job, title, 99);
            sscanf(title, "%s", exch_type);
            if (exch_type[0] != 'X' || exch_type[1] != 'C' && job->taskid == 0) {
              fprintf(file.out,"HYBRID Hamiltonian requires XC potential input\n");
              MPI_Finalize();
              exit(1);
             }

                   if (!strcmp(exch_type, "XC_HYB_GGA_XC_B3PW91")) {
                 job->xc_typ[i]   = 401;      
                 sprintf(x_functional,"%s","B3PW91");
                 sprintf(c_functional,"%s","B3PW91");
                }
              else if (!strcmp(exch_type, "XC_HYB_GGA_XC_B3LYP")) {
                 job->xc_typ[i]   = 402;      
                 sprintf(x_functional,"%s","B3LYP");
                 sprintf(c_functional,"%s","B3LYP");
                }
              else if (!strcmp(exch_type, "XC_HYB_GGA_XC_B3P86")) {
                 job->xc_typ[i]   = 403;      
                }
              else if (!strcmp(exch_type, "XC_HYB_GGA_XC_O3LYP")) {
                 job->xc_typ[i]   = 404;      
                }
              else if (!strcmp(exch_type, "XC_HYB_GGA_XC_PBEH")) {
                 job->xc_typ[i]   = 406;      
                }
              else if (!strcmp(exch_type, "XC_HYB_GGA_XC_B97")) {
                 job->xc_typ[i]   = 407;      
                }
              else if (!strcmp(exch_type, "XC_HYB_GGA_XC_B97_1")) {
                 job->xc_typ[i]   = 408;      
                }
              else if (!strcmp(exch_type, "XC_HYB_GGA_XC_B97_2")) {
                 job->xc_typ[i]   = 410;      
                }
              else if (!strcmp(exch_type, "XC_HYB_GGA_XC_X3LYP")) {
                 job->xc_typ[i]   = 411;      
                }
              else if (!strcmp(exch_type, "XC_HYB_GGA_XC_B1WC")) {
                 job->xc_typ[i]   = 412;      
                }
              else if (!strcmp(exch_type, "XC_HYB_GGA_XC_B97_K")) {
                 job->xc_typ[i]   = 413;      
                }
              else if (!strcmp(exch_type, "XC_HYB_GGA_XC_B97_3")) {
                 job->xc_typ[i]   = 414;      
                }
              else if (!strcmp(exch_type, "XC_HYB_GGA_XC_mPW3PW")) {
                 job->xc_typ[i]   = 415;      
                }
              else if (!strcmp(exch_type, "XC_HYB_GGA_XC_B1LYP")) {
                 job->xc_typ[i]   = 416;      
                }
              else if (!strcmp(exch_type, "XC_HYB_GGA_XC_B1PW91")) {
                 job->xc_typ[i]   = 417;      
                }
              else if (!strcmp(exch_type, "XC_HYB_GGA_XC_mPW1PW")) {
                 job->xc_typ[i]   = 418;      
                }
              else if (!strcmp(exch_type, "XC_HYB_GGA_XC_mPW3LYP")) {
                 job->xc_typ[i]   = 419;      
                }
              else if (!strcmp(exch_type, "XC_HYB_GGA_XC_SB98_1a")) {
                 job->xc_typ[i]   = 420;      
                }
              else if (!strcmp(exch_type, "XC_HYB_GGA_XC_SB98_1b")) {
                 job->xc_typ[i]   = 421;      
                }
              else if (!strcmp(exch_type, "XC_HYB_GGA_XC_SB98_1c")) {
                 job->xc_typ[i]   = 422;      
                }
              else if (!strcmp(exch_type, "XC_HYB_GGA_XC_SB98_2a")) {
                 job->xc_typ[i]   = 423;      
                }
              else if (!strcmp(exch_type, "XC_HYB_GGA_XC_SB98_2b")) {
                 job->xc_typ[i]   = 424;      
                }
              else if (!strcmp(exch_type, "XC_HYB_GGA_XC_SB98_2c")) {
                 job->xc_typ[i]   = 425;      
                }
              else if (!strcmp(exch_type, "XC_MGGA_X_LTA")) {
                 job->xc_typ[i]   = 201;      
                }
              else if (!strcmp(exch_type, "XC_MGGA_X_TPSS")) {
                 job->xc_typ[i]   = 202;      
                }
              else if (!strcmp(exch_type, "XC_MGGA_X_M06L")) {
                 job->xc_typ[i]   = 203;      
                }
              else if (!strcmp(exch_type, "XC_MGGA_X_GVT4")) {
                 job->xc_typ[i]   = 204;      
                }
              else if (!strcmp(exch_type, "XC_MGGA_X_TAU_HCTH")) {
                 job->xc_typ[i]   = 205;      
                }
              else if (!strcmp(exch_type, "XC_MGGA_X_BR89")) {
                 job->xc_typ[i]   = 206;      
                }
              else if (!strcmp(exch_type, "XC_MGGA_X_BJ06")) {
                 job->xc_typ[i]   = 207;      
                }
              else if (!strcmp(exch_type, "XC_MGGA_X_TB09")) {
                 job->xc_typ[i]   = 208;      
                }
              else if (!strcmp(exch_type, "XC_MGGA_X_RPP09")) {
                 job->xc_typ[i]   = 209;      
                }
              else if (!strcmp(exch_type, "XC_MGGA_C_TPSS")) {
                 job->xc_typ[i]   = 231;      
                }
              else if (!strcmp(exch_type, "XC_MGGA_C_VSXC")) {
                 job->xc_typ[i]   = 231;      
                }
              else if (!strcmp(exch_type, "XC_LCA_OMC")) {
                 job->xc_typ[i]   = 301;      
                }
              else if (!strcmp(exch_type, "XC_LCA_LCH")) {
                 job->xc_typ[i]   = 302;      
                }
              else {
             //    xc_func_type func;
             //    if(xc_func_init(&func, job->xc_typ[i], 1) != 0){
                 if (job->taskid == 0) {
                   fprintf(file.out, "Functional %s not found\n",exch_type);
                   fprintf(file.out,"Choices are HARTREE_FOCK or at www.tddft.org/programs/octopus/wiki/index.php/Libxc:manual\n");
                  }
                   MPI_Finalize();
                   exit(0);
             //     }
             //    xc_func_end(&func);
                }
               } // close loop on i

             } // close HYBRID input

            if (!strcmp(jobname1, "DFT_GRID")) {
              read_line(file.job, title, 99);
              sscanf(title, "%s", jobname1);
              if (!strcmp(jobname1, "STANDARD")) {
              job->xc_grd = 0;
              job->xc_lmx = 13;
              job->xc_rad = 64;
              sprintf(dft_grid,"%s","STANDARD");
             } 
              else if (!strcmp(jobname1, "LARGE")) {
              job->xc_grd = 1;
              job->xc_lmx = 13;
              job->xc_rad = 96;
              sprintf(dft_grid,"%s","LARGE");
             } 
              else if (!strcmp(jobname1, "XLARGE")) {
              job->xc_grd = 2;
              job->xc_lmx = 16;
              job->xc_rad = 96;
              sprintf(dft_grid,"%s","XLARGE");
             } 
              else if (!strcmp(jobname1, "XXLARGE")) {
              job->xc_grd = 3;
              job->xc_lmx = 19;
              job->xc_rad = 96;
              sprintf(dft_grid,"%s","XXLARGE");
             } 
              else {
              if (job->taskid == 0) 
              fprintf(file.out,"ERROR: DFT_GRID must be one of STANDARD, LARGE, XLARGE, XXLARGE\n");
              MPI_Finalize();
              exit(0);
             } 
            } 

        if (!strcmp(jobname1, "UHF")) {
           read_line(file.job, title, 99);
           sscanf(title, "%d", &job->spin_tot);
           job->spin_dim = 2;
           job->spin_fac = 1;
           job->spin_pol = 1;
           job->spin_polarisation = 1;
           sprintf(hamiltonian,"%s","UHF");
           sprintf(spin_pol,"%s","ON");
          }


        if (!strcmp(jobname1, "ATOM_SPIN")) {
            int spin_tot = 0;
            read_line(file.job, title, 99);
            sscanf(title, "%d", &number_of_magnetic_ions);
            for (i = 0; i < number_of_magnetic_ions; i++) {
              read_line(file.job, title, 99);
              sscanf(title, "%d %d", &ion, &spin);
              atoms->magnetic[ion - 1] = 2;
              atoms->spin[ion - 1] = spin;
              spin_tot += spin;
             }
            if (spin_tot != job->spin_tot) {
              fprintf(file.out,"Sum of spins of ions %3d does not equal total spin specified %3d\n",spin_tot, job->spin_tot);
              MPI_Finalize();
              exit(0);
             }
            }

      if (!strcmp(jobname1, "RESTART")) {
      //read_line(file.job, title, 99);
      //sscanf(title, "%s", guess);
      //if (guess[0] == 'A')
      //job->guess_type = 0;
      //else if (guess[0] == 'R')
      job->guess_type = 1;
      //if (guess[0] != 'A' && guess[0] != 'R') {
      //fprintf(file.out,"ERROR: GUESS %s must be Atoms (new calculation) or Restart\n",guess);
      //MPI_Finalize();
      //exit(1);
     }

      if (!strcmp(jobname1, "INTEXIST")) {
      //read_line(file.job, title, 99);
      //sscanf(title, "%s", intexist);
      //if (intexist[0] == 'I')
      job->int_exist = 1;      // if SCF integrals exit on disk then do not recalculate
     }

      if (!strcmp(jobname1, "MAX_CYCLES")) {
      read_line(file.job, title, 99);
      sscanf(title, "%d", &job->max_cycle);
      if (job->max_cycle < 1 || job->max_cycle > 999) {
      fprintf(file.out,"ERROR: MAX_CYCLES must range between 1 and 999\n");
      MPI_Finalize();
      exit(1);
     }
    }

      if (!strcmp(jobname1, "SCF_TOLERANCE")) {
      read_line(file.job, title, 99);
      sscanf(title, "%lf", &job->scf_tol);
      if (job->scf_tol < 1 || job->scf_tol > 20) {
      fprintf(file.out,"ERROR: SCF_TOLERANCE must range between 1 and 20\n");
      MPI_Finalize();
      exit(1);
     }
      job->scf_tol = pow(10.00,-job->scf_tol);
     }

      if (!strcmp(jobname1, "FMIXING")) {
      read_line(file.job, title, 99);
      sscanf(title, "%lf", &job->fock_mixing);
      if (job->fock_mixing < 0.00 || job->fock_mixing > 0.99) {
      fprintf(file.out,"ERROR: MAX_CYCLES must range between 0.00 and 0.99\n");
      MPI_Finalize();
      exit(1);
     }
    }

      if (!strcmp(jobname1, "PULAY_ORDER")) {
      read_line(file.job, title, 99);
      sscanf(title, "%d", &job->mixing_order);
      if (job->mixing_order < 1 || job->mixing_order > 20) {
      fprintf(file.out,"ERROR: PULAY ORDER must range between 1 and 20\n");
      MPI_Finalize();
      exit(1);
     }
     }

      if (!strcmp(jobname1, "CANONICAL")) {
      job->scf_trans = 0;
      if (job->taskid == 0) {
      fprintf(file.out,"CANONICAL TRANSFORMATION OF BASIS SELECTED\n");
     }
    }

      if (!strcmp(jobname1, "DIRECT_SCF")) {
      job->scf_direct = 1;
      sprintf(direct,"ON");
     }

      if (!strcmp(jobname1, "DIRECT_SCF_OFF")) {
      job->scf_direct = 0;
      sprintf(direct,"OFF");
     }

      if (!strcmp(jobname1, "COU_OFF")) {
      job->scf_coulomb = 0;
      if (job->taskid == 0)
      fprintf(file.out,"2E COULOMB INTEGRALS SWITCHED OFF\n");
     }

      if (!strcmp(jobname1, "EXC_OFF")) {
      job->scf_exchange = 0;
      if (job->taskid == 0)
      fprintf(file.out,"2E EXCHANGE INTEGRALS SWITCHED OFF\n");
     }

      if (!strcmp(jobname1, "DENSITY_FITTING")) {
      job->scf_denfit = 1;                               // not available for SCF at present
      for (i = 0; i < atoms_ax->number_of_atoms_in_unit_cell; i++) {
        if (atoms_ax->basis_set[i] == -1) {
        if (job->taskid == 0)
        fprintf(file.out,"NO AUXILLARY BASIS FOR ATOM %3d\n",i);
        MPI_Finalize();
        exit(1);
       }
      }
      fprintf(file.out,"2E INTEGRALS CALCULATED BY DENSITY FITTING\n");
     }

  } while (strcmp(jobname1, "END"));

 if (hamiltonian[0] == 0)  sprintf(hamiltonian,"%s","RHF");
 if (x_functional[0] == 0) sprintf(x_functional,"%s","HF");
 if (c_functional[0] == 0) sprintf(c_functional,"%s","NONE");
 if (dft_grid[0] == 0)     sprintf(dft_grid,"%s","STANDARD");
 if (spin_pol[0] == 0)     sprintf(spin_pol,"%s","OFF");
 if (guess[0] == 0)        sprintf(guess,"%s","ATOM DENSITY");
 if (direct[0] == 0)       sprintf(direct,"%s","OFF");
 if (integrals[0] == 0)    sprintf(integrals,"%s","EXACT");

 if (job->taskid == 0) {

 fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
 fprintf(file.out,"| SCF HAMILTONIAN   %7s | EXCHANGE %14s | CORRELATION %11s | SPIN POLARIZATION   %3s |\n", \
 hamiltonian,x_functional, c_functional, spin_pol);
 fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
 if (dft_grid[0] != 'N') {
 fprintf(file.out,"| DFT GRID         %8s |                         |                         |                         |\n", \
 dft_grid);
 fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
 }
 fprintf(file.out,"| MAX CYCLES           %4d | SCF TOLERANCE   %-7.1e | MIXING PARAMETER   %4.1lf | PULAY MIXING ORDER  %3d |\n",\
 job->max_cycle, job->scf_tol, job->fock_mixing, job->mixing_order);
 fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
 fprintf(file.out,"| INTEGRALS       %9s | DIRECT SCF          %3s | GUESS     %13s | MONKHORST-PACK %8s |\n", \
 integrals, direct, guess, monkhorst_pack);
 fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
 if (crystal->type[0] == 'C' || crystal->type[0] == 'S' || crystal->type[0] == 'P') { 
 fprintf(file.out,\
 "| COULOMB ITOL    %4.1f %4.1f | EXC ITOL %4.1f %4.1f %4.1f | MXR %2d R cutoff %4.1f A  | MXG %2d G cutoff %4.1f 1/A|\n", \
 job->itol1,job->itol2,job->itol3,job->itol4,job->itol5,job->mxr,R->cutoff * bohr_to_AA,job->mxg,G->cutoff / bohr_to_AA);
 fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
 fprintf(file.out,\
 "| RMAX %6d RMARG  %6d | EWALD %5d LAST %5d  | G MAX %4d G LAST %4d  | SGS %s PMS %s TRS %s KSS %s |\n", \
    R->max_vector,R->margin_vector,R->last_vector,R->last_ewald_vector,G->max_vector,G->last_vector,sgs,pms,trs,kss); //}
 fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n\n");
 }

 } // if (job->taskid == 0)

  switch (crystal->type[0]) {
   case 'C':
   case 'S':
   case 'P':

     break ;

   case 'M':

     scf_molecule(&fermi,atoms,atoms_ax,atom_p,shells,gaussians,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);

     break ;
  }

 free_fermi(&fermi,job);

 }

    // *****JOB: GW CALCULATION ***************************************************************************************

    else if (!strcmp(jobname, "GW")) {

      int_value[0] = 1;
      int_value[1] = atoms->number_of_sh_bfns_in_unit_cell;
      double_value[0] = lower_energy_limit;
      double_value[1] = upper_energy_limit;

      // Initialise char arrays

      double_range[0] = 0;
      int_range[0] = 0;
      fermi_homo[0] = 0; 
      fermi_bands[0] = 0;
      fermi_bands1[0] = 0;
      nspec[0] = 0;
      spectrum[0] = 0;
      ntrans[0] = 0;
      monkhorst_pack[0] = 0;

      sprintf(double_range, "%12s", "ENERGY RANGE");
      sprintf(int_range, "%11s", "ENERGY BAND");

      // Set job defaults

      fermi_default(&fermi,crystal,atoms,job,file);
      job->self_plot = 0;
      job->bse_ham = 2;
      job->bse_tda = 1;
      job->npoints  = 1000;    // Number of samples of self-energy in energy range
      energy_range[0] = -50.0; // Energy range for self-energy calculation
      energy_range[1] =  50.0;
      
      // Read directory for MPI File I/O

      read_line(file.job, title, 99);
      sscanf(title, "%s", file.scf_eigvec);

      // Read options for GW calculation

    do {
 
       read_line(file.job, title, 99);
       sscanf(title, "%s", jobname1);   
  
      if (!strcmp(jobname1, "MONKHORST_PACK")) {   // read in Monkhorst-Pack net for periodic system
   
        if (crystal->type[0] != 'M') {
          read_line(file.job, title, 99);
          sscanf(title, "%d %d %d", &fermi.is[0], &fermi.is[1], &fermi.is[2]);  // periodic systems only
          knet_check(fermi.is,jobname,crystal,file);
          //fprintf(file.out,"IS is %3d %3d %3d\n",fermi.is[0],fermi.is[1],fermi.is[2]);
         }
      
	else { 
          read_line(file.job, title, 99);
          sscanf(title, "%d %d %d", &fermi.is[0], &fermi.is[1], &fermi.is[2]);  // ignore Monkhorst-Pack input for molecules
	  fermi.is[0] = 1; fermi.is[1] =  1; fermi.is[2] = 1; 
	 }

    }

      if (!strcmp(jobname1, "RPA_VECTORS")) {   // use a limited number of CASIDA eigenvectors in GW calculation

      read_line(file.job, title, 99);
      sscanf(title, "%d", &job->bse_lim);

    }

       if (!strcmp(jobname1, "MO_RANGE")) {

       read_line(file.job, title, 99);
       sscanf(title, "%d %d", &fermi.bands[0], &fermi.bands[1]);

     }
   
       if (!strcmp(jobname1, "SPECTRUM")) {

       read_line(file.job, title, 99);
       sscanf(title, "%lf %lf %d", &energy_range[0], &energy_range[1], &job->npoints);
       read_line(file.job, title, 99);
       sscanf(title, "%d", &job->nspectra);
       job->self_plot = 1;

       if (job->nspectra > 10) {
       if (job->taskid == 0)
       fprintf(file.out,"NUMBER OF SELF-ENERGY PLOTS %3d MUST NOT EXCEED 10\n", job->nspectra);
       MPI_Finalize();
       exit(1);
      }
       read_line(file.job, title, 99);
       sscanf(title, "%d %d %d %d %d %d %d %d %d %d", \
       &self[0],&self[1], &self[2], &self[3], &self[4], &self[5], &self[6], &self[7], &self[8], &self[9]);

       for (i = 0; i < job->nspectra; i++) {
         if ((self[i] < fermi.bands[0] || self[i] > fermi.bands[1]) && job->spin_polarisation == 0) {
         if (job->taskid == 0) 
         fprintf(file.out,"RANGE OF BANDS FOR SELF-ENERGY PLOT MUST LIE WITHIN MO RANGE %d %d\n",fermi.bands[0],fermi.bands[1]);
         MPI_Finalize();
         exit(0);
        }

         else if ((self[i] < fermi.bands[0] || self[i] > fermi.bands[1] || self[i] < fermi.bands[2] || self[i] > fermi.bands[3]) \
         && job->spin_polarisation == 1) { 
         if (job->taskid == 0) 
         fprintf(file.out,"RANGE OF BANDS FOR SELF-ENERGY PLOT MUST LIE WITHIN MO RANGE %d %d and %d %d\n", \
         fermi.bands[0],fermi.bands[1],fermi.bands[2],fermi.bands[3]);
         MPI_Finalize();
         exit(0);
        }

         fermi.plot_bands[i] = self[i];
       }

     }
  
     } while (strcmp(jobname1, "END"));

      // Allocate fermi array and calculate nbands, etc.

      fermi.nkunique = 1;
      allocate_fermi(&fermi,atoms,job,file);
      fermi.occupied[0] = fermi.homo[0] - fermi.bands[0] + 1;
      nbands = fermi.bands[1] - fermi.bands[0] + 1;
      nocc = fermi.occupied[0];
      nvir = nbands - fermi.occupied[0];
      ntransitions = nocc * nvir;

      // Check input

      knet_check(fermi.is,jobname,crystal,file);
      int_range_check(fermi.bands,int_value,int_range,jobname,file);
      double_range_check(energy_range,double_value,double_range,jobname,file);
      
      // Arrays for string output

      if (fermi_homo[0] == 0)  sprintf(fermi_homo, "       %d",fermi.homo[0]);
      if (fermi_bands[0] == 0) sprintf(fermi_bands,"  %d - %d",fermi.bands[0],fermi.bands[1]);
      if (job->spin_polarisation == 1 && fermi_bands1[0] == 0) sprintf(fermi_bands1,"  %d - %d",fermi.bands[2],fermi.bands[3]);
      if (nspec[0] == 0)       sprintf(nspec,"      %d",job->nspectra);
      if (spectrum[0] == 0)    sprintf(spectrum,"%4.1lf - %4.1lf",energy_range[0],energy_range[1]);
      if (ntrans[0] == 0)      sprintf(ntrans,"   %d",ntransitions);

      if (job->taskid == 0) {
   
      fprintf(file.out,"\n\n===========================================================================================================\n");
      if (crystal->type[0] == 'M') {
      fprintf(file.out,"| GW CALCULATION            | HOMO LEVEL    %8s | MO     %16s | NTRANS     %12s |\n", \
      fermi_homo, fermi_bands, ntrans);
      if (job->self_plot == 1) {
      fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
      fprintf(file.out,"| SELF-ENERGY PLOTTED       | NBANDS         %8s | E RANGE  %14s |                         |\n", \
      nspec,spectrum);
     }
    }
      if (crystal->type[0] == 'C' || crystal->type[0] == 'S' || crystal->type[0] == 'P') {
      fprintf(file.out,"| GW CALCULATION            | FERMI LEVEL    %8s | BANDS  %16s | MONKHORST-PACK %8s |\n", \
      fermi_homo,fermi_bands, monkhorst_pack);
      //fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
     }
      fprintf(file.out,"===========================================================================================================\n");
      fflush(file.out);

  }

      // Convert input to atomic units

      job->energy_range[0] = energy_range[0] / au_to_eV;
      job->energy_range[1] = energy_range[1] / au_to_eV;

        switch (crystal->type[0]) {

    case 'C':
    case 'S':
    case 'P':

    break;

    case 'M':

       gw_molecule(&fermi,atoms,atom_p,shells,gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);

    break;

   } // close switch

       free_fermi(&fermi,job);
   }

    // *****JOB: RPA CALCULATION ***************************************************************************************

    else if (!strcmp(jobname, "RPA")) {

      int_value[0] = 1;
      int_value[1] = atoms->number_of_sh_bfns_in_unit_cell;

      // Initialise char arrays

      fermi_homo[0] = 0; 
      fermi_bands[0] = 0;
      fermi_bands1[0] = 0;
      ntrans[0] = 0;
      monkhorst_pack[0] = 0;
      int_range[0] = 0;

      sprintf(int_range, "%11s", "ENERGY BAND");

      // Set job defaults

      fermi_default(&fermi,crystal,atoms,job,file);
      job->bse_int = 2;
      job->bse_ham = 2;
      job->bse_tda = 1;
      
      // Read directory for MPI File I/O

   if (crystal->type[0] != 'M') {
      read_line(file.job, title, 99);
      sscanf(title, "%d %d %d", &fermi.is[0], &fermi.is[1], &fermi.is[2]);
      knet_check(fermi.is,jobname,crystal,file);
      //fprintf(file.out,"IS is %3d %3d %3d\n",fermi.is[0],fermi.is[1],fermi.is[2]);
     }
   else { fermi.is[0] = 1; fermi.is[1] = 1; fermi.is[2] = 1; }

      read_line(file.job, title, 99);
      sscanf(title, "%s", file.scf_eigvec);
  
      // Read options for RPA calculation

    do {
 
       read_line(file.job, title, 99);
       sscanf(title, "%s", jobname1);
       //printf("%s\n",jobname1);
  
      if (!strcmp(jobname1, "RPA_VECTORS")) {
      read_line(file.job, title, 99);
      sscanf(title, "%d", &job->bse_lim);
      if (job->taskid == 0) {
     }
    }

       if (!strcmp(jobname1, "MO_RANGE")) {
       read_line(file.job, title, 99);
       sscanf(title, "%d %d", &fermi.bands[0], &fermi.bands[1]);
     }
   
     } while (strcmp(jobname1, "END"));

      // Allocate fermi array and calculate nbands, etc.

      fermi.nkunique = 1;
      allocate_fermi(&fermi,atoms,job,file);
      fermi.occupied[0] = fermi.homo[0] - fermi.bands[0] + 1;
      nbands = fermi.bands[1] - fermi.bands[0] + 1;
      nocc = fermi.occupied[0];
      nvir = nbands - fermi.occupied[0];
      ntransitions = nocc * nvir;

      // Check input

      knet_check(fermi.is,jobname,crystal,file);
      int_range_check(fermi.bands,int_value,int_range,jobname,file);

      // Arrays for string output

      if (fermi_homo[0] == 0)  sprintf(fermi_homo, "       %d",fermi.homo[0]);
      if (fermi_bands[0] == 0) sprintf(fermi_bands,"  %d - %d",fermi.bands[0],fermi.bands[1]);
      //if (job->spin_polarisation == 1 && fermi_bands1[0] == 0) sprintf(fermi_bands1,"  %d - %d",fermi.bands[2],fermi.bands[3]);
      if (ntrans[0] == 0)      sprintf(ntrans,"   %d",ntransitions);

   if (job->taskid == 0) {

   fprintf(file.out,"\n\n===========================================================================================================\n");
   if (crystal->type[0] == 'M') {
   fprintf(file.out,"| RPA CALCULATION           | HOMO LEVEL    %8s | MO     %16s | NTRANS     %12s |\n", \
   fermi_homo, fermi_bands, ntrans);
  }
   if (crystal->type[0] == 'C' || crystal->type[0] == 'S' || crystal->type[0] == 'P') {
   fprintf(file.out,"| GW CALCULATION            | FERMI LEVEL    %8s | BANDS  %16s | MONKHORST-PACK %8s |\n", \
   fermi_homo,fermi_bands, monkhorst_pack);
   //fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
  }
   fprintf(file.out,"===========================================================================================================\n");

   fflush(file.out);

  }
      // Convert input to atomic units

        switch (crystal->type[0]) {

    case 'C':

    break;

    case 'M':

       rpa_molecule(&fermi,atoms,atom_p,shells,gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);

    break;

   } // close switch

       free_fermi(&fermi,job);

    }

    // *****JOB: BETHE-SALPETER EQUATION HAMILTONIAN DIAGONALISATION AND OPTICAL SPECTRUM *****************************

    else if (!strcmp(jobname, "BSE")) {

      double_value[0] = lower_energy_limit;
      double_value[1] = upper_energy_limit;
      int_value[0] = 1;
      int_value[1] = atoms->number_of_sh_bfns_in_unit_cell;

      // Initialise char arrays

      hamiltonian[0] = 0;
      fermi_homo[0] = 0; 
      fermi_bands[0] = 0;
      fermi_bands1[0] = 0;
      monkhorst_pack[0] = 0;
      field[0] = 0;
      tamm_dancoff[0] = 0;
      spin_state[0] = 0;
      scissor_shift[0] = 0;
      scale_factor[0] = 0;
      ntrans[0] = 0;
      linear_algebra[0] = 0;
      int_range[0] = 0;
      double_range[0] = 0;

      sprintf(int_range, "%11s", "ENERGY BAND");
      sprintf(double_range, "%11s", "ENERGY RANGE");

      // Set job defaults

      job->type = 1;

      job->guess_type = 1;
      job->vectors    = 2;
      job->values     = 2;

      job->kpoints    = 0;

      job->bse_ham    = 1;     // default is BSE 
      job->bse_tda    = 0;     // default is full BSE matrix 
      job->bse_spin   = 0;     // default is spin singlet in BSE and similar calculations
      job->bse_spk    = 0;     // default is not to use SCALAPACK in BSE and TDHF calculations
      job->bse_int    = 0;     // default is integrals need to be calculated and written to disk
      job->bse_cou    = 0;     // default is do not calculate Coulomb energy of occupied states in MO RANGE
      job->bse_exc    = 0;     // default is do not calculate exchange energy of occupied states in MO RANGE
      job->bse_lim    = 0;     // default is include all transitions implied by MO RANGE in BSE calculation

      fermi_default(&fermi,crystal,atoms,job,file);

      job->npoints = 400;                                 // Default range for optical spectrum is 0 - 10 eV with 200 points
      energy_range[0] = 0.01;
      energy_range[1] = 20.0;

      job->field_dirs = 3;                                // Default is fields along Cartesian x, y and z
      job->e_field[0].comp1 = k_one;
      job->e_field[0].comp2 = k_zero;
      job->e_field[0].comp3 = k_zero;
      job->e_field[1].comp1 = k_zero;
      job->e_field[1].comp2 = k_one;
      job->e_field[1].comp3 = k_zero;
      job->e_field[2].comp1 = k_zero;
      job->e_field[2].comp2 = k_zero;
      job->e_field[2].comp3 = k_one;

      job->scissor   = 0.00 ; // default virtual state shift for periodic TDHF
      job->scalefac  = k_one; // default scaling factor for electron-hole attraction matrix elements in periodic TDHF

      // Read directory for MPI File I/O

      read_line(file.job, title, 99);
      sscanf(title, "%s", file.scf_eigvec);
  
      // Read options for BSE calculation

   do {

      read_line(file.job, title, 99);
      sscanf(title, "%s", jobname1);

      if (!strcmp(jobname1, "COULOMB_TEST")) {
      job->bse_cou = 1;
     }

      if (!strcmp(jobname1, "EXCHANGE_TEST")) {
      job->bse_exc = 1;
     }

      if (!strcmp(jobname1, "SCISSOR")) {
      read_line(file.job, title, 99);
      sscanf(title, "%lf", &job->scissor);
     }

      if (!strcmp(jobname1, "SCALEFACTOR")) {
      read_line(file.job, title, 99);
      sscanf(title, "%lf", &job->scalefac);
     }

      if (!strcmp(jobname1, "TDA")) {
      job->bse_tda = 1;
    }

      if (!strcmp(jobname1, "TRIPLET")) {
      job->bse_spin = 1;
    }

      if (!strcmp(jobname1, "TDHF")) {
      job->bse_ham = 0;
    }

      if (!strcmp(jobname1, "BSE_VECTORS")) {
      read_line(file.job, title, 99);
      sscanf(title, "%d", &job->bse_lim);
    }

      if (!strcmp(jobname1, "MONKHORST_PACK")) {   // read in Monkhorst-Pack net for periodic system
   
        if (crystal->type[0] != 'M') {
          read_line(file.job, title, 99);
          sscanf(title, "%d %d %d", &fermi.is[0], &fermi.is[1], &fermi.is[2]);  // periodic systems only
          knet_check(fermi.is,jobname,crystal,file);
         }
      
	else { 
          read_line(file.job, title, 99);
          sscanf(title, "%d %d %d", &fermi.is[0], &fermi.is[1], &fermi.is[2]);  // ignore Monkhorst-Pack input for molecules
	  fermi.is[0] = 1; fermi.is[1] =  1; fermi.is[2] = 1; 
	 }

    }

      if (!strcmp(jobname1, "FERMI_LEVEL")) {

      read_line(file.job, title, 99);
      if (job->spin_polarisation == 0) {
      sscanf(title, "%d", &fermi.homo[0]);
     }
      else if (job->spin_polarisation == 1) {
      sscanf(title, "%d %d", &fermi.homo[0],&fermi.homo[1]);
     }
    }

      if (!strcmp(jobname1, "MO_RANGE")) {
      read_line(file.job, title, 99);
      sscanf(title, "%d %d", &fermi.bands[0], &fermi.bands[1]);
      if (job->spin_polarisation == 1) {
      read_line(file.job, title, 99);
      sscanf(title, "%d %d", &fermi.bands[2], &fermi.bands[3]);
     }
    }

      if (!strcmp(jobname1, "SPECTRUM")) {
      read_line(file.job, title, 99);
      sscanf(title, "%d %lf %lf", &job->npoints, &energy_range[0], &energy_range[1]);
    }

      if (!strcmp(jobname1, "LINEWIDTH")) {
      read_line(file.job, title, 99);
      sscanf(title, "%lf", &job->linewidth);
    }

      if (!strcmp(jobname1, "FIELD")) {
      read_line(file.job, title, 99);
      sscanf(title, "%d", &job->field_dirs);
      if (job->taskid == 0 && job->field_dirs > 3) {
      fprintf(file.out,"Number of field directions %d exceeds maximum of 3\n",job->field_dirs);
      MPI_Finalize();
      exit(1);
     }
      VECTOR_DOUBLE fld[job->field_dirs];
      for (i = 0; i < job->field_dirs; i++) {
      read_line(file.job, title, 99);
      sscanf(title, "%lf %lf %lf", &fld[i].comp1,  &fld[i].comp2, &fld[i].comp3);
      double fld_tmp = sqrt(double_vec_dot(&fld[i],&fld[i]));
      fld[i].comp1 /= fld_tmp;
      fld[i].comp2 /= fld_tmp;
      fld[i].comp3 /= fld_tmp;
      job->e_field[i].comp1 = fld[i].comp1;
      job->e_field[i].comp2 = fld[i].comp2;
      job->e_field[i].comp3 = fld[i].comp3;
     }
    }

      if (!strcmp(jobname1, "FRAGMENT_DIPOLE")) {
      for (j = 0; j < numfrag; j++) natoms[j] = 0;
      for (j = 0; j < max_atoms; j++) nat1[j] = 0;
      read_line(file.job, title, 99);
      sscanf(title, "%d", &numfrag);
      read_line(file.job, title, 99);
      sscanf(title, "%d %d", &natoms[0], &natoms[1]); // change with numfrag
      for (i = 0; i < numfrag; i++) {
      read_line(file.job, title, 99);
      // change with max_atoms
      sscanf(title, "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d", \
      &nat1[0],&nat1[1], &nat1[2], &nat1[3], &nat1[4], &nat1[5], &nat1[6], &nat1[7], &nat1[8], &nat1[9], \
      &nat1[10], &nat1[11], &nat1[12], &nat1[13], &nat1[14], &nat1[15], &nat1[16], &nat1[17], &nat1[18], &nat1[19]); 
      if (natoms[i] > 20) {
      read_line(file.job, title, 99);
      sscanf(title, "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d", \
      &nat1[20], &nat1[21], &nat1[22], &nat1[23], &nat1[24], &nat1[25], &nat1[26], &nat1[27], &nat1[28], &nat1[29], \
      &nat1[30], &nat1[31], &nat1[32], &nat1[33], &nat1[34], &nat1[35], &nat1[36], &nat1[37], &nat1[38], &nat1[39]); 
     }
      for (j = 0; j < natoms[i]; j++) {
      if (nat1[j] < 1 || nat1[j] > atoms->number_of_atoms_in_unit_cell) {
      if (job->taskid == 0) {
      fprintf(file.out,"ERROR: Atom [%3d] in fragment [%2d] in BSE must lie beween 1 and number of atoms in unit cell\n",nat1[j],i);
     }
      MPI_Finalize();
      exit(0);
     }
      nat[j][i] = nat1[j] - 1;
    }
   }
  }

    } while (strcmp(jobname1, "END"));

      // Allocate fermi array and calculate nbands, etc.

      fermi.nkunique = 1;
      allocate_fermi(&fermi,atoms,job,file);
      fermi.occupied[0] = fermi.homo[0] - fermi.bands[0] + 1;
      nbands = fermi.bands[1] - fermi.bands[0] + 1;
      nocc = fermi.occupied[0];
      nvir = nbands - fermi.occupied[0];
      ntransitions = nocc * nvir;

      // Check input

      knet_check(fermi.is,jobname,crystal,file);

      int_range_check(fermi.bands,int_value,int_range,jobname,file);

      double_range_check(energy_range,double_value,double_range,jobname,file);

      if (job->linewidth <= 0.0049 || job->linewidth >= 1.0001) {
      if (job->taskid == 0)
      fprintf(file.out,"\nLINEWIDTH %9.2e eV OUT OF RANGE: CHOOSE VALUE BETWEEN 0.005 AND 1.0 eV\n",job->linewidth);
      MPI_Finalize();
      exit(1);
     }

      // Arrays for string output

      sprintf(field," %d  %4.1lf %4.1lf %4.1lf    %4.1lf %4.1lf %4.1lf      %4.1lf %4.1lf %4.1lf",job->field_dirs, \
      job->e_field[0].comp1, job->e_field[0].comp2, job->e_field[0].comp3, \
      job->e_field[1].comp1, job->e_field[1].comp2, job->e_field[1].comp3, \
      job->e_field[2].comp1, job->e_field[2].comp2, job->e_field[2].comp3);

      sprintf(spectrum,"%4.1lf - %4.1lf",energy_range[0],energy_range[1]);

      if (job->spin_polarisation == 0) {
      sprintf(fermi_homo, "     %3d",fermi.homo[0]);
      sprintf(fermi_bands,"  %d - %d",fermi.bands[0],fermi.bands[1]);
     }
      else if (job->spin_polarisation == 1) {
      sprintf(fermi_homo, "%3d %3d",fermi.homo[0],fermi.homo[1]);
      sprintf(fermi_bands,"  %d - %d %d - %d",fermi.bands[0],fermi.bands[1],fermi.bands[2],fermi.bands[3]);
     }

      sprintf(scale_factor,"%8.2f",job->scalefac);

      sprintf(scissor_shift,"%8.2f",job->scissor);

      sprintf(linear_algebra,"%s","SCALAPACK");

      if (job->bse_ham == 0) 
      sprintf(hamiltonian,"%s","TDHF");

      else if (job->bse_ham == 1) 
      sprintf(hamiltonian,"%s","BSE");

      if (job->bse_tda == 0) 
      sprintf(tamm_dancoff,"%s","OFF");

      else if (job->bse_tda == 1) 
      sprintf(tamm_dancoff,"%s","ON");

      if (job->bse_spin == 0) 
      sprintf(spin_state,"%s","SINGLET");

      else if (job->bse_spin == 1) 
      sprintf(spin_state,"%s","TRIPLET");

      sprintf(monkhorst_pack,"%d %d %d",fermi.is[0], fermi.is[1], fermi.is[2]);

      // Convert input to atomic units

      job->scissor /= au_to_eV;
      job->linewidth /= au_to_eV;
      job->energy_range[0] = energy_range[0] / au_to_eV;
      job->energy_range[1] = energy_range[1] / au_to_eV;

   if (job->taskid == 0) {
   fprintf(file.out,"\n\n===========================================================================================================\n");
   if (crystal->type[0] == 'M') {
   fprintf(file.out,"| BSE CALCULATION           | HOMO LEVEL     %8s | RANGE  %16s | MONKHORST-PACK %8s |\n", \
   fermi_homo,fermi_bands, monkhorst_pack);
   fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
   fprintf(file.out,"| HAMILTONIAN      %8s | TAMM-DANCOFF   %8s | SPIN STATE     %8s | LINEAR ALG.   %9s |\n", \
   hamiltonian,tamm_dancoff,spin_state,linear_algebra);
  }
   if (crystal->type[0] == 'C' || crystal->type[0] == 'S' || crystal->type[0] == 'P') {
   fprintf(file.out,"| BSE CALCULATION           | FERMI LEVEL    %8s | RANGE  %16s | MONKHORST-PACK %8s |\n", \
   fermi_homo,fermi_bands, monkhorst_pack);
   fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
   fprintf(file.out,"| TDHF HAMILTONIAN          | VIRT SHIFT     %8s | SCALEFACTOR  %10s |                         | \n", \
   scissor_shift,scale_factor);
  }
   fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
   fprintf(file.out,"| RANGE (eV) %14s | FIELDS   %66s |\n", \
   spectrum,field);
   fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
   fflush(file.out);
  }

        switch (crystal->type[0]) {
        case 'C':

        //bse_crystal1(&fermi,atoms,atom_p,&numfrag,natoms,nat,shells,gaussians,atoms_ax,shells_ax,gaussians_ax,crystal, \
        symmetry,R,R_tables,G,job,file);

        //optical_spectrum_crystal2(&fermi,atoms,atom_p,&numfrag,natoms,nat,shells,gaussians,atoms_ax,shells_ax,gaussians_ax, \
        crystal,symmetry,R,R_tables,G,job,file);

        case 'S':
        case 'P':

        break;

        case 'M':

        bse_molecule(&fermi,atoms,atom_p,&numfrag,natoms,nat,shells,gaussians,atoms_ax,shells_ax,gaussians_ax,crystal, \
        symmetry,R,R_tables,G,job,file);

        optical_spectrum_molecule(&fermi,atoms,atom_p,&numfrag,natoms,nat,shells,gaussians,atoms_ax,shells_ax,gaussians_ax, \
        crystal,symmetry,R,R_tables,G,job,file);

        break;

       }

      free_fermi(&fermi,job);

    }

    // *****JOB: BETHE-SALPETER EQUATION HAMILTONIAN DIAGONALISATION AND OPTICAL SPECTRUM *****************************

    else if (!strcmp(jobname, "OPTICAL_SPECTRUM")) {

      double_value[0] = lower_energy_limit;
      double_value[1] = upper_energy_limit;
      int_value[0] = 1;
      int_value[1] = atoms->number_of_sh_bfns_in_unit_cell;

      // Initialise char arrays

      hamiltonian[0] = 0;
      fermi_homo[0] = 0; 
      fermi_bands[0] = 0;
      fermi_bands1[0] = 0;

      // Set job defaults

      //job->vectors    = 2;
      //job->values     = 2;

      //job->kpoints    = 0;

      job->bse_ham    = 1;     // default is BSE 
      job->bse_tda    = 0;     // default is full BSE matrix 
      job->bse_spin   = 0;     // default is spin singlet in BSE and similar calculations
      job->bse_spk    = 1;     // default is to use SCALAPACK in BSE and TDHF calculations
      job->bse_cou    = 0;     // default is do not calculate Coulomb energy of occupied states in MO RANGE
      job->bse_exc    = 0;     // default is do not calculate exchange energy of occupied states in MO RANGE
      job->bse_lim    = 0;     // default is include all transitions implied by MO RANGE in BSE calculation

      fermi_default(&fermi,crystal,atoms,job,file);

      job->npoints = 400;                                 // Default range for optical spectrum is 0 - 10 eV with 200 points
      energy_range[0] = 0.01;
      energy_range[1] = 20.0;

      job->field_dirs = 3;                                // Default is fields along Cartesian x, y and z
      job->e_field[0].comp1 = k_one;
      job->e_field[0].comp2 = k_zero;
      job->e_field[0].comp3 = k_zero;
      job->e_field[1].comp1 = k_zero;
      job->e_field[1].comp2 = k_one;
      job->e_field[1].comp3 = k_zero;
      job->e_field[2].comp1 = k_zero;
      job->e_field[2].comp2 = k_zero;

      // Read directory for MPI File I/O

      read_line(file.job, title, 99);
      sscanf(title, "%s", file.scf_eigvec);
  
      // Read options for OPTICAL SPECTRUM calculation

   do {

      read_line(file.job, title, 99);
      sscanf(title, "%s", jobname1);

      if (!strcmp(jobname1, "TDA")) {
      job->bse_tda = 1;
    }

      if (!strcmp(jobname1, "TRIPLET")) {
      job->bse_spin = 1;
    }

      if (!strcmp(jobname1, "TDHF")) {
      job->bse_ham = 0;
    }

      if (!strcmp(jobname1, "BSE_VECTORS")) {
      read_line(file.job, title, 99);
      sscanf(title, "%d", &job->bse_lim);
    }

      if (!strcmp(jobname1, "MONKHORST_PACK")) {   // read in Monkhorst-Pack net for periodic system
   
        if (crystal->type[0] != 'M') {
          read_line(file.job, title, 99);
          sscanf(title, "%d %d %d", &fermi.is[0], &fermi.is[1], &fermi.is[2]);  // periodic systems only
          knet_check(fermi.is,jobname,crystal,file);
         }
      
	else { 
          read_line(file.job, title, 99);
          sscanf(title, "%d %d %d", &fermi.is[0], &fermi.is[1], &fermi.is[2]);  // ignore Monkhorst-Pack input for molecules
	  fermi.is[0] = 1; fermi.is[1] =  1; fermi.is[2] = 1; 
	 }

    }

      if (!strcmp(jobname1, "FERMI_LEVEL")) {

      read_line(file.job, title, 99);
      if (job->spin_polarisation == 0) {
      sscanf(title, "%d", &fermi.homo[0]);
     }
      else if (job->spin_polarisation == 1) {
      sscanf(title, "%d %d", &fermi.homo[0],&fermi.homo[1]);
     }
    }

      if (!strcmp(jobname1, "MO_RANGE")) {
      read_line(file.job, title, 99);
      sscanf(title, "%d %d", &fermi.bands[0], &fermi.bands[1]);
      if (job->spin_polarisation == 1) {
      read_line(file.job, title, 99);
      sscanf(title, "%d %d", &fermi.bands[2], &fermi.bands[3]);
     }
    }

      if (!strcmp(jobname1, "SPECTRUM")) {
      read_line(file.job, title, 99);
      sscanf(title, "%d %lf %lf", &job->npoints, &energy_range[0], &energy_range[1]);
    }

      if (!strcmp(jobname1, "LINEWIDTH")) {
      read_line(file.job, title, 99);
      sscanf(title, "%lf", &job->linewidth);
    }

      if (!strcmp(jobname1, "FIELD")) {
      read_line(file.job, title, 99);
      sscanf(title, "%d", &job->field_dirs);
      if (job->taskid == 0 && job->field_dirs > 3) {
      fprintf(file.out,"Number of field directions %d exceeds maximum of 3\n",job->field_dirs);
      MPI_Finalize();
      exit(1);
     }
      VECTOR_DOUBLE fld[job->field_dirs];
      for (i = 0; i < job->field_dirs; i++) {
      read_line(file.job, title, 99);
      sscanf(title, "%lf %lf %lf", &fld[i].comp1,  &fld[i].comp2, &fld[i].comp3);
      double fld_tmp = sqrt(double_vec_dot(&fld[i],&fld[i]));
      fld[i].comp1 /= fld_tmp;
      fld[i].comp2 /= fld_tmp;
      fld[i].comp3 /= fld_tmp;
      job->e_field[i].comp1 = fld[i].comp1;
      job->e_field[i].comp2 = fld[i].comp2;
      job->e_field[i].comp3 = fld[i].comp3;
     }
    }

    } while (strcmp(jobname1, "END"));

      // Allocate fermi array and calculate nbands, etc.

      fermi.nkunique = 1;
      allocate_fermi(&fermi,atoms,job,file);
      fermi.occupied[0] = fermi.homo[0] - fermi.bands[0] + 1;
      nbands = fermi.bands[1] - fermi.bands[0] + 1;
      nocc = fermi.occupied[0];
      nvir = nbands - fermi.occupied[0];
      ntransitions = nocc * nvir;

      // Check input

      // Arrays for string output

      sprintf(field," %d  %4.1lf %4.1lf %4.1lf    %4.1lf %4.1lf %4.1lf      %4.1lf %4.1lf %4.1lf",job->field_dirs, \
      job->e_field[0].comp1, job->e_field[0].comp2, job->e_field[0].comp3, \
      job->e_field[1].comp1, job->e_field[1].comp2, job->e_field[1].comp3, \
      job->e_field[2].comp1, job->e_field[2].comp2, job->e_field[2].comp3);

      sprintf(spectrum,"%4.1lf - %4.1lf",energy_range[0],energy_range[1]);

      if (job->spin_polarisation == 0) {
      sprintf(fermi_homo, "     %3d",fermi.homo[0]);
      sprintf(fermi_bands,"  %d - %d",fermi.bands[0],fermi.bands[1]);
     }
      else if (job->spin_polarisation == 1) {
      sprintf(fermi_homo, "%3d %3d",fermi.homo[0],fermi.homo[1]);
      sprintf(fermi_bands,"  %d - %d %d - %d",fermi.bands[0],fermi.bands[1],fermi.bands[2],fermi.bands[3]);
     }

      sprintf(linear_algebra,"%s","SCALAPACK");

      if (job->bse_ham == 0) 
      sprintf(hamiltonian,"%s","TDHF");

      else if (job->bse_ham == 1) 
      sprintf(hamiltonian,"%s","BSE");

      if (job->bse_tda == 0) 
      sprintf(tamm_dancoff,"%s","OFF");

      else if (job->bse_tda == 1) 
      sprintf(tamm_dancoff,"%s","ON");

      if (job->bse_spin == 0) 
      sprintf(spin_state,"%s","SINGLET");

      else if (job->bse_spin == 1) 
      sprintf(spin_state,"%s","TRIPLET");

      sprintf(monkhorst_pack,"%d %d %d",fermi.is[0], fermi.is[1], fermi.is[2]);

      // Convert input to atomic units

      job->scissor /= au_to_eV;
      job->linewidth /= au_to_eV;
      job->energy_range[0] = energy_range[0] / au_to_eV;
      job->energy_range[1] = energy_range[1] / au_to_eV;

   if (job->taskid == 0) {
   fprintf(file.out,"\n\n===========================================================================================================\n");
   if (crystal->type[0] == 'M') {
   fprintf(file.out,"| OPTICAL SPECTRUM          | HOMO LEVEL     %8s | RANGE  %16s | MONKHORST-PACK %8s |\n", \
   fermi_homo,fermi_bands, monkhorst_pack);
   fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
   fprintf(file.out,"| HAMILTONIAN      %8s | TAMM-DANCOFF   %8s | SPIN STATE     %8s | LINEAR ALG.   %9s |\n", \
   hamiltonian,tamm_dancoff,spin_state,linear_algebra);
  }
   if (crystal->type[0] == 'C' || crystal->type[0] == 'S' || crystal->type[0] == 'P') {
   fprintf(file.out,"| OPTICAL SPECTRUM          | FERMI LEVEL    %8s | RANGE  %16s | MONKHORST-PACK %8s |\n", \
   fermi_homo,fermi_bands, monkhorst_pack);
   fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
   fprintf(file.out,"| TDHF HAMILTONIAN          | VIRT SHIFT     %8s | SCALEFACTOR  %10s |                         | \n", \
   scissor_shift,scale_factor);
  }
   fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
   fprintf(file.out,"| RANGE (eV) %14s | FIELDS   %66s |\n", \
   spectrum,field);
   fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
   fflush(file.out);
  }

        switch (crystal->type[0]) {
        case 'C':

        //optical_spectrum_crystal2(&fermi,atoms,atom_p,&numfrag,natoms,nat,shells,gaussians,atoms_ax,shells_ax,gaussians_ax, \
        crystal,symmetry,R,R_tables,G,job,file);

        case 'S':
        case 'P':

        break;

        case 'M':

        optical_spectrum_molecule(&fermi,atoms,atom_p,&numfrag,natoms,nat,shells,gaussians,atoms_ax,shells_ax,gaussians_ax, \
        crystal,symmetry,R,R_tables,G,job,file);

        break;

       }

      free_fermi(&fermi,job);

    }

    // *****JOB: BETHE-SALPETER EQUATION WAVEFUNCTION PLOTTING ********************************************************

    else if (!strcmp(jobname, "BSE_PLOT") || !strcmp(jobname, "BSE_TRANS_DEN") || !strcmp(jobname, "ELECTRON_HOLE_PLOT")) {

      // Initialise char arrays

      int bse_state[2];
      int grid_type, plot_dim = 0;
      int grid_par[3];
      int num_points = 0;
      double grid_sum[job->spin_dim];
      VECTOR_DOUBLE points[4],points1[10];

      plot_type[0] = 0;
      fermi_homo[0] = 0;
      fermi_bands[0] = 0;

      // Set job defaults

      fermi_default(&fermi,crystal,atoms,job,file);

      // Read directory for MPI File I/O

      read_line(file.job, title, 99);
      sscanf(title, "%s", file.scf_eigvec);

      // Read integer dimensions of plotting grid

      read_line(file.job, title, 99);
      sscanf(title, "%d %d %d", &grid_par[0], &grid_par[1], &grid_par[2]);

      for (i = 0; i < 3; i++) {
        if (grid_par[i] > 1)
          plot_dim++;
      }

      // Read vectors defining GRIDVEC plotting grid 

      read_line(file.job, title, 99); 
      sscanf(title, "%8s", title2);

      if (!strcmp(title2, "GRIDVEC")) {

        grid_type = 0;

        //
        // Insert the points which define the plotting region
        // points[0] == ORIGIN
        // axes may be non-orthogonal
        // For the plane (grid_par[2]=1) the points[4] will be ignored
        // 1D     points[3] will be ignored
        //

        for (i = 0; i <= plot_dim; i++) {
          read_line(file.job, title, 99);
          sscanf(title, "%lf %lf %lf", &points[i].comp1, &points[i].comp2, &points[i].comp3);
          points[i].comp1 /= bohr_to_AA;
          points[i].comp2 /= bohr_to_AA;
          points[i].comp3 /= bohr_to_AA;
        }
      } // end of GRIDVEC input

      // Read options for BSE_PLOT or ELECTRON_HOLE_PLOT calculations

   do {

      read_line(file.job, title, 99);
      sscanf(title, "%s", jobname1);

      if (!strcmp(jobname1, "MONKHORST_PACK_NET")) {
      read_line(file.job, title, 99);
      switch (crystal->type[0]) {
      case 'C':
      sscanf(title, "%d %d %d", &fermi.is[0], &fermi.is[1], &fermi.is[2]);
      if (fermi.is[0] < 1 || fermi.is[0] > 32 || fermi.is[1] < 1 || fermi.is[1] > 32 || fermi.is[2] < 1 || fermi.is[2] > 32) {
      fprintf(file.out,"ERROR: Monkhorst-Pack Net integers [%2d,%2d,%2d] must range between 1 and 32\n",\
      fermi.is[0],fermi.is[1],fermi.is[2]);
      MPI_Finalize();
      exit(1);
     }
      break;
      case 'S':
      sscanf(title, "%d %d", &fermi.is[0], &fermi.is[1]);
      if (fermi.is[0] < 1 || fermi.is[0] > 32 || fermi.is[1] < 1 || fermi.is[1] > 32) {
      fprintf(file.out,"ERROR: Monkhorst-Pack Net integers [%2d,%2d] must range between 1 and 32\n",fermi.is[0],fermi.is[1]);
      MPI_Finalize();
      exit(1);
     }
      break;
      case 'P':
      sscanf(title, "%d", &fermi.is[2]);
      if (fermi.is[2] < 1 || fermi.is[2] > 32) {
      fprintf(file.out,"ERROR: Monkhorst-Pack Net integer [%2d] must range between 1 and 32\n",fermi.is[2]);
      MPI_Finalize();
      exit(1);
     }
      break;
      case 'M':
      read_line(file.job, title, 99);
      sscanf(title, "%d %d %d", &fermi.is[0], &fermi.is[1], &fermi.is[2]);  // ignore Monkhorst-Pack input for molecules
      break;
     }
    }

      if (!strcmp(jobname1, "FERMI_LEVEL")) {

      read_line(file.job, title, 99);
      if (job->spin_polarisation == 0) {
      sscanf(title, "%d", &fermi.homo[0]);
     }
      else if (job->spin_polarisation == 1) {
      sscanf(title, "%d %d", &fermi.homo[0],&fermi.homo[1]);
     }
    }

      if (!strcmp(jobname1, "MO_RANGE")) {
      read_line(file.job, title, 99);
      sscanf(title, "%d %d", &fermi.bands[0], &fermi.bands[1]);
      if (job->spin_polarisation == 1)
      sprintf(fermi_bands,"  %d - %d",fermi.bands[2],fermi.bands[3]);
    }

      if (!strcmp(jobname1, "BSE_RANGE")) {
      read_line(file.job, title, 99);
      sscanf(title, "%d %d", &bse_state[0], &bse_state[1]);
    }

      if (!strcmp(jobname1, "NUMPOINTS")) {
      read_line(file.job, title, 99);
      sscanf(title, "%d", &num_points);
      if (num_points > 10) {
      if (job->taskid == 0) fprintf(file.out,"ERROR: NUMBER OF POINTS %3d MUST NOT EXCEED 10\n",num_points);
      MPI_Finalize();
      exit(0);
     }
      for (i = 0; i < num_points; i++) {
        read_line(file.job, title, 99);
        sscanf(title, "%lf %lf %lf", &points1[i].comp1, &points1[i].comp2, &points1[i].comp3);
        points1[i].comp1 /= bohr_to_AA;
        points1[i].comp2 /= bohr_to_AA;
        points1[i].comp3 /= bohr_to_AA;
      }
    }

    } while (strcmp(jobname1, "END"));

      // Allocate fermi array and calculate nbands, etc.

      fermi.nkunique = 1;
      allocate_fermi(&fermi,atoms,job,file);
      fermi.occupied[0] = fermi.homo[0] - fermi.bands[0] + 1;
      nbands = fermi.bands[1] - fermi.bands[0] + 1;
      nocc = fermi.occupied[0];
      nvir = nbands - fermi.occupied[0];
      ntransitions = nocc * nvir;

      // Check input

      // Arrays for string output

      if (!strcmp(jobname, "BSE_PLOT")) sprintf(plot_type,"%8s","BSE WFN");
      if (!strcmp(jobname, "ELECTRON_HOLE_PLOT")) sprintf(plot_type,"%8s","SCF WFN");

      if (job->spin_polarisation == 0) {
      sprintf(fermi_homo,"      %d",fermi.homo[0]);
      sprintf(fermi_bands,"  %d - %d",fermi.bands[0],fermi.bands[1]);
     }

      if (job->spin_polarisation == 1) {
      sprintf(fermi_homo,"  %d - %d",fermi.homo[0],fermi.homo[1]);
      sprintf(fermi_bands,"%d - %d %d - %d",fermi.bands[0],fermi.bands[1],fermi.bands[2],fermi.bands[3]);
     }

      sprintf(monkhorst_pack,"%d %d %d",fermi.is[0], fermi.is[1], fermi.is[2]);

      // Convert input to atomic units

    if (job->taskid == 0) {
    fprintf(file.out,"===========================================================================================================\n");
    if (!strcmp(jobname, "BSE_PLOT"))
    fprintf(file.out,"|                                CORRELATED ELECTRON-HOLE ISOSURFACE PLOT                                 |\n");
    if (!strcmp(jobname, "ELECTRON_HOLE_PLOT"))
    fprintf(file.out,"|               ELECTRON AND HOLE WAVE FUNCTION AND TRANSITION CHARGE DENSITY ISOSURFACE PLOT             |\n");
    fprintf(file.out,"===========================================================================================================\n");

    //fprintf(file.out,"| %-99s     |\n", title1);
    //fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    fprintf(file.out,"| HOMO LEVEL     %8s   | MO RANGE %9s      |                                                   |\n", 
    fermi_homo, fermi_bands);
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    fprintf(file.out,"| TYPE OF PLOT   %8s   | GRID TYPE %8s      | GRID PARAM     %2d %2d %2d |                         |\n", 
    plot_type, title2, grid_par[0], grid_par[1], grid_par[2]);
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    //fprintf(file.out,"| SHRINKING FAC %2d %2d %2d  | NUMBER OF K POINTS %3d  | NUMBER OF PLOTS %2d      | BAND RANGE %4d TO %4d |\n",\
    is_plot[0], is_plot[1], is_plot[2], nkpoints, nplots, vectors[0], vectors[1]);
    //fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    fprintf(file.out,"| GRID ORIGIN                        %10.4lf %10.4lf %10.4lf                                     |\n", 
    points[0].comp1 *bohr_to_AA, points[0].comp2 * bohr_to_AA, points[0].comp3 * bohr_to_AA);
    for (i = 1; i < 4; i++) 
    fprintf(file.out,"| GRID VECTOR %3d                    %10.4lf %10.4lf %10.4lf                                     |\n", 
    i , points[i].comp1 * bohr_to_AA, points[i].comp2 * bohr_to_AA, points[i].comp3 * bohr_to_AA);
    //fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    //fprintf(file.out,"|                         K POINTS IN SHRINKING FACTOR AND ANGS^-1 UNITS                                  |\n");
    //fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    //for (i = 0; i < nkpoints; i++) {
    //j = kplot[i].comp1 * is_plot[1] * is_plot[2] + kplot[i].comp2 * is_plot[2] + kplot[i].comp3;
    //fprintf(file.out,"|%3d            %3d %3d %3d          %10.4lf %10.4lf %10.4lf                                     |\n", \
    i + 1, kplot1[i].comp1, kplot1[i].comp2, kplot1[i].comp3, knet.cart[j].comp1 * bohr_to_AA, knet.cart[j].comp2 * bohr_to_AA,\
    knet.cart[j].comp3 * bohr_to_AA);
   //}
    fprintf(file.out,"===========================================================================================================\n");
   }

    if (!strcmp(jobname, "BSE_PLOT")) {
      if (crystal->type[0] == 'M')
      plot_correlated_electron_hole(grid_par,points,num_points,points1,bse_state,&fermi,atoms,atom_p,shells,gaussians,crystal,\
      symmetry,R,R_tables,job,file);
      else if (crystal->type[0] == 'C' || crystal->type[0] == 'S' || crystal->type[0] == 'P') {
      ;
      //plot_correlated_electron_hole_crystal(grid_sum,bse_state,grid_par,points,num_points,points1,&fermi,atoms,atom_p,shells,\
      gaussians,crystal,R,R_tables,G,symmetry,job,file);
     }
    }

    else if (!strcmp(jobname, "ELECTRON_HOLE_PLOT")) {
      if (crystal->type[0] == 'M' && job->taskid == 0)
      plot_electron_hole_molecule(grid_par,points,&fermi,atoms,atom_p,shells,gaussians,crystal,symmetry,R,R_tables,\
      job,file);
      else if (crystal->type[0] == 'C' || crystal->type[0] == 'S' || crystal->type[0] == 'P') {
      ;
      //plot_electron_hole_crystal(grid_sum,bse_state,grid_par,points,num_points,points1,&fermi,atoms,atom_p,shells,gaussians,\
      crystal,R,R_tables,G,symmetry,job,file);
     }
    }

  free_fermi(&fermi,job);

  }

    /*
    / *****JOB: WAVEFUNCTION_PLOT ***************************************************

    else if (!strcmp(jobname, "WAVEFUNCTION_PLOT")) {

      int grid_type;
      int grid_par[3];
      int atom_label[4];
      int plot_dim = 0;
      int num_sym = symmetry->number_of_operators;
      KPOINT_TRAN knet;
      VECTOR_DOUBLE points[4];

      job->guess_type = 1;
      job->vectors    = 2;
      job->values     = 0;
      job->density    = 1;
      job->max_cycle  = 1;
      job->kpoints    = 1;

      //
      // Insert number of points along each vector.
      //  For 2D plot set grid_par[2]=1,
      //      1D plot     grid_par[1]=1 AND grid_par[2]=1
      //

      read_line(file.job, title, 99);
      sscanf(title, "%s", file.directory1);
      //sscanf(title, "%s", file.scf_eigvec);


      read_line(file.job, title, 99);
      sscanf(title, "%d %d %d", &grid_par[0], &grid_par[1], &grid_par[2]);

      for (i = 0; i < 3; i++) {
        if (grid_par[i] > 1)
          plot_dim++;
      }

      read_line(file.job, title, 99); // type of grid
      sscanf(title, "%8s", title2);

      if (!strcmp(title2, "GRIDVEC")) {

        grid_type = 0;

        //
        // Insert the points which define the plotting region
        // points[0] == ORIGIN
        // axes may be non-orthogonal
        // For the plane (grid_par[2]=1) the points[4] will be ignored
        // 1D     points[3] will be ignored
        //

        for (i = 0; i <= plot_dim; i++) {
          read_line(file.job, title, 99);
          sscanf(title, "%lf %lf %lf", &points[i].comp1, &points[i].comp2, &points[i].comp3);
          points[i].comp1 /= bohr_to_AA;
          points[i].comp2 /= bohr_to_AA;
          points[i].comp3 /= bohr_to_AA;
        }
      } // end of GRIDVEC input

      else if (!strcmp(title2, "GRIDATOM")) {

        grid_type = 1;

        //
        // Insert the labelling numbers of atoms which define the plotting region.
        // atom_label[0] == ORIGIN, atom_label[1], atom_label[2], atom_label[3] define x, y, z axes
        // axes may be non-orthogonal
        //

        for (i = 0; i <= plot_dim; i++) {
          read_line(file.job, title, 99);
          sscanf(title, "%d", &atom_label[i]);
          points[i].comp1 = atoms->cell_vector[i].comp1;
          points[i].comp2 = atoms->cell_vector[i].comp2;
          points[i].comp3 = atoms->cell_vector[i].comp3;
          fprintf(file.out, "GRID VERTEX AT ATOM %d LOCATED AT %lf %lf %lf\n", atom_label[i], points[i].comp1,
          points[i].comp2, points[i].comp3);
        }

      } // end of GRIDATOM input

      int nplots, is_plot[3];
      int vectors[2], nkpoints,ksize;

      read_line(file.job, title, 99); // number of plots to perform
      sscanf(title, "%d", &nplots);

      read_line(file.job, title, 99); // type of plot
      sscanf(title, "%8s", title3);

      read_line(file.job, title, 99); // title of plot
      sscanf(title, "%99s", title1);

      read_line(file.job, title, 99); // shrinking factor for plots, number of k points, first and last vector, spin
      sscanf(title, "%d %d %d", &is_plot[0], &is_plot[1], &is_plot[2]);

      read_line(file.job, title, 99); // shrinking factor for plots, number of k points, first and last vector, spin
      sscanf(title, "%d %d %d",&nkpoints, &vectors[0], &vectors[1]);

      if (nkpoints > 10) {
      if (job->taskid == 0)
      fprintf(file.out, "NUMBER OF K POINTS %3d MUST BE 10 OR FEWER\n",nkpoints);
      MPI_Finalize();
      exit(1);
     }

      count_k_points(&knet,is_plot,crystal,symmetry,job,file);
      allocate_k_points(&knet,crystal,job,file);
      if (job->C09 == 1)
      read_XCBD_crystal_09(&knet,crystal,job,file);
      generate_k_points(&knet,is_plot,crystal,symmetry,job,file);

      VECTOR_INT kplot[nkpoints], kplot1[nkpoints];

      for (i = 0; i < nkpoints; i++) {
        read_line(file.job, title, 99); // k point coordinates to plot
        sscanf(title, "%d %d %d", &kplot1[i].comp1, &kplot1[i].comp2, &kplot1[i].comp3);
        kplot[i].comp1 = kplot1[i].comp1;
        kplot[i].comp2 = kplot1[i].comp2;
        kplot[i].comp3 = kplot1[i].comp3;
        if (kplot[i].comp1 < 0 || kplot[i].comp2 < 0 || kplot[i].comp3 < 0) {
        if (job->taskid == 0) 
        fprintf(file.out,"K POINTS FOR PLOTTING MUST HAVE POSITIVE COMPONENT VALUES\n");
        MPI_Finalize();
        exit(1);
       }
        if (kplot[i].comp1 / is_plot[0] > 0)
        kplot[i].comp1 -= (kplot[i].comp1 / is_plot[0]) * is_plot[0];
        if (kplot[i].comp2 / is_plot[1] > 0)
        kplot[i].comp2 -= (kplot[i].comp2 / is_plot[1]) * is_plot[1];
        if (kplot[i].comp3 / is_plot[2] > 0)
        kplot[i].comp3 -= (kplot[i].comp3 / is_plot[2]) * is_plot[2];
      }

    if (job->taskid == 0) {
    fprintf(file.out,"===========================================================================================================\n");
    fprintf(file.out,"|                                  WAVE FUNCTION MODULUS ISOSURFACE PLOT                                  |\n");
    fprintf(file.out,"===========================================================================================================\n");
    fprintf(file.out,"| %-99s     |\n", title1);
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    fprintf(file.out,"| TYPE OF PLOT %8s   | GRID TYPE %8s        | GRID PARAM %3d %3d %3d  |                         |\n", 
    title3, title2, grid_par[0], grid_par[1], grid_par[2]);
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    fprintf(file.out,"| SHRINKING FAC %2d %2d %2d  | NUMBER OF K POINTS %3d  | NUMBER OF PLOTS %2d      | BAND RANGE %4d TO %4d |\n", 
    is_plot[0], is_plot[1], is_plot[2], nkpoints, nplots, vectors[0], vectors[1]);
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    fprintf(file.out,"| GRID ORIGIN                        %10.4lf %10.4lf %10.4lf                                     |\n", 
    points[0].comp1 *bohr_to_AA, points[0].comp2 * bohr_to_AA, points[0].comp3 * bohr_to_AA);
    for (i = 1; i < 4; i++) 
    fprintf(file.out,"| GRID VECTOR %3d                    %10.4lf %10.4lf %10.4lf                                     |\n", 
    i , points[i].comp1 * bohr_to_AA, points[i].comp2 * bohr_to_AA, points[i].comp3 * bohr_to_AA);
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    fprintf(file.out,"|                         K POINTS IN SHRINKING FACTOR AND ANGS^-1 UNITS                                  |\n");
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    for (i = 0; i < nkpoints; i++) {
    j = kplot[i].comp1 * is_plot[1] * is_plot[2] + kplot[i].comp2 * is_plot[2] + kplot[i].comp3;
    fprintf(file.out,"|%3d            %3d %3d %3d          %10.4lf %10.4lf %10.4lf                                     |\n", \
    i + 1, kplot1[i].comp1, kplot1[i].comp2, kplot1[i].comp3, knet.cart[j].comp1 * bohr_to_AA, knet.cart[j].comp2 * bohr_to_AA, 
    knet.cart[j].comp3 * bohr_to_AA);
   }
    fprintf(file.out,"===========================================================================================================\n");
   }

      switch (title3[0]) {

        case 'W': // wavefunction plot

        //  if (plot_dim < 3)
            ////eigenvec_plot(is_plot, vectors, nkpoints, &knet, eigvec1, kplot, grid_par, points, atoms, shells, \
            gaussians, crystal, symmetry, job, file);
        //  if (plot_dim == 3) {
            eigenvec_isosurface(is_plot, vectors, nkpoints, &knet, kplot, grid_par, points, atoms, atom_p, shells, \
            gaussians, crystal, R, symmetry, job, file);

       //    }
          break;

      } // close switch (title2

        free_k_points(&knet,job);

    }

    // *****JOB: DENSITY_PLOT ***************************************************

    else if (!strcmp(jobname, "DENSITY_PLOT")) {

      int grid_type;
      int grid_par[3];
      int atom_label[4];
      int is_plot[3];
      int vectors[2 * job->spin_dim], ksize;
      double grid_sum[job->spin_dim];
      KPOINT_TRAN knet;
      VECTOR_DOUBLE points[4];

      job->guess_type = 1;
      job->vectors    = 2;
      job->values     = 0;
      job->density    = 1;
      job->max_cycle  = 1;
      job->kpoints    = 1;

      for (i = 0; i < job->spin_dim; i++) 
      grid_sum[i] = k_zero;

      //
      // Insert number of points along each vector.
      //  For 2D plot set grid_par[2]=1,
      //      1D plot     grid_par[1]=1 AND grid_par[2]=1
      //

      read_line(file.job, title, 99); // title of plot
      sscanf(title, "%99s", title1);

      read_line(file.job, title, 99);
      sscanf(title, "%s", file.scf_eigvec);

      read_line(file.job, title, 99);
      sscanf(title, "%d %d %d", &grid_par[0], &grid_par[1], &grid_par[2]);

      read_line(file.job, title, 99); // type of grid
      sscanf(title, "%8s", title2);

      if (!strcmp(title2, "GRIDVEC")) {

        grid_type = 0;

        //
        // Insert the points which define the plotting region
        // points[0] == ORIGIN
        // axes may be non-orthogonal
        // For the plane (grid_par[2]=1) the points[4] will be ignored
        // 1D     points[3] will be ignored
        //

        for (i = 0; i < 4; i++) {
          read_line(file.job, title, 99);
          sscanf(title, "%lf %lf %lf", &points[i].comp1, &points[i].comp2, &points[i].comp3);
          points[i].comp1 /= bohr_to_AA;
          points[i].comp2 /= bohr_to_AA;
          points[i].comp3 /= bohr_to_AA;
        }
      } // end of GRIDVEC input

      else if (!strcmp(title2, "GRIDATOM")) {

        grid_type = 1;

        //
        // Insert the labelling numbers of atoms which define the plotting region.
        // atom_label[0] == ORIGIN, atom_label[1], atom_label[2], atom_label[3] define x, y, z axes
        // axes may be non-orthogonal
        //

        for (i = 0; i < 4; i++) {
          read_line(file.job, title, 99);
          sscanf(title, "%d", &atom_label[i]);
          points[i].comp1 = atoms->cell_vector[i].comp1;
          points[i].comp2 = atoms->cell_vector[i].comp2;
          points[i].comp3 = atoms->cell_vector[i].comp3;
          fprintf(file.out, "GRID VERTEX AT ATOM %d LOCATED AT %lf %lf %lf\n", atom_label[i], points[i].comp1,
          points[i].comp2, points[i].comp3);
        }

      } // end of GRIDATOM input

      read_line(file.job, title, 99); // shrinking factor for plots
      if (job->spin_polarisation == 0)
      sscanf(title, "%d %d %d %d %d", &is_plot[0], &is_plot[1], &is_plot[2], &vectors[0], &vectors[1]);
      else if (job->spin_polarisation == 1)
      sscanf(title, "%d %d %d %d %d %d %d", &is_plot[0], &is_plot[1], &is_plot[2], &vectors[0], &vectors[1], &vectors[2], &vectors[3]);

      count_k_points(&knet,is_plot,crystal,symmetry,job,file);
      allocate_k_points(&knet,crystal,job,file);
      if (job->C09 == 1)
      read_XCBD_crystal_09(&knet,crystal,job,file);
      generate_k_points(&knet,is_plot,crystal,symmetry,job,file);

    if (job->taskid == 0) {
    fprintf(file.out,"===========================================================================================================\n");
    fprintf(file.out,"|                                     CHARGE DENSITY ISOSURFACE PLOT                                      |\n");
    fprintf(file.out,"===========================================================================================================\n");
    fprintf(file.out,"| %-99s     |\n", title1);
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    fprintf(file.out,"| TYPE OF PLOT            | GRID TYPE %8s        | GRID PARAM %3d %3d %3d  |                         |\n", 
    title2, grid_par[0], grid_par[1], grid_par[2]);
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    fprintf(file.out,"| SHRINKING FAC %2d %2d %2d  | NUMBER OF K POINTS        | NUMBER OF PLOTS         | BAND RANGE %4d TO %4d |\n", 
    is_plot[0], is_plot[1], is_plot[2], vectors[0], vectors[1]);
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    fprintf(file.out,"| GRID ORIGIN                        %10.4lf %10.4lf %10.4lf                                     |\n", 
    points[0].comp1 *bohr_to_AA, points[0].comp2 * bohr_to_AA, points[0].comp3 * bohr_to_AA);
    for (i = 1; i < 4; i++) 
    fprintf(file.out,"| GRID VECTOR %3d                    %10.4lf %10.4lf %10.4lf                                     |\n", 
    i , points[i].comp1 * bohr_to_AA, points[i].comp2 * bohr_to_AA, points[i].comp3 * bohr_to_AA);
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
   }

     density_isosurface1(grid_sum, is_plot, vectors, grid_par, points, &knet, atoms, atom_p, shells, gaussians, crystal, R, R_tables, symmetry, job, file);
      //density_isosurface(grid_sum, is_plot, vectors, grid_par, points, &knet, atoms, atom_p, shells, gaussians, crystal, R, R_tables, symmetry, job, file);

    if (job->taskid == 0) {
    if (job->spin_dim == 1) 
    fprintf(file.out,"| INTEGRATED DENSITY %8.2lf  |         |          |  |\n", 
    grid_sum[0]);
    else if (job->spin_dim == 2) 
    fprintf(file.out,"| INTEGRATED DENSITY %8.2lf %8.2lf |         |          |  |\n", 
    grid_sum[0], grid_sum[1]);
    fprintf(file.out,"===========================================================================================================\n");
   }

        free_k_points(&knet,job);

    }
    */

  } while (strcmp(jobname, "END"));

  if (job->taskid == 0) 
  fprintf(file.out, " Function startjob run successfully \n");

  return 0;

}

void fermi_default(FERMI *fermi, CRYSTAL *crystal, ATOM *atoms, JOB_PARAM *job, FILES file)

{
      
int i;
char fermi_homo[9], fermi_bands[9], monkhorst_pack[9];

        switch (crystal->type[0]) {

    case 'C':
      fermi->is[0] = 4;
      fermi->is[1] = 4;
      fermi->is[2] = 4;
      break;

    case 'S':
      fermi->is[0] = 4;
      fermi->is[1] = 4;
      fermi->is[2] = 1;
      break;

    case 'P':
      fermi->is[0] = 1;
      fermi->is[1] = 1;
      fermi->is[2] = 4;
      break;

    case 'M':
      fermi->is[0] = 1;
      fermi->is[1] = 1;
      fermi->is[2] = 1;
      break;

    } // close switch
     
      int electron_count = 0;
      for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++)   // Default is uncharged, unpolarised system
      electron_count += atoms->atomic_number[i];
      fermi->homo[0] = electron_count / 2;
      fermi->homo[1] = electron_count / 2;

      if (fermi->homo[0] > 4) {
      fermi->bands[0] = fermi->homo[0] - 4;                    // Default is 4 active occupied and 4 active unoccupied orbitals
      fermi->bands[1] = fermi->homo[0] + 4;
      fermi->bands[2] = fermi->homo[1] - 4;                    // Default is 4 active occupied and 4 active unoccupied orbitals
      fermi->bands[3] = fermi->homo[1] + 4;
      if (fermi->bands[1] > atoms->number_of_sh_bfns_in_unit_cell) {
      fermi->bands[1] = atoms->number_of_sh_bfns_in_unit_cell;
      fermi->bands[3] = atoms->number_of_sh_bfns_in_unit_cell; }
     }
      else if (fermi->homo[0] <= 4) {
      fermi->bands[0] = 0;           
      fermi->bands[1] = fermi->homo[0] + 1;           
      fermi->bands[2] = 0;           
      fermi->bands[3] = fermi->homo[1] + 1;           
     }
      else {
      if (job->taskid == 0)
      fprintf(file.out,"Error in default values for fermi->bands\n");
      MPI_Finalize();
      exit(1);
    }

}

