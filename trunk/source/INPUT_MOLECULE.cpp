
/*
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <fstream>
#include <xc.h>
#include "ALLOCATE_MEMORY.h"
#include "MATRIX_UTIL.h"
#include "SCF.h"
#include "SCF_CRYSTAL.h"
#include "BETHE_SALPETER.h"
#include "BETHE_SALPETER1.h"
#include "GW.h"
#include "GW_CRYSTAL.h"
#include "TDHF_CRYSTAL.h"
#include "GW_BSE_SCALAPACK.h"
#include "GW_3D.h"
#include "ANALYSIS.h"
#include "PLOTTING.h"
#include "ROTATION_OPERATORS.h"
#include "KPOINTS.h"
#include "CRYSTAL09.h"
#include "INTEGRALS_TEST.h"
*/
#include <mpi.h>
#include <cstring>
#include "mycomplex.h"
#include "conversion_factors.h"
#include "myconstants.h"
#include "mylogical.h"
#include "USER_DATA.h"
#include "TOOLS.h"
#include "LIMITS.h"
#include "ALLOCATE_MEMORY_MOLECULE.h"
#include "SCF_MOLECULE.h"
#include "OPTICAL_SPECTRUM_MOLECULE.h"
#include "GW_BSE_MOLECULE.h"
#include "ERRORS.h"
#include "INPUT_MOLECULE.h"

using namespace std;

int startjob(ATOM *atoms, ATOM *atoms_ax, ATOM_TRAN *atom_p, SHELL *shells, GAUSSIAN *gaussians, SHELL *shells_ax, GAUSSIAN *gaussians_ax, CRYSTAL *crystal, SYMMETRY *symmetry, REAL_LATTICE *R, REAL_LATTICE_TABLES *R_tables, RECIPROCAL_LATTICE *G, JOB_PARAM *job, FILES file)

{

  int i, j, k;

  // set default values for job parameters

  job->itol1         =  7.0; // real space cutoff for product of two Gaussians in Bohr
  job->itol2         =  7.0; // real space cutoff for product of two Gaussians in Bohr
  job->itol3         =  7.0; // real space cutoff for product of two Gaussians in Bohr
  job->itol4         = 16.0; // real space cutoff for product of two Gaussians in Bohr
  job->itol5         = 22.0; // real space cutoff for product of two Gaussians in Bohr

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

  job->kpoints = 0;
  job->values  = 0;

  job->fix_occ = 0;             //!< Occupancy for density matrix calculation determined by Fermi level

  job->linewidth = 0.01 / au_to_eV; // default linewidth for plotting 0.01 eV

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
      job->gw_spk = 0;          // default is not to use SCALAPACK for GW
      job->gw_int = 0;          // default is integrals need to be calculated and written to disk
      job->gw = 0;              // default is not to calculate GW correction for BSE 

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
  char *double_range;
  char *int_range;
  double double_value[2];
  int  int_value[2];

  do {

    read_line(file.job, title, 99);
    sscanf(title, "%s", jobname);
    //printf("%s\n",jobname); fflush(stdout);

    // *****JOB: SPIN POLARISATION ***************************************************

    if (!strcmp(jobname, "SPIN_POLARISED") || !strcmp(jobname, "UHF")) {

      read_line(file.job, title, 99);
      sscanf(title, "%d %d", &job->spin_polarisation,&job->spin_orb);
      printf("%d %d \n",job->spin_orb,job->spin_polarisation); fflush(stdout);

       switch (job->spin_polarisation) {
         case TRUE:
           job->spin_dim = 2;
           job->spin_fac = 1;
           job->spin_pol = 1;
            if (job->taskid == 0) 
             fprintf(file.out, "SPIN POLARISED CALCULATION\n\n");
            if (job->taskid == 0 && job->spin_orb == 1) 
             fprintf(file.out, "SPIN ORBIT COUPLING ON\n\n");
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
      char hamiltonian[7], spin_pol[4], x_functional[15], c_functional[12], dft_grid[9], integrals[10], direct[4];
      char guess[13], intexist[13], monkhorst_pack[9];
      double w[2];
      FERMI fermi;

      // Set Job Defaults

      //sprintf(file.scf_eigvec,"."); // Default for MPI file IO is current directory
      read_line(file.job, title, 99);
      sscanf(title, "%s", file.scf_eigvec);

      job->mpi_io = 0;          // Default is to use local disk
      job->mpp = 0;             // Default is to use cores on maximum 1 node 

      job->type = 0;

      job->vectors = 2;
      job->values  = 2;
      job->density = 2;
      job->kpoints = 0;

      job->xc_hfx  = 1;          // Default is Hartree-Fock
      job->xc_num  = 0;          // Default is no DFT functionals
      job->xc_typ[0] = -1;          
      job->xc_typ[1] = -1;          
      sprintf(hamiltonian,"%s","HF");
      sprintf(x_functional,"%s","HF");
      sprintf(c_functional,"%s","None");
      sprintf(dft_grid,"%s","None");

      job->spin_dim = 1;         // Default is no spin polarisation
      job->spin_fac = 2;
      job->spin_pol = 0;
      job->spin_orb = 0;
      sprintf(spin_pol,"%s","Off");
      for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++) {
         atoms->magnetic[i] = 1;
         atoms->spin[i]     = 0;
        }

      job->guess_type = 0;       // Default is new calculation from atomic densities
      sprintf(guess,"%s","Atom Density");

      job->max_cycle  = 50;      // Default is 50 cycles
      job->scf_tol    = 1.0e-06; // Default is 1e-06

      if (atoms->number_of_sh_bfns_in_unit_cell <= 3000) {
      job->scf_direct = 0;  
      sprintf(direct,"%s","Off");
     }
      else if (atoms->number_of_sh_bfns_in_unit_cell >  3000) {
      job->scf_direct = 1;       // Default is On for more than 3000 basis functions per unit cell
      sprintf(direct,"%s","On");
     }

      job->coul_int = 0;
      sprintf(integrals,"%s","Exact");

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
      sprintf(monkhorst_pack,"%s","Off");
      //sprintf(monkhorst_pack,"%d %d %d",fermi.is[0], fermi.is[1], fermi.is[2]);
      break;

    } // close switch

      // read input file

   do {

      read_line(file.job, title, 99);
      sscanf(title, "%s", jobname1);

        if (!strcmp(jobname1, "ITOL")) {
          read_line(file.job, title, 99);
          sscanf(title, "%lf %lf %lf %lf %lf", &job->itol1, &job->itol2, &job->itol3, &job->itol4, &job->itol5);
         }

        if (!strcmp(jobname1, "MPI_IO")) {
          job->mpi_io = 1; 
         }

        if (!strcmp(jobname1, "MPP")) {
          job->mpp = 1; 
         }

        if (!strcmp(jobname1, "DIIS_ON")) {
          job->diis = 1; 
         }

        if (!strcmp(jobname1, "HF")) {
          job->xc_num = 0; 
          job->xc_hfx    =  1;          
          sprintf(hamiltonian,"%s","HFT");
          sprintf(x_functional,"%s","HF");
          sprintf(c_functional,"%s","None");
         }

        else if (!strcmp(jobname1, "DFT")) {

          job->xc_grd = 0;           // Default DFT grid is STANDARD
          job->xc_lmx = 13;
          job->xc_rad = 64;
          sprintf(dft_grid,"%s","Standard");
          sprintf(hamiltonian,"%s","DFT");

            read_line(file.job, title, 99);
            sscanf(title, "%d", &job->xc_num);
            if (job->xc_num < 0 || job->xc_num > 2 && job->taskid == 0) {
              fprintf(file.out,"ERROR: DFT Hamiltonian requires number of XC potentials\n");
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

	      /*
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
                 xc_func_type func;
                 if(xc_func_init(&func, job->xc_typ[i], 1) != 0){
                 if (job->taskid == 0) {
                   fprintf(file.out, "Functional %s not found\n",exch_type);
                   fprintf(file.out,"Choices are HARTREE_FOCK or at www.tddft.org/programs/octopus/wiki/index.php/Libxc:manual\n");
                  }
                   MPI_Finalize();
                   exit(0);
                  }
                 xc_func_end(&func);
                }
	      */
               } // close loop on i
	/*

        else if (!strcmp(jobname1, "HYBRID")) {
          job->xc_grd = 0;           // Default DFT grid is STANDARD
          job->xc_lmx = 13;
          job->xc_rad = 64;
          job->xc_hfx = 1;
          sprintf(dft_grid,"%s","Standard");
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
                 xc_func_type func;
                 if(xc_func_init(&func, job->xc_typ[i], 1) != 0){
                 if (job->taskid == 0) {
                   fprintf(file.out, "Functional %s not found\n",exch_type);
                   fprintf(file.out,"Choices are HARTREE_FOCK or at www.tddft.org/programs/octopus/wiki/index.php/Libxc:manual\n");
                  }
                   MPI_Finalize();
                   exit(0);
                  }
                 xc_func_end(&func);
                }
               } // close loop on i
              }
  */

            if (!strcmp(jobname1, "DFT_GRID")) {
              read_line(file.job, title, 99);
              sscanf(title, "%s", jobname1);
              if (!strcmp(jobname1, "STANDARD")) {
              job->xc_grd = 0;
              job->xc_lmx = 13;
              job->xc_rad = 64;
              sprintf(dft_grid,"%s","Standard");
             } 
              else if (!strcmp(jobname1, "LARGE")) {
              job->xc_grd = 1;
              job->xc_lmx = 13;
              job->xc_rad = 96;
              sprintf(dft_grid,"%s","Large");
             } 
              else if (!strcmp(jobname1, "XLARGE")) {
              job->xc_grd = 2;
              job->xc_lmx = 16;
              job->xc_rad = 96;
              sprintf(dft_grid,"%s","XLarge");
             } 
              else if (!strcmp(jobname1, "XXLARGE")) {
              job->xc_grd = 3;
              job->xc_lmx = 19;
              job->xc_rad = 96;
              sprintf(dft_grid,"%s","XXLarge");
             } 
              else {
              if (job->taskid == 0) 
              fprintf(file.out,"ERROR: DFT_GRID must be one of STANDARD, LARGE, XLARGE, XXLARGE\n");
              MPI_Finalize();
              exit(0);
             } 
            } 

   } // close DFT input

        if (!strcmp(jobname1, "SPIN_POLARISED") || !strcmp(jobname1, "UHF")) {
           read_line(file.job, title, 99);
           sscanf(title, "%d", &job->spin_tot);
           job->spin_dim = 2;
           job->spin_fac = 1;
           job->spin_pol = 1;
           job->spin_polarisation = 1;
           sprintf(spin_pol,"%s","On");
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

      if (!strcmp(jobname1, "GUESS")) {
      read_line(file.job, title, 99);
      sscanf(title, "%s", guess);
      if (guess[0] == 'A')
      job->guess_type = 0;
      else if (guess[0] == 'R')
      job->guess_type = 1;
      if (guess[0] != 'A' && guess[0] != 'R') {
      fprintf(file.out,"ERROR: GUESS %s must be Atoms (new calculation) or Restart\n",guess);
      MPI_Finalize();
      exit(1);
     }
      read_line(file.job, title, 99);
      sscanf(title, "%s", intexist);
      if (intexist[0] == 'I')
      job->int_exist = 1;
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

      if (!strcmp(jobname1, "DIRECT_SCF_ON")) {
      job->scf_direct = 1;
      sprintf(direct,"On");
     }

      if (!strcmp(jobname1, "DIRECT_SCF_OFF")) {
      job->scf_direct = 0;
      sprintf(direct,"Off");
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

      if (!strcmp(jobname1, "DENSITY_FITTING_ON")) {
      job->scf_denfit = 1;
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
 //if (crystal->type[0] == 'C' || crystal->type[0] == 'S' || crystal->type[0] == 'P') { 
 fprintf(file.out,\
 "| COULOMB ITOL    %4.1f %4.1f | EXC ITOL %4.1f %4.1f %4.1f | MXR %2d R cutoff %4.1f A  | MXG %2d G cutoff %4.1f 1/A|\n", \
 job->itol1,job->itol2,job->itol3,job->itol4,job->itol5,job->mxr,R->cutoff * bohr_to_AA,job->mxg,G->cutoff / bohr_to_AA);
 fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
 fprintf(file.out,\
 "| RMAX %6d RMARG  %6d | EWALD %5d LAST %5d  | G MAX %4d G LAST %4d  | SGS %s PMS %s TRS %s KSS %s |\n", \
    R->max_vector,R->margin_vector,R->last_vector,R->last_ewald_vector,G->max_vector,G->last_vector,sgs,pms,trs,kss); //}
 //else if (crystal->type[0] == 'M') fprintf(file.out,\
 "| COULOMB ITOL    %4.1f %4.1f | EXC ITOL %4.1f           |                         |                         |\n", \
 job->itol1,job->itol2,job->itol3);
 fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n\n");

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
          //} while (strcmp(jobname, "END"));

/*

    // *****JOB: PRINT EIGENVECTORS AND/OR EIGENVALUES *****************************

    else if (!strcmp(jobname, "EXCITON_VECTORS_VALUES")) {

      int shr[3];
      KPOINT_TRAN knet;
      VECTOR_INT k1;
      FERMI fermi;

      int_value[0] = 1;
      int_value[1] = atoms->number_of_sh_bfns_in_unit_cell;
      //int_range = "ENERGY BAND";
      char int_range[12] = "ENERGY BAND";

      job->guess_type = 1;
      //job->vectors    = 0;
      job->vectors    = 2;
      job->mpp        = 0;
      job->values     = 2;
      job->density    = 1;
      job->max_cycle  = 1;
      job->kpoints    = 1;
      job->fix_occ    = 1;

      read_line(file.job, title, 99);
      sscanf(title, "%d %d %d %d %d", &shr[0], &shr[1], &shr[2], &fermi.bands[0], &fermi.bands[1]);

      read_line(file.job, title, 99);
      sscanf(title, "%s", file.scf_eigvec);

      knet_check(shr,jobname,crystal,file);
      int_range_check(fermi.bands,int_value,int_range,jobname,file);

        if (job->taskid == 0) {
        fprintf(file.out,"===========================================================================================================\n");
        fprintf(file.out,"|                               EXCITON EIGENVECTOR AND EIGENVALUE CALCULATION                            |\n");
        fprintf(file.out,"===========================================================================================================\n");
        fprintf(file.out,"| MONKHORST-PACK   %2d %2d %2d | BAND RANGE %4d to %4d |                                                   |\n", \
        shr[0], shr[1], shr[2], fermi.bands[0], fermi.bands[1]);
        fprintf(file.out,"===========================================================================================================\n");
       }

      int nbands = fermi.bands[1] - fermi.bands[0] + 1;
      count_k_points(&knet,shr,crystal,symmetry,job,file);
      allocate_k_points(&knet,crystal,job,file);
      generate_k_points(&knet,shr,crystal,symmetry,job,file);
      fermi.npoints = knet.unique;
      fermi.knet_list = (VECTOR_KNET *) malloc(knet.unique * sizeof(VECTOR_KNET));
      for (i = 0; i < knet.unique; i++) {
      k1 = decompose_k_point(shr, knet.ibz[i], crystal, job, file);
      getcart(&k1, &fermi.knet_list[i].cart, shr, crystal);
      //fprintf(file.out,"%3d %3d %3d %10.4lf %10.4lf %10.4lf\n",k1.comp1,k1.comp2,k1.comp3,fermi.knet_list[i].cart.comp1, fermi.knet_list[i].cart.comp2,\
      fermi.knet_list[i].cart.comp3);
     }

      crystal_scf(&fermi,atoms,atoms_ax,atom_p,shells,gaussians,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
      //CHANGES2015crystal_scf(&fermi, atoms, atom_p, shells, gaussians, crystal, symmetry, R, R_tables, G, job, file);

      free_fermi(&fermi,job);
      free_k_points(&knet,job);

    }
    */

    // *****JOB: GW CALCULATION ***************************************************************************************

    else if (!strcmp(jobname, "GW")) {

      FERMI fermi;
      int nspectra, npoints, self[10];
      double energy_range[2];
      char fermi_homo[9], fermi_bands[9], spectrum[14], monkhorst_pack[9];
      double_value[0] = lower_energy_limit;
      double_value[1] = upper_energy_limit;
      //double_range="ENERGY RANGE";
      char double_range[13] = "ENERGY RANGE";
      //int_range = "ENERGY BAND";
      char int_range[12] = "ENERGY BAND";
      int_value[0] = 1;
      int_value[1] = atoms->number_of_sh_bfns_in_unit_cell;

      // Defaults

      fermi_default(&fermi,crystal,atoms,job,file);
      job->self_plot = 0;
      job->bse_ham = 2;
      job->bse_tda = 1;
      
      read_line(file.job, title, 99);
      sscanf(title, "%s", file.scf_eigvec);

   if (crystal->type[0] != 'M') {
      read_line(file.job, title, 99);
      sscanf(title, "%d %d %d", &fermi.is[0], &fermi.is[1], &fermi.is[2]);
      ////knet_check(fermi.is,jobname,crystal,file);
      //fprintf(file.out,"IS is %3d %3d %3d\n",fermi.is[0],fermi.is[1],fermi.is[2]);
     }
   else { fermi.is[0] = 1; fermi.is[1] = 1; fermi.is[2] = 1; }

      //read_line(file.job, title, 99);
      //sscanf(title, "%s", file.scf_eigvec);
 
    do {
 
       read_line(file.job, title, 99);
       sscanf(title, "%s", jobname1);
  
      if (!strcmp(jobname1, "INTEGRALS_EXIST")) {
      job->gw_int = 1;
     }

      //if (!strcmp(jobname1, "IN_CORE")) {
      //job->gw_int = 2;
     //}

      if (!strcmp(jobname1, "SCALAPACK")) {
      job->gw_spk = 1;
      job->bse_ham = 2;
     }

      if (!strcmp(jobname1, "RPA_VECTORS")) {
      read_line(file.job, title, 99);
      sscanf(title, "%d", &job->rpa_lim);
      if (job->taskid == 0) {
     }
    }
       if (!strcmp(jobname1, "MO_RANGE")) {
       read_line(file.job, title, 99);
       sscanf(title, "%d %d", &fermi.bands[0], &fermi.bands[1]);
       ////sprintf(fermi_homo, "       %d",fermi.homo[0]);
       ////sprintf(fermi_bands,"  %d - %d",fermi.bands[0],fermi.bands[1]);
       ////if (job->spin_polarisation == 1)
       ////sprintf(fermi_bands,"  %d - %d",fermi.bands[2],fermi.bands[3]);
     }
   
       if (!strcmp(jobname1, "SPECTRUM")) {
       read_line(file.job, title, 99);
       sscanf(title, "%lf %lf %d", &energy_range[0], &energy_range[1], &npoints);
       read_line(file.job, title, 99);
       sscanf(title, "%d", &nspectra);
       ////sprintf(spectrum,"%4.1lf - %4.1lf",energy_range[0],energy_range[1]);
       job->self_plot = 1;
       if (nspectra > 10) {
       if (job->taskid == 0)
       fprintf(file.out,"NUMBER OF SELF-ENERGY PLOTS %3d MUST NOT EXCEED 10\n", nspectra);
       MPI_Finalize();
       exit(1);
      }
       read_line(file.job, title, 99);
       sscanf(title, "%d %d %d %d %d %d %d %d %d %d", \
       &self[0],&self[1], &self[2], &self[3], &self[4], &self[5], &self[6], &self[7], &self[8], &self[9]);
       for (i = 0; i < nspectra; i++) {
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

      ////knet_check(fermi.is,jobname,crystal,file);
      ////int_range_check(fermi.bands,int_value,int_range,jobname,file);

      job->npoints  = npoints;
      job->nspectra = nspectra;
      job->energy_range[0] = energy_range[0] / au_to_eV;
      job->energy_range[1] = energy_range[1] / au_to_eV;

      ////double_range_check(energy_range,double_value,double_range,jobname,file);

        switch (crystal->type[0]) {

    case 'C':

      //temporary until shifted to BSE input

      energy_range[0] = 0.01;
      energy_range[1] = 85.0;
      job->energy_range[0] = energy_range[0] / au_to_eV;
      job->energy_range[1] = energy_range[1] / au_to_eV;
      job->npoints  = 1000;
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

       //gw_crystal(&fermi,atoms,atom_p,shells,gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
       //GW3D(&fermi,atoms,atom_p,shells,gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);

    break;

    case 'M':

       //if (job->gw_spk == 0)
       //GW(&fermi,atoms,atom_p,shells,gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
       //else if (job->gw_spk == 1)
       gw_molecule(&fermi,atoms,atom_p,shells,gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);

    break;

   } // close switch

       free_fermi(&fermi,job);
    }

    // *****JOB: RPA CALCULATION ***************************************************************************************

    else if (!strcmp(jobname, "RPA")) {

      FERMI fermi;
      char fermi_homo[9], fermi_bands[9], monkhorst_pack[9];
      char int_range[12] = "ENERGY BAND";
      int_value[0] = 1;
      int_value[1] = atoms->number_of_sh_bfns_in_unit_cell;

      // Defaults

      job->bse_int = 2;
      job->bse_ham = 2;
      job->bse_tda = 1;
      fermi_default(&fermi,crystal,atoms,job,file);
      
   if (crystal->type[0] != 'M') {
      read_line(file.job, title, 99);
      sscanf(title, "%d %d %d", &fermi.is[0], &fermi.is[1], &fermi.is[2]);
      //knet_check(fermi.is,jobname,crystal,file);
      fprintf(file.out,"IS is %3d %3d %3d\n",fermi.is[0],fermi.is[1],fermi.is[2]);
     }
   else { fermi.is[0] = 1; fermi.is[1] = 1; fermi.is[2] = 1; }

      read_line(file.job, title, 99);
      sscanf(title, "%s", file.scf_eigvec);
  
    do {
 
       read_line(file.job, title, 99);
       sscanf(title, "%s", jobname1);
       //printf("%s\n",jobname1);
  
      //if (!strcmp(jobname1, "INTEGRALS_EXIST")) {
      //job->gw_int = 1;
     //}

      //if (!strcmp(jobname1, "IN_CORE")) {
      //job->gw_int = 2;
     //}

      //if (!strcmp(jobname1, "SCALAPACK")) {
      //job->gw_spk = 1;
      //job->bse_ham = 2;
     //}

      if (!strcmp(jobname1, "RPA_VECTORS")) {
      read_line(file.job, title, 99);
      sscanf(title, "%d", &job->rpa_lim);
      if (job->taskid == 0) {
     }
    }

       if (!strcmp(jobname1, "MO_RANGE")) {
       read_line(file.job, title, 99);
       sscanf(title, "%d %d", &fermi.bands[0], &fermi.bands[1]);
       printf("%3d %3d\n",fermi.bands[0],fermi.bands[1]);
       ////sprintf(fermi_homo, "       %d",fermi.homo[0]);
       ////sprintf(fermi_bands,"  %d - %d",fermi.bands[0],fermi.bands[1]);
       ////if (job->spin_polarisation == 1)
       ////sprintf(fermi_bands,"  %d - %d",fermi.bands[2],fermi.bands[3]);
     }
   
     } while (strcmp(jobname1, "END"));

      knet_check(fermi.is,jobname,crystal,file);
      int_range_check(fermi.bands,int_value,int_range,jobname,file);

        switch (crystal->type[0]) {

    case 'C':

    break;

    case 'M':

       rpa_molecule1(&fermi,atoms,atom_p,shells,gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);

    break;

   } // close switch

       free_fermi(&fermi,job);

    }

    // *****JOB: BETHE-SALPETER EQUATION HAMILTONIAN DIAGONALISATION AND OPTICAL SPECTRUM *****************************

    else if (!strcmp(jobname, "BETHE_SALPETER") || !strcmp(jobname, "BSE")) {

      int npoints;
      double energy_range[2];
      char optical_type = 'O';
      char hamiltonian[9], monkhorst_pack[9], mpi_io[80], spectrum[14], fermi_homo[9], fermi_bands[17], field[80];
      char tamm_dancoff[9],spin_state[9],linear_algebra[10];
      char scissor_shift[9],scale_factor[9];
      int numfrag = 2, max_atoms = 40;
      int natoms[numfrag], nat1[max_atoms];
      int nat[max_atoms][2];
      FERMI fermi;
      double_value[0] = lower_energy_limit;
      double_value[1] = upper_energy_limit;
      //double_range="ENERGY RANGE";
      char double_range[13] = "ENERGY RANGE";
      //int_range = "ENERGY BAND";
      char int_range[12] = "ENERGY BAND";
      int_value[0] = 1;
      int_value[1] = atoms->number_of_sh_bfns_in_unit_cell;

      // Defaults

      read_line(file.job, title, 99);
      sscanf(title, "%s", file.scf_eigvec);
  
      //read_line(file.job, title, 99);
      //sscanf(title, "%s", file.directory1);
      //printf("%s\n",file.directory1);

      //sprintf(file.directory1,"."); // Default for MPI file IO is current directory

      job->mpi_io = 0;          // Default is to use local disk

      job->type = 1;

      job->guess_type = 1;
      job->vectors    = 2;
      job->values     = 2;

      job->kpoints    = 0;

      job->bse_ham  = 1;     // default is BSE 
      job->bse_tda  = 0;     // default is full BSE matrix 
      job->bse_spin = 0;     // default is spin singlet in BSE and similar calculations
      job->bse_spk  = 0;     // default is not to use SCALAPACK in BSE and TDHF calculations
      job->bse_int  = 0;     // default is integrals need to be calculated and written to disk
      job->gw_spk   = 0;     // default is not to use SCALAPACK in GW calculations
      job->gw_int   = 0;     // default is integrals need to be calculated and written to disk

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
     
      int electron_count = 0;
      for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++)   // Default is uncharged, unpolarised system
      electron_count += atoms->atomic_number[i];
      fermi.homo[0] = electron_count / 2;
//HERE      fermi.homo = electron_count / 2;
      fermi.homo[1] = electron_count / 2;

      //if (fermi.homo[0] > 4) {
//HERE      fermi.bands[0] = fermi.homo[0] - 4;                    // Default is 4 active occupied and 4 active unoccupied orbitals
      //fermi.bands[1] = fermi.homo[0] + 4;
      //fermi.bands[2] = fermi.homo[1] - 4;                    // Default is 4 active occupied and 4 active unoccupied orbitals
      //fermi.bands[3] = fermi.homo[1] + 4;
      if (fermi.homo[0] > 4) {
      fermi.bands[0] = fermi.homo[0] - 4;                    // Default is 4 active occupied and 4 active unoccupied orbitals
      fermi.bands[1] = fermi.homo[0] + 4;
      fermi.bands[2] = fermi.homo[1] - 4;                    // Default is 4 active occupied and 4 active unoccupied orbitals
      fermi.bands[3] = fermi.homo[1] + 4;
      if (fermi.bands[1] > atoms->number_of_sh_bfns_in_unit_cell) {
      fermi.bands[1] = atoms->number_of_sh_bfns_in_unit_cell;
      fermi.bands[3] = atoms->number_of_sh_bfns_in_unit_cell; }
     }
      //else if (fermi.homo <= 4) {
      //fermi.bands[0] = 0;           
      //fermi.bands[1] = fermi.homo + 1;           
      //fermi.bands[2] = 0;           
      //fermi.bands[3] = fermi.homo + 1;           
     //}
      else if (fermi.homo[0] <= 4) {
      fermi.bands[0] = 0;           
      fermi.bands[1] = fermi.homo[0] + 1;           
      fermi.bands[2] = 0;           
      fermi.bands[3] = fermi.homo[1] + 1;           
     }
      else {
      if (job->taskid == 0)
      fprintf(file.out,"Error in default values for fermi.bands in BSE calculation input\n");
      MPI_Finalize();
      exit(1);
    }

        switch (crystal->type[0]) {

    case 'C':
    case 'S':
    case 'P':
      //sprintf(fermi_homo, "       %d",fermi.homo[0]);
      //HERE sprintf(fermi_homo, "       %d",fermi.homo);
      //sprintf(fermi_bands,"  %d - %d",fermi.bands[0],fermi.bands[1]);
      //if (job->spin_polarisation == 1)
      //sprintf(fermi_bands,"  %d - %d",fermi.bands[2],fermi.bands[3]);
      break;
    case 'M':
      //sprintf(fermi_homo, "       %d",fermi.homo);
      //sprintf(fermi_bands,"  %d - %d",fermi.bands[0],fermi.bands[1]);
      //HERE sprintf(fermi_homo, "     %3d",fermi.homo);
      //sprintf(fermi_homo, "     %3d",fermi.homo[0]);
      //sprintf(fermi_bands," %3d-%3d",fermi.bands[0],fermi.bands[1]);
      //if (job->spin_polarisation == 1)
      //sprintf(fermi_bands,"  %d - %d",fermi.bands[2],fermi.bands[3]);
      break;

    } // close switch

      npoints = 400;                                      // Default range for optical spectrum is 0 - 10 eV with 200 points
      energy_range[0] = 0.01;
      energy_range[1] = 20.0;
      sprintf(spectrum,"%4.1lf - %4.1lf",energy_range[0],energy_range[1]);

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
      sprintf(field," %d     %4.1lf %4.1lf %4.1lf    %4.1lf %4.1lf %4.1lf   %4.1lf %4.1lf %4.1lf",job->field_dirs, \
      job->e_field[0].comp1, job->e_field[0].comp2, job->e_field[0].comp3, \
      job->e_field[1].comp1, job->e_field[1].comp2, job->e_field[1].comp3, \
      job->e_field[2].comp1, job->e_field[2].comp2, job->e_field[2].comp3);

      job->scissor   = 0.00 ; // default virtual state shift for periodic TDHF
      job->scalefac  = k_one; // default scaling factor for electron-hole attraction matrix elements in periodic TDHF


   do {

      read_line(file.job, title, 99);
      sscanf(title, "%s", jobname1);

      if (!strcmp(jobname1, "COULOMB_TEST")) {
      job->bse_cou = 1;
     }

      if (!strcmp(jobname1, "EXCHANGE_TEST")) {
      job->bse_exc = 1;
     }

      if (!strcmp(jobname1, "IN_CORE")) {
      job->bse_int = 2;
     }

      if (!strcmp(jobname1, "INTEGRALS_EXIST")) {
      job->bse_int = 1;
     }

      if (!strcmp(jobname1, "SCALAPACK")) {
      job->bse_spk = 1;
     }

      if (!strcmp(jobname1, "SCISSOR")) {
      read_line(file.job, title, 99);
      sscanf(title, "%lf", &job->scissor);
      //printf("sciss %lf %s\n",job->scissor,title);
     }

      if (!strcmp(jobname1, "SCALEFACTOR")) {
      read_line(file.job, title, 99);
      sscanf(title, "%lf", &job->scalefac);
      //printf("scale %lf\n",job->scalefac);
     }

      if (!strcmp(jobname1, "TDA")) {
      job->bse_tda = 1;
      if (job->taskid == 0) {
     }
    }

      if (!strcmp(jobname1, "TRIPLET")) {
      job->bse_spin = 1;
      if (job->taskid == 0) {
     }
    }

      if (!strcmp(jobname1, "TDHF")) {
      job->bse_ham = 0;
      if (job->taskid == 0) {
     }
    }

      if (!strcmp(jobname1, "BSE_VECTORS")) {
      read_line(file.job, title, 99);
      sscanf(title, "%d", &job->bse_lim);
      if (job->taskid == 0) {
     }
    }

      if (!strcmp(jobname1, "MPI_IO")) {
        job->mpi_io = 1; 
       }

      if (!strcmp(jobname1, "MONKHORST_PACK_NET")) {
      read_line(file.job, title, 99);
      sscanf(title, "%d %d %d", &fermi.is[0], &fermi.is[1], &fermi.is[2]);
      sprintf(monkhorst_pack,"%d %d %d",fermi.is[0], fermi.is[1], fermi.is[2]);
      if (crystal->type[0] == 'M' && job->taskid == 0) {
      fprintf(file.out,"Monkhorst-Pack net keyword not relevant for molecules\n");
      MPI_Finalize();
      exit(1);
     }
    }

      if (!strcmp(jobname1, "FERMI_LEVEL")) {
      read_line(file.job, title, 99);
      if (job->spin_polarisation == 0) {
      sscanf(title, "%d", &fermi.homo[0]);
      sprintf(fermi_homo, "     %3d",fermi.homo[0]);
     }
      else if (job->spin_polarisation == 1) {
      sscanf(title, "%d %d", &fermi.homo[0],&fermi.homo[1]);
      sprintf(fermi_homo, "%3d %3d",fermi.homo[0],fermi.homo[1]);
     }
    }

      if (!strcmp(jobname1, "MO_RANGE")) {
      read_line(file.job, title, 99);
      sscanf(title, "%d %d", &fermi.bands[0], &fermi.bands[1]);
      //sprintf(fermi_homo, "     %3d",fermi.homo[0]);
      //HERE sprintf(fermi_homo, "     %3d",fermi.homo);
      //sprintf(fermi_homo, "       %d",fermi.homo);
      //sprintf(fermi_bands,"  %d - %d",fermi.bands[0],fermi.bands[1]);
      if (job->spin_polarisation == 1) {
      read_line(file.job, title, 99);
      sscanf(title, "%d %d", &fermi.bands[2], &fermi.bands[3]);
      //sprintf(fermi_bands,"  %d - %d",fermi.bands[2],fermi.bands[3]);
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

      if (!strcmp(jobname1, "SPECTRUM")) {
      read_line(file.job, title, 99);
      sscanf(title, "%d %lf %lf", &npoints, &energy_range[0], &energy_range[1]);
      sprintf(spectrum,"%4.1lf - %4.1lf",energy_range[0],energy_range[1]);
    }

      if (!strcmp(jobname1, "LINEWIDTH")) {
      read_line(file.job, title, 99);
      sscanf(title, "%lf", &job->linewidth);
      if (job->linewidth <= 0.0049 || job->linewidth >= 1.0001) {
      if (job->taskid == 0)
      fprintf(file.out,"\nLINEWIDTH %9.2e eV OUT OF RANGE: CHOOSE VALUE BETWEEN 0.005 AND 1.0 eV\n",job->linewidth);
      MPI_Finalize();
      exit(1);
     }
      job->linewidth /= au_to_eV;
    }

      if (!strcmp(jobname1, "DENSITY_FITTING_ON")) {
      job->bse_denfit = 1;
      if (job->taskid == 0) {
      fprintf(file.out,"USING DENSITY FITTED INTEGRALS IN BSE CALCULATION\n");
    }
    }

      if (!strcmp(jobname1, "SCREENED_EXCHANGE_ON")) {
      job->bse_denfit = 1;
      job->bse_screxc = 1;
      if (job->taskid == 0) {
      fprintf(file.out,"USING SCREENED EXCHANGE VIA DENSITY FITTED INTEGRALS IN BSE CALCULATION\n");
    }
    }

      if (!strcmp(jobname1, "FIELD")) {
      read_line(file.job, title, 99);
      sscanf(title, "%d", &job->field_dirs);
      if (job->taskid == 0 && job->field_dirs > 3) {
      fprintf(file.out,"Number of field directions %d exceeds maximum value (3)\n",job->field_dirs);
      MPI_Finalize();
      exit(1);
     }
      VECTOR_DOUBLE fld[job->field_dirs];
      for (i = 0; i < job->field_dirs; i++) {
      read_line(file.job, title, 99);
      sscanf(title, "%lf %lf %lf", &fld[i].comp1,  &fld[i].comp2, &fld[i].comp3);
      ////double fld_tmp = sqrt(double_vec_dot(&fld[i],&fld[i]));
      ////fld[i].comp1 /= fld_tmp;
      ////fld[i].comp2 /= fld_tmp;
      ////fld[i].comp3 /= fld_tmp;
      job->e_field[i].comp1 = fld[i].comp1;
      job->e_field[i].comp2 = fld[i].comp2;
      job->e_field[i].comp3 = fld[i].comp3;
      sprintf(field," %d  %4.1lf %4.1lf %4.1lf    %4.1lf %4.1lf %4.1lf      %4.1lf %4.1lf %4.1lf",job->field_dirs, \
      job->e_field[0].comp1, job->e_field[0].comp2, job->e_field[0].comp3, \
      job->e_field[1].comp1, job->e_field[1].comp2, job->e_field[1].comp3, \
      job->e_field[2].comp1, job->e_field[2].comp2, job->e_field[2].comp3);
     }
    }

    } while (strcmp(jobname1, "END"));

      if (job->spin_polarisation == 0) {
      sprintf(fermi_homo, "     %3d",fermi.homo[0]);
      sprintf(fermi_bands,"  %d - %d",fermi.bands[0],fermi.bands[1]);
     }
      else if (job->spin_polarisation == 1) {
      sprintf(fermi_homo, "%3d %3d",fermi.homo[0],fermi.homo[1]);
      sprintf(fermi_bands,"  %d - %d %d - %d",fermi.bands[0],fermi.bands[1],fermi.bands[2],fermi.bands[3]);
     }
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
      if (job->bse_spk == 0) 
      sprintf(linear_algebra,"%s","LAPACK");
      else if (job->bse_spk == 1) 
      sprintf(linear_algebra,"%s","SCALAPACK");

      sprintf(scale_factor,"%8.2f",job->scalefac);
      sprintf(scissor_shift,"%8.2f",job->scissor);
      //printf("%s %f %s %f \n",scale_factor, job->scalefac, scissor_shift, job->scissor);
      job->scissor /= au_to_eV;

      ////knet_check(fermi.is,jobname,crystal,file);
      ////int_range_check(fermi.bands,int_value,int_range,jobname,file);

      job->npoints = npoints;
      job->energy_range[0] = energy_range[0] / au_to_eV;
      job->energy_range[1] = energy_range[1] / au_to_eV;

      ////double_range_check(energy_range,double_value,double_range,jobname,file);

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
        ////bse_crystal1(&fermi,atoms,atom_p,&numfrag,natoms,nat,shells,gaussians,atoms_ax,shells_ax,gaussians_ax,crystal, \
        symmetry,R,R_tables,G,job,file);
        ////optical_spectrum_crystal2(&fermi,atoms,atom_p,&numfrag,natoms,nat,shells,gaussians,atoms_ax,shells_ax,gaussians_ax, \
        crystal,symmetry,R,R_tables,G,job,file);
        case 'S':
        case 'P':
        //bethe_salpeter_crystal(&fermi, atoms, atom_p, shells, gaussians, crystal, symmetry, R, R_tables, G, job, file);
        break;
        case 'M':
        double *temp3;
        //casida_molecule(&fermi, atoms, atom_p, shells, gaussians, atoms_ax, shells_ax, gaussians_ax, crystal, symmetry, \
        R, R_tables, G, job, file);
        //if (job->bse_spk == 0)
        //BSE_molecule(&fermi,atoms,atom_p,&numfrag,natoms,nat,shells,gaussians,atoms_ax,shells_ax,gaussians_ax,crystal,symmetry,R, \
        R_tables,G,job,file);
        //else if (job->bse_spk == 1) {
        bse_molecule(&fermi,atoms,atom_p,&numfrag,natoms,nat,shells,gaussians,atoms_ax,shells_ax,gaussians_ax,crystal, \
        symmetry,R,R_tables,G,job,file);
        optical_spectrum_molecule(&fermi,atoms,atom_p,&numfrag,natoms,nat,shells,gaussians,atoms_ax,shells_ax,gaussians_ax, \
        crystal,symmetry,R,R_tables,G,job,file);
        //BSE_molecule_spk(&fermi,atoms,atom_p,&numfrag,natoms,nat,shells,gaussians,atoms_ax,shells_ax,gaussians_ax,crystal, \
        symmetry,R,R_tables,G,job,file);
       }
        break;

       ////}

      ////free_fermi(&fermi,job);

    }

    /*

    // *****JOB: BETHE-SALPETER EQUATION WAVEFUNCTION PLOTTING ********************************************************

    else if (!strcmp(jobname, "BSE_PLOT") || !strcmp(jobname, "BSE_TRANS_DEN") || !strcmp(jobname, "ELECTRON_HOLE_PLOT")) {

      char hamiltonian[9], monkhorst_pack[9], mpi_io[80], spectrum[14], fermi_homo[9], fermi_bands[9], field[80];
      char title2[15];
      FERMI fermi;
      int electron_count = 0;
      int bse_state[2];
      int grid_type, plot_dim = 0;
      int grid_par[3];
      int num_points = 0;
      double grid_sum[job->spin_dim];
      VECTOR_DOUBLE points[4],points1[10];
      for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++)   // Default is uncharged, unpolarised system
      electron_count += atoms->atomic_number[i];
      //HEREfermi.homo = electron_count / 2;
      fermi.homo[0] = electron_count / 2;
      fermi.homo[1] = electron_count / 2;

      fermi.bands[0] = 4;
      fermi.bands[1] = 6;
      fermi.bands[2] = 4;
      fermi.bands[3] = 6;

      read_line(file.job, title, 99);
      sscanf(title, "%s", file.bse_eigvec);
      //printf("eigenvector directory: %s\n",file.bse_eigvec);

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

   do {

      read_line(file.job, title, 99);
      sscanf(title, "%s", jobname1);

      if (!strcmp(jobname1, "MONKHORST_PACK_NET")) {
      read_line(file.job, title, 99);
      switch (crystal->type[0]) {
      case 'C':
      sscanf(title, "%d %d %d", &fermi.is[0], &fermi.is[1], &fermi.is[2]);
      sprintf(monkhorst_pack,"%d %d %d",fermi.is[0], fermi.is[1], fermi.is[2]);
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
      sprintf(monkhorst_pack,"%d %d %d",fermi.is[0], fermi.is[1], fermi.is[2]);
      fprintf(file.out,"WARNING: MONKHORST-PACK NET Off for molecules\n");
      break;
     }
    }

      if (!strcmp(jobname1, "MO_RANGE")) {
      read_line(file.job, title, 99);
      sscanf(title, "%d %d", &fermi.bands[0], &fermi.bands[1]);
      sprintf(fermi_homo, "       %d",fermi.homo);
      sprintf(fermi_bands,"  %d - %d",fermi.bands[0],fermi.bands[1]);
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

    if (!strcmp(jobname, "BSE_PLOT")) {
      if (crystal->type[0] == 'M' && job->taskid == 0)
      plot_correlated_electron_hole(grid_par,points,num_points,points1,bse_state,&fermi,atoms,atom_p,shells,gaussians,crystal,\
      symmetry,R,R_tables,job,file);
      else if (crystal->type[0] == 'C' || crystal->type[0] == 'S' || crystal->type[0] == 'P') {
      plot_correlated_electron_hole_crystal(grid_sum,bse_state,grid_par,points,num_points,points1,&fermi,atoms,atom_p,shells,\
      gaussians,crystal,R,R_tables,G,symmetry,job,file);
     }
    }

    else if (!strcmp(jobname, "BSE_TRANS_DEN")) {
      if (crystal->type[0] == 'M' && job->taskid == 0)
      plot_bse_transition_density(grid_par,points,points1,bse_state,&fermi,atoms,atom_p,shells,gaussians,crystal,symmetry,R,\
      R_tables,job,file);
      else if (crystal->type[0] == 'C' || crystal->type[0] == 'S' || crystal->type[0] == 'P') {
      fprintf(file.out,"BSE_TRANS_DEN only available for molecules at present\n");
      MPI_Finalize();
      exit(0);
      //plot_bse_transition_density_crystal(grid_sum,bse_state,grid_par,points,num_points,points1,&fermi,atoms,atom_p,shells,\
      gaussians,crystal,R,R_tables,G,symmetry,job,file);
     }
    }

    else if (!strcmp(jobname, "ELECTRON_HOLE_PLOT")) {
      if (crystal->type[0] == 'M' && job->taskid == 0)
      plot_electron_hole_transition_density(grid_par,points,&fermi,atoms,atom_p,shells,gaussians,crystal,symmetry,R,R_tables,\
      job,file);
      else if (crystal->type[0] == 'C' || crystal->type[0] == 'S' || crystal->type[0] == 'P') {
      plot_electron_hole_transition_density_crystal(grid_sum,bse_state,grid_par,points,num_points,points1,&fermi,atoms,atom_p,\
      shells,gaussians,crystal,R,R_tables,G,symmetry,job,file);
      //fprintf(file.out,"ELECTRON_HOLE_PLOT only available for molecules at present\n");
      //MPI_Finalize();
      //exit(0);
     }
    }
      if (job->taskid == 0) printf("%s\n",jobname);

  //free_fermi(&fermi,job);

  }

    // *****JOB: BAND STRUCTURE ****************************************************

    else if (!strcmp(jobname, "BAND_STRUCTURE")) {

      int number_of_lines, number_of_points, npoints, segment_points[10], shr[3], count;
      double segment_length[10], length = k_zero;
      VECTOR_INT is1[10], is2[10];
      VECTOR_DOUBLE k1, k2, k3;
      FERMI fermi;
      //char fermi_homo[9];

      job->guess_type = 1;
      job->vectors    = 0;
      job->values     = 2;
      job->density    = 1;
      job->max_cycle  = 1;
      job->kpoints    = 1;

      int_value[0] = 1;
      int_value[1] = atoms->number_of_sh_bfns_in_unit_cell;
      //int_range = "ENERGY BAND";
      char int_range[12] = "ENERGY BAND";

      if (job->taskid == 0) 
      fprintf(file.out, "BAND STRUCTURE CALCULATION\n");

      read_line(file.job, title, 99);
      sscanf(title, "%d %d", &number_of_lines, &number_of_points);
      if (number_of_lines > 9) {
      if (job->taskid == 0) 
      fprintf(file.out,"Number of lines in band structure must be 9 or fewer\n");
      MPI_Finalize();
      exit(0);
     }
      read_line(file.job, title, 99);
      sscanf(title, "%d %d %d %d %d", &shr[0], &shr[1], &shr[2], &fermi.bands[0], &fermi.bands[1]);
 
      if (job->taskid == 0) {
      fprintf(file.out, "SHRINKING FACTOR %d %d %d\tLOWER BAND LIMIT %d\tUPPER BAND LIMIT %d\n\n", \
      shr[0], shr[1], shr[2], fermi.bands[0], fermi.bands[1]);
      //fprintf(file.out,"%d k-points in Reciprocal Space\n",number_of_points);
      //fprintf(file.out,"%d line Segments in Reciprocal Space\n",number_of_lines);
     }

      for (i = 0; i < number_of_lines; i++) {
      read_line(file.job, title, 99);
      sscanf(title, "%d %d %d %d %d %d", &is1[i].comp1, &is1[i].comp2, &is1[i].comp3, &is2[i].comp1, &is2[i].comp2, &is2[i].comp3);
      //fprintf(file.out, "\t%d %d %d\t %d %d %d\n", is1[i].comp1, is1[i].comp2, is1[i].comp3, \
      is2[i].comp1, is2[i].comp2, is2[i].comp3);
      getcart(&is1[i], &k1, shr, crystal);
      getcart(&is2[i], &k2, shr, crystal);
      double_vec_diff(&k1,&k2,&k3);
      segment_length[i] = sqrt(double_vec_dot(&k3,&k3));
      length += segment_length[i];
     }
      //fprintf(file.out,"\n");

      read_line(file.job, title, 99);
      sscanf(title, "%s", file.scf_eigvec);


      //knet_check(fermi.is,jobname,crystal,file);
      //int_range_check(fermi.bands,int_value,int_range,jobname,file);

      npoints = 0;
      segment_points[0] = 1;
      for (i = 0; i < number_of_lines; i++) {
      npoints += int ((number_of_points * segment_length[i]) / length);
      segment_points[i + 1] = npoints;
      //fprintf(file.out,"Number of segment points%d\n\n",segment_points[i+1]);
     }
      fermi.npoints = npoints;
      //fprintf(file.out,"Number of points used %d\n\n",npoints);

      fermi.knet_list = (VECTOR_KNET *) malloc(npoints * sizeof(VECTOR_KNET));

      count = 0;
      for (i = 0; i < number_of_lines; i++) {
      npoints = int ((number_of_points * segment_length[i]) / length);
      getcart(&is1[i], &k1, shr, crystal);
      getcart(&is2[i], &k2, shr, crystal);
      double_vec_diff(&k1,&k2,&k3);
      //fprintf(file.out,"%10.4lf %10.4lf %10.4lf\n",fermi.knet_list[0].cart.comp1, fermi.knet_list[0].cart.comp2,\
      fermi.knet_list[0].cart.comp3);
      for (j = 0; j < npoints; j++) {
      fermi.knet_list[count].cart.comp1 = k1.comp1 + (double) j * k3.comp1 / (double) npoints;
      fermi.knet_list[count].cart.comp2 = k1.comp2 + (double) j * k3.comp2 / (double) npoints;
      fermi.knet_list[count].cart.comp3 = k1.comp3 + (double) j * k3.comp3 / (double) npoints;
      //fprintf(file.out,"%10.4lf %10.4lf %10.4lf\n",fermi.knet_list[count].cart.comp1, fermi.knet_list[count].cart.comp2,\
      fermi.knet_list[count].cart.comp3);
      count++;
     }
     //fprintf(file.out,"\n");
    }

      scf_crystal(&fermi,atoms,atoms_ax,atom_p,shells,gaussians,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
      //crystal_scf(&fermi,atoms,atoms_ax,atom_p,shells,gaussians,shells_ax,gaussians_ax,crystal,symmetry,R,R_tables,G,job,file);
      //crystal_scf1(&fermi, atoms, atom_p, atom_i, shells, gaussians, crystal, symmetry, R, R_tables, G, job, file);
      //CHANGES2015crystal_scf(&fermi, atoms, atom_p, shells, gaussians, crystal, symmetry, R, R_tables, G, job, file);

      char buf1[110], xx[10] = "/evalfile";
      double *eigenvalues, emin = 9999.0, emax = -9999.0;
      int value_size  = (fermi.bands[1] - fermi.bands[0] + 1);
      int dim1 = job->spin_dim * fermi.npoints * value_size;
      FILE *band_structure_0, *band_structure_1, *band_gnu;
      AllocateDoubleArray(&eigenvalues,&dim1,job);
      strcpy(buf1,file.scf_eigvec);
      strcat(buf1,xx);

      MPI_File eh;
      MPI_File_open(MPI_COMM_WORLD,buf1,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&eh) ;
      MPI_File_seek(eh, 0, MPI_SEEK_SET) ;
      MPI_File_read(eh, eigenvalues, job->spin_dim * fermi.npoints * value_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
      MPI_File_close(&eh);

      if (job->taskid == 0) {
      band_structure_0 = fopen("band_structure_0.dat", "w");
      if (band_structure_0 == NULL) { fprintf(file.out, "cannot open file band_structure_0.dat\n"); MPI_Finalize(); exit(1); }
      for (i = 0; i < value_size; i++) {
      count = 1;
      for (j = 0; j < fermi.npoints; j++) {
      //printf("%3d %d %d %d %d %d %10.4lf\n",count,i,value_size,fermi.npoints,j,i + value_size * j,eigenvalues[i + 2 * value_size * j]);
      fprintf(band_structure_0,"%3d %10.4lf\n",count, (eigenvalues[i + job->spin_dim * value_size * j] - job->fermi_energy) * au_to_eV);
      if (eigenvalues[i + j * job->spin_dim * value_size] * au_to_eV > emax) {
      emax = eigenvalues[i + j * job->spin_dim * value_size] * au_to_eV;
     }
      if (eigenvalues[i + j * job->spin_dim * value_size] * au_to_eV < emin) {
      emin = eigenvalues[i + j * job->spin_dim * value_size] * au_to_eV;
     }
      count++;
     }
     fprintf(band_structure_0,"\n");
    }
     fflush(band_structure_0);

      if (job->spin_polarisation == 1) {
      band_structure_1 = fopen("band_structure_1.dat", "w");
      if (job->taskid == 0 && band_structure_1 == NULL) { fprintf(file.out, "cannot open file band_structure_1.dat\n"); 
      MPI_Finalize(); exit(1); }
      for (i = 0; i < value_size; i++) {
      count = 1;
      for (j = 0; j < fermi.npoints; j++) {
      fprintf(band_structure_1,"%3d %10.4lf\n",count, (eigenvalues[i + value_size + job->spin_dim * value_size * j] - job->fermi_energy) * \
      au_to_eV);
      //printf("%3d %d %d %d %d %d %10.4lf\n",count, i, value_size, fermi.npoints, j, \
      i + value_size + job->spin_dim * value_size * j,eigenvalues[i + value_size + job->spin_dim * value_size * j]);
      if (eigenvalues[i + value_size + job->spin_dim * value_size * j] * au_to_eV > emax) {
      emax = eigenvalues[i + value_size + job->spin_dim * value_size * j] * au_to_eV;
     }
      if (eigenvalues[i + value_size + job->spin_dim * value_size * j] * au_to_eV < emin) {
      emin = eigenvalues[i + value_size + job->spin_dim * value_size * j] * au_to_eV;
     }
      count++;
     }
     fprintf(band_structure_1,"\n");
    }
     fflush(band_structure_1);
   }
   emax -= job->fermi_energy * au_to_eV;
   emin -= job->fermi_energy * au_to_eV;
 }

     if (job->taskid == 0) {
     band_gnu = fopen("band_structure.gnu", "w");
     if (band_gnu == NULL) { fprintf(file.out, "cannot open file band.gnu\n"); MPI_Finalize(); exit(1); }
     fprintf(band_gnu, "set term x11 enhanced font 'Helvetica,20'\n");
     fprintf(band_gnu, "set grid xtics noytics lt -1 lw 1.0\n");
     fprintf(band_gnu, "set border 31 lt -1 lw 1.5\n");
     fprintf(band_gnu, "set nokey\n");
     fprintf(band_gnu, "set size ratio 1.8\n");
     fprintf(band_gnu, "set xrange[0:%d]\n",fermi.npoints - 1);
     fprintf(band_gnu, "set yrange[%lf:%lf]\n",floor(emin/5)*5,ceil(emax/5)*5);
     fprintf(band_gnu, "set ylabel 'Energy (eV)'\n");
     fprintf(band_gnu, "set tics in\n");
     fprintf(band_gnu, "set xtics border mirror norotate (");
     for (i = 0; i < number_of_lines; i++) {
     fprintf(band_gnu, "'X' %d, ",segment_points[i]);
     //CHANGES2014fprintf(band_gnu, "'X' %d, ",segment_points[i] - 1);
    }
     fprintf(band_gnu, "'X' %d)\n",segment_points[number_of_lines]);
     fprintf(band_gnu, "set ytics border mirror norotate\n");
     fprintf(band_gnu, "f(x) = 0.00\n");
     fprintf(band_gnu, "plot 'band_structure_0.dat' using 1:2 w l lw 2.0, f(x) lw 2.0\n");
     fprintf(band_gnu, "set output 'band_structure_0.eps'\n");
     fprintf(band_gnu, "set term post portrait enhanced monochrome\n");
     fprintf(band_gnu, "plot 'band_structure_0.dat' using 1:2 w l lw 2.0, f(x) lw 2.0\n");
     fprintf(band_gnu, "pause -1\n");
     if (job->spin_polarisation == 1) {
     fprintf(band_gnu, "plot 'band_structure_1.dat' using 1:2 w l lw 2.0, f(x) lw 2.0\n");
     fprintf(band_gnu, "set output 'band_structure_1.eps'\n");
     fprintf(band_gnu, "set term post portrait enhanced monochrome\n");
     fprintf(band_gnu, "plot 'band_structure_1.dat' using 1:2 w l lw 2.0, f(x) lw 2.0\n");
     fprintf(band_gnu, "pause -1\n");
    }
   }
 
      free(fermi.knet_list);

      DestroyDoubleArray(&eigenvalues,&dim1,job);
      free_fermi(&fermi,job);

    }

    // *****JOB: BULK BAND PROJECTION ****************************************************

    else if (!strcmp(jobname, "BULK_BAND_PROJECTION")) {

      int number_of_lines, number_of_points, npoints, segment_points[10], shr[3], count;
      int map, count1;
      double dum1, dum2, dum3;
      double segment_length[10], length = k_zero;
      double *eigenvalues;
      char buf1[110], xx[10] = "/evalfile";
      VECTOR_INT is1[10], is2[10], is3, is4, is5;
      VECTOR_DOUBLE k1, k2, k3, k4, k5, k6;
      FERMI fermi;
      FILE *band_project_0, *band_project_1;

      band_project_0 = fopen("band_project_0.dat", "w");
      if (band_project_0 == NULL) { fprintf(file.out, "cannot open file band_project_0.dat\n"); MPI_Finalize(); exit(1); }
      band_project_1 = fopen("band_project_1.dat", "w");
      if (band_project_1 == NULL) { fprintf(file.out, "cannot open file band_project_1.dat\n"); MPI_Finalize(); exit(1); }

      job->guess_type = 1;
      job->vectors    = 0;
      job->values     = 2;
      job->density    = 1;
      job->max_cycle  = 1;
      job->kpoints    = 1;

      int_value[0] = 1;
      int_value[1] = atoms->number_of_sh_bfns_in_unit_cell;
      //int_range = "ENERGY BAND";
      char int_range[12] = "ENERGY BAND";

      read_line(file.job, title, 99);
      sscanf(title, "%s", file.scf_eigvec);

      read_line(file.job, title, 99); // title of plot
      sscanf(title, "%99s", title3);

      read_line(file.job, title, 99); // normal to the slab plane
      sscanf(title, "%d %d %d",&is3.comp1, &is3.comp2, &is3.comp3);

      read_line(file.job, title, 99); // normal to the slab plane
      sscanf(title, "%d %d %d",&is4.comp1, &is4.comp2, &is4.comp3);

      read_line(file.job, title, 99); // normal to the slab plane
      sscanf(title, "%d %d %d",&is5.comp1, &is5.comp2, &is5.comp3);

      read_line(file.job, title, 99);
      sscanf(title, "%d %d", &number_of_lines, &number_of_points);
      if (number_of_lines > 9) {
      if (job->taskid == 0) 
      fprintf(file.out,"Number of lines in band structure must be 9 or fewer\n");
      MPI_Finalize();
      exit(0);
     }

      read_line(file.job, title, 99);
      sscanf(title, "%d %d %d %d %d", &shr[0], &shr[1], &shr[2], &fermi.bands[0], &fermi.bands[1]);
 
      for (i = 0; i < number_of_lines; i++) {
      read_line(file.job, title, 99);
      sscanf(title, "%d %d %d %d %d %d", &is1[i].comp1, &is1[i].comp2, &is1[i].comp3, &is2[i].comp1, &is2[i].comp2, &is2[i].comp3);
      //fprintf(file.out, "\t%d %d %d\t %d %d %d\n", is1[i].comp1, is1[i].comp2, is1[i].comp3, is2[i].comp1, is2[i].comp2, is2[i].comp3);
      getcart(&is1[i], &k1, shr, crystal);
      getcart(&is2[i], &k2, shr, crystal);
      double_vec_diff(&k1,&k2,&k3);
      segment_length[i] = sqrt(double_vec_dot(&k3,&k3));
      length += segment_length[i];
     }
      fprintf(file.out,"\n");

      //knet_check(fermi.is,jobname,crystal,file);
      //int_range_check(fermi.bands,int_value,int_range,jobname,file);

      npoints = 0;
      segment_points[0] = 1;
      for (i = 0; i < number_of_lines; i++) {
      npoints += int ((number_of_points * segment_length[i]) / length);
      segment_points[i + 1] = npoints;
      //fprintf(file.out,"Number of segment points%d\n\n",segment_points[i+1]);
     }

      fermi.npoints = 0;
      for (map = 2; map < 3 ; map++) {
      for (i = 0; i < number_of_lines; i++) {
      npoints = int ((number_of_points * segment_length[i]) / length);
      fermi.npoints += npoints * npoints;
     }
     }

      int value_size  = (fermi.bands[1] - fermi.bands[0] + 1);
      int dim1 = job->spin_dim * fermi.npoints * value_size;
      AllocateDoubleArray(&eigenvalues,&dim1,job);
      strcpy(buf1,file.scf_eigvec);
      strcat(buf1,xx);

      fermi.knet_list = (VECTOR_KNET *) malloc(fermi.npoints * sizeof(VECTOR_KNET));

      getcart(&is3, &k4, shr, crystal);
      fprintf(file.out,"k4 %10.4lf %10.4lf %10.4lf\n",k4.comp1,k4.comp2,k4.comp3);

      getcart(&is4, &k5, shr, crystal);
      fprintf(file.out,"k5 %10.4lf %10.4lf %10.4lf\n",k5.comp1,k5.comp2,k5.comp3);

      getcart(&is5, &k6, shr, crystal);
      fprintf(file.out,"k6 %10.4lf %10.4lf %10.4lf\n",k6.comp1,k6.comp2,k6.comp3);

      for (map = 0; map < 2 ; map++) {
      count = 0;
      for (i = 0; i < number_of_lines; i++) {
      npoints = int ((number_of_points * segment_length[i]) / length);
      getcart(&is1[i], &k1, shr, crystal);
      getcart(&is2[i], &k2, shr, crystal);
      fprintf(file.out,"k1 %3d %3d %3d   %10.4lf %10.4lf %10.4lf\n",is1[0],is1[1],is1[2],k1.comp1,k1.comp2,k1.comp3);
      fprintf(file.out,"k2 %3d %3d %3d   %10.4lf %10.4lf %10.4lf\n",is1[0],is1[1],is1[2],k2.comp1,k2.comp2,k2.comp3);
      k2.comp1 += map * k5.comp1 ;
      k2.comp2 += map * k5.comp2 ;
      k2.comp3 += map * k5.comp3 ;
      k1.comp1 += map * k5.comp1 ;
      k1.comp2 += map * k5.comp2 ;
      k1.comp3 += map * k5.comp3 ;
      double_vec_diff(&k1,&k2,&k3);
      count1 = 0;
      for (j = 0; j < npoints; j++) {
      dum1 = k1.comp1 + (double) j * k3.comp1 / (double) npoints;
      dum2 = k1.comp2 + (double) j * k3.comp2 / (double) npoints;
      dum3 = k1.comp3 + (double) j * k3.comp3 / (double) npoints;
      for (k = 0; k < npoints ; k++) {
      fermi.knet_list[count1 + k].cart.comp1 = dum1 + (double) k * k4.comp1 / (double) npoints;
      fermi.knet_list[count1 + k].cart.comp2 = dum2 + (double) k * k4.comp2 / (double) npoints;
      fermi.knet_list[count1 + k].cart.comp3 = dum3 + (double) k * k4.comp3 / (double) npoints;
      fprintf(file.out,"fermi.knet_list[count1].cart.comp1 %3d %3d %3d %3d %10.4f %10.4lf %10.4lf\n",map,i,j,k, fermi.knet_list[count1 + k].cart.comp1,\
      fermi.knet_list[count1 + k].cart.comp2,fermi.knet_list[count1 + k].cart.comp3);
      }
      count1 += npoints;
     }
    }
   }

      crystal_scf(&fermi, atoms, atoms_ax, atom_p, shells, gaussians, shells_ax, gaussians_ax, crystal, symmetry, R, R_tables, G, job, file);
      //CHANGES2015crystal_scf(&fermi, atoms, atom_p, shells, gaussians, crystal, symmetry, R, R_tables, G, job, file);

      MPI_File eh;
      MPI_File_open(MPI_COMM_WORLD,buf1,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&eh) ;
      MPI_File_seek(eh, 0, MPI_SEEK_SET) ;
      MPI_File_read(eh, eigenvalues, job->spin_dim * fermi.npoints * value_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
      MPI_File_close(&eh);

      for (map = 0; map < 2 ; map++) {
      count = 0;
      for (i = 0; i < number_of_lines; i++) {
      npoints = int ((number_of_points * segment_length[i]) / length);
      for (j = 0; j < npoints; j++) {
      for (k = 0; k < npoints * value_size ; k++) {
      fprintf(band_project_0,"%3d %10.4lf\n",j, (eigenvalues[count] - job->fermi_energy) * au_to_eV);
      count++;
      }
     }
    }
   }

    if (job->taskid == 0) {
    fprintf(file.out,"===========================================================================================================\n");
    fprintf(file.out,"|                                         BULK BAND PROJECTION                                            |\n");
    fprintf(file.out,"===========================================================================================================\n");
    fprintf(file.out,"| %-99s     |\n", title3);
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    fprintf(file.out,"| SHRINKING FAC %2d %2d %2d  | NUMBER OF LINES %3d  | NUMBER OF POINTS %2d     | BAND RANGE %4d TO %4d |\n", 
    shr[0], shr[1], shr[2], number_of_lines, npoints, fermi.bands[0], fermi.bands[1]);
    fprintf(file.out,"===========================================================================================================\n");
   }

      free(fermi.knet_list);
      //free_fermi(&fermi,job);
      DestroyDoubleArray(&eigenvalues,&dim1,job);

    }

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

    // *****JOB: BAND PLOT ***************************************************

    else if (!strcmp(jobname, "BAND_PLOT")) {

      int nproj;
      int nat1[20], natoms[5];
      int is_plot[3];
      int vectors[2], nkpoints;
      double atom_proj[atoms->number_of_sh_bfns_in_unit_cell][5];
      KPOINT_TRAN knet;

      job->guess_type = 1;
      job->vectors    = 2;
      job->values     = 0;
      job->density    = 1;
      job->max_cycle  = 1;
      job->kpoints    = 1;

      read_line(file.job, title, 99);
      sscanf(title, "%s", file.directory1);

      read_line(file.job, title, 99); // title of plot
      sscanf(title, "%99s", title3);

      read_line(file.job, title, 99); // shrinking factor for plots
      sscanf(title, "%d %d %d", &is_plot[0], &is_plot[1], &is_plot[2]);

      read_line(file.job, title, 99); // number of k points, first and last vector
      sscanf(title, "%d %d %d",&nkpoints, &vectors[0], &vectors[1]);

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

       read_line(file.job, title, 99); // number of projections
       sscanf(title, "%d",&nproj);

       if (nproj < 1 || nproj > 4) {
       if (job->taskid == 0) {
       fprintf(file.out,"ERROR: Number of projections in BAND_PLOT %3d must lie beween 1 and 4\n",nproj);
      }
       MPI_Finalize();
       exit(1);
      }

       natoms[0] = 0;
       natoms[1] = 0;
       natoms[2] = 0;
       natoms[3] = 0;
       read_line(file.job, title, 99); // number of atoms to project onto, first projection is unit operator
       sscanf(title, "%d %d %d %d", &natoms[0], &natoms[1], &natoms[2], &natoms[3]);
       int nat[nproj][20];
       for (i = 0; i < nproj; i++) {
       for (j = 0; j < 20; j++) {
       nat[i][j] = 0;
      }
      }

       for (j = 0; j < atoms->number_of_sh_bfns_in_unit_cell; j++) 
       atom_proj[j][0] = k_one;
       for (i = 0; i < nproj; i++) {
       for (j = 0; j < atoms->number_of_sh_bfns_in_unit_cell; j++) 
       atom_proj[j][i + 1] = k_zero;
      }

       for (i = 0; i < nproj; i++) {
       read_line(file.job, title, 99); // list of atoms to project onto
       sscanf(title, "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d", 
       &nat1[0],&nat1[1], &nat1[2], &nat1[3], &nat1[4], &nat1[5], &nat1[6], &nat1[7], &nat1[8], &nat1[9],
       &nat1[10], &nat1[11], &nat1[12], &nat1[13], &nat1[14], &nat1[15], &nat1[16], &nat1[17], &nat1[18], &nat1[19]);
       for (j = 0; j < 20; j++) 
       nat[i][j] = nat1[j];
       for (j = 0; j < natoms[i]; j++) {
         for (k = 0; k < atoms->bfnnumb_sh[nat1[j] - 1]; k++) {
           atom_proj[atoms->bfnposn_sh[nat1[j] - 1] + k][i + 1] = k_one;
          }
         }
        }

    if (job->taskid == 0) {
    fprintf(file.out,"===========================================================================================================\n");
    fprintf(file.out,"|                                      ATOM PROJECTED BAND STRUCTURE                                      |\n");
    fprintf(file.out,"===========================================================================================================\n");
    fprintf(file.out,"| %-99s     |\n", title3);
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    fprintf(file.out,"| SHRINKING FAC %3d %3d %3d | NUMBER OF K POINTS %3d  |NUMBER OF PROJECTIONS %2d | BAND RANGE %4d TO %4d |\n", 
    is_plot[0], is_plot[1], is_plot[2], nkpoints, nproj, vectors[0], vectors[1]);
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    fprintf(file.out,"|                                      ATOMS FOR PROJECTION                                               |\n");
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    for (i = 0; i < nproj; i++) {
    fprintf(file.out,"| PROJECTION %3d:",i + 1);
    for (j = 0; j < 20; j++) {
    if (j < natoms[i]) fprintf(file.out," %3d",nat[i][j]);    
    else fprintf(file.out,"    ");
   }
    fprintf(file.out,"         |\n");
   }
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    fprintf(file.out,"|                         K POINTS IN SHRINKING FACTOR AND ANGS^-1 UNITS                                  |\n");
    fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
    for (i = 0; i < nkpoints; i++) {
    j = kplot[i].comp1 * is_plot[1] * is_plot[2] + kplot[i].comp2 * is_plot[2] + kplot[i].comp3;
    fprintf(file.out,"|%3d            %3d %3d %3d          %10.4lf %10.4lf %10.4lf                                     |\n", \
    i + 1, kplot1[i].comp1, kplot1[i].comp2, kplot1[i].comp3,knet.cart[j].comp1 * bohr_to_AA,knet.cart[j].comp2 * bohr_to_AA, \
    knet.cart[j].comp3 * bohr_to_AA);
   }
    fprintf(file.out,"===========================================================================================================\n");
   }

    band_plot(is_plot,vectors,nkpoints,nproj,atom_proj,&knet,kplot,kplot1,atom_p,atoms,shells,gaussians,crystal,R,R_tables,G,symmetry,job,file);

    free_k_points(&knet,job);

    }

    // *****JOB: BULK BAND PROJECTION ***************************************************

    else if (!strcmp(jobname, "BULK_BAND_PROJECTION")) {

      KPOINT_TRAN knet;

      job->guess_type = 1;
      job->vectors    = 2;
      job->values     = 0;
      job->density    = 1;
      job->max_cycle  = 1;
      job->kpoints    = 1;

      read_line(file.job, title, 99);
      sscanf(title, "%s", file.directory1);

      int is_plot[3], is_proj[3];
      int vectors[2], ksize;

      read_line(file.job, title, 99); // title of plot
      sscanf(title, "%99s", title3);
      if (job->taskid == 0)
        fprintf(file.out, "TITLE OF PLOT %8s\n", title3);

      read_line(file.job, title, 99); // shrinking factor for plots
      sscanf(title, "%d %d %d", &is_plot[0], &is_plot[1], &is_plot[2]);

      read_line(file.job, title, 99); // direction in k-space for projection in oblique coordinates
      sscanf(title, "%d %d %d",&is_proj[0], &is_proj[1], &is_proj[2]);

      read_line(file.job, title, 99); // first and last vector
      sscanf(title, "%d %d",&vectors[0], &vectors[1]);

      if (job->taskid == 0)
        fprintf(file.out, "SHRINKING FACTOR %d %d %d\nBANDS FROM %d TO %d\n", is_plot[0],
            is_plot[1], is_plot[2], vectors[0], vectors[1]);

      //knet_size(&ksize,is_plot,crystal);
      //count_k_points(num_sym,is_plot,&nkunique,crystal,symmetry,job,file);
      //count_k_points(is_plot,&nkunique,&ksize,crystal,symmetry,job,file);
      //allocate_k_points(&knet,nkunique,ksize,crystal,job,file);
      //read_XCBD_crystal_09(nkunique,&knet,crystal,job,file);
      // change 6generate_k_points(&knet,num_sym,is_plot,nkunique,crystal,symmetry,job,file);
      count_k_points(&knet,is_plot,crystal,symmetry,job,file);
      allocate_k_points(&knet,crystal,job,file);
      if (job->C09 == 1)
      read_XCBD_crystal_09(&knet,crystal,job,file);
      generate_k_points(&knet,is_plot,crystal,symmetry,job,file);

      //int i;
      int nbands = vectors[1] - vectors[0] + 1;
      int dim1 = job->spin_dim * nbands * knet.unique;
      //int dim1 = job->spin_dim * nbands * nkunique;
      double *eigval;
      AllocateDoubleArray(&eigval,&dim1,job);
      ResetDoubleArray(eigval,&dim1);

      MPI_File gh;

      char buf[110];
      char xx[10] = "/evalfile";
      strcpy(buf,file.directory1);
      strcat(buf,xx);

      MPI_File_open(MPI_COMM_WORLD,buf,MPI_MODE_RDWR | MPI_MODE_CREATE,MPI_INFO_NULL,&gh) ;
      MPI_File_seek(gh, 0, MPI_SEEK_SET) ;
      MPI_File_read(gh, eigval, job->spin_dim * knet.unique * nbands, MPI_DOUBLE, MPI_STATUS_IGNORE) ;
      //MPI_File_read(gh, eigval, job->spin_dim * nkunique * nbands, MPI_DOUBLE, MPI_STATUS_IGNORE) ;

      if (job->taskid == 0 && job->verbosity > 1) {
        fprintf(file.out,"EIGENVALUES FOR PLOTTING\n");
        for (k = 0; k < knet.unique; k++) {
        //for (k = 0; k < nkunique; k++) {
          for (i = 0; i < nbands; i++)     {
            fprintf(file.out,"%3d %3d   %10.4lf\n",k,i, eigval[k * nbands + i]);
           }
          }
         fflush(file.out);
        }

        ////bulk_band_projection(is_plot,is_proj,vectors,&knet,eigval,atom_p,atom_i,atoms,shells,gaussians,crystal,R,R_tables,G,symmetry,job,file);

        free_k_points(&knet,job);
        DestroyDoubleArray(&eigval,&dim1,job);
        free(eigval);
        MPI_File_close(&gh);

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
      sprintf(monkhorst_pack,"%d %d %d",fermi->is[0], fermi->is[1], fermi->is[2]);
      break;

    case 'S':
      fermi->is[0] = 4;
      fermi->is[1] = 4;
      fermi->is[2] = 1;
      sprintf(monkhorst_pack,"%d %d %d",fermi->is[0], fermi->is[1], fermi->is[2]);
      break;

    case 'P':
      fermi->is[0] = 1;
      fermi->is[1] = 1;
      fermi->is[2] = 4;
      sprintf(monkhorst_pack,"%d %d %d",fermi->is[0], fermi->is[1], fermi->is[2]);
      break;

    case 'M':
      fermi->is[0] = 1;
      fermi->is[1] = 1;
      fermi->is[2] = 1;
      sprintf(monkhorst_pack,"%s","Off");
      break;

    } // close switch
     
      int electron_count = 0;
      for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++)   // Default is uncharged, unpolarised system
      electron_count += atoms->atomic_number[i];
      fermi->homo[0] = electron_count / 2;
      fermi->homo[1] = electron_count / 2;
      //HERE fermi->homo = electron_count / 2;

      //HERE if (fermi->homo > 4) {
      //fermi->bands[0] = fermi->homo - 4;                    // Default is 4 active occupied and 4 active unoccupied orbitals
      //fermi->bands[1] = fermi->homo + 4;
      //fermi->bands[2] = fermi->homo - 4;                    // Default is 4 active occupied and 4 active unoccupied orbitals
      //fermi->bands[3] = fermi->homo + 4;
      if (fermi->homo[0] > 4) {
      fermi->bands[0] = fermi->homo[0] - 4;                    // Default is 4 active occupied and 4 active unoccupied orbitals
      fermi->bands[1] = fermi->homo[0] + 4;
      fermi->bands[2] = fermi->homo[1] - 4;                    // Default is 4 active occupied and 4 active unoccupied orbitals
      fermi->bands[3] = fermi->homo[1] + 4;
      if (fermi->bands[1] > atoms->number_of_sh_bfns_in_unit_cell) {
      fermi->bands[1] = atoms->number_of_sh_bfns_in_unit_cell;
      fermi->bands[3] = atoms->number_of_sh_bfns_in_unit_cell; }
     }
      //else if (fermi->homo <= 4) {
      //fermi->bands[0] = 0;           
      //fermi->bands[1] = fermi->homo + 1;           
      //fermi->bands[2] = 0;           
      //fermi->bands[3] = fermi->homo + 1;           
     //}
      else if (fermi->homo[0] <= 4) {
      fermi->bands[0] = 0;           
      fermi->bands[1] = fermi->homo[0] + 1;           
      fermi->bands[2] = 0;           
      fermi->bands[3] = fermi->homo[1] + 1;           
     }
      else {
      if (job->taskid == 0)
      fprintf(file.out,"Error in default values for fermi->bands in BSE calculation input\n");
      MPI_Finalize();
      exit(1);
    }

/*

        switch (crystal->type[0]) {

    case 'C':
    case 'S':
    case 'P':
      sprintf(fermi_homo, "       %d",fermi->homo[0]);
      sprintf(fermi_bands,"  %d - %d",fermi->bands[0],fermi->bands[1]);
      if (job->spin_polarisation == 1)
      sprintf(fermi_bands,"  %d - %d",fermi->bands[2],fermi->bands[3]);
      break;
    case 'M':
      sprintf(fermi_homo, "       %d",fermi->homo[0]);
      sprintf(fermi_bands,"  %d - %d",fermi->bands[0],fermi->bands[1]);
      if (job->spin_polarisation == 1)
      sprintf(fermi_bands,"  %d - %d",fermi->bands[2],fermi->bands[3]);
      break;

    } // close switch
*/
}
