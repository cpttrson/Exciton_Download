#include <mpi.h>
#include <cstring>
#include "conversion_factors.h"
#include "myconstants.h"
#include "mylogical.h"
#include "USER_DATA.h"
#include "MATRIX_UTIL.h"
#include "SYMMETRY.h"
#include "TOOLS.h"
#include "SETUP_ATOMS.h"

using namespace std;

  // ******************************************************************************************
  // *  Count unique atomic positions                                                         *
  // ******************************************************************************************

void count_unique_atoms(ATOM *atoms, JOB_PARAM *job, FILES file)

{

  char title[150];

  read_line(file.in, title, 150);
  sscanf(title, "%3d", &atoms->number_of_unique_atoms);
  if (job->taskid == 0 && job->verbosity > 1)
  fprintf(file.out, "%5d UNIQUE IONS\n", atoms->number_of_unique_atoms);

}

void read_unique_atoms(ATOM *atoms, ATOMTYPE *iratom, CRYSTAL *crystal, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // *  Read unique atomic positions                                                          *
  // ******************************************************************************************

  int i;
  char title[150];
  int line_number;
  VECTOR_DOUBLE fractional_coord[atoms->number_of_unique_atoms];

  line_number = 0;

  switch (crystal->type[0]) {

    case 'C':

      for (i = 0; i < atoms->number_of_unique_atoms; i++) {
        read_line(file.in, title, 150);
        line_number++;
        sscanf(title, "%d %lf %lf %lf", &iratom[i].atomic_number, &fractional_coord[i].comp1,
        &fractional_coord[i].comp2, &fractional_coord[i].comp3);
        if (job->taskid == 0 && job->verbosity > 1)
        fprintf(file.out,"%3d %20.15lf %20.15lf %20.15lf\n",iratom[i].atomic_number,fractional_coord[i].comp1,\
        fractional_coord[i].comp2,fractional_coord[i].comp3) ;
       }

      for (i = 0; i < atoms->number_of_unique_atoms; i++) {
        iratom[i].posn.comp1 = fractional_coord[i].comp1 * crystal->conventional_cell[0].comp1 + fractional_coord[i].comp2
        * crystal->conventional_cell[1].comp1 + fractional_coord[i].comp3 * crystal->conventional_cell[2].comp1;
        iratom[i].posn.comp2 = fractional_coord[i].comp1 * crystal->conventional_cell[0].comp2 + fractional_coord[i].comp2
        * crystal->conventional_cell[1].comp2 + fractional_coord[i].comp3 * crystal->conventional_cell[2].comp2;
        iratom[i].posn.comp3 = fractional_coord[i].comp1 * crystal->conventional_cell[0].comp3 + fractional_coord[i].comp2
        * crystal->conventional_cell[1].comp3 + fractional_coord[i].comp3 * crystal->conventional_cell[2].comp3;
        if (job->taskid == 0 && job->verbosity > 1)
        fprintf(file.out,"%3d %20.15lf %20.15lf %20.15lf\n",iratom[i].atomic_number,iratom[i].posn.comp1/bohr_to_AA,\
        iratom[i].posn.comp2/bohr_to_AA,iratom[i].posn.comp3/bohr_to_AA) ;
       }

      break;

    case 'S':

      for (i = 0; i < atoms->number_of_unique_atoms; i++) {
        read_line(file.in, title, 150);
        line_number++;
        sscanf(title, "%d %lf %lf %lf", &iratom[i].atomic_number, &fractional_coord[i].comp1,
        &fractional_coord[i].comp2, &iratom[i].posn.comp3);
        iratom[i].posn.comp3 /= bohr_to_AA;
        if (job->taskid == 0 && job->verbosity > 1)
        fprintf(file.out,"%3d %20.15lf %20.15lf %20.15lf\n",iratom[i].atomic_number,fractional_coord[i].comp1,\
        fractional_coord[i].comp2,iratom[i].posn.comp3) ;
       }

      for (i = 0; i < atoms->number_of_unique_atoms; i++) {
        iratom[i].posn.comp1 = fractional_coord[i].comp1 * crystal->primitive_cell[0].comp1 + fractional_coord[i].comp2
        * crystal->primitive_cell[1].comp1;
        iratom[i].posn.comp2 = fractional_coord[i].comp1 * crystal->primitive_cell[0].comp2 + fractional_coord[i].comp2
        * crystal->primitive_cell[1].comp2;
        if (job->taskid == 0 && job->verbosity > 1)
        fprintf(file.out,"%3d %20.15lf %20.15lf %20.15lf\n",iratom[i].atomic_number,iratom[i].posn.comp1/bohr_to_AA,\
        iratom[i].posn.comp2/bohr_to_AA,iratom[i].posn.comp3) ;
       }

      break;

    case 'P':

      for (i = 0; i < atoms->number_of_unique_atoms; i++) {
        read_line(file.in, title, 150);
        line_number++;
        sscanf(title, "%d %lf %lf %lf", &iratom[i].atomic_number, &iratom[i].posn.comp1, &iratom[i].posn.comp2,
        &fractional_coord[i].comp3);
        iratom[i].posn.comp1 /= bohr_to_AA;
        iratom[i].posn.comp2 /= bohr_to_AA;
        if (fabs(fractional_coord[i].comp3) > k_one) {
        if (job->taskid == 0 && job->verbosity > 1)
        fprintf(file.out,"Third coordinate %10.5e for polymer ion %4d must be a fractional coordinate in range [0,1] \n",fractional_coord[i].comp3,i+1);
        MPI_Finalize();
        exit(1);
       } 
      }

      for (i = 0; i < atoms->number_of_unique_atoms; i++) {
        iratom[i].posn.comp3 = fractional_coord[i].comp3 * crystal->primitive_cell[0].comp3;
        if (job->taskid == 0 && job->verbosity > 1)
        fprintf(file.out,"Unique atom coords %3d %20.15lf %20.15lf %20.15lf\n",iratom[i].atomic_number,iratom[i].posn.comp1,\
        iratom[i].posn.comp2,iratom[i].posn.comp3) ;
      }

      break;

    case 'M':

      for (i = 0; i < atoms->number_of_unique_atoms; i++) {
        read_line(file.in, title, 150);
        line_number++;
        sscanf(title, "%d %lf %lf %lf", &iratom[i].atomic_number, &iratom[i].posn.comp1, &iratom[i].posn.comp2, &iratom[i].posn.comp3);
        iratom[i].posn.comp1 /= bohr_to_AA;
        iratom[i].posn.comp2 /= bohr_to_AA;
        iratom[i].posn.comp3 /= bohr_to_AA;
        if (job->taskid == 0 && job->verbosity > 1)
        fprintf(file.out,"%3d %20.15lf %20.15lf %20.15lf\n",iratom[i].atomic_number,iratom[i].posn.comp1,iratom[i].posn.comp2,iratom[i].posn.comp3);
      }

      break;

  } // close switch (crystal->type[0]) {

  read_line(file.in, title, 150); // reading number of atoms
  line_number++;

}

void count_all_atoms(CRYSTAL *crystal, ATOM *atoms, ATOMTYPE *iratom, SYMMETRY *symmetry, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Count the number of all atoms in the unit cell for memory allocation                   *
  // ******************************************************************************************

  int mxatm;
  int atm;
  int i, k;
  double *p_irr;
  mxatm = atoms->number_of_unique_atoms * symmetry->number_of_operators;
  VECTOR_DOUBLE cell_vector_irr[atoms->number_of_unique_atoms];
  VECTOR_DOUBLE cell_vector_tmp[mxatm];
  VECTOR_DOUBLE cell_vec1, O_tmp, Rvec_tmp[2];

  for (i = 0; i < mxatm; i++) {
    cell_vector_tmp[i].comp1 = 9999.9;
    cell_vector_tmp[i].comp2 = 9999.9;
    cell_vector_tmp[i].comp3 = 9999.9;
   }

   atoms->number_of_atoms_in_unit_cell = 0 ;

         switch (crystal->type[0]) {

           case 'C':
           case 'S':
           case 'P':

  // ******************************************************************************************
  // * Map unique ion positions to zeroth cell                                                *
  // ******************************************************************************************

              for (i = 0; i < atoms->number_of_unique_atoms; i++)  
                map_to_wigner(crystal,&iratom[i].posn, &cell_vector_irr[i], &Rvec_tmp[0]);

              if (job->taskid == 0 && job->verbosity > 1) {
              fprintf(file.out, "  S             S VECTOR               WIGNER MAPPED S\n");
              for (i = 0; i < atoms->number_of_unique_atoms; i++)  {
                fprintf(file.out, "%3d %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf %10.6lf\n", i, iratom[i].posn.comp1,
                iratom[i].posn.comp2, iratom[i].posn.comp3, cell_vector_irr[i].comp1, cell_vector_irr[i].comp2, cell_vector_irr[i].comp3);
               }
              fprintf(file.out, "\n");
             }

  // ******************************************************************************************
  // * Generate equivalent ions (CSPM) and map ion positions to zeroth cell (CSP)             *
  // ******************************************************************************************

                  atoms->number_of_atoms_in_unit_cell = 0;
                  for (atm = 0; atm < atoms->number_of_unique_atoms; atm++) {
           	    for (k = 0; k < symmetry->number_of_operators; k++) {
                      p_irr = symmetry->irr + k * 9;
                      rotate_vector3(p_irr, &cell_vector_irr[atm], &Rvec_tmp[1]);
                      Rvec_tmp[1].comp1 += symmetry->taur[k].comp1;
                      Rvec_tmp[1].comp2 += symmetry->taur[k].comp2;
                      Rvec_tmp[1].comp3 += symmetry->taur[k].comp3;
                      map_to_wigner(crystal,&Rvec_tmp[1], &cell_vec1, &O_tmp);
                      for (i = 0; i <= atoms->number_of_atoms_in_unit_cell; i++) {
	                if (check_vec(&cell_vec1, &cell_vector_tmp[i]) == 1)
                          break;
	                  if (check_vec(&cell_vec1, &cell_vector_tmp[i]) == 0 && i == atoms->number_of_atoms_in_unit_cell) {
             	            cell_vector_tmp[i].comp1 = cell_vec1.comp1;
	                    cell_vector_tmp[i].comp2 = cell_vec1.comp2;
             	            cell_vector_tmp[i].comp3 = cell_vec1.comp3;
                            atoms->number_of_atoms_in_unit_cell++;
                            break;
                           }
                          } // close loop over i
                         } // close loop over k
                        } // close loop over atm

                       break;

           case 'M':

                  for (atm = 0; atm < atoms->number_of_unique_atoms; atm++) {
           	       for (k = 0; k < symmetry->number_of_operators; k++) {
                         p_irr = symmetry->irr + k * 9;
                         rotate_vector3(p_irr, &iratom[atm].posn, &cell_vec1);
                         for (i = 0; i <= atoms->number_of_atoms_in_unit_cell; i++) {
	                   if (check_vec(&cell_vec1, &cell_vector_tmp[i]) == 1)
                           break;
	                   if (check_vec(&cell_vec1, &cell_vector_tmp[i]) == 0 && i == atoms->number_of_atoms_in_unit_cell) {
             	             cell_vector_tmp[i].comp1 = cell_vec1.comp1;
	                     cell_vector_tmp[i].comp2 = cell_vec1.comp2;
             	             cell_vector_tmp[i].comp3 = cell_vec1.comp3;
                             atoms->number_of_atoms_in_unit_cell++;
                             break;
                            }
                           } // close loop over i
                          } // close loop over k
                         } // close loop over atm
        
                       break;

        } // end switch

}

void generate_all_atoms(CRYSTAL *crystal, REAL_LATTICE *R, SYMMETRY *symmetry, ATOMTYPE *iratom, ATOMTYPE *basis, ATOM *atoms, ATOM_TRAN *atom_p, ATOM_TRAN *atom_i, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Generate full set of atomic positions and their transformations under symmops          *
  // ******************************************************************************************

  int i, j, k;
  int atm, bas, count, count1[atoms->number_of_atoms_in_unit_cell];
  double *p_irr;
  VECTOR_DOUBLE cell_vector_irr[atoms->number_of_unique_atoms];
  VECTOR_DOUBLE Rvec_tmp[2], O_tmp, cell_vec1;

  for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++) {
    atoms->cell_vector[i].comp1 = 99999.9;
    atoms->cell_vector[i].comp2 = 99999.9;
    atoms->cell_vector[i].comp3 = 99999.9;
    atoms->basis_set[i] = -1;
    count1[i] = 0;
  }

  // ******************************************************************************************
  // * Zero elements of ATOM structure                                                        *
  // ******************************************************************************************

  for (i=0;i<atoms->number_of_atoms_in_unit_cell;i++) {
  atoms->cell_vector[i].comp1=99999.9 ; atoms->cell_vector[i].comp2=99999.9 ; atoms->cell_vector[i].comp3=99999.9 ; }
  for (i = 0; i < atoms->number_of_atoms_in_unit_cell * symmetry->number_of_operators; i++) {
  atom_p->K[i] = -1;
  atom_p->O[i] = -1;
  atom_i->K[i] = -1;
  atom_i->O[i] = -1;
 }

    atoms->number_of_atoms_in_unit_cell = 0;
    atoms->number_of_shells_in_unit_cell = 0 ;
    atoms->number_of_gaussians_in_unit_cell = 0 ;
    atoms->number_of_sh_shells_in_unit_cell = 0 ;
    atoms->number_of_sh_gaussians_in_unit_cell = 0 ;
    atoms->number_of_electrons_in_unit_cell = 0 ;

    for (i = 0; i < 2 * atoms->number_of_atoms_in_unit_cell; i++) 
    atoms->pop[i] = k_zero;

    switch (crystal->type[0]) {

      case 'C':
      case 'S':
      case 'P':

  // ******************************************************************************************
  // * Transformation of CSP atom positions under symmops ATOMS->cell_vector array            *
  // ******************************************************************************************

    for (i = 0; i < atoms->number_of_unique_atoms; i++) 
      map_to_wigner(crystal,&iratom[i].posn, &cell_vector_irr[i], &Rvec_tmp[0]);

     for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++) count1[i] = 0 ;
       for (atm = 0; atm < atoms->number_of_unique_atoms; atm++) {
         for (k = 0; k < symmetry->number_of_operators; k++) {
           p_irr = symmetry->irr + k * 9;
           rotate_vector3(p_irr, &cell_vector_irr[atm], &Rvec_tmp[1]);
           Rvec_tmp[1].comp1 += symmetry->taur[k].comp1;
           Rvec_tmp[1].comp2 += symmetry->taur[k].comp2;
           Rvec_tmp[1].comp3 += symmetry->taur[k].comp3;
           //fprintf(file.out,"Rvec_tmp[1]%3d%3d%10.4lf%10.4lf%10.4lf %10.4lf%10.4lf%10.4lf %10.4lf%10.4lf%10.4lf\n",\
           k,atoms->number_of_atoms_in_unit_cell,Rvec_tmp[1].comp1,Rvec_tmp[1].comp2,Rvec_tmp[1].comp3,\
           cell_vec1.comp1,cell_vec1.comp2,cell_vec1.comp3,Rvec_tmp[0].comp1,Rvec_tmp[0].comp2,Rvec_tmp[0].comp3) ; fflush(file.out);
           map_to_wigner(crystal,&Rvec_tmp[1], &cell_vec1, &O_tmp);
           for (i = 0; i <= atoms->number_of_atoms_in_unit_cell; i++) {
           //fprintf(file.out,"%3d %3d\n",check_vec(&cell_vec1,&atoms->cell_vector[i]),i) ; fflush(file.out);
           if (check_vec(&cell_vec1, &atoms->cell_vector[i]) == 1)
           break;
           if (check_vec(&cell_vec1, &atoms->cell_vector[i]) == 0 && i == atoms->number_of_atoms_in_unit_cell) {
           atoms->cell_vector[i].comp1 = cell_vec1.comp1;
           atoms->cell_vector[i].comp2 = cell_vec1.comp2;
           atoms->cell_vector[i].comp3 = cell_vec1.comp3;
           atoms->uniq[i] = atm;
           count1[atm]++;
           for (bas = 0; bas < atoms->number_of_basis_sets; bas++) {
             //fprintf(file.out,"at %d %d %d %d %d\n",\
             atoms->number_of_basis_sets,atm,iratom[atm].atomic_number,bas,basis[bas].atomic_number);
             if (iratom[atm].atomic_number == basis[bas].atomic_number % 300) {
             //fprintf(file.out,"At %d %d\n",iratom[atm].atomic_number,basis[bas].atomic_number);
             atoms->atomic_number[i] = iratom[atm].atomic_number % 200;
             atoms->number_of_shells_in_unit_cell    += basis[bas].number_of_shells;
             atoms->number_of_sh_shells_in_unit_cell += basis[bas].number_of_shells;
             for (j=0;j<basis[bas].number_of_shells;j++) {
               atoms->number_of_gaussians_in_unit_cell += basis[bas].shell[j].number_of_Gaussians;
               atoms->number_of_sh_gaussians_in_unit_cell += basis[bas].shell[j].number_of_Gaussians;
               if (basis[bas].shell[j].shell_type == 1) {
               atoms->number_of_sh_gaussians_in_unit_cell += basis[bas].shell[j].number_of_Gaussians;
               atoms->number_of_sh_shells_in_unit_cell += 1;
              }
             }
            }
           }
          atoms->number_of_atoms_in_unit_cell++;
          break;
         }
        } // close loop over i
       } // close loop over k
//        if (job->taskid == 0)
//        fprintf(file.out,"|                                                                                                         |\n");
       } // close loop over atm
        if (job->taskid == 0 && job->verbosity > 1)
        fprintf(file.out,"number of atoms %d number of shells %d number of sh shells %d number of gaussians %d sh gaussians %d\n",\
        atoms->number_of_atoms_in_unit_cell,atoms->number_of_shells_in_unit_cell, atoms->number_of_sh_shells_in_unit_cell, \
        atoms->number_of_gaussians_in_unit_cell, atoms->number_of_sh_gaussians_in_unit_cell);
        //fprintf(file.out,"%5d ION POSITIONS IN FRACTIONAL COORDINATES\n",number_of_atoms_in_unit_cell) ;
        //for (atm=0;atm<number_of_unique_atoms;atm++) {
        //fprintf(file.out,"%3d %20.15lf %20.15lf %20.15lf\n",atomirr[i].atomic_number,fractional_coord[i].comp1,\
        //fractional_coord[i].comp2,fractional_coord[i].comp3) ; }
        //fprintf(file.out, "\n") ;
        //}

  // ******************************************************************************************
  // * Transformation of CSP atom positions under symmops ATOMS->atom_p array                 *
  // ******************************************************************************************

    count = 0;
    if (job->taskid == 0 && job->verbosity > 1)
    fprintf(file.out, "TRANSFORMATION OF CSP ATOMIC POSITIONS BY SYMMETRY OPERATIONS\n");
    for (atm = 0; atm < atoms->number_of_unique_atoms; atm++) {
      for (i = count; i < count + count1[atm]; i++) {
        for (k = 0; k < symmetry->number_of_operators; k++) {
          p_irr = symmetry->irr + k * 9;
          rotate_vector3(p_irr, &atoms->cell_vector[i], &Rvec_tmp[1]);
          Rvec_tmp[1].comp1 += symmetry->taur[k].comp1;
          Rvec_tmp[1].comp2 += symmetry->taur[k].comp2;
          Rvec_tmp[1].comp3 += symmetry->taur[k].comp3;
          map_to_wigner(crystal,&Rvec_tmp[1], &cell_vec1, &O_tmp);
          //fprintf(file.out,"Rvec %3d %3d %10.4lf %10.4lf %10.4lf  %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf\n",\
          i,k,Rvec_tmp[1].comp1,Rvec_tmp[1].comp2,Rvec_tmp[1].comp3,cell_vec1.comp1,cell_vec1.comp2,cell_vec1.comp3,\
          O_tmp.comp1,O_tmp.comp2,O_tmp.comp3) ;
          for (j = 0; j < R->max_vector; j++) {
            if (check_vec(&O_tmp, &R->vec_ai[j]) == 1) {
            atom_p->O[i * symmetry->number_of_operators + k] = j;
            break;
           }
          }
         for (j = 0; j < atoms->number_of_atoms_in_unit_cell; j++) {
           if (check_vec(&cell_vec1, &atoms->cell_vector[j]) == 1) {
           atom_p->K[i * symmetry->number_of_operators + k] = j;
           atom_p->P[i * symmetry->number_of_operators + k] = 1; // change in site
           if (j == i)
           atom_p->P[i * symmetry->number_of_operators + k] = 0; // no change in site
          }
         } // close loop over j
        } // close loop over k
       } // close loop over i
      atom_p->numb[atm] = count1[atm];
      atom_p->posn[atm] = count; // CHANGES2016
      if (job->taskid == 0 && job->verbosity > 1) {
      for (i = count; i < count + count1[atm]; i++) {
        for (k = 0; k < symmetry->number_of_operators; k++) {
          fprintf(file.out, "atom_p[%2d][%2d] elem O P numb posn %3d %3d %3d %3d %3d \n", \
          i, k, atom_p->K[i * symmetry->number_of_operators + k], atom_p->O[i * symmetry->number_of_operators + k], \
          atom_p->P[i * symmetry->number_of_operators + k],atom_p->numb[atm],atom_p->posn[atm]);
         }
        }
       }
      for (i = count; i < count + count1[atm]; i++) {
        for (k = 0; k < symmetry->number_of_operators; k++) {
          if (atom_p->K[i * symmetry->number_of_operators + k] == -1 || atom_p->O[i * symmetry->number_of_operators + k] == -1) {
          if (job->taskid == 0) fprintf(file.out,"Transformation table atoms->atom_p incorrectly initialised\n");
          MPI_Finalize();
          exit(1);
         }
        }
       }
       count += count1[atm];
      } // close loop over atm

  // ******************************************************************************************
  // * Transformation of CSP atom positions under symmop inverses ATOMS->atom_i array         *
  // ******************************************************************************************

    count = 0;
    if (job->taskid == 0 && job->verbosity > 1)
    fprintf(file.out, "TRANSFORMATION OF CSP ATOMIC POSITIONS BY SYMMETRY OPERATIONS\n");
    for (atm = 0; atm < atoms->number_of_unique_atoms; atm++) {
      for (i = count; i < count + count1[atm]; i++) {
        for (k = 0; k < symmetry->number_of_operators; k++) {
          p_irr = symmetry->irr + symmetry->inverse[k] * 9;
          Rvec_tmp[0].comp1 = atoms->cell_vector[i].comp1 - symmetry->taur[k].comp1;
          Rvec_tmp[0].comp2 = atoms->cell_vector[i].comp2 - symmetry->taur[k].comp2;
          Rvec_tmp[0].comp3 = atoms->cell_vector[i].comp3 - symmetry->taur[k].comp3;
          rotate_vector3(p_irr, &Rvec_tmp[0], &Rvec_tmp[1]);
          //fprintf(file.out,"Rvec_tmp %3d %3d %3d %10.4lf %10.4lf %10.4lf\n",atm,i,k, \
          Rvec_tmp[1].comp1,Rvec_tmp[1].comp2,Rvec_tmp[1].comp3) ;
          map_to_wigner(crystal,&Rvec_tmp[1], &cell_vec1, &O_tmp);
          //fprintf(file.out,"tmp %3d %3d %3d %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf\n",atm,i,k,\
          cell_vec1.comp1,cell_vec1.comp2,cell_vec1.comp3,O_tmp.comp1,O_tmp.comp2,O_tmp.comp3) ;
          for (j = 0; j < R->max_vector; j++) {
            if (check_vec(&O_tmp, &R->vec_ai[j]) == 1) {
              atom_i->O[i * symmetry->number_of_operators + k] = j;
              break;
             }
            }
          for (j = 0; j < atoms->number_of_atoms_in_unit_cell; j++) {
            if (check_vec(&cell_vec1, &atoms->cell_vector[j]) == 1) {
            atom_i->K[i * symmetry->number_of_operators + k] = j;
            atom_i->P[i * symmetry->number_of_operators + k] = 1; // change in site
            if (j == i)
            atom_i->P[i * symmetry->number_of_operators + k] = 0; // no change in site
           } // no change in site
          } // close loop over j
         } // close loop over k
        } // close loop over i
        atom_i->numb[atm] = count1[atm];
        if (job->verbosity > 1) {
        for (i = count; i < count + count1[atm]; i++) {
          for (k = 0; k < symmetry->number_of_operators; k++) {
            //fprintf(file.out, "atom_i[%2d][%2d] elem O P numb %3d %3d %3d %3d \n", i, k, atom_i[i].elem[k],\
            atom_i[i].O[k], atom_i[i].P[k], atom_i[i].numb);
            if (job->taskid == 0 && job->verbosity > 1)
            fprintf(file.out,"atom_i[%2d][%2d] elem O P numb %3d %3d %3d %3d \n",i,k,\
            atom_i->K[i * symmetry->number_of_operators + k],atom_i->O[i * symmetry->number_of_operators + k], \
            atom_i->P[i * symmetry->number_of_operators + k],atom_i->numb[atm]);
            }
           }
          if (job->taskid == 0 && job->verbosity > 1)
          fprintf(file.out, "\n");
         }
         for (i = count; i < count + count1[atm]; i++) {
           for (k = 0; k < symmetry->number_of_operators; k++) {
             if (atom_i->K[i * symmetry->number_of_operators + k] == -1 || atom_i->O[i * symmetry->number_of_operators + k] == -1) {
             if (job->taskid == 0) fprintf(file.out,"Transformation table atoms->atom_i incorrectly initialised\n");
             MPI_Finalize();
             exit(1);
            }
           }
          }
         count += count1[atm];
        } // close loop over atm

           break;

           case 'M':

  // ******************************************************************************************
  // * Transformation of M atom positions under symmops ATOMS->cell_vector array              *
  // ******************************************************************************************

    for (i=0;i<atoms->number_of_atoms_in_unit_cell;i++) count1[i] = 0 ;
      for (atm = 0; atm < atoms->number_of_unique_atoms; atm++) {
        for (k = 0; k < symmetry->number_of_operators; k++) {
          p_irr = symmetry->irr + k * 9;
          rotate_vector3(p_irr, &iratom[atm].posn, &cell_vec1);
          //fprintf(file.out,"cell_vector_irr %3d %3d %3d %10.4lf %10.4lf %10.4lf\n",atm,i,k,cell_vec1.comp1,cell_vec1.comp2,cell_vec1.comp3) ;
          for (i = 0; i <= atoms->number_of_atoms_in_unit_cell; i++) {
            //fprintf(file.out,"%3d %3d\n",check_vec(&cell_vec1,&atoms->cell_vector[i]),i) ;
            if (check_vec(&cell_vec1, &atoms->cell_vector[i]) == 1)
            break;
            if (check_vec(&cell_vec1, &atoms->cell_vector[i]) == 0 && i == atoms->number_of_atoms_in_unit_cell) {
            atoms->cell_vector[i].comp1 = cell_vec1.comp1;
            atoms->cell_vector[i].comp2 = cell_vec1.comp2;
            atoms->cell_vector[i].comp3 = cell_vec1.comp3;
            atoms->uniq[i] = atm;
            //fprintf(file.out,"%3d %3d %3d\n",i,atm,atoms->uniq[i]);
            count1[atm]++;
            for (bas=0;bas<atoms->number_of_basis_sets;bas++) {
              //fprintf(file.out,"at %d %d %d %d %d\n",atoms->number_of_basis_sets,atm,iratom[atm].atomic_number,bas,basis[bas].atomic_number);
              if (iratom[atm].atomic_number == basis[bas].atomic_number % 300) {
              atoms->atomic_number[i] = iratom[atm].atomic_number;
              atoms->number_of_shells_in_unit_cell    += basis[bas].number_of_shells;
              atoms->number_of_sh_shells_in_unit_cell += basis[bas].number_of_shells;
              //fprintf(file.out,"at no %d %d %d\n", atoms->atomic_number[i],atoms->number_of_shells_in_unit_cell,atoms->number_of_sh_shells_in_unit_cell);
              for (j=0;j<basis[bas].number_of_shells;j++) {
                atoms->number_of_gaussians_in_unit_cell += basis[bas].shell[j].number_of_Gaussians;
                atoms->number_of_sh_gaussians_in_unit_cell += basis[bas].shell[j].number_of_Gaussians;
                if (basis[bas].shell[j].shell_type == 1) {
	        atoms->number_of_sh_gaussians_in_unit_cell += basis[bas].shell[j].number_of_Gaussians;
                atoms->number_of_sh_shells_in_unit_cell += 1;
               }
             //fprintf(file.out,"gn %d %d %d\n", atoms->atomic_number[i],atoms->number_of_sh_gaussians_in_unit_cell,atoms->number_of_gaussians_in_unit_cell);
              }
             }
            }
           atoms->number_of_atoms_in_unit_cell++;
           break;
          }
         } // close loop over i
        } // close loop over k
       } // close loop over atm

  // ******************************************************************************************
  // * Transformation of M atom positions under symmops ATOMS->atom_p array                   *
  // ******************************************************************************************

    count = 0;
    if (job->taskid == 0 && job->verbosity > 1)
    fprintf(file.out, "TRANSFORMATION OF M ATOMIC POSITIONS BY SYMMETRY OPERATIONS\n");
    for (atm = 0; atm < atoms->number_of_unique_atoms; atm++) {
      for (i = count; i < count + count1[atm]; i++) {
        for (k = 0; k < symmetry->number_of_operators; k++) {
          p_irr = symmetry->irr + k * 9;
          rotate_vector3(p_irr, &atoms->cell_vector[i], &cell_vec1);
          //fprintf(file.out,"%3d %3d %3d  %10.4lf %10.4lf %10.4lf   %10.4lf %10.4lf %10.4lf\n",atm,i,k,atoms->cell_vector[i].comp1,\
          atoms->cell_vector[i].comp2,atoms->cell_vector[i].comp3,cell_vec1.comp1,cell_vec1.comp2,cell_vec1.comp3); fflush(file.out);
          for (j = 0; j < atoms->number_of_atoms_in_unit_cell; j++) {
          if (check_vec(&cell_vec1, &atoms->cell_vector[j]) == 1) {
          atom_p->K[i * symmetry->number_of_operators + k] = j;
          atom_p->P[i * symmetry->number_of_operators + k] = 1; // change in site
          atom_p->O[i * symmetry->number_of_operators + k] = 0;
          if (j == i)
          atom_p->P[i * symmetry->number_of_operators + k] = 0; // no change in site
         }
        } // close loop over j
       } // close loop over k
      } // close loop over i

     atom_p->numb[atm] = count1[atm];
     atom_p->posn[atm] = count;
     if (job->verbosity > 1) {
     for (i = count; i < count + count1[atm]; i++) {
       for (k = 0; k < symmetry->number_of_operators; k++) {
         if (job->taskid == 0 && job->verbosity > 1)
         fprintf(file.out, "atom_p[%2d][%2d] elem P numb posn %3d %3d %3d %3d \n", i, k, \
         atom_p->K[i * symmetry->number_of_operators + k],atom_p->P[i * symmetry->number_of_operators + k], \
         atom_p->numb[atm],atom_p->posn[atm]);
        }
       }
      if (job->taskid == 0 && job->verbosity > 1)
      fprintf(file.out, "\n");
     }
     for (i = count; i < count + count1[atm]; i++) {
       for (k = 0; k < symmetry->number_of_operators; k++) {
         if (atom_p->K[i * symmetry->number_of_operators + k] == -1 || atom_p->O[i * symmetry->number_of_operators + k] == -1) {
         if (job->taskid == 0) fprintf(file.out,"Transformation table atoms->atom_p incorrectly initialised\n");
         MPI_Finalize();
         exit(1);
        }
       }
      }
     count += count1[atm];
    } // close loop over atm

  // ******************************************************************************************
  // * Transformation of M atom positions under inverses of symmops ATOMS->atom_i array       *
  // ******************************************************************************************

    count = 0;
    if (job->taskid == 0 && job->verbosity > 1)
    fprintf(file.out, "TRANSFORMATION OF M ATOMIC POSITIONS BY INVERSE SYMMETRY OPERATIONS\n");
    for (atm = 0; atm < atoms->number_of_unique_atoms; atm++) {
      for (i = count; i < count + count1[atm]; i++) {
        for (k = 0; k < symmetry->number_of_operators; k++) {
          p_irr = symmetry->irr + symmetry->inverse[k] * 9;
          rotate_vector3(p_irr, &atoms->cell_vector[i], &cell_vec1);
          for (j = 0; j < atoms->number_of_atoms_in_unit_cell; j++) {
          if (check_vec(&cell_vec1, &atoms->cell_vector[j]) == 1) {
          atom_i->K[i * symmetry->number_of_operators + k] = j;
          atom_i->P[i * symmetry->number_of_operators + k] = 1; // change in site
          atom_i->O[i * symmetry->number_of_operators + k] = 0;
          if (j == i)
          atom_i->P[i * symmetry->number_of_operators + k] = 0; // no change in site
         } // no change in site
        } // close loop over j
       } // close loop over k
      } // close loop over i
      atom_i->numb[atm] = count1[atm];
      if (job->verbosity > 1) {
        for (i = count; i < count + count1[atm]; i++) {
          for (k = 0; k < symmetry->number_of_operators; k++) {
            if (job->taskid == 0 && job->verbosity > 1)
            fprintf(file.out, "atom_i[%2d][%2d] elem P numb %3d %3d %3d \n", i, k, atom_i->K[i * symmetry->number_of_operators + k],
            atom_i->P[i * symmetry->number_of_operators + k], atom_i->numb[atm]);
           }
          }
         if (job->taskid == 0 && job->verbosity > 1)
         fprintf(file.out, "\n");
        }
         for (i = count; i < count + count1[atm]; i++) {
           for (k = 0; k < symmetry->number_of_operators; k++) {
             if (atom_i->K[i * symmetry->number_of_operators + k] == -1 || atom_i->O[i * symmetry->number_of_operators + k] == -1) {
             if (job->taskid == 0) fprintf(file.out,"Transformation table atoms->atom_i incorrectly initialised\n");
             MPI_Finalize();
             exit(1);
            }
           }
          }
         count += count1[atm];
        } // close loop over atm

       break;

         } // close switch

      if (job->taskid == 0) {

      fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
      fprintf(file.out,"|  ATOM   |  ATOMIC NUMBER  |                                 COORDINATES                                 |\n");
      fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
      for (i=0;i<atoms->number_of_atoms_in_unit_cell;i++) {
      fprintf(file.out,"|    %3d  |     %3d         |     %20.15lf  %20.15lf  %20.15lf        |\n", i + 1, \
      atoms->atomic_number[i],atoms->cell_vector[i].comp1 * bohr_to_AA, atoms->cell_vector[i].comp2 * bohr_to_AA, atoms->cell_vector[i].comp3 * bohr_to_AA);
     }
      fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n\n");

    }

}

void generate_ATOM_TRAN(CRYSTAL *crystal, REAL_LATTICE *R, SYMMETRY *symmetry, ATOM *atoms, ATOM_TRAN *atom_p, JOB_PARAM *job, FILES file)

{

int i, j, k;
int atm, count;
double *p_irr;
VECTOR_DOUBLE Rvec_tmp[2], O_tmp, cell_vec1;

    for (i = 0; i < atoms->number_of_unique_atoms; i++) 
      atom_p->numb[i] = 0; 
    for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++) { 
      atm = atoms->uniq[i];
      atom_p->numb[atm]++; 
     }

    switch (crystal->type[0]) {

      case 'C':
      case 'S':
      case 'P':

  // ******************************************************************************************
  // * Transformation of CSP atom positions under symmops ATOMS->atom_p array                 *
  // ******************************************************************************************

    count = 0;
    if (job->taskid == 0 && job->verbosity > 1)
    fprintf(file.out, "TRANSFORMATION OF CSP ATOMIC POSITIONS BY SYMMETRY OPERATIONS\n");
    for (atm = 0; atm < atoms->number_of_unique_atoms; atm++) {
      //fprintf(file.out,"%3d %3d %3d\n",atm,atom_p->numb[atm],symmetry->number_of_operators);
      for (i = count; i < count + atom_p->numb[atm]; i++) {
      //for (i = count; i < count + count1[atm]; i++) {
        for (k = 0; k < symmetry->number_of_operators; k++) {
          p_irr = symmetry->irr + k * 9;
          rotate_vector3(p_irr, &atoms->cell_vector[i], &Rvec_tmp[1]);
          Rvec_tmp[1].comp1 += symmetry->taur[k].comp1;
          Rvec_tmp[1].comp2 += symmetry->taur[k].comp2;
          Rvec_tmp[1].comp3 += symmetry->taur[k].comp3;
          map_to_wigner(crystal,&Rvec_tmp[1], &cell_vec1, &O_tmp);
          //fprintf(file.out,"Rvec %3d %3d %10.4lf %10.4lf %10.4lf   %10.4lf %10.4lf %10.4lf   %10.4lf %10.4lf %10.4lf\n", \
          i,k,Rvec_tmp[1].comp1, Rvec_tmp[1].comp2,Rvec_tmp[1].comp3,cell_vec1.comp1,cell_vec1.comp2,cell_vec1.comp3, \
          O_tmp.comp1,O_tmp.comp2,O_tmp.comp3) ;
          fflush(file.out);
          for (j = 0; j < R->max_vector; j++) {
            if (check_vec(&O_tmp, &R->vec_ai[j]) == 1) {
            atom_p->O[i * symmetry->number_of_operators + k] = j;
            break;
           }
          }
         for (j = 0; j < atoms->number_of_atoms_in_unit_cell; j++) {
           if (check_vec(&cell_vec1, &atoms->cell_vector[j]) == 1) {
           atom_p->K[i * symmetry->number_of_operators + k] = j;
           atom_p->P[i * symmetry->number_of_operators + k] = 1; // change in site
           if (j == i)
           atom_p->P[i * symmetry->number_of_operators + k] = 0; // no change in site
          }
         } // close loop over j
        } // close loop over k
       } // close loop over i
      ////atom_p->numb[atm] = count1[atm];
      atom_p->posn[atm] = count; // CHANGES2016
      if (job->taskid == 0 && job->verbosity > 1) {
      ////for (i = count; i < count + count1[atm]; i++) {
      for (i = count; i < count + atom_p->numb[atm]; i++) {
        for (k = 0; k < symmetry->number_of_operators; k++) {
          fprintf(file.out, "atom_p[%2d][%2d] elem O P numb posn %3d %3d %3d %3d %3d \n", i, k, \
          atom_p->K[i * symmetry->number_of_operators + k], atom_p->O[i * symmetry->number_of_operators + k], \
          atom_p->P[i * symmetry->number_of_operators + k],atom_p->numb[atm],atom_p->posn[atm]);
         }
        }
       }
      ////for (i = count; i < count + count1[atm]; i++) {
      for (i = count; i < count + atom_p->numb[atm]; i++) {
        for (k = 0; k < symmetry->number_of_operators; k++) {
          if (atom_p->K[i * symmetry->number_of_operators + k] == -1 || atom_p->O[i * symmetry->number_of_operators + k] == -1) {
          if (job->taskid == 0) fprintf(file.out,"Transformation table atoms->atom_p incorrectly initialised\n");
          MPI_Finalize();
          exit(1);
         }
        }
       }
       //count += count1[atm];
       count += atom_p->numb[atm];
      } // close loop over atm

           break;

           case 'M':

  // ******************************************************************************************
  // * Transformation of M atom positions under symmops ATOMS->atom_p array                   *
  // ******************************************************************************************

    count = 0;
    if (job->taskid == 0 && job->verbosity > 1)
    fprintf(file.out, "TRANSFORMATION OF M ATOMIC POSITIONS BY SYMMETRY OPERATIONS\n");
    for (atm = 0; atm < atoms->number_of_unique_atoms; atm++) {
      for (i = count; i < count + atom_p->numb[atm]; i++) {
        for (k = 0; k < symmetry->number_of_operators; k++) {
          p_irr = symmetry->irr + k * 9;
          rotate_vector3(p_irr, &atoms->cell_vector[i], &cell_vec1);
          //fprintf(file.out,"%3d %3d %3d  %10.4lf %10.4lf %10.4lf  %10.4lf %10.4lf %10.4lf\n",atm,i,k,atoms->cell_vector[i].comp1, \
          atoms->cell_vector[i].comp2,atoms->cell_vector[i].comp3,cell_vec1.comp1,cell_vec1.comp2,cell_vec1.comp3); fflush(file.out);
          for (j = 0; j < atoms->number_of_atoms_in_unit_cell; j++) {
          if (check_vec(&cell_vec1, &atoms->cell_vector[j]) == 1) {
          atom_p->K[i * symmetry->number_of_operators + k] = j;
          atom_p->P[i * symmetry->number_of_operators + k] = 1; // change in site
          atom_p->O[i * symmetry->number_of_operators + k] = 0;
          if (j == i)
          atom_p->P[i * symmetry->number_of_operators + k] = 0; // no change in site
         }
        } // close loop over j
       } // close loop over k
      } // close loop over i

     atom_p->numb[atm] = atom_p->numb[atm];
     atom_p->posn[atm] = count;
     if (job->verbosity > 1) {
     for (i = count; i < count + atom_p->numb[atm]; i++) {
       for (k = 0; k < symmetry->number_of_operators; k++) {
         if (job->taskid == 0 && job->verbosity > 1)
          fprintf(file.out, "atom_p[%2d][%2d] elem O P numb posn %3d %3d %3d %3d %3d \n", i, k, \
          atom_p->K[i * symmetry->number_of_operators + k], atom_p->O[i * symmetry->number_of_operators + k], \
          atom_p->P[i * symmetry->number_of_operators + k],atom_p->numb[atm],atom_p->posn[atm]);
        }
       }
      if (job->taskid == 0 && job->verbosity > 1)
      fprintf(file.out, "\n");
     }
     for (i = count; i < count + atom_p->numb[atm]; i++) {
       for (k = 0; k < symmetry->number_of_operators; k++) {
         if (atom_p->K[i * symmetry->number_of_operators + k] == -1 || atom_p->O[i * symmetry->number_of_operators + k] == -1) {
         if (job->taskid == 0) fprintf(file.out,"Transformation table atoms->atom_p incorrectly initialised\n");
         MPI_Finalize();
         exit(1);
        }
       }
      }
     count += atom_p->numb[atm];
    } // close loop over atm

       break;

         } // close switch

}

void count_all_crystal_atoms(ATOM *atoms, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Count all atoms in cell in Crystal09 OUTPUT file                                       *
  // ******************************************************************************************

  FILE *output;
  char *title;

      output = fopen("OUTPUT", "r");
      if (output == NULL) {
      if (job->taskid == 0)
      fprintf(file.out, "Cannot open Crystal09 file OUTPUT for reading in count_all_crystal_atoms");
      MPI_Finalize();
      exit(1);
     }
      for (;;) {
      	if ((search_file(output, "UNIT")) != TRUE) {
        if (job->taskid == 0)
        fprintf(file.out, "ERROR: Keyword UNIT not found in Crystal09 file OUTPUT\n");
        break;
       }
	if (fscanf(output, "%5s", title) != 1) {
	fprintf(file.out,"String not found\n");
	MPI_Finalize();
	exit(1);
	}
	//fscanf(output, "%5s", &title);
       	if (!strcmp(title, "CELL:")) {
        break;
       }
      }
      if (fscanf(output, "%d\n", &atoms->number_of_atoms_in_unit_cell) != 1) {
	fprintf(file.out,"String not found\n");
	MPI_Finalize();
	exit(1);
	}
 
      if (job->taskid == 0 && job->verbosity > 1)
      fprintf(file.out, "Count_all_crystal_atoms: Number of atoms in unit cell = %4d %4d\n", \
      atoms->number_of_atoms_in_unit_cell, atoms->number_of_unique_atoms) ;

      fclose(output);

}

void generate_all_crystal_atoms(CRYSTAL *crystal, REAL_LATTICE *R, SYMMETRY *symmetry, ATOMTYPE *basis, ATOM *atoms, ATOM_TRAN *atom_p, ATOM_TRAN *atom_i, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // * Generate all atoms in cell in Crystal09 OUTPUT file                                    *
  // ******************************************************************************************

  int i, j, k;
  int atm, bas, count, count1[atoms->number_of_atoms_in_unit_cell];
  double *p_irr;
  FILE *output;
  char *title;
  char TF[8];
  char *line;
  VECTOR_DOUBLE Rvec_tmp[2], O_tmp, cell_vec1;
  VECTOR_DOUBLE O_tmp_crys[atoms->number_of_atoms_in_unit_cell], cell_vec2;

  // ******************************************************************************************
  // * Zero elements of ATOM structure                                                        *
  // ******************************************************************************************

  for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++) 
  count1[i] = 0;
  atoms->number_of_shells_in_unit_cell = 0 ;
  atoms->number_of_gaussians_in_unit_cell = 0 ;
  atoms->number_of_sh_shells_in_unit_cell = 0 ;
  atoms->number_of_sh_gaussians_in_unit_cell = 0 ;
  atoms->number_of_electrons_in_unit_cell = 0 ;

 for (i = 0; i < 2 * atoms->number_of_atoms_in_unit_cell; i++) 
 atoms->pop[i] = k_zero;
    
      output = fopen("OUTPUT", "r");
      if (output == NULL) {
      if (job->taskid == 0)
      fprintf(file.out, "Cannot open Crystal09 file OUTPUT for reading in generate_all_crystal_atoms");
      MPI_Finalize();
      exit(1);
     }
      for (;;) {
      	if ((search_file(output, "UNIT")) != TRUE) {
        if (job->taskid == 0)
        fprintf(file.out, "ERROR: Keyword UNIT not found in Crystal09 file OUTPUT\n");
        break;
       }
	if (fscanf(output, "%5s", title) != 1) {
	fprintf(file.out,"String not found\n");
	MPI_Finalize();
	exit(1);
	}
	//fscanf(output, "%5s", &title);
       	if (!strcmp(title, "CELL:")) {
        break;
       }
      }
      if (fscanf(output, "%d\n", &atoms->number_of_atoms_in_unit_cell) != 1) {
	fprintf(file.out,"String not found\n");
	MPI_Finalize();
	exit(1);
	}
      //fprintf(file.out, "Number of atoms in unit cell = %4d %4d\n", atoms->number_of_atoms_in_unit_cell, atoms->number_of_unique_atoms) ;
      read_line(output, line, 149);
      read_line(output, line, 149);    

  // ******************************************************************************************
  // * Count all CSPM atoms, shells and gaussians and assign atomic numbers                   *
  // ******************************************************************************************

      atm = -1;
      for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++) {
        read_line(output, line, 149);
        sscanf(line, "%d %s %d", &j, TF, &atoms->atomic_number[i]);
        atoms->atomic_number[i] = atoms->atomic_number[i] % 200;
  	if (TF[0] == 'T') atm++;
        //fprintf(file.out,"%d %d atomic number %d\n",atm,i,atoms->atomic_number[i]); fflush(file.out);
        count1[atm]++;
        for (bas=0;bas<atoms->number_of_basis_sets;bas++) {
          if (atoms->atomic_number[i] == basis[bas].atomic_number % 300) {
          atoms->number_of_shells_in_unit_cell    += basis[bas].number_of_shells;
          atoms->number_of_sh_shells_in_unit_cell += basis[bas].number_of_shells;
          //fprintf(file.out,"Number of shells %d %d on atom %d %d %d\n",basis[bas].number_of_shells,iratom[atm].atomic_number,i,atm,bas);
          for (j=0;j<basis[bas].number_of_shells;j++) {
            atoms->number_of_gaussians_in_unit_cell += basis[bas].shell[j].number_of_Gaussians;
            atoms->number_of_sh_gaussians_in_unit_cell += basis[bas].shell[j].number_of_Gaussians;
            if (basis[bas].shell[j].shell_type == 1) {
            atoms->number_of_sh_gaussians_in_unit_cell += basis[bas].shell[j].number_of_Gaussians;
            atoms->number_of_sh_shells_in_unit_cell += 1;
           }
          }
         }
        }
       }

       if (job->taskid == 0 && job->verbosity > 1)
       fprintf(file.out,"number of atoms %d number of shells %d number of sh shells %d number of gaussians %d number of sh gaussians %d\n",\
        atoms->number_of_atoms_in_unit_cell,atoms->number_of_shells_in_unit_cell, atoms->number_of_sh_shells_in_unit_cell, \
        atoms->number_of_gaussians_in_unit_cell, atoms->number_of_sh_gaussians_in_unit_cell);

  // ******************************************************************************************
  // * Search for absolute atomic coordinates in Crystal09 OUTPUT file                        *
  // ******************************************************************************************

      rewind(output);
      for (;;) {
        if ((search_file(output, "ATOM")) != TRUE) {
        fprintf(file.out, "ERROR: Keyword ATOM not found in Crystal09 file OUTPUT\n");
        MPI_Finalize();
        exit(1);
        break;
       }
        if (fscanf(output, "%12s", line) != 1) {
	fprintf(file.out,"String not found\n");
	MPI_Finalize();
	exit(1);
	}
        //Scan for line containing X(ANGSTROM) and break loop
        if (!strcmp(line, "X(ANGSTROM)")) {
        if (job->taskid == 0 && job->verbosity > 1)
        fprintf(file.out, " ");
	break;
       }
      }
       read_line(output, line, 149);
       read_line(output, line, 149);

  // ******************************************************************************************
  // * Read absolute atomic coordinates from Crystal09 OUTPUT file                            *
  // ******************************************************************************************

      for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++) {
        read_line(output, line, 149);
        sscanf(line, "%d %d %s %le %le %le", &j, &atoms->atomic_number[i], title,
        &atoms->cell_vector[i].comp1, &atoms->cell_vector[i].comp2,\
        &atoms->cell_vector[i].comp3);
        atoms->atomic_number[i] = atoms->atomic_number[i] % 200;
        //if (job->taskid == 0 job->verbosity > 1)
        //fprintf(file.out, "%5d %5d             %15.10lf %15.10lf %15.10lf\n",i + 1, atoms->atomic_number[i],\
        atoms->cell_vector[i].comp1, atoms->cell_vector[i].comp2, atoms->cell_vector[i].comp3);
        atoms->cell_vector[i].comp1 /= bohr_to_AA;
        atoms->cell_vector[i].comp2 /= bohr_to_AA;
        atoms->cell_vector[i].comp3 /= bohr_to_AA;
       }

      fclose(output);

      if (job->taskid == 0) {
      fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
      fprintf(file.out,"|  ATOM   |  ATOMIC NUMBER  |                                 COORDINATES                                 |\n");
      fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
      for (i=0;i<atoms->number_of_atoms_in_unit_cell;i++) {
      fprintf(file.out,"|    %3d  |     %3d         |     %20.15lf  %20.15lf  %20.15lf        |\n", i + 1, \
      atoms->atomic_number[i],atoms->cell_vector[i].comp1 * bohr_to_AA, atoms->cell_vector[i].comp2 * bohr_to_AA, atoms->cell_vector[i].comp3 * bohr_to_AA);
     }
      fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n\n");
    }

        for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++) {
        map_to_wigner(crystal,&atoms->cell_vector[i], &cell_vec1, &O_tmp_crys[i]);
        //fprintf(file.out,"avec %3d %3d %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf\n",i,k,\
        atoms->cell_vector[i].comp1,atoms->cell_vector[i].comp2,atoms->cell_vector[i].comp3,\
        cell_vec1.comp1,cell_vec1.comp2,cell_vec1.comp3,O_tmp_crys[i].comp1,O_tmp_crys[i].comp2,O_tmp_crys[i].comp3) ; 
       }

  // ******************************************************************************************
  // * Transformation of CSP atom positions under symmops ATOMS->atom_p array                 *
  // ******************************************************************************************

  count = 0;
  if (job->taskid == 0 && job->verbosity > 1)
  fprintf(file.out, "TRANSFORMATION OF ATOMIC POSITIONS BY SYMMETRY OPERATIONS\n");
  for (atm = 0; atm < atoms->number_of_unique_atoms; atm++) {
    for (i = count; i < count + count1[atm]; i++) {
      for (k = 0; k < symmetry->number_of_operators; k++) {
        p_irr = symmetry->irr + k * 9;
        rotate_vector3(p_irr, &atoms->cell_vector[i], &Rvec_tmp[1]);
        Rvec_tmp[1].comp1 += symmetry->taur[k].comp1;
        Rvec_tmp[1].comp2 += symmetry->taur[k].comp2;
        Rvec_tmp[1].comp3 += symmetry->taur[k].comp3;
        map_to_wigner(crystal,&Rvec_tmp[1], &cell_vec1, &O_tmp);
        //fprintf(file.out,"Rvec %3d %3d %10.4lf %10.4lf %10.4lf taur %10.4lf %10.4lf %10.4lf\n",i,k,Rvec_tmp[1].comp1,\
        Rvec_tmp[1].comp2,Rvec_tmp[1].comp3,symmetry->taur[k].comp1,symmetry->taur[k].comp2,symmetry->taur[k].comp3) ;
        //fprintf(file.out,"Cvec %3d %3d %10.4lf %10.4lf %10.4lf otmp %10.4lf %10.4lf %10.4lf\n",i,k,cell_vec1.comp1,\
        cell_vec1.comp2,cell_vec1.comp3,O_tmp.comp1,O_tmp.comp2,O_tmp.comp3) ;
        for (j = 0; j < atoms->number_of_atoms_in_unit_cell; j++) {
          cell_vec2.comp1 = atoms->cell_vector[j].comp1 + O_tmp_crys[j].comp1;
          cell_vec2.comp2 = atoms->cell_vector[j].comp2 + O_tmp_crys[j].comp2;
          cell_vec2.comp3 = atoms->cell_vector[j].comp3 + O_tmp_crys[j].comp3;
          //fprintf(file.out,"2vec %3d %3d %10.4lf %10.4lf %10.4lf atom %10.4lf %10.4lf %10.4lf\n",j,k,cell_vec2.comp1,\
          cell_vec2.comp2,cell_vec2.comp3,cell_vec1.comp1,cell_vec1.comp2,cell_vec1.comp3) ;
          if (check_vec(&cell_vec2, &cell_vec1) == 1) {
          //if (check_vec(&cell_vec1, &atoms->cell_vector[j]) == 1) {
            atom_p->K[i * symmetry->number_of_operators + k] = j;
            atom_p->P[i * symmetry->number_of_operators + k] = 1; // change in site
            if (j == i)
            atom_p->P[i * symmetry->number_of_operators + k] = 0; // no change in site
            O_tmp.comp1 -= O_tmp_crys[j].comp1;
            O_tmp.comp2 -= O_tmp_crys[j].comp2;
            O_tmp.comp3 -= O_tmp_crys[j].comp3;
            //fprintf(file.out,"Otmp %3d %3d                                       %10.4lf %10.4lf %10.4lf\n",j,k,\
            O_tmp.comp1,O_tmp.comp2,O_tmp.comp3) ;
            //fprintf(file.out, "\n");
            break;
           }
          } // close loop over j

        for (j = 0; j < R->max_vector; j++) {
          if (check_vec(&O_tmp, &R->vec_ai[j]) == 1) {
          atom_p->O[i * symmetry->number_of_operators + k] = j;
          break;
         } 
         if (j == R->max_vector && job->taskid == 0) fprintf(file.out,"Vector not found in generate_all_crystal_atoms\n");
        }
       } // close loop over k
      } // close loop over i

    atom_p->numb[atm] = count1[atm];
    if (job->taskid == 0 && job->verbosity > 1) {
    for (i = count; i < count + count1[atm]; i++) {
      for (k = 0; k < symmetry->number_of_operators; k++) {
        fprintf(file.out, "atom_p[%2d][%2d] elem O P numb %3d %3d %3d %3d \n", i, k, \
        atom_p->K[i * symmetry->number_of_operators + k],atom_p->O[i * symmetry->number_of_operators + k], atom_p->P[i * symmetry->number_of_operators + k],\
        atom_p->numb[atm]);
       }
      }
      if (job->taskid == 0 && job->verbosity > 1)
      fprintf(file.out, "\n");

      for (i = count; i < count + count1[atm]; i++) {
        for (k = 0; k < symmetry->number_of_operators; k++) {
          if (atom_p->K[i * symmetry->number_of_operators + k] == -1 || atom_p->O[i * symmetry->number_of_operators + k] == -1) {
          if (job->taskid == 0) fprintf(file.out,"Transformation table atoms->atom_p incorrectly initialised\n");
          MPI_Finalize();
          exit(1);
         }
        }
       }
      }
     count += count1[atm];
    } // close loop over atm

  // ******************************************************************************************
  // * Transformation of CSP atom positions under symmop inverses ATOMS->atom_i array         *
  // ******************************************************************************************

  count = 0;
  if (job->taskid == 0 && job->verbosity > 1)
    fprintf(file.out, "TRANSFORMATION OF ATOMIC POSITIONS BY SYMMETRY OPERATIONS\n");
  for (atm = 0; atm < atoms->number_of_unique_atoms; atm++) {
    for (i = count; i < count + count1[atm]; i++) {
      for (k = 0; k < symmetry->number_of_operators; k++) {
        p_irr = symmetry->irr + symmetry->inverse[k] * 9;
        Rvec_tmp[0].comp1 = atoms->cell_vector[i].comp1 - symmetry->taur[k].comp1;
        Rvec_tmp[0].comp2 = atoms->cell_vector[i].comp2 - symmetry->taur[k].comp2;
        Rvec_tmp[0].comp3 = atoms->cell_vector[i].comp3 - symmetry->taur[k].comp3;
        rotate_vector3(p_irr, &atoms->cell_vector[i], &Rvec_tmp[0]);
        //fprintf(file.out,"Rvec_tmp %3d %3d %3d %10.4lf %10.4lf %10.4lf\n",atm,i,k,Rvec_tmp[1].comp1,Rvec_tmp[1].comp2,Rvec_tmp[1].comp3);
        map_to_wigner(crystal,&Rvec_tmp[1], &cell_vec1, &O_tmp);
        //fprintf(file.out,"Rvec %3d %3d %10.4lf %10.4lf %10.4lf taur %10.4lf %10.4lf %10.4lf\n",i,k,Rvec_tmp[1].comp1,\
        Rvec_tmp[1].comp2,Rvec_tmp[1].comp3,symmetry->taur[k].comp1,symmetry->taur[k].comp2,symmetry->taur[k].comp3) ;
        //fprintf(file.out,"Cvec %3d %3d %10.4lf %10.4lf %10.4lf otmp %10.4lf %10.4lf %10.4lf\n",i,k,cell_vec1.comp1,\
        cell_vec1.comp2,cell_vec1.comp3,O_tmp.comp1,O_tmp.comp2,O_tmp.comp3) ;
        //fprintf(file.out,"tmp %3d %3d %3d %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf %10.4lf\n",atm,i,k,\
        cell_vec1.comp1,cell_vec1.comp2,cell_vec1.comp3,O_tmp.comp1,O_tmp.comp2,O_tmp.comp3) ;
        for (j = 0; j < atoms->number_of_atoms_in_unit_cell; j++) {
          cell_vec2.comp1 = atoms->cell_vector[j].comp1 + O_tmp_crys[j].comp1;
          cell_vec2.comp2 = atoms->cell_vector[j].comp2 + O_tmp_crys[j].comp2;
          cell_vec2.comp3 = atoms->cell_vector[j].comp3 + O_tmp_crys[j].comp3;
          //fprintf(file.out,"2vec %3d %3d %10.4lf %10.4lf %10.4lf atom %10.4lf %10.4lf %10.4lf\n",j,k,cell_vec2.comp1,\
          cell_vec2.comp2,cell_vec2.comp3,atoms->cell_vector[j].comp1,atoms->cell_vector[j].comp2,atoms->cell_vector[j].comp3) ;
          if (check_vec(&cell_vec2, &cell_vec1) == 1) {
          atom_i->K[i * symmetry->number_of_operators + k] = j;
          atom_i->P[i * symmetry->number_of_operators + k] = 1; // change in site
          if (j == i)
          atom_i->P[i * symmetry->number_of_operators + k] = 0; // no change in site
          O_tmp.comp1 -= O_tmp_crys[j].comp1;
          O_tmp.comp2 -= O_tmp_crys[j].comp2;
          O_tmp.comp3 -= O_tmp_crys[j].comp3;
          //fprintf(file.out,"Otmp %3d %3d                                       %10.4lf %10.4lf %10.4lf\n",j,k,\
          O_tmp.comp1,O_tmp.comp2,O_tmp.comp3) ;
          //fprintf(file.out, "\n");
          break;
         } // no change in site
        } // close loop over j

        for (j = 0; j < R->max_vector; j++) {
          if (check_vec(&O_tmp, &R->vec_ai[j]) == 1) {
          atom_i->O[i * symmetry->number_of_operators + k] = j;
          break;
         }
         if (j == R->max_vector && job->taskid == 0) fprintf(file.out,"Vector not found in generate_all_crystal_atoms\n");
        }
       } // close loop over k
      } // close loop over i

    atom_i->numb[atm] = count1[atm];
    if (job->taskid == 0 && job->verbosity > 1) {
    for (i = count; i < count + count1[atm]; i++) {
      for (k = 0; k < symmetry->number_of_operators; k++) {
        fprintf(file.out, "atom_i[%2d][%2d] elem O P numb %3d %3d %3d %3d \n", i, k, \
        atom_i->K[i * symmetry->number_of_operators + k],atom_i->O[i * symmetry->number_of_operators + k], atom_i->P[i * symmetry->number_of_operators + k],\
        atom_i->numb[atm]);
       }
      }
      if (job->taskid == 0 && job->verbosity > 1)
      fprintf(file.out, "\n");
     }

      for (i = count; i < count + count1[atm]; i++) {
        for (k = 0; k < symmetry->number_of_operators; k++) {
          if (atom_i->K[i * symmetry->number_of_operators + k] == -1 || atom_i->O[i * symmetry->number_of_operators + k] == -1) {
          if (job->taskid == 0) fprintf(file.out,"Transformation table atoms->atom_i incorrectly initialised\n");
          MPI_Finalize();
          exit(1);
         }
        }
       }
      count += count1[atm];
     } // close loop over atm

}

void count_basis_sets(ATOM *atoms, int flag, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // *  Count number of basis sets                                                            *
  // ******************************************************************************************

  char title[150], junk;
  int j, k;
  ATOMTYPE irtmp;
  int atomic_number;

 if (job->taskid == 0 && job->verbosity > 1) {
  fprintf(file.out, "BASIS SET INFORMATION\n");
  fprintf(file.out, "  ATOMIC NO.    NO. SHELLS    SHELL TYPE    NO. GAUSSIANS\n");
  fflush(file.out);
 }

  atoms->number_of_basis_sets = 0;

  do {
    read_line(file.in, title, 150); // reading atomic number and number of shells
    sscanf(title, "%d %d", &irtmp.atomic_number, &irtmp.number_of_shells);
    for (j = 0; j < irtmp.number_of_shells; j++) {
      read_line(file.in, title, 150);
      sscanf(title, "%c %d %d %lf", &junk, &irtmp.shell[j].shell_type, &irtmp.shell[j].number_of_Gaussians,
      &irtmp.shell[j].number_of_electrons);
      if (job->taskid == 0 && job->verbosity > 1)
      fprintf(file.out, "    %2d             %2d           %2d              %2d\n", irtmp.atomic_number,
      irtmp.number_of_shells, irtmp.shell[j].shell_type, irtmp.shell[j].number_of_Gaussians);

      switch (irtmp.shell[j].shell_type) {

        case 0:
        case 1:
        case 2:
        case 3:

          break;

        case 4:

            if (job->l_max < 3) job->l_max = 3;

          break;

        case 5:

            if (job->l_max < 4) job->l_max = 4;

          break;

      } // end switch

      switch (irtmp.shell[j].shell_type) {

        case 0:
        case 2:
        case 3:
        case 4:
        case 5:

          for (k = 0; k < irtmp.shell[j].number_of_Gaussians; k++) {
            read_line(file.in, title, 150);
            sscanf(title, "%lf %lf", &irtmp.shell[j].exponent[k], &irtmp.shell[j].coefc[k]);
            if (job->taskid == 0 && job->verbosity > 1)
            fprintf(file.out, "%10.5lf %10.5lf\n", irtmp.shell[j].exponent[k], irtmp.shell[j].coefc[k]);
            if (fabs(irtmp.shell[j].exponent[k]) < 1.0e-10 || fabs(irtmp.shell[j].coefc[k]) < 1.0e-10) {
            if (job->taskid == 0) {
            fprintf(file.out,"Read zero exponent or expansion coefficient: check for D+01 format in INPUT\n");
            fprintf(file.out, "%10.5lf %10.5lf\n", irtmp.shell[j].exponent[k], irtmp.shell[j].coefc[k]);
           }
            MPI_Finalize();
            exit(1);
           }
          }

          break;

        case 1:

          for (k = 0; k < irtmp.shell[j].number_of_Gaussians; k++) {
            read_line(file.in, title, 150);
            sscanf(title, "%lf %lf %lf", &irtmp.shell[j].exponent[k], &irtmp.shell[j].coefc[k],
            &irtmp.shell[j].coefp[k]);
            if (job->taskid == 0 && job->verbosity > 1)
            fprintf(file.out, "%10.5lf %10.5lf %10.5lf\n", irtmp.shell[j].exponent[k], irtmp.shell[j].coefc[k],
            irtmp.shell[j].coefp[k]);
            if (fabs(irtmp.shell[j].exponent[k]) < 1.0e-10 || fabs(irtmp.shell[j].coefc[k]) < 1.0e-10 || \
            fabs(irtmp.shell[j].coefp[k]) < 1.0e-10) {
            if (job->taskid == 0 && job->taskid == 0)
            fprintf(file.out,"Read zero exponent or expansion coefficient: check for D+01 format in INPUT\n");
            fprintf(file.out, "%10.5lf %10.5lf %10.5lf\n", irtmp.shell[j].exponent[k], irtmp.shell[j].coefc[k],\
            irtmp.shell[j].coefp[k]);
            MPI_Finalize();
            exit(1);
           }
          }

          break;

        default:

          fprintf(file.out, "ERROR: Basis set type unrecognised\n");
          fclose(file.in);
          if (job->taskid == 0)
          fclose(file.out);
          MPI_Finalize();
          exit(1);

      } // end switch
    } // end loop on shells

    atomic_number = irtmp.atomic_number;
    if      (flag == 0 && atomic_number <  300 && atomic_number != 99)  (atoms->number_of_basis_sets)++;   // Wave function basis set
    else if (flag == 1 && atomic_number >= 300)                         (atoms->number_of_basis_sets)++;   // Auxillary basis set

  } while (atomic_number != 99);

    if (job->taskid == 0)
    fprintf(file.out, "\n");

}

void read_basis_sets(ATOM *atoms, int flag, ATOMTYPE *iratom, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // *  Read and normalize basis sets from INPUT                                              *
  // ******************************************************************************************

  char title[150], junk;
  int i, j, k, l;
  int atomic_number;

  if (job->taskid == 0 && job->verbosity > 1) {
  fprintf(file.out, "\n                    BASIS SET INFORMATION\n");
  fprintf(file.out, "  ATOMIC NO.    NO. SHELLS    SHELL TYPE    NO. GAUSSIANS\n");
 }

  i = 0;

  do {
    read_line(file.in, title, 150); // reading atomic number and number of shells
    sscanf(title, "%3d %3d", &iratom[i].atomic_number, &iratom[i].number_of_shells);
    atomic_number = iratom[i].atomic_number;
    if (iratom[i].atomic_number == 99)

      break;

    for (j = 0; j < iratom[i].number_of_shells; j++) {
      read_line(file.in, title, 150);
      sscanf(title, "%c %d %d %lf", &junk, &iratom[i].shell[j].shell_type, &iratom[i].shell[j].number_of_Gaussians,
      &iratom[i].shell[j].number_of_electrons);

      if ((flag == 0 && iratom[i].atomic_number < 300 && iratom[i].atomic_number != 99) && (job->taskid == 0 && job->verbosity > 1))
      fprintf(file.out, "    %3d             %3d           %3d              %3d\n", iratom[i].atomic_number,
      iratom[i].number_of_shells, iratom[i].shell[j].shell_type, iratom[i].shell[j].number_of_Gaussians);
      else if ((flag == 1 && iratom[i].atomic_number >= 300) && (job->taskid == 0 && job->verbosity > 1))
      fprintf(file.out, "    %3d             %3d           %3d              %3d\n", iratom[i].atomic_number,
      iratom[i].number_of_shells, iratom[i].shell[j].shell_type, iratom[i].shell[j].number_of_Gaussians);

      switch (iratom[i].shell[j].shell_type) {

        case 0:
        case 2:
        case 3:
        case 4:
        case 5:

    for (k = 0; k < iratom[i].shell[j].number_of_Gaussians; k++) {
      read_line(file.in, title, 150);
      sscanf(title, "%lf %lf", &iratom[i].shell[j].exponent[k], &iratom[i].shell[j].coefc[k]);
      if ((flag == 0 && iratom[i].atomic_number < 300 && iratom[i].atomic_number != 99) && (job->taskid == 0 && job->verbosity > 1))
      fprintf(file.out, "%16.6lf %16.6lf\n", iratom[i].shell[j].exponent[k], iratom[i].shell[j].coefc[k]);
      else if ((flag == 1 && iratom[i].atomic_number >= 300) && (job->taskid == 0 && job->verbosity > 1))
      fprintf(file.out, "%16.6lf %16.6lf\n", iratom[i].shell[j].exponent[k], iratom[i].shell[j].coefc[k]);
     }

          break;

        case 1:

    for (k = 0; k < iratom[i].shell[j].number_of_Gaussians; k++) {
      read_line(file.in, title, 150);
      sscanf(title, "%lf %lf %lf", &iratom[i].shell[j].exponent[k], &iratom[i].shell[j].coefc[k],
      &iratom[i].shell[j].coefp[k]);
      if ((flag == 0 && iratom[i].atomic_number < 300 && iratom[i].atomic_number != 99) && (job->taskid == 0 && job->verbosity > 1))
      fprintf(file.out, "%16.6lf %16.6lf %16.6lf\n", iratom[i].shell[j].exponent[k], iratom[i].shell[j].coefc[k],
      iratom[i].shell[j].coefp[k]);
      else if ((flag == 1 && iratom[i].atomic_number >= 300) && (job->taskid == 0 && job->verbosity > 1))
      fprintf(file.out, "%16.5lf %16.5lf %16.5lf\n", iratom[i].shell[j].exponent[k], iratom[i].shell[j].coefc[k],
      iratom[i].shell[j].coefp[k]);
     }

          break;

        default:

          fprintf(file.out, "ERROR: Basis set type unrecognised\n");
          fclose(file.in);
          if (job->taskid == 0)
          fclose(file.out);
          MPI_Finalize();
          exit(1);

      } // end switch
     }

    if      (flag == 0 && iratom[i].atomic_number <  300 && iratom[i].atomic_number != 99) i++;
    else if (flag == 1 && iratom[i].atomic_number >= 300) i++;

  } while (atomic_number != 99);

    if (job->taskid == 0)
    fprintf(file.out, "\n");

    if (job->taskid == 0 && job->verbosity > 1)
    fprintf(file.out, "Number of basis sets read: %3d flag %3d (0 = wfn basis, 1 = aux basis)\n", \
    atoms->number_of_basis_sets, flag); 

  // ******************************************************************************************
  // *  Normalise basis sets                                                                  *
  // ******************************************************************************************

  for (i = 0; i < atoms->number_of_basis_sets; i++) {
    for (j = 0; j < iratom[i].number_of_shells; j++) {
      for (k = 0; k < iratom[i].shell[j].number_of_Gaussians; k++) {

        switch (iratom[i].shell[j].shell_type) {

          case 0:

            iratom[i].shell[j].coefc[k] *= pow(2.0 * iratom[i].shell[j].exponent[k] / pi, 0.75);

            break;

          case 1:

            iratom[i].shell[j].coefc[k] *= pow(2.0 * iratom[i].shell[j].exponent[k] / pi, 0.75);
            iratom[i].shell[j].coefp[k] *= pow(128.0 * iratom[i].shell[j].exponent[k] / pi / pi / pi, 0.25)
            * iratom[i].shell[j].exponent[k];

            break;

          case 2:

            iratom[i].shell[j].coefc[k] *= pow(128.0 * iratom[i].shell[j].exponent[k] / pi / pi / pi, 0.25)
            * iratom[i].shell[j].exponent[k];

            break;

          case 3:

            iratom[i].shell[j].coefc[k] *= pow(2048.0 / pi / pi / pi, 0.25) * iratom[i].shell[j].exponent[k] * \
            iratom[i].shell[j].exponent[k] * pow(iratom[i].shell[j].exponent[k], -0.25);

            break;

          case 4:

            iratom[i].shell[j].coefc[k] *= pow(32768.0 * iratom[i].shell[j].exponent[k] / pi / pi / pi, 0.25) *
            iratom[i].shell[j].exponent[k] * iratom[i].shell[j].exponent[k];

            break;

          case 5:

            iratom[i].shell[j].coefc[k] *= pow(524288.0 / three * iratom[i].shell[j].exponent[k] / pi / pi / pi, 0.25) *
            iratom[i].shell[j].exponent[k] * iratom[i].shell[j].exponent[k] * sqrt(iratom[i].shell[j].exponent[k]);

            break;

        } // end switch

       }

      double norm = k_zero, norm1 = k_zero, pa;

      for (k = 0; k < iratom[i].shell[j].number_of_Gaussians; k++) {
        for (l = 0; l < iratom[i].shell[j].number_of_Gaussians; l++) {
          pa = iratom[i].shell[j].exponent[k] + iratom[i].shell[j].exponent[l];
          //fprintf(file.out,"pa %3d %3d shell type  %3d  %10.4lf\n",k,l,iratom[i].shell[j].shell_type,pa) ;
          //fprintf(file.out,"%lf %lf\n",iratom[i].shell[j].coefc[k],iratom[i].shell[j].coefc[l]) ;
          switch (iratom[i].shell[j].shell_type) {

            case 0:

              norm += iratom[i].shell[j].coefc[k] * iratom[i].shell[j].coefc[l] * pow(pi / pa, 1.5);
              //fprintf(file.out,"%lf %lf\n",iratom[i].shell[j].coefc[k],iratom[i].shell[j].coefc[l]) ;

              break;

            case 1:

              norm += iratom[i].shell[j].coefc[k] * iratom[i].shell[j].coefc[l] * pow(pi / pa, 1.5);
              norm1 += iratom[i].shell[j].coefp[k] * iratom[i].shell[j].coefp[l] * pow(pi / pa, 1.5) / 2.0 / pa;
              //fprintf(file.out,"pa %10.4lf\n",pa) ;
              //fprintf(file.out,"%lf %lf\n",iratom[i].shell[j].coefc[k],iratom[i].shell[j].coefc[l]) ;
              //fprintf(file.out,"%lf %lf\n",iratom[i].shell[j].coefp[k],iratom[i].shell[j].coefp[l]) ;

              break;

            case 2:

              norm += iratom[i].shell[j].coefc[k] * iratom[i].shell[j].coefc[l] * pow(pi / pa, 1.5) / 2.0 / pa;

              break;

            case 3:

              norm += iratom[i].shell[j].coefc[k] * iratom[i].shell[j].coefc[l] * pow(pi / pa, 1.5) / 4.0 / pa / pa;

              break;

            case 4:

              norm += iratom[i].shell[j].coefc[k] * iratom[i].shell[j].coefc[l] * pow(pi / pa, 1.5) / 8.0 / pa / pa / pa;
              //fprintf(file.out,"%lf %lf\n",iratom[i].shell[j].coefc[k],iratom[i].shell[j].coefc[l]) ;

              break;

            case 5:

              norm += iratom[i].shell[j].coefc[k] * iratom[i].shell[j].coefc[l] * rtthree / sixteen * pow(pi / pa, 1.5) / \
              pa / pa / pa / pa;
              //fprintf(file.out,"%lf %lf\n",iratom[i].shell[j].coefc[k],iratom[i].shell[j].coefc[l]) ;

              break;

          } // end switch
         }
        }

      for (k = 0; k < iratom[i].shell[j].number_of_Gaussians; k++) {

        switch (iratom[i].shell[j].shell_type) {

          case 0:
          case 2:
          case 3:
          case 4:
          case 5:

            iratom[i].shell[j].coefc[k] /= sqrt(norm);
            //fprintf(file.out,"Normalised %10.5lf %10.5lf %10.5lf\n",norm,iratom[i].shell[j].exponent[k],iratom[i].shell[j].coefc[k]) ;

            break;

          case 1:

            iratom[i].shell[j].coefc[k] /= sqrt(norm);
            iratom[i].shell[j].coefp[k] /= sqrt(norm1);
            //fprintf(file.out, "Normalised %10.5lf %10.5lf %10.5lf\n", iratom[i].shell[j].exponent[k],iratom[i].shell[j].coefc[k],\
            iratom[i].shell[j].coefp[k]) ;

            break;

        } // end switch

       }
      }
     }

}

void generate_runtime_arrays(ATOM *atoms, ATOM_TRAN *atom_p, ATOMTYPE *iratom, SHELL *shell, GAUSSIAN *gaussian, JOB_PARAM *job, FILES file)

{

  // ******************************************************************************************
  // *  Generate arrays in ATOM, SHELL and GAUSSIAN structures                                *
  // ******************************************************************************************

  int i, j, k, l;
  int shell_index;
  int bfn_index;
  int shl_index;
  int g, g_index;
  double min_expo, min_coef;

  atoms->number_of_bfns_in_unit_cell = 0;
  atoms->number_of_electrons_in_unit_cell = k_zero;
  bfn_index = 0;
  shl_index = 0;
  for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++) {
  min_expo = 999999.0;
    atoms->bfnposn[i] = bfn_index;
    atoms->shlposn[i] = shl_index;
    for (j = 0; j < atoms->number_of_basis_sets; j++) {
      //fprintf(file.out, "runtime %3d %3d %3d %3d %3d %3d\n", \
      i, j, iratom[j].atomic_number, iratom[j].atomic_number % 200, atoms->atomic_number[i], atoms->atomic_number[i] % 200);
      if (atoms->atomic_number[i] == iratom[j].atomic_number % 300 || atoms->atomic_number[i] == iratom[j].atomic_number % 200) {
        atoms->nshel[i] = iratom[j].number_of_shells;
         shell_index = 0;
          g_index = 0;
          for (l = 0; l < j; l++) {
           shell_index += iratom[l].number_of_shells;
            for (k = 0; k < iratom[l].number_of_shells; k++) {
             g_index += iratom[l].shell[k].number_of_Gaussians; 
           }
           }
            atoms->shelposn[i] = shell_index;
            atoms->gausposn[i] = g_index;
            atoms->basis_set[i] = j;
          if (job->taskid == 0 && job->verbosity > 2)
          fprintf(file.out, "shell.posn[%3d] %3d basis.posn[%3d] %3d gaussian.posn[%3d] %3d nshel[%3d] %3d\n", 
          i, shl_index, i, bfn_index, i, g_index, i, atoms->nshel[i]);
        for (k = 0; k < iratom[j].number_of_shells; k++) {
          atoms->number_of_electrons_in_unit_cell += iratom[j].shell[k].number_of_electrons;
          shell->ng[shell_index] = iratom[j].shell[k].number_of_Gaussians;
          shell->nele[shell_index] = iratom[j].shell[k].number_of_electrons;
          if (job->taskid == 0 && job->verbosity > 2)
          fprintf(file.out, "atom in cell %3d basis set %3d shell no. %3d %3d %9.2e\n", i,j,k,iratom[j].shell[k].shell_type,\
          iratom[j].shell[k].number_of_electrons) ;
          switch (iratom[j].shell[k].shell_type) {
            case 0:
             shell->imax[shell_index] = 0;
             shell->type[shell_index] = 1;
             shell->type1[shell_index] = 1;
             shell->ord[shell_index] = 0;
             shell->opp[shell_index] = 0;
             //fprintf(file.out,"shells->type1[%2d] %2d\n",shell_index,shell->type1[shell_index]);
              for (g = 0; g < iratom[j].shell[k].number_of_Gaussians; g++) {
                gaussian->expo[g_index] = iratom[j].shell[k].exponent[g];
                gaussian->sc[g_index] = iratom[j].shell[k].coefc[g];
                gaussian->pc[g_index] = k_zero;
                gaussian->dc[g_index] = k_zero;
                gaussian->fc[g_index] = k_zero;
                gaussian->gc[g_index] = k_zero;
                //fprintf(file.out, "g_index %d  %lf %lf \n", g_index, gaussian->sc[g_index],gaussian->pc[g_index]);
                if (gaussian->expo[g_index] < min_expo) {
                min_expo = gaussian->expo[g_index];
                min_coef = pow(2.0 * gaussian->expo[g_index] / pi, 0.75);
               }
                //fprintf(file.out, "%d  %lf %lf %lf %lf\n", g_index, gaussian->expo[g_index], gaussian->sc[g_index],\
                gaussian->pc[g_index], gaussian->dc[g_index]) ;
                g_index++;
              }
              bfn_index++;
              shl_index++;
              atoms->number_of_bfns_in_unit_cell += 1;
              //fprintf(file.out,"shell bfnposn %3d %3d\n",shl_index,bfn_index) ;
              //fprintf(file.out,"shell bfnposn %3d %3d %e %e\n",shell_index,bfn_index,shell->min_coef[shell_index],\
              shell->min_expo[shell_index]) ;
              break;
            case 1:
             shell->imax[shell_index] = 1;
             shell->type[shell_index] = 4;
             shell->type1[shell_index] = 4;
             shell->ord[shell_index] = 1;
             shell->opp[shell_index] = 0;
             //fprintf(file.out,"shells->type1[%2d] %2d\n",shell_index,shell->type1[shell_index]);
              for (g = 0; g < iratom[j].shell[k].number_of_Gaussians; g++) {
                gaussian->expo[g_index] = iratom[j].shell[k].exponent[g];
                gaussian->sc[g_index] = iratom[j].shell[k].coefc[g];
                gaussian->pc[g_index] = iratom[j].shell[k].coefp[g];
                gaussian->dc[g_index] = k_zero;
                gaussian->fc[g_index] = k_zero;
                gaussian->gc[g_index] = k_zero;
                //fprintf(file.out, "g_index %d  %lf %lf \n", g_index, gaussian->sc[g_index],gaussian->pc[g_index]);
                if (gaussian->expo[g_index] < min_expo) {
                min_expo = gaussian->expo[g_index];
                min_coef = pow(2.0 * gaussian->expo[g_index] / pi, 0.75);
               }
                //fprintf(file.out, "%d  %lf %lf %lf %lf\n", g_index, gaussian->expo[g_index], gaussian->sc[g_index],\
                gaussian->pc[g_index], gaussian->dc[g_index]) ;
                g_index++;
              }
              bfn_index += 4;
              shl_index++;
              atoms->number_of_bfns_in_unit_cell += 4;
              //fprintf(file.out,"shell bfnposn %3d %3d\n",shl_index,bfn_index) ;
              //fprintf(file.out,"shell bfnposn g_index %3d %3d %3d\n",shell_index,bfn_index,g_index) ;
              //fprintf(file.out,"shell bfnposn %3d %3d %e %e\n",shell_index,bfn_index,shell->min_coef[shell_index],\
              shell->min_expo[shell_index]) ;
              break;
            case 2:
             shell->imax[shell_index] = 1;
             shell->type[shell_index] = 3;
             shell->type1[shell_index] = 3;
             shell->ord[shell_index] = 2;
             shell->opp[shell_index] = 0;
             //fprintf(file.out,"shells->type1[%2d] %2d\n",shell_index,shell->type1[shell_index]);
              for (g = 0; g < iratom[j].shell[k].number_of_Gaussians; g++) {
                gaussian->expo[g_index] = iratom[j].shell[k].exponent[g];
                gaussian->sc[g_index] = k_zero;
                gaussian->pc[g_index] = iratom[j].shell[k].coefc[g];
                gaussian->dc[g_index] = k_zero;
                gaussian->fc[g_index] = k_zero;
                gaussian->gc[g_index] = k_zero;
                if (gaussian->expo[g_index] < min_expo) {
                min_expo = gaussian->expo[g_index];
                min_coef = pow(2.0 * gaussian->expo[g_index] / pi, 0.75);
               }
                //fprintf(file.out, "%d  %lf %lf %lf %lf\n", g_index, gaussian->expo[g_index], gaussian->sc[g_index],\
                gaussian->pc[g_index], gaussian->dc[g_index]) ;
                g_index++;
              }
              bfn_index += 3;
              shl_index++;
              atoms->number_of_bfns_in_unit_cell += 3;
              //fprintf(file.out,"shell bfnposn %3d %3d\n",shl_index,bfn_index) ;
              break;
            case 3:
             shell->imax[shell_index] = 2;
             shell->type[shell_index] = 6;
             shell->type1[shell_index] = 5;
             shell->ord[shell_index] = 3;
             shell->opp[shell_index] = 4;
             //fprintf(file.out,"shells->type1[%2d] %2d\n",shell_index,shell->type1[shell_index]);
              for (g = 0; g < iratom[j].shell[k].number_of_Gaussians; g++) {
                gaussian->expo[g_index] = iratom[j].shell[k].exponent[g];
                gaussian->sc[g_index] = k_zero;
                gaussian->pc[g_index] = k_zero;
                gaussian->dc[g_index] = iratom[j].shell[k].coefc[g];
                gaussian->fc[g_index] = k_zero;
                gaussian->gc[g_index] = k_zero;
                //fprintf(file.out, "g_index %d  %lf %lf \n", g_index, gaussian->dc[g_index],gaussian->pc[g_index]);
                if (gaussian->expo[g_index] < min_expo) {
                min_expo = gaussian->expo[g_index];
                min_coef = pow(2.0 * gaussian->expo[g_index] / pi, 0.75);
               }
                //fprintf(file.out, "%d  %lf %lf %lf %lf\n", g_index, gaussian->expo[g_index], gaussian->sc[g_index],\
                gaussian->pc[g_index], gaussian->dc[g_index]) ;
                g_index++;
              }
              bfn_index += 6;
              shl_index++;
              atoms->number_of_bfns_in_unit_cell += 6;
              //fprintf(file.out,"shell bfnposn %3d %3d\n",shl_index,bfn_index) ;
              break;
          case 4:
           shell->imax[shell_index] = 3;
           shell->type[shell_index] = 10;
           shell->type1[shell_index] = 7;
           shell->ord[shell_index] = 4;
           shell->opp[shell_index] = 12;
           //fprintf(file.out,"shells->type1[%2d] %2d\n",shell_index,shell->type1[shell_index]);
              for (g = 0; g < iratom[j].shell[k].number_of_Gaussians; g++) {
                gaussian->expo[g_index] = iratom[j].shell[k].exponent[g];
                gaussian->sc[g_index] = k_zero;
                gaussian->pc[g_index] = k_zero;
                gaussian->dc[g_index] = k_zero;
                gaussian->fc[g_index] = iratom[j].shell[k].coefc[g];
                gaussian->gc[g_index] = k_zero;
                if (gaussian->expo[g_index] < min_expo) {
                min_expo = gaussian->expo[g_index];
                min_coef = pow(2.0 * gaussian->expo[g_index] / pi, 0.75);
               }
                //fprintf(file.out, "%d  %lf %lf %lf %lf\n", g_index, gaussian->expo[g_index], gaussian->sc[g_index],\
                gaussian->pc[g_index], gaussian->fc[g_index]) ;
                g_index++;
             }
               bfn_index += 10;
               shl_index++;
               atoms->number_of_bfns_in_unit_cell += 10;
               //fprintf(file.out,"shell bfnposn %3d %3d\n",shl_index,bfn_index) ;
           break;
          case 5:
           shell->imax[shell_index] = 4;
           shell->type[shell_index] = 15;
           shell->type1[shell_index] = 9;
           shell->ord[shell_index] = 5;
           shell->opp[shell_index] = 28;
           //fprintf(file.out,"shells->type1[%2d] %2d\n",shell_index,shell->type1[shell_index]);
              for (g = 0; g < iratom[j].shell[k].number_of_Gaussians; g++) {
                gaussian->expo[g_index] = iratom[j].shell[k].exponent[g];
                gaussian->sc[g_index] = k_zero;
                gaussian->pc[g_index] = k_zero;
                gaussian->dc[g_index] = k_zero;
                gaussian->fc[g_index] = k_zero;
                gaussian->gc[g_index] = iratom[j].shell[k].coefc[g];
                if (gaussian->expo[g_index] < min_expo) {
                min_expo = gaussian->expo[g_index];
                min_coef = pow(2.0 * gaussian->expo[g_index] / pi, 0.75);
               }
                g_index++;
             }
               bfn_index += 15;
               shl_index++;
               atoms->number_of_bfns_in_unit_cell += 15;
               //fprintf(file.out,"shell bfnposn %3d %3d\n",shl_index,bfn_index) ;
           break;
          } // close switch
          shell_index++;
          //fprintf(file.out,"shell index %d %d\n",shell_index,atoms->number_of_shells_in_unit_cell);
        } // close loop over k
      } // close if (at..
    } // close loop over j
    atoms->bfnnumb[i] = bfn_index - atoms->bfnposn[i]; // number of basis functions on ith atom
    atom_p->expo[i] = min_expo;
    atom_p->coef[i] = min_coef;
     if (job->taskid == 0 && job->verbosity > 1)
     fprintf(file.out,"Number of cart basis functions %4d on atom %4d bfnposn %4d\n",atoms->bfnnumb[i],i,atoms->bfnposn[i]); 
  } // close loop over i
  if (job->taskid == 0 && job->verbosity > 1)
  fprintf(file.out,"\n");

int  bfn_index_sh = 0;
int  shl_index_sh = 0;
int  shell_index_sh;
int g_index_sh;
atoms->number_of_sh_bfns_in_unit_cell = 0;

  for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++) {

    for (j = 0; j < atoms->number_of_basis_sets; j++) {
     //fprintf(file.out, "%3d%3d%3d%3d\n", i, j, iratom[j].atomic_number, atoms->atomic_number[i]);
      if (atoms->atomic_number[i] == iratom[j].atomic_number % 300 || atoms->atomic_number[i] == iratom[j].atomic_number % 200) {
         shell_index = 0;
         shell_index_sh = 0;
          g_index = 0;
          g_index_sh = 0;
          for (l = 0; l < j; l++) {
           shell_index += iratom[l].number_of_shells;
           shell_index_sh += iratom[l].number_of_shells;
            for (k = 0; k < iratom[l].number_of_shells; k++) {
             g_index += iratom[l].shell[k].number_of_Gaussians; 
             g_index_sh += iratom[l].shell[k].number_of_Gaussians; 
             if (iratom[l].shell[k].shell_type == 1) {
              shell_index_sh++;
              g_index_sh += iratom[l].shell[k].number_of_Gaussians;
             }
           }
           }

        //fprintf(file.out,"indices for atom %d atomic number %d  %d %d %d %d\n",i,atoms->atomic_number[i],\
        shell_index,shell_index_sh,g_index,g_index_sh);

          atoms->bfnposn_sh[i] = bfn_index_sh;
          atoms->shlposn_sh[i] = shl_index_sh;
          atoms->nshel_sh[i] = atoms->nshel[i];
          atoms->shelposn_sh[i] = shell_index_sh;
          atoms->gausposn_sh[i] = g_index_sh;

    for (l = 0; l < atoms->nshel[i]; l++) {
         if (iratom[j].shell[l].shell_type == 1) 
          atoms->nshel_sh[i]++;
        }

        if (job->taskid == 0 && job->verbosity > 2)
          fprintf(file.out, "shell.posn[%2d] %3d basis.posn[%2d] %3d gaussian.posn[%2d] %3d nshel[%2d] %3d\n", 
          i, atoms->shelposn_sh[i], i, atoms->bfnposn_sh[i], i, atoms->gausposn_sh[i], i, atoms->nshel_sh[i]);

    for (k = 0; k < atoms->nshel[i]; k++) {

        shell->min_expo_sh[shell_index_sh] = 999999.0;

        switch (shell->type[shell_index]) {
          case 1:
           shell->imax_sh[shell_index_sh] = 0;
           shell->type_sh[shell_index_sh] = 1;
           shell->type1_sh[shell_index_sh] = 1;
           shell->cart[shell_index_sh] = 1;
           shell->shar[shell_index_sh] = 1;
           shell->ord_sh[shell_index_sh] = 0;
           shell->ng_sh[shell_index_sh] = shell->ng[shell_index];;
           shell->opp_sh[shell_index_sh] = 0;
           //fprintf(file.out,"shells->type1_sh[%2d] %2d\n",shell_index_sh,shell->type1_sh[shell_index_sh]);
             for (g = 0; g < shell->ng_sh[shell_index_sh]; g++) {
                gaussian->expo_sh[g_index_sh] = gaussian->expo[g_index];
                gaussian->c_sh[g_index_sh] = gaussian->sc[g_index];
                //fprintf(file.out,"shell type %d %d   %d %d gauss %d %d    %e %e   %e  %e\n",shell->type[shell_index],\
                shell->type_sh[shell_index_sh],shell_index,shell_index_sh,g_index,g_index_sh,gaussian->expo[g_index],\
                gaussian->expo_sh[g_index_sh],gaussian->sc[g_index],gaussian->c_sh[g_index_sh]);
                if (gaussian->expo_sh[g_index_sh] < shell->min_expo_sh[shell_index_sh]) {
                shell->min_expo_sh[shell_index_sh] = gaussian->expo_sh[g_index_sh];
                shell->min_coef_sh[shell_index_sh] = pow(2.0 * gaussian->expo_sh[g_index_sh] / pi, 0.75);
                //fprintf(file.out,"shell bfnposn %3d %3d %e %e\n",shell_index_sh,bfn_index_sh,shell->min_coef_sh[shell_index_sh], \
                shell->min_expo_sh[shell_index_sh]) ;
               }
                g_index++;
                g_index_sh++;
             }
               bfn_index_sh += 1;
               shl_index_sh++;
               atoms->number_of_sh_bfns_in_unit_cell += 1;
           break;
          case 3:
           shell->imax_sh[shell_index_sh] = 1;
           shell->type_sh[shell_index_sh] = 3;
           shell->type1_sh[shell_index_sh] = 3;
           shell->cart[shell_index_sh] = 3;
           shell->shar[shell_index_sh] = 3;
           shell->ord_sh[shell_index_sh] = 2;
           shell->ng_sh[shell_index_sh] = shell->ng[shell_index];
           shell->opp_sh[shell_index_sh] = 0;
           //fprintf(file.out,"shells->type1_sh[%2d] %2d\n",shell_index_sh,shell->type1_sh[shell_index_sh]);
             for (g = 0; g < shell->ng_sh[shell_index_sh]; g++) {
                gaussian->expo_sh[g_index_sh] = gaussian->expo[g_index];
                gaussian->c_sh[g_index_sh] = gaussian->pc[g_index];
                if (gaussian->expo_sh[g_index_sh] < shell->min_expo_sh[shell_index_sh]) {
                shell->min_expo_sh[shell_index_sh] = gaussian->expo_sh[g_index_sh];
                shell->min_coef_sh[shell_index_sh] = pow(2.0 * gaussian->expo_sh[g_index_sh] / pi, 0.75);
                //fprintf(file.out,"shell bfnposn %3d %3d %e %e\n",shell_index_sh,bfn_index_sh,shell->min_coef_sh[shell_index_sh], \
                shell->min_expo_sh[shell_index_sh]) ;
               }
                g_index++;
                g_index_sh++;
             }
               bfn_index_sh += 3;
               shl_index_sh++;
               atoms->number_of_sh_bfns_in_unit_cell += 3;
           break;
          case 4:
           shell->imax_sh[shell_index_sh] = 0;
           shell->type_sh[shell_index_sh] = 1;
           shell->type1_sh[shell_index_sh] = 1;
           shell->cart[shell_index_sh] = 1;
           shell->shar[shell_index_sh] = 1;
           shell->ord_sh[shell_index_sh] = 0;
           shell->opp_sh[shell_index_sh] = 0;
           shell->ng_sh[shell_index_sh] = shell->ng[shell_index];
           //fprintf(file.out,"shells->type1_sh[%2d] %2d\n",shell_index_sh,shell->type1_sh[shell_index_sh]);
             for (g = 0; g < shell->ng_sh[shell_index_sh]; g++) {
                gaussian->expo_sh[g_index_sh] = gaussian->expo[g_index];
                gaussian->c_sh[g_index_sh] = gaussian->sc[g_index];
                //fprintf(file.out,"shell type %d %d   %d %d gauss %d %d    %e %e   %e  %e\n",shell->type[shell_index],\
                shell->type_sh[shell_index_sh],shell_index,shell_index_sh,g_index,g_index_sh,gaussian->expo[g_index],\
                gaussian->expo_sh[g_index_sh],gaussian->sc[g_index],gaussian->c_sh[g_index_sh]);
                if (gaussian->expo_sh[g_index_sh] <= shell->min_expo_sh[shell_index_sh]) {
                shell->min_expo_sh[shell_index_sh] = gaussian->expo_sh[g_index_sh];
                shell->min_coef_sh[shell_index_sh] = pow(2.0 * gaussian->expo_sh[g_index_sh] / pi, 0.75);
                //fprintf(file.out,"shell bfnposn %3d %3d %e %e\n",shell_index_sh,bfn_index_sh,shell->min_coef_sh[shell_index_sh], \
                shell->min_expo_sh[shell_index_sh]) ;
               }
                g_index++;
                g_index_sh++;
             }
               bfn_index_sh += 1;
               shl_index_sh++;
               atoms->number_of_sh_bfns_in_unit_cell += 1;
           shell_index_sh++;

           shell->min_expo_sh[shell_index_sh] = 999999.0;

           shell->imax_sh[shell_index_sh] = 1;
           shell->type_sh[shell_index_sh] = 3;
           shell->type1_sh[shell_index_sh] = 3;
           shell->cart[shell_index_sh] = 3;
           shell->shar[shell_index_sh] = 3;
           shell->ord_sh[shell_index_sh] = 2;
           shell->ng_sh[shell_index_sh] = shell->ng[shell_index];;
           shell->opp_sh[shell_index_sh] = 0;
           g_index -= shell->ng_sh[shell_index_sh];
           //fprintf(file.out,"shells->type1_sh[%2d] %2d\n",shell_index_sh,shell->type1_sh[shell_index_sh]);
             for (g = 0; g < shell->ng_sh[shell_index_sh]; g++) {
                gaussian->expo_sh[g_index_sh] = gaussian->expo[g_index];
                gaussian->c_sh[g_index_sh] = gaussian->pc[g_index];
                //fprintf(file.out,"shell type %d %d   %d %d gauss %d %d    %e %e   %e  %e\n",shell->type[shell_index],\
                shell->type_sh[shell_index_sh],shell_index,shell_index_sh,g_index,g_index_sh,gaussian->expo[g_index],\
                gaussian->expo_sh[g_index_sh],gaussian->pc[g_index],gaussian->c_sh[g_index_sh]);
                if (gaussian->expo_sh[g_index_sh] < shell->min_expo_sh[shell_index_sh]) {
                shell->min_expo_sh[shell_index_sh] = gaussian->expo_sh[g_index_sh];
                shell->min_coef_sh[shell_index_sh] = pow(2.0 * gaussian->expo_sh[g_index_sh] / pi, 0.75);
                //fprintf(file.out,"shell bfnposn %3d %3d %e %e\n",shell_index_sh,bfn_index_sh,shell->min_coef_sh[shell_index_sh], \
                shell->min_expo_sh[shell_index_sh]) ;
               }
                g_index++;
                g_index_sh++;
             }
               bfn_index_sh += 3;
               shl_index_sh++;
               atoms->number_of_sh_bfns_in_unit_cell += 3;
           break;
          case 6:
           shell->imax_sh[shell_index_sh] = 2;
           shell->type_sh[shell_index_sh] = 5;
           shell->type1_sh[shell_index_sh] = 6;
           shell->cart[shell_index_sh] = 6;
           shell->shar[shell_index_sh] = 5;
           shell->ord_sh[shell_index_sh] = 3;
           shell->opp_sh[shell_index_sh] = 4;
           shell->ng_sh[shell_index_sh] = shell->ng[shell_index];
           //fprintf(file.out,"shells->type1_sh[%2d] %2d\n",shell_index_sh,shell->type1_sh[shell_index_sh]);
             for (g = 0; g < shell->ng_sh[shell_index_sh]; g++) {
                gaussian->expo_sh[g_index_sh] = gaussian->expo[g_index];
                gaussian->c_sh[g_index_sh] = gaussian->dc[g_index];
                //fprintf(file.out,"shell type %d %d   %d %d gauss %d %d    %e %e   %e  %e\n",shell->type[shell_index],\
                shell->type_sh[shell_index_sh],shell_index,shell_index_sh,g_index,g_index_sh,gaussian->expo[g_index],\
                gaussian->expo_sh[g_index_sh],gaussian->dc[g_index],gaussian->c_sh[g_index_sh]);
                if (gaussian->expo_sh[g_index_sh] < shell->min_expo_sh[shell_index_sh]) {
                shell->min_expo_sh[shell_index_sh] = gaussian->expo_sh[g_index_sh];
                shell->min_coef_sh[shell_index_sh] = pow(2.0 * gaussian->expo_sh[g_index_sh] / pi, 0.75);
               }
                g_index++;
                g_index_sh++;
             }
               bfn_index_sh += 5;
               shl_index_sh++;
               atoms->number_of_sh_bfns_in_unit_cell += 5;
           break;
          case 10:
           shell->imax_sh[shell_index_sh] = 3;
           shell->type_sh[shell_index_sh] = 7;
           shell->type1_sh[shell_index_sh] = 10;
           shell->cart[shell_index_sh] = 10;
           shell->shar[shell_index_sh] = 7;
           shell->ord_sh[shell_index_sh] = 4;
           shell->ng_sh[shell_index_sh] = shell->ng[shell_index];
           shell->opp_sh[shell_index_sh] = 12;
           //fprintf(file.out,"shells->type1_sh[%2d] %2d\n",shell_index_sh,shell->type1_sh[shell_index_sh]);
             for (g = 0; g < shell->ng_sh[shell_index_sh]; g++) {
                gaussian->expo_sh[g_index_sh] = gaussian->expo[g_index];
                gaussian->c_sh[g_index_sh] = gaussian->fc[g_index];
                //fprintf(file.out,"shell type %d %d   %d %d gauss %d %d    %e %e   %e  %e\n",shell->type[shell_index],\
                shell->type_sh[shell_index_sh],shell_index,shell_index_sh,g_index,g_index_sh,gaussian->expo[g_index],\
                gaussian->expo_sh[g_index_sh],gaussian->fc[g_index],gaussian->c_sh[g_index_sh]);
                if (gaussian->expo_sh[g_index_sh] < shell->min_expo_sh[shell_index_sh]) {
                shell->min_expo_sh[shell_index_sh] = gaussian->expo_sh[g_index_sh];
                shell->min_coef_sh[shell_index_sh] = pow(2.0 * gaussian->expo_sh[g_index_sh] / pi, 0.75);
               }
                g_index++;
                g_index_sh++;
             }
               bfn_index_sh += 7;
               shl_index_sh++;
               atoms->number_of_sh_bfns_in_unit_cell += 7;
           break;
          case 15:
           shell->imax_sh[shell_index_sh] = 4;
           shell->type_sh[shell_index_sh] = 9;
           shell->type1_sh[shell_index_sh] = 15;
           shell->cart[shell_index_sh] = 15;
           shell->shar[shell_index_sh] = 9;
           shell->ord_sh[shell_index_sh] = 5;
           shell->ng_sh[shell_index_sh] = shell->ng[shell_index];
           shell->opp_sh[shell_index_sh] = 28;
           //fprintf(file.out,"shells->type1_sh[%2d] %2d\n",shell_index_sh,shell->type1_sh[shell_index_sh]);
             for (g = 0; g < shell->ng_sh[shell_index_sh]; g++) {
                gaussian->expo_sh[g_index_sh] = gaussian->expo[g_index];
                gaussian->c_sh[g_index_sh] = gaussian->gc[g_index];
                if (gaussian->expo_sh[g_index_sh] < shell->min_expo_sh[shell_index_sh]) {
                shell->min_expo_sh[shell_index_sh] = gaussian->expo_sh[g_index_sh];
                shell->min_coef_sh[shell_index_sh] = pow(2.0 * gaussian->expo_sh[g_index_sh] / pi, 0.75);
               }
                g_index++;
                g_index_sh++;
             }
               bfn_index_sh += 9;
               shl_index_sh++;
               atoms->number_of_sh_bfns_in_unit_cell += 9;
           break;
      }

           shell_index++;
           shell_index_sh++;

        } // close loop over k

      } // close if (at..

    } // close loop over j

    atoms->bfnnumb_sh[i] = bfn_index_sh - atoms->bfnposn_sh[i]; // number of basis functions on ith atom
    if (job->taskid == 0 && job->verbosity > 1)
    fprintf(file.out,"Number of shar basis functions %4d on atom %4d bfnposn_sh %4d\n",atoms->bfnnumb_sh[i],i,atoms->bfnposn_sh[i]);

  }

      if (job->taskid == 0 && atoms->number_of_basis_sets > 0) {

      fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
      fprintf(file.out,"| BASIS SETS                | BASIS FUNCTIONS    %4d | UNIQUE SHELLS      %4d | UNIQUE GAUSSIANS   %4d |\n", \
      atoms->number_of_sh_bfns_in_unit_cell,atoms->number_of_sh_shells_in_unit_cell, atoms->number_of_sh_gaussians_in_unit_cell);
      fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
      int taken[atoms->number_of_basis_sets];
      char shell_type[2];
      for (i = 0; i < atoms->number_of_basis_sets; i++) 
      taken[i] = -1;
      for (i = 0; i < atoms->number_of_atoms_in_unit_cell; i++) {
        for (j = 0; j < atoms->number_of_basis_sets; j++) {
          if ((atoms->atomic_number[i] == iratom[j].atomic_number % 300 || atoms->atomic_number[i] == iratom[j].atomic_number % 200)
              && taken[j] == -1)  {
            shell_index_sh = atoms->shelposn_sh[i];
            g_index_sh     = atoms->gausposn_sh[i];
            for (k = 0; k < atoms->nshel_sh[i]; k++) {
              switch (shell->type_sh[k + atoms->shelposn_sh[i]]) {
                case 1:
                sprintf(shell_type,"s");
                break;
                case 3:
                sprintf(shell_type,"p");
                break;
                case 5:
                sprintf(shell_type,"d");
                break;
                case 7:
                sprintf(shell_type,"f");
                break;
                case 9:
                sprintf(shell_type,"g");
                break;
               }
              for (g = 0; g < shell->ng_sh[shell_index_sh]; g++) {
                if (k == 0 && g == 0)
                fprintf(file.out,"| ATOMIC NUMBER %3d BASIS   | %1s  %11.4lf %8.4lf                                                     |\n", \
                atoms->atomic_number[i],shell_type,gaussian->expo_sh[g_index_sh],gaussian->c_sh[g_index_sh]);
                else
                fprintf(file.out,"|                           | %1s  %11.4lf %8.4lf                                                     |\n", \
                shell_type,gaussian->expo_sh[g_index_sh],gaussian->c_sh[g_index_sh]);
                g_index_sh++;
                if (k < atoms->nshel_sh[i] - 1 && g == shell->ng_sh[shell_index_sh] - 1)
                fprintf(file.out,"|                           |                                                                             |\n");
               }
              shell_index_sh++;
             }
            taken[j] = 0;
            fprintf(file.out,"-----------------------------------------------------------------------------------------------------------\n");
           }
          }
         }
        fprintf(file.out,"\n");

    }

    for (i = 0; i < 5; i++) {
      for (j = 0; j < 15; j++) {
        for (k = 0; k < 3; k++) {
        shell->tuv[i][j][k] = 0;
        }
      }
    }

    
    shell->tuv[1][0][0] = 1; // l = 1, x, x
    shell->tuv[1][1][1] = 1; // l = 1, y, y
    shell->tuv[1][2][2] = 1; // l = 1, z, z

    shell->tuv[2][0][0] = 2; // l = 2, xx, xx
    shell->tuv[2][1][0] = 1; // l = 2, xy, x
    shell->tuv[2][2][0] = 1; // l = 2, xz, x
    shell->tuv[2][1][1] = 1; // l = 2, xy, y
    shell->tuv[2][3][1] = 2; // l = 2, yy, yy
    shell->tuv[2][4][1] = 1; // l = 2, yz, y
    shell->tuv[2][2][2] = 1; // l = 2, xz, z
    shell->tuv[2][4][2] = 1; // l = 2, yz, z
    shell->tuv[2][5][2] = 2; // l = 2, zz, zz

    shell->tuv[3][0][0] = 3; // l = 3, xxx, xxx
    shell->tuv[3][1][0] = 2; // l = 3, xxy, xx
    shell->tuv[3][2][0] = 2; // l = 3, xxz, xx
    shell->tuv[3][3][0] = 1; // l = 3, xyy, x
    shell->tuv[3][4][0] = 1; // l = 3, xyz, x
    shell->tuv[3][5][0] = 1; // l = 3, xzz, x
    shell->tuv[3][1][1] = 1; // l = 3, xxy, y
    shell->tuv[3][3][1] = 2; // l = 3, xyy, yy
    shell->tuv[3][4][1] = 1; // l = 3, xyz, y
    shell->tuv[3][6][1] = 3; // l = 3, yyy, yyy
    shell->tuv[3][7][1] = 2; // l = 3, yyz, yy
    shell->tuv[3][8][1] = 1; // l = 3, yzz, y
    shell->tuv[3][2][2] = 1; // l = 3, xxz, z
    shell->tuv[3][4][2] = 1; // l = 3, xyz, z
    shell->tuv[3][5][2] = 2; // l = 3, xzz, zz
    shell->tuv[3][7][2] = 1; // l = 3, yyz, z
    shell->tuv[3][8][2] = 2; // l = 3, yzz, zz
    shell->tuv[3][9][2] = 3; // l = 3, zzz, zzz

    shell->tuv[4][0][0]  = 4; // l = 4, xxxx, xxxx
    shell->tuv[4][1][0]  = 3; // l = 4, xxxy, xxx
    shell->tuv[4][2][0]  = 3; // l = 4, xxxz, xxx
    shell->tuv[4][3][0]  = 2; // l = 4, xxyy, xx
    shell->tuv[4][4][0]  = 2; // l = 4, xxyz, xx
    shell->tuv[4][5][0]  = 2; // l = 4, xxzz, xx
    shell->tuv[4][6][0]  = 1; // l = 4, xyyy, x
    shell->tuv[4][7][0]  = 1; // l = 4, xyyz, x
    shell->tuv[4][8][0]  = 1; // l = 4, xyzz, x
    shell->tuv[4][9][0]  = 1; // l = 4, xzzz, x
    shell->tuv[4][1][1]  = 1; // l = 4, xxxy, y
    shell->tuv[4][3][1]  = 2; // l = 4, xxyy, yy
    shell->tuv[4][4][1]  = 1; // l = 4, xxyz, y
    shell->tuv[4][6][1]  = 3; // l = 4, xyyy, yyy
    shell->tuv[4][7][1]  = 2; // l = 4, xyyz, yy
    shell->tuv[4][8][1]  = 1; // l = 4, xyzz, y
    shell->tuv[4][10][1] = 4; // l = 4, yyyy, yyyy
    shell->tuv[4][11][1] = 3; // l = 4, yyyz, yyy
    shell->tuv[4][12][1] = 2; // l = 4, yyzz, yy
    shell->tuv[4][13][1] = 1; // l = 4, yzzz, y
    shell->tuv[4][2][2]  = 1; // l = 4, xxxz, z
    shell->tuv[4][4][2]  = 1; // l = 4, xxyz, z
    shell->tuv[4][5][2]  = 2; // l = 4, xxzz, zz
    shell->tuv[4][7][2]  = 1; // l = 4, xyyz, z
    shell->tuv[4][8][2]  = 2; // l = 4, xyzz, zz
    shell->tuv[4][9][2]  = 3; // l = 4, xzzz, zzz
    shell->tuv[4][11][2] = 1; // l = 4, yyyz, z
    shell->tuv[4][12][2] = 2; // l = 4, yyzz, zz
    shell->tuv[4][13][2] = 3; // l = 4, yzzz, zzz
    shell->tuv[4][14][2] = 4; // l = 4, zzzz, zzzz

    shell->ind_i[0] = 0; // s, sp or p
    shell->ind_i[1] = 1;
    shell->ind_i[2] = 2;
    shell->ind_i[3] = 3;

    shell->ind_i[4] = 0; // d
    shell->ind_i[5] = 3;
    shell->ind_i[6] = 5;
    shell->ind_i[7] = 2;
    shell->ind_i[8] = 4;
    shell->ind_i[9] = 0;
    shell->ind_i[10] = 3;
    shell->ind_i[11] = 1;

    shell->ind_i[12] = 1; // f 
    shell->ind_i[13] = 6;  
    shell->ind_i[14] = 4;  
    shell->ind_i[15] = 1;  
    shell->ind_i[16] = 6;  
    shell->ind_i[17] = 8;  
    shell->ind_i[18] = 2;  
    shell->ind_i[19] = 7;  
    shell->ind_i[20] = 9;  
    shell->ind_i[21] = 0;  
    shell->ind_i[22] = 3;  
    shell->ind_i[23] = 5;  
    shell->ind_i[24] = 2;  
    shell->ind_i[25] = 7;  
    shell->ind_i[26] = 0;  
    shell->ind_i[27] = 3;  

    shell->ind_i[28] = 1; // g
    shell->ind_i[29] = 6; 
    shell->ind_i[30] = 4; 
    shell->ind_i[31] = 11; 
    shell->ind_i[32] = 1; 
    shell->ind_i[33] = 6; 
    shell->ind_i[34] = 8; 
    shell->ind_i[35] = 4; 
    shell->ind_i[36] = 11; 
    shell->ind_i[37] = 13; 
    shell->ind_i[38] = 0; 
    shell->ind_i[39] = 3; 
    shell->ind_i[40] = 5; 
    shell->ind_i[41] = 10; 
    shell->ind_i[42] = 12; 
    shell->ind_i[43] = 14; 
    shell->ind_i[44] = 2; 
    shell->ind_i[45] = 7; 
    shell->ind_i[46] = 9; 
    shell->ind_i[47] = 0; 
    shell->ind_i[48] = 5; 
    shell->ind_i[49] = 10; 
    shell->ind_i[50] = 12; 
    shell->ind_i[51] = 2; 
    shell->ind_i[52] = 7; 
    shell->ind_i[53] = 0; 
    shell->ind_i[54] = 3; 
    shell->ind_i[55] = 10; 

    shell->ind_j[0] = 0; // s, sp or p
    shell->ind_j[1] = 1;
    shell->ind_j[2] = 2;
    shell->ind_j[3] = 3;

    shell->ind_j[4] = 0; // d
    shell->ind_j[5] = 0;
    shell->ind_j[6] = 0;
    shell->ind_j[7] = 1;
    shell->ind_j[8] = 2;
    shell->ind_j[9] = 3;
    shell->ind_j[10] = 3;
    shell->ind_j[11] = 4;

    shell->ind_j[12] = 0; // f 
    shell->ind_j[13] = 0; 
    shell->ind_j[14] = 1; 
    shell->ind_j[15] = 2; 
    shell->ind_j[16] = 2; 
    shell->ind_j[17] = 2; 
    shell->ind_j[18] = 3; 
    shell->ind_j[19] = 3; 
    shell->ind_j[20] = 3; 
    shell->ind_j[21] = 4; 
    shell->ind_j[22] = 4; 
    shell->ind_j[23] = 4; 
    shell->ind_j[24] = 5; 
    shell->ind_j[25] = 5; 
    shell->ind_j[26] = 6; 
    shell->ind_j[27] = 6; 

    shell->ind_j[28] = 0; // g
    shell->ind_j[29] = 0; 
    shell->ind_j[30] = 1; 
    shell->ind_j[31] = 1; 
    shell->ind_j[32] = 2; 
    shell->ind_j[33] = 2; 
    shell->ind_j[34] = 2; 
    shell->ind_j[35] = 3; 
    shell->ind_j[36] = 3; 
    shell->ind_j[37] = 3; 
    shell->ind_j[38] = 4; 
    shell->ind_j[39] = 4; 
    shell->ind_j[40] = 4; 
    shell->ind_j[41] = 4; 
    shell->ind_j[42] = 4; 
    shell->ind_j[43] = 4; 
    shell->ind_j[44] = 5; 
    shell->ind_j[45] = 5; 
    shell->ind_j[46] = 5; 
    shell->ind_j[47] = 6; 
    shell->ind_j[48] = 6; 
    shell->ind_j[49] = 6; 
    shell->ind_j[50] = 6; 
    shell->ind_j[51] = 7; 
    shell->ind_j[52] = 7; 
    shell->ind_j[53] = 8; 
    shell->ind_j[54] = 8; 
    shell->ind_j[55] = 8; 

    shell->num_ij[0] = 1;
    shell->num_ij[1] = 4;
    shell->num_ij[2] = 3;
    shell->num_ij[3] = 8;
    shell->num_ij[4] = 16;
    shell->num_ij[5] = 28;

    shell->rot[0] = k_one; // s, sp or p
    shell->rot[1] = k_one;
    shell->rot[2] = k_one;
    shell->rot[3] = k_one;

    shell->rot[4] = -k_one / two / sqrt(three); // d
    shell->rot[5] = -k_one / two / sqrt(three);
    shell->rot[6] = +k_one / sqrt(three);
    shell->rot[7] = +k_one;
    shell->rot[8] = +k_one;
    shell->rot[9] = +k_one / two;
    shell->rot[10] = -k_one / two;
    shell->rot[11] = +k_one;

    shell->rot[12] = +three/two/sqrt(six);     // f 3,-3
    shell->rot[13] = -k_one/two/sqrt(six);  
    shell->rot[14] = +k_one;                   // f 3,-2
    shell->rot[15] = -k_one/two/sqrt(ten);     // f 3,-1
    shell->rot[16] = -k_one/two/sqrt(ten); 
    shell->rot[17] = +four/two/sqrt(ten); 
    shell->rot[18] = -three/two/sqrt(fifteen); // f 3,0
    shell->rot[19] = -three/two/sqrt(fifteen); 
    shell->rot[20] = +k_one/sqrt(fifteen); 
    shell->rot[21] = -k_one/two/sqrt(ten);     // f 3,1
    shell->rot[22] = -k_one/two/sqrt(ten); 
    shell->rot[23] = +four/two/sqrt(ten); 
    shell->rot[24] = +k_one/two;               // f 3,2
    shell->rot[25] = -k_one/two;          
    shell->rot[26] = +k_one/two/sqrt(six);     // f 3,3
    shell->rot[27] = -three/two/sqrt(six);

    shell->rot[28] = +k_one/sqrt(four * rtthree);         // g 4,-4
    shell->rot[29] = -k_one/sqrt(four * rtthree);     
    shell->rot[30] = +three/sqrt(eight * rtthree);        // g 4,-3
    shell->rot[31] = -k_one/sqrt(eight * rtthree); 
    shell->rot[32] = -k_one/sqrt(twenty_eight * rtthree); // g 4,-2
    shell->rot[33] = -k_one/sqrt(twenty_eight * rtthree);  
    shell->rot[34] = +six/sqrt(twenty_eight * rtthree);  
    shell->rot[35] = -three/sqrt(fifty_six * rtthree);    // g 4,-1
    shell->rot[36] = -three/sqrt(fifty_six * rtthree);  
    shell->rot[37] = +four/sqrt(fifty_six * rtthree);  
    shell->rot[38] = +three/sqrt(2240.0 * rtthree);       // g 4,0
    shell->rot[39] = +six/sqrt(2240.0 * rtthree);   
    shell->rot[40] = -twenty_four/sqrt(2240.0 * rtthree);  
    shell->rot[41] = +three/sqrt(2240.0 * rtthree); 
    shell->rot[42] = -twenty_four/sqrt(2240.0 * rtthree);  
    shell->rot[43] = +eight/sqrt(2240.0 * rtthree); 
    shell->rot[44] = -three/sqrt(fifty_six * rtthree);    // g 4,1
    shell->rot[45] = -three/sqrt(fifty_six * rtthree);    
    shell->rot[46] = +four/sqrt(fifty_six * rtthree);    
    shell->rot[47] = -k_one/sqrt(112.0 * rtthree);        // g 4,2
    shell->rot[48] = +six/sqrt(112.0 * rtthree);    
    shell->rot[49] = +k_one/sqrt(112.0 * rtthree);  
    shell->rot[50] = -six/sqrt(112.0 * rtthree);    
    shell->rot[51] = +k_one/sqrt(eight * rtthree);        // g 4,3
    shell->rot[52] = -three/sqrt(eight * rtthree);  
    shell->rot[53] = +k_one/eight/sqrt(rtthree);          // g 4,4
    shell->rot[54] = -six/eight/sqrt(rtthree);
    shell->rot[55] = +k_one/eight/sqrt(rtthree);

}
