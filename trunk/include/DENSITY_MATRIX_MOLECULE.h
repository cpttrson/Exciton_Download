
  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 C. H. Patterson                             *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

void initial_density_matrix(double*, double*, PAIR_TRAN*, FERMI*, ATOM*, ATOM_TRAN*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void read_density_matrix(FERMI*,double**,int*,int*,ATOM*,JOB_PARAM*,FILES);

void read_JOB_PARAM_array(FERMI*, int*, double*, JOB_PARAM*, FILES);

void write_JOB_PARAM_array(FERMI*, int*, double*, JOB_PARAM*, FILES);

void reduced_density_matrix_molecule(int*, double*, ComplexMatrix*, FERMI*, double*, KPOINT_TRAN*, int*, REAL_LATTICE*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void reduced_density_matrix_molecule_mpp(int, int*, double*, ComplexMatrix*, FERMI*, double*, KPOINT_TRAN*, int*, REAL_LATTICE*, PAIR_TRAN*, ATOM_TRAN*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void expand_density_matrix(double*, double*, PAIR_TRAN*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

void print_density_matrix_molecule(int, double*, PAIR_TRAN*, ATOM*, JOB_PARAM*, FILES);

void atom_shell_populations2(INT_1E*, double*, PAIR_TRAN*, ATOM*, SHELL*, SYMMETRY*, JOB_PARAM*, FILES);

double print_total_population2(ATOM*, SHELL*, JOB_PARAM*, FILES);

void print_atom_populations2(ATOM*, SHELL*, JOB_PARAM*, FILES);

void print_shell_populations2(ATOM*, SHELL*, JOB_PARAM*, FILES);

void init_mix(DoubleMatrix*, DoubleMatrix*, DoubleMatrix*, double*, double*);

void mix_density(DoubleMatrix*, DoubleMatrix*, DoubleMatrix*, double*, JOB_PARAM*, FILES);

