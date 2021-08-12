
  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 C. H. Patterson                             *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

void generate_integral_buffers_molecule_ija(double*, double*, int*, int*, FERMI*, ATOM*, ATOM_TRAN*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*,CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void contract_integrals_molecule_ija(int*, MPI_File, double*, FERMI*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void generate_coulomb_matrix_inverse(DoubleMatrix*, ATOM_TRAN*, ATOM*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void contract_coulomb_integrals(int*, int*, DoubleMatrix*, double*, double*, FERMI*, ATOM_TRAN*, ATOM*, SHELL*, GAUSSIAN*, ATOM*, SHELL*, GAUSSIAN*, CRYSTAL*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, RECIPROCAL_LATTICE*, JOB_PARAM*, FILES);

void initialise_ring_rotate(int, int*, int*, int*, JOB_PARAM*, FILES);

void ring_rotate_once(double**, int, int, int*, int*, int*, JOB_PARAM*, FILES);

void ring_transfer(double*, double*, int, int, int, JOB_PARAM*, FILES);

