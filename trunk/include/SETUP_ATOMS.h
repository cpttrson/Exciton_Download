
  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 C. H. Patterson                             *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

void count_unique_atoms(ATOM*, JOB_PARAM*, FILES);

void read_unique_atoms(ATOM*, ATOMTYPE*, CRYSTAL*, JOB_PARAM*, FILES);

void count_all_atoms(CRYSTAL*, ATOM*, ATOMTYPE*, SYMMETRY*, JOB_PARAM*, FILES);

void generate_all_atoms(CRYSTAL*, REAL_LATTICE*, SYMMETRY*, ATOMTYPE*, ATOMTYPE*, ATOM*, ATOM_TRAN*, ATOM_TRAN*, JOB_PARAM*, FILES);

void generate_ATOM_TRAN(CRYSTAL*, REAL_LATTICE*, SYMMETRY*, ATOM*, ATOM_TRAN*, JOB_PARAM*, FILES);

void count_all_crystal_atoms(ATOM*, JOB_PARAM*, FILES);

void generate_all_crystal_atoms(CRYSTAL*, REAL_LATTICE*, SYMMETRY*, ATOMTYPE*, ATOM*, ATOM_TRAN*, ATOM_TRAN*, JOB_PARAM*, FILES);

void count_basis_sets(ATOM*, int, JOB_PARAM*, FILES);

void read_basis_sets(ATOM*, int, ATOMTYPE*, JOB_PARAM*, FILES);

void generate_runtime_arrays(ATOM*, ATOM_TRAN*, ATOMTYPE*, SHELL*, GAUSSIAN*, JOB_PARAM*, FILES);

