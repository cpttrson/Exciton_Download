
  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 C. H. Patterson                             *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

void count_pairs4(PAIR_TRAN*, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE*, JOB_PARAM*, FILES);

void select_pairs4(PAIR_TRAN*, ATOM*, REAL_LATTICE*, JOB_PARAM*, FILES);

void count_density_pairs4(PAIR_TRAN*, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE *R, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void count_range_selected_pairs(PAIR_TRAN*, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE *R, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void generate_density_pairs4(PAIR_TRAN*, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void generate_range_selected_pairs(PAIR_TRAN*, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void count_triples1(TRIPLE_TRAN*, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void generate_triples1(TRIPLE_TRAN*, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void count_triples1_reversed(int, TRIPLE_TRAN*, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void generate_triples1_reversed(int, TRIPLE_TRAN*, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void generate_molecule_quads(PAIR_TRAN*, QUAD_TRAN* ,int, int, int, int, int, int, int, int, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void generate_c_quads(PAIR_TRAN*, QUAD_TRAN*, int, int, int, int, int, int, int, int, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void generate_e_quads(PAIR_TRAN*, QUAD_TRAN*, int, int, int, int, int, int, int, int, ATOM*, ATOM_TRAN*, SYMMETRY*, REAL_LATTICE_TABLES*, JOB_PARAM*, FILES);

void print_pairs(PAIR_TRAN*, ATOM*, REAL_LATTICE*, JOB_PARAM*, FILES);

