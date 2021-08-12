
  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 C. H. Patterson                             *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

void read_symmetry_group(CRYSTAL*, SYMMETRY*, int, JOB_PARAM*, FILES);

void generate_cartesian_symmetry_operators(CRYSTAL*, SYMMETRY*, JOB_PARAM*, FILES);

void reorder_symmetry_operators_by_class(SYMMETRY*, JOB_PARAM*, FILES);

void generate_operator_inverses(SYMMETRY*, JOB_PARAM*, FILES);

void generate_group_multiplication_table(SYMMETRY*, int, FILES);

void generate_group_conjugacy_classes(SYMMETRY*, int, FILES);

void generate_Dirac_characters(SYMMETRY*, JOB_PARAM*, FILES);

void generate_Dirac_characters_complex(SYMMETRY*, JOB_PARAM*, FILES);

void generate_permutation_group_table(SYMMETRY*, JOB_PARAM*, FILES);

void print_symmetry_operators(SYMMETRY*, JOB_PARAM*, FILES);

void count_little_k_group_operators(int, SYMMETRY*, SYMMETRY*, CRYSTAL*, KPOINT_TRAN*, FERMI*, JOB_PARAM*, FILES);

void generate_little_k_group_operstors(int, SYMMETRY*, SYMMETRY*, CRYSTAL*, KPOINT_TRAN*, FERMI*, JOB_PARAM*, FILES);

void generate_little_k_group(int, SYMMETRY*, FERMI*, KPOINT_TRAN*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

