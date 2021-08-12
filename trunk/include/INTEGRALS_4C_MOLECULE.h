
  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 C. H. Patterson                             *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

void integrals_molecule_ijkl(INTEGRAL_LIST*, int, int, int, int, int*, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, SYMMETRY*, JOB_PARAM*, FILES);

void integrals_molecule_ijij(double*, int, int, REAL_LATTICE*, ATOM*, SHELL*, GAUSSIAN*, SYMMETRY*, CRYSTAL*, JOB_PARAM*, FILES);

void pack_write_molecule_2e_integrals(INTEGRAL_LIST*, QUAD_TRAN*, FILE*, JOB_PARAM*, FILES);

void read_unpack_molecule_2e_integrals(INTEGRAL_LIST*, QUAD_TRAN*, FILE*, JOB_PARAM*, FILES);

