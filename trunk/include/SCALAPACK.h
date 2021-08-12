
  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 C. H. Patterson                             *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

void initialise_spk_grid(int*, int*, int*, JOB_PARAM*, FILES);

void initialise_spk_grid_crystal(int*, FERMI*, int, int*, int*, int*, JOB_PARAM*, FILES);

void setup_procs(int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, FERMI*, FILES, JOB_PARAM*);

void block_cyclic_to_linear(int*, int*, int*, double*, char*, JOB_PARAM*, FILES);

void block_cyclic_to_linear_limit(int*, int*, int*, int, double*, char*, JOB_PARAM*, FILES);

void block_cyclic_to_linear_limit_complex(int*, int*, int*, int, Complex*, char*, JOB_PARAM*, FILES);

void block_cyclic_zero_triangle(char*, int*, int*, int*, double*);

