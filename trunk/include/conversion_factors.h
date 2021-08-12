
  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 C. H. Patterson                             *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

/* This file define conversion factors to CGS units    *
 *--- All conversion factors are multiplicative, e.g.
 *--- 3eV * eV_to_erg = 3*(1.60219d-12) = 4.806d-12 erg
 *--- Atomic units of energy are Hartree.
 *--- Factors derived from units in IOP diary...!   */

#define ev_to_erg  (double)(1.60217733E-12)
#define ev_to_freq (double)(1.519215)
#define au_to_eV   (double)(27.211396130)
#define bohr_to_AA (double)(0.5291772083)

