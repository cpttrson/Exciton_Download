
  // ******************************************************************************************
  //                                                                                          *
  //                           Copyright (C) 2021 C. H. Patterson                             *
  //                                                                                          *
  //  This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.    *
  //  If a copy of the MPL was not distributed with this file, you can obtain one at          *
  //  http://mozilla.org/MPL/2.0/.                                                            *
  //                                                                                          *
  // ******************************************************************************************

void mcmurchie_davidson_ij(double*,int,int,int,int,int,double*,double*,double*,double*,double*,double*,double*,SHELL*,JOB_PARAM*,FILES);

void mcmurchie_davidson_ij_complex(Complex*,Complex*,int,int,int,int,int,double*,double*,double*,double*,double*,double*,SHELL*,JOB_PARAM*,FILES);

void mcmurchie_davidson_ija(double*, int, int, int, int, int, int, int, int, double*, double*, double*, double*, double*, double*, double*, SHELL*, SHELL*, JOB_PARAM*, FILES);

void mcmurchie_davidson_ija_complex(Complex*, Complex*, int, int, int, int, int, int, int, int, double*, double*, double*, double*, double*, double*, SHELL*, SHELL*, JOB_PARAM*, FILES);

void mcmurchie_davidson_ijkl(double*, int, int, int, int,  double*, double*, double*, double*, double*, double*, double*, SHELL*, JOB_PARAM*, FILES);

void mcmurchie_davidson_ijij(double*, int, int, int, int, int, double*, double*, double*, double*, SHELL*, JOB_PARAM*, FILES);

void mcmurchie_davidson_ijij_complex(Complex*, int, int, int, int, int, double*, double*, double*, Complex*, SHELL*, JOB_PARAM*, FILES);

