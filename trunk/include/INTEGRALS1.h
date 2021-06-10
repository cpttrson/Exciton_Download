#ifndef INTEGRALS1H
#define INTEGRALS1H

void fock_element_1e(INT_1E *, int, PAIR_TRAN *, int [], REAL_LATTICE *R, RECIPROCAL_LATTICE *G, ATOM *, SHELL *, GAUSSIAN *, CRYSTAL *, JOB_PARAM *, FILES);

double cosfactor(int, double);

double e(int, int, int, double, double, double);

#endif
