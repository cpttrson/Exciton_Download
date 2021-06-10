#ifndef LIMITSH
#define LIMITSH

#define lower_energy_limit             -5.0E+03
#define upper_energy_limit             +5.0E+03
#define initial_fermi_level_guess      -1.0E+04
//#define integral_rejection_threshold    1.0E-14
//#define integral_rejection_threshold    1.0E-12
#define integral_rejection_threshold    1.0E-10
//#define Fock_matrix_rejection_threshold 1.0E-11
#define Fock_matrix_rejection_threshold 1.0E-13
#define real_space_cutoff               1.0E-18
#define reciprocal_space_cutoff         1.0E-16
//#define real_space_cutoff               1.0E-24
//#define reciprocal_space_cutoff         1.0E-22
#define slab_real_space_z_period        2.0E+01
#define rod_real_space_x_y_period       2.0E+01

#endif
