#ifndef MYCONSTANTSH
#define MYCONSTANTSH
/* Fundamental constants in MKS units

 hbar = Planck constant /2 / pi
 m0 = electron rest mass in kg
 el = electronic charge in coulomb
 a0 = Bohr radius in m
 */

#define pi   (double)(3.14159265358979323846)

#define rtpi (double)(1.77245385090551588191)

#define pi32 (double)(5.5683279968317078454)

#define deg_rad (double)(pi/180.000000000000)

#define hbar (double)(1.054571628E-34)

#define m0 (double)(9.10938215E-31)

#define epsilon_0 (double)(8.854187817E-12)

#define el (double)(1.602176487E-19)

#define a0 (double)(5.2917720859E-11)

//#define zero           (double) 0.00000000000000000000
//#define one            (double) 1.00000000000000000000
#define k_tol          (double) 0.000000000000001000000000
//#define k_tol          (double) 0.00000000001000000000
//#define ktol           (double) 0.00000000010000000000
#define k_zero         (double) 0.00000000000000000000
#define k_one          (double) 1.00000000000000000000
#define two            (double) 2.00000000000000000000
#define three          (double) 3.00000000000000000000
#define four           (double) 4.00000000000000000000
#define five           (double) 5.00000000000000000000
#define six            (double) 6.00000000000000000000
#define seven          (double) 7.00000000000000000000
#define eight          (double) 8.00000000000000000000
#define nine           (double) 9.00000000000000000000
#define ten            (double) 10.0000000000000000000
#define twelve         (double) 12.0000000000000000000
#define fourteen       (double) 14.0000000000000000000
#define fifteen        (double) 15.0000000000000000000
#define sixteen        (double) 16.0000000000000000000
#define eighteen       (double) 18.0000000000000000000
#define twenty_two     (double) 22.0000000000000000000
#define twenty_four    (double) 24.0000000000000000000
#define twenty_eight   (double) 28.0000000000000000000
#define fifty_six      (double) 56.0000000000000000000
#define half           (double) 0.50000000000000000000
#define third          (double) 0.33333333333333333333
#define fourth         (double) 0.25000000000000000000
#define fifth          (double) 0.20000000000000000000
#define sixth          (double) 0.16666666666666666667
#define quarter        (double) 0.25000000000000000000
#define two_thirds     (double) 0.66666666666666666667
#define five_sixths    (double) 0.83333333333333333333
#define three_quarters (double) 0.75000000000000000000
#define rttwo          (double) 1.41421356237309504880
#define rtthree        (double) 1.73205080756887729352

/* BZ integration */

#define INT_FCC (double)(13.897645562159250)

#define INT_SC (double)(15.672495234738540)

#define INT_BCC (double)(43.198066515915080)

#endif
