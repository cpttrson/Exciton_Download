#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>
#include <ctime>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <xc.h>
#include <sys/types.h>

#include "mylogical.h"
#include "mycomplex.h"
#include "conversion_factors.h"
#include "USER_DATA.h"
#include "MATRIX_UTIL.h"
#include "PRINT_UTIL.h"
#include "LINEAR_ALGEBRA_UTIL.h"
#include "myconstants.h"
#include "ALLOCATE_MEMORY.h"
#include "TOOLS.h"
#include "INTEGRALS1.h"
#include "LEBEDEV_LAIKOV.h"
//#include "PSEUDOPOTENTIAL.h"
#include "INCOMPLETE_GAMMA.h"
#include "SCF_ATOM.h"

using namespace std;

void orbital_grid(double *grid_x, double *grid_y, double *grid_z, double *grid_r, double *grid_basis, int *atm, int nradial, int nsphere, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, JOB_PARAM *job, FILES file)

{

  int i, j, i4;
  int index_i, sheli, shelposi, gausposi;
  int count, count1;
  double rsqrd, sum[atoms->nshel_sh[*atm]];
 
   count  = 0;
   count1 = 0;

   for (i = 0; i < nradial; i++) {
     rsqrd = grid_r[i] * grid_r[i];
     shelposi = atoms->shelposn_sh[*atm];
     gausposi = atoms->gausposn_sh[*atm];
       for (index_i = 0; index_i < atoms->nshel_sh[*atm]; index_i++) {
         sheli = shells->type_sh[shelposi + index_i];
         sum[index_i] = k_zero;
         for (i4 = 0; i4 < shells->ng_sh[shelposi + index_i]; i4++) {
           sum[index_i] += gaussians->c_sh[gausposi + i4] * exp(-gaussians->expo_sh[gausposi + i4] * rsqrd);
           //fprintf(file.out,"i %d %d %d %d %d %d %lf %lf %lf %lf\n",i,i4,shelposi,gausposi,index_i,sheli,gaussians->c_sh[gausposi + i4],\
           gaussians->expo_sh[gausposi + i4],rsqrd,gaussians->c_sh[gausposi + i4]*exp(-gaussians->expo_sh[gausposi + i4]*rsqrd));
          }
           gausposi += shells->ng_sh[shelposi + index_i];
         }
          for (j = 0; j < nsphere; j++) {
            shelposi = atoms->shelposn_sh[*atm];
            for (index_i = 0; index_i < atoms->nshel_sh[*atm]; index_i++) {
              sheli = shells->type_sh[shelposi + index_i];
              switch (sheli) {
                case 1:
                  grid_basis[count] = sum[index_i];
                  //fprintf(file.out,"grid_basis %d %d %d %lf\n",count,i,j,grid_basis[count]);
                  count++;
                  break;
                case 3:
                  grid_basis[count]     = grid_x[count1] * sum[index_i];
                  grid_basis[count + 1] = grid_y[count1] * sum[index_i];
                  grid_basis[count + 2] = grid_z[count1] * sum[index_i];
                  //fprintf(file.out,"grip_basis %d %d %d %20.15lf\n",count,i,j,grid_basis[count]);
                  //fprintf(file.out,"grip_basis %d %d %d %20.15lf\n",count+1,i,j,grid_basis[count+1]);
                  //fprintf(file.out,"grip_basis %d %d %d %20.15lf\n",count+2,i,j,grid_basis[count+2]);
                  count += 3;
                  break;
                case 5:
                  grid_basis[count] = (grid_z[count1] * grid_z[count1] - \
                 (grid_x[count1] * grid_x[count1] + grid_y[count1] * grid_y[count1]) / two) / sqrt(three) * sum[index_i];
                  grid_basis[count + 1] =  grid_x[count1] * grid_z[count1] * sum[index_i];
                  grid_basis[count + 2] =  grid_y[count1] * grid_z[count1] * sum[index_i];
                  grid_basis[count + 3] = (grid_x[count1] * grid_x[count1] - grid_y[count1] * grid_y[count1]) / two * sum[index_i];
                  grid_basis[count + 4] =  grid_x[count1] * grid_y[count1] * sum[index_i];
                  count += 5;
                  break;
                case 7:
                  if (job->taskid == 0)
                  fprintf(file.out,"f coefficients not entered yet\n");
                  exit(1);
                  break;
              } // close switch
             }
            count1++;
           } // close loop over j
          } // close loop over i

}

void orbital_gradient_grid(double *grid_x, double *grid_y, double *grid_z, double *grid_r, VECTOR_DOUBLE *grad_basis, int *atm, int nradial, int nsphere, ATOM *atoms, SHELL *shells, GAUSSIAN *gaussians, JOB_PARAM *job, FILES file)

{

  int i, j, i4;
  int index_i, sheli, shelposi, gausposi;
  int count, count1;
  double temp, rsqrd, sum[atoms->nshel_sh[*atm]], sum_alpha[atoms->nshel_sh[*atm]];
 
   count  = 0;
   count1 = 0;

   for (i = 0; i < nradial; i++) {
     rsqrd = grid_r[i] * grid_r[i];
     shelposi = atoms->shelposn_sh[*atm];
     gausposi = atoms->gausposn_sh[*atm];
       for (index_i = 0; index_i < atoms->nshel_sh[*atm]; index_i++) {
         sheli = shells->type_sh[shelposi + index_i];
         sum[index_i]       = k_zero;
         sum_alpha[index_i] = k_zero;
         for (i4 = 0; i4 < shells->ng_sh[index_i]; i4++) {
           temp = gaussians->c_sh[gausposi + i4] * exp(-gaussians->expo_sh[gausposi + i4] * rsqrd);
           sum[index_i]       += temp;
           sum_alpha[index_i] += -two * gaussians->expo_sh[gausposi + i4] * temp;
           //fprintf(file.out,"i %d %d %d %d %d %d %lf %lf %lf %lf\n",i,i4,shelposi,gausposi,index_i,sheli,gaussians->c_sh[gausposi + i4],\
           gaussians->expo_sh[gausposi + i4],rsqrd,gaussians->c_sh[gausposi + i4]*exp(-gaussians->expo_sh[gausposi + i4]*rsqrd));
          }
           gausposi += shells->ng_sh[shelposi + index_i];
         }
          for (j = 0; j < nsphere; j++) {
            shelposi = atoms->shelposn_sh[*atm];
            for (index_i = 0; index_i < atoms->nshel_sh[*atm]; index_i++) {
              sheli = shells->type_sh[shelposi + index_i];
              switch (sheli) {
                case 1:
                  grad_basis[count].comp1 = grid_x[count1] * sum_alpha[index_i];
                  grad_basis[count].comp2 = grid_y[count1] * sum_alpha[index_i];
                  grad_basis[count].comp3 = grid_z[count1] * sum_alpha[index_i];
                  //fprintf(file.out,"grad_basis %d %d %d %lf\n",count,i,j,grad_basis[count].comp1);
                  count++;
                  break;
                case 3:
                  grad_basis[count    ].comp1 = grid_x[count1] * grid_x[count1] * sum_alpha[index_i] + sum[index_i];
                  grad_basis[count + 1].comp1 = grid_y[count1] * grid_x[count1] * sum_alpha[index_i];
                  grad_basis[count + 2].comp1 = grid_z[count1] * grid_x[count1] * sum_alpha[index_i];
                  grad_basis[count    ].comp2 = grid_x[count1] * grid_y[count1] * sum_alpha[index_i];
                  grad_basis[count + 1].comp2 = grid_y[count1] * grid_y[count1] * sum_alpha[index_i] + sum[index_i];
                  grad_basis[count + 2].comp2 = grid_z[count1] * grid_y[count1] * sum_alpha[index_i];
                  grad_basis[count    ].comp3 = grid_x[count1] * grid_z[count1] * sum_alpha[index_i];
                  grad_basis[count + 1].comp3 = grid_y[count1] * grid_z[count1] * sum_alpha[index_i];
                  grad_basis[count + 2].comp3 = grid_z[count1] * grid_z[count1] * sum_alpha[index_i] + sum[index_i];
                  //fprintf(file.out,"grap_basis %d %d %d %20.15lf\n",count,i,j,grad_basis[count].comp1);
                  count += 3;
                  break;
                case 5:
                  grad_basis[count    ].comp1 = (grid_z[count1] * grid_z[count1] * grid_x[count1] * sum_alpha[index_i] - \
                 (grid_x[count1] * grid_x[count1] * grid_x[count1] * sum_alpha[index_i] + two * grid_x[count1] * sum[index_i] + \
                  grid_y[count1] * grid_y[count1] * grid_x[count1] * sum_alpha[index_i]) / two) / sqrt(three);
                  grad_basis[count + 1].comp1 = grid_x[count1] * grid_z[count1] * sum_alpha[index_i] + grid_z[count1] * sum_alpha[index_i];
                  grad_basis[count + 2].comp1 = grid_y[count1] * grid_z[count1] * sum_alpha[index_i];
                  grad_basis[count + 3].comp1 = (grid_x[count1] * grid_x[count1] * grid_x[count1] * sum_alpha[index_i] + \
                  two * grid_x[count1] * sum[index_i] - grid_y[count1] * grid_y[count1] * grid_x[count1] * sum_alpha[index_i]) / two;
                  grad_basis[count + 4].comp1 = grid_x[count1] * grid_y[count1] * sum_alpha[index_i] + grid_y[count1] * sum[index_i];
                  grad_basis[count    ].comp2 = (grid_z[count1] * grid_z[count1] * grid_y[count1] * sum_alpha[index_i] - \
                 (grid_x[count1] * grid_x[count1] * grid_y[count1] * sum_alpha[index_i] + \
                  grid_y[count1] * grid_y[count1] * grid_y[count1] * sum_alpha[index_i] + \
                  two * grid_y[count1] * sum[index_i]) / two) / sqrt(three);
                  grad_basis[count + 1].comp2 = grid_x[count1] * grid_z[count1] * sum_alpha[index_i];
                  grad_basis[count + 2].comp2 = grid_y[count1] * grid_z[count1] * sum_alpha[index_i] + grid_y[count1] * sum[index_i];
                  grad_basis[count + 3].comp2 = (grid_x[count1] * grid_x[count1] * grid_x[count1] * sum_alpha[index_i] - \
                  grid_y[count1] * grid_y[count1] * grid_y[count1] * sum_alpha[index_i] - two * grid_y[count1] * sum[index_i]) / two;
                  grad_basis[count + 4].comp2 = grid_x[count1] * grid_y[count1] * sum_alpha[index_i] + grid_x[count1] * sum[index_i];
                  grad_basis[count    ].comp3 = (grid_z[count1] * grid_z[count1] * grid_z[count1] * sum_alpha[index_i] + \
                  two * grid_z[count1] * sum[index_i] - (grid_x[count1] * grid_x[count1] * grid_z[count1] * sum_alpha[index_i] + \
                  grid_y[count1] * grid_y[count1] * grid_z[count1] * sum_alpha[index_i]) / two) / sqrt(three);
                  grad_basis[count + 1].comp3 = grid_x[count1] * grid_z[count1] * sum_alpha[index_i] + grid_x[count1] * sum[index_i];
                  grad_basis[count + 2].comp3 = grid_y[count1] * grid_z[count1] * sum_alpha[index_i] + grid_y[count1] * sum[index_i];
                  grad_basis[count + 3].comp3 = (grid_x[count1] * grid_x[count1] * grid_z[count1] * sum_alpha[index_i] - \
                  grid_y[count1] * grid_y[count1] * grid_z[count1] * sum_alpha[index_i]) / two;
                  grad_basis[count + 4].comp3 = grid_x[count1] * grid_y[count1] * sum_alpha[index_i];
                  count += 5;
                  break;
                case 7:
                  if (job->taskid == 0)
                  fprintf(file.out,"f coefficients not entered yet\n");
                  exit(1);
                  break;
              } // close switch
             }
            count1++;
           } // close loop over j
          } // close loop over i

}

void Gauss_Chebyshev_weights_abscissa(int nradial, double *y, double *radial_weight, JOB_PARAM *job, FILES file)

{

int i;

    switch (nradial) {

       case 1:
        y[0] = -k_one;
        *radial_weight =  k_one;
        break;

       case 4:
        static const double x4[] = {
        9.23879532511286756101e-01,    3.82683432365089771723e-01};
        static const double A4 =       7.85398163397448309628e-01;

        for(i = 0;i < nradial / 2;i++) {
        y[i] = -x4[i];
        y[nradial - i - 1] = x4[i];
        *radial_weight =  A4;
        }
        break;

       case 16:
        static const double x16[] = {
        9.95184726672196886231e-01,    9.56940335732208864931e-01,
        8.81921264348355029715e-01,    7.73010453362736960797e-01,
        6.34393284163645498203e-01,    4.71396736825997648545e-01,
        2.90284677254462367645e-01,    9.80171403295606019957e-02};
        static const double A16 =      1.96349540849362077407e-01;
        for(i = 0;i < nradial / 2;i++) {
        y[i] = -x16[i];
        y[nradial - i - 1] = x16[i];
        *radial_weight =  A16;
        }
        break;

       case 24:
        static const double x24[] = {
        9.97858923238603506725e-01,    9.80785280403230449119e-01,
        9.46930129495105664243e-01,    8.96872741532688303910e-01,
        8.31469612302545237081e-01,    7.51839807478977396426e-01,
        6.59345815100068868419e-01,    5.55570233019602224757e-01,
        4.42288690219001281995e-01,    3.21439465303161580700e-01,
        1.95090322016128267843e-01,    6.54031292301430668151e-02};
        static const double A24 =      1.30899693899574718267e-01;
        for(i = 0;i < nradial / 2;i++) {
        y[i] = -x24[i];
        y[nradial - i - 1] = x24[i];
        *radial_weight =  A24;
        }
        break;

       case 32:
        static const double x32[] = {
        9.98795456205172392701e-01,    9.89176509964780973456e-01,
        9.70031253194543992616e-01,    9.41544065183020778391e-01,
        9.03989293123443331582e-01,    8.57728610000272069893e-01,
        8.03207531480644909799e-01,    7.40951125354959091193e-01,
        6.71558954847018400619e-01,    5.95699304492433343462e-01,
        5.14102744193221726607e-01,    4.27555093430282094315e-01,
        3.36889853392220050703e-01,    2.42980179903263889945e-01,
        1.46730474455361751659e-01,    4.90676743274180142536e-02};
        static const double A32 =      9.81747704246810387035e-02;
        for(i = 0;i < nradial / 2;i++) {
        y[i] = -x32[i];
        y[nradial - i - 1] = x32[i];
        *radial_weight =  A32;
        }
        break;

       case 48:
        static const double x48[] = {
        9.99464587476365644422e-01,    9.95184726672196886231e-01,
        9.86643332084879004723e-01,    9.73876979277333648167e-01,
        9.56940335732208864931e-01,    9.35905926757325700295e-01,
        9.10863824921175818552e-01,    8.81921264348355029715e-01,
        8.49202181526578887653e-01,    8.12846684591615216604e-01,
        7.73010453362736960797e-01,    7.29864072697835657329e-01,
        6.83592302022871280522e-01,    6.34393284163645498203e-01,
        5.82477696867802149179e-01,    5.28067850650367995877e-01,
        4.71396736825997648545e-01,    4.12707029804394737052e-01,
        3.52250047921233506521e-01,    2.90284677254462367645e-01,
        2.27076263034373207589e-01,    1.62895473394588739485e-01,
        9.80171403295606019957e-02,    3.27190828217761420651e-02};
        static const double A48 =      6.54498469497873591334e-02;
        for(i = 0;i < nradial / 2;i++) {
        y[i] = -x48[i];
        y[nradial - i - 1] = x48[i];
        *radial_weight =  A48;
        }
        break;

       case 64:
        static const double x64[] = {
        9.99698818696204220097e-01,    9.97290456678690216132e-01,
        9.92479534598709998180e-01,    9.85277642388941244766e-01,
        9.75702130038528544446e-01,    9.63776065795439866677e-01,
        9.49528180593036667215e-01,    9.32992798834738887737e-01,
        9.14209755703530654630e-01,    8.93224301195515320339e-01,
        8.70086991108711418636e-01,    8.44853565249707073252e-01,
        8.17584813151583696480e-01,    7.88346427626606262036e-01,
        7.57208846506484547589e-01,    7.24247082951466920962e-01,
        6.89540544737066924622e-01,    6.53172842953776764080e-01,
        6.15231590580626845491e-01,    5.75808191417845300763e-01,
        5.34997619887097210678e-01,    4.92898192229784036869e-01,
        4.49611329654606600042e-01,    4.05241314004989870918e-01,
        3.59895036534988148786e-01,    3.13681740398891476651e-01,
        2.66712757474898386336e-01,    2.19101240156869797227e-01,
        1.70961888760301226360e-01,    1.22410675199216198496e-01,
        7.35645635996674235297e-02,    2.45412285229122880321e-02};
        static const double A64 =      4.90873852123405193518e-02;
        for(i = 0;i < nradial / 2;i++) {
        y[i] = -x64[i];
        y[nradial - i - 1] = x64[i];
        *radial_weight =  A64;
        }
        break;

       case 96:
        static const double x96[] = {
        9.99866137909561782863e-01,    9.98795456205172392701e-01,
        9.96655239309180324928e-01,    9.93447779019444395529e-01,
        9.89176509964780973456e-01,    9.83846005927077416097e-01,
        9.77461974943571863372e-01,    9.70031253194543992616e-01,
        9.61561797682961947144e-01,    9.52062677713924257114e-01,
        9.41544065183020778391e-01,    9.30017223684012117049e-01,
        9.17494496447491307935e-01,    9.03989293123443331582e-01,
        8.89516075421856035267e-01,    8.74090341626758851525e-01,
        8.57728610000272069893e-01,    8.40448401094438021037e-01,
        8.22268218989775107855e-01,    8.03207531480644909799e-01,
        7.83286749228650365376e-01,    7.62527203906388096352e-01,
        7.40951125354959091193e-01,    7.18581617779698057173e-01,
        6.95442635009611651129e-01,    6.71558954847018400619e-01,
        6.46956152534857365421e-01,    6.21660573370077408053e-01,
        5.95699304492433343462e-01,    5.69100145878898230586e-01,
        5.41891580574751716148e-01,    5.14102744193221726607e-01,
        4.85763393716340056256e-01,    4.56903875630420676552e-01,
        4.27555093430282094315e-01,    3.97748474527011052041e-01,
        3.67515936594703565403e-01,    3.36889853392220050703e-01,
        3.05903020096553462752e-01,    2.74588618184932341487e-01,
        2.42980179903263889945e-01,    2.11111552358965165921e-01,
        1.79016861276632682038e-01,    1.46730474455361751659e-01,
        1.14286964966846398116e-01,    8.17210741336682237475e-02,
        4.90676743274180142536e-02,    1.63617316264867816422e-02};
        static const double A96 =      3.27249234748936795667e-02;
        for(i = 0; i < nradial / 2; i++) {
        y[i] = -x96[i];
        y[nradial - i - 1] = x96[i];
        *radial_weight =  A96;
        }
        break;
 
        default:
        if (job->taskid == 0)
        fprintf(file.out,"Number of radial intervals %d in Gauss_Chebyshev_weights_abscissa not allowed.\nAllowed values 16, 24, 32, 48, 64, 96\n",nradial);
        MPI_Finalize();
        exit(1);

       } // close switch

}
void Gauss_Legendre_weights_abscissa(int nradial, double *y, double *radial_weight, JOB_PARAM *job, FILES file)

{

int i;

       // source mymathlib.webtrellis.net/c_source/quadrature/gauss_legendre

    switch (nradial) {

       case 2:
        static const double x2[] = {
        5.77350269189625764507e-01};
        static const double A2[] = {
        1.00000000000000000000e+00};

        for(i = 0;i < nradial / 2;i++) {
        y[i] =                           -x2[i];
        y[nradial - i - 1] =              x2[i];
        radial_weight[i] =                A2[i];
        radial_weight[nradial - i - 1] =  A2[i];
        }
        break;

       case 4:
        static const double x4[] = {
      //8.61136311594052575248e-01,  3.39981043584856264792e-01};
        8.611363115940525752239e-01, 3.399810435848562648027e-01};
        static const double A4[] = {
      //3.47854845137453857383e-01,  6.52145154862546142644e-01};
        3.478548451374538573731e-01, 6.5214515486254614263e-01};


        for(i = 0;i < nradial / 2;i++) {
        y[i] = -x4[i];
        y[nradial - i - 1] = x4[i];
        radial_weight[i] =  A4[i];
        radial_weight[nradial - i - 1] =  A4[i];
        }
        break;

       case 8:
        static const double x8[] = {
      /*9.60289856497536231661e-01,  7.96666477413626739567e-01,
        5.25532409916328985830e-01,  1.83434642495649804936e-01}; */
        9.602898564975362316836e-01, 7.966664774136267395916e-01,
        5.255324099163289858177e-01, 1.834346424956498049395e-01};
        static const double A8[] = {
      /*1.01228536290376259154e-01, 2.22381034453374470546e-01,
        3.13706645877887287338e-01, 3.62683783378361982976e-01};*/
        1.01228536290376259153e-01, 2.22381034453374470544e-01,
        3.13706645877887287338e-01, 3.62683783378361982965e-01};

        for(i = 0;i < nradial / 2;i++) {
        y[i] = -x8[i];
        y[nradial - i - 1] = x8[i];
        radial_weight[i] =  A8[i];
        radial_weight[nradial - i - 1] =  A8[i];
        }
        break;

       case 16:
        static const double x16[] = {
        9.89400934991649932601e-01, 9.44575023073232576090e-01,
        8.65631202387831743866e-01, 7.55404408355003033891e-01,
        6.17876244402643748452e-01, 4.58016777657227386350e-01,
        2.81603550779258913231e-01, 9.50125098376374401877e-02};
        static const double A16[] = {
        2.71524594117540948514e-02, 6.22535239386478928628e-02,
        9.51585116824927848073e-02, 1.24628971255533872056e-01,
        1.49595988816576732080e-01, 1.69156519395002538183e-01,
        1.82603415044923588872e-01,1.89450610455068496287e-01};

        for(i = 0;i < nradial / 2;i++) {
        y[i] = -x16[i];
        y[nradial - i - 1] = x16[i];
        radial_weight[i] =  A16[i];
        radial_weight[nradial - i - 1] =  A16[i];
        }
        break;

       case 24:
        static const double x24[] = {
        9.95187219997021360195e-01, 9.74728555971309498199e-01,
        9.38274552002732758539e-01, 8.86415527004401034190e-01,
        8.20001985973902921981e-01, 7.40124191578554364260e-01,
        6.48093651936975569268e-01, 5.45421471388839535649e-01,
        4.33793507626045138478e-01, 3.15042679696163374398e-01,
        1.91118867473616309153e-01, 6.40568928626056260827e-02};
        static const double A24[] = {
        1.23412297999871995469e-02, 2.85313886289336631809e-02,
        4.42774388174198061695e-02, 5.92985849154367807461e-02,
        7.33464814110803057346e-02, 8.61901615319532759152e-02,
        9.76186521041138882720e-02, 1.07444270115965634785e-01,
        1.15505668053725601353e-01, 1.21670472927803391202e-01,
        1.25837456346828296117e-01, 1.27938195346752156976e-01};

        for(i = 0;i < nradial / 2;i++) {
        y[i] =                           -x24[i];
        y[nradial - i - 1] =              x24[i];
        radial_weight[i] =                A24[i];
        radial_weight[nradial - i - 1] =  A24[i];
        }
        break;

       case 32:
        static const double x32[] = {
        9.97263861849481563534e-01, 9.85611511545268335400e-01,
        9.64762255587506430761e-01, 9.34906075937739689159e-01,
        8.96321155766052123971e-01, 8.49367613732569970160e-01,
        7.94483795967942406965e-01, 7.32182118740289680412e-01,
        6.63044266930215200960e-01, 5.87715757240762329066e-01,
        5.06899908932229390044e-01, 4.21351276130635345353e-01,
        3.31868602282127649782e-01, 2.39287362252137074544e-01,
        1.44471961582796493484e-01, 4.83076656877383162364e-02};
        static const double A32[] = {
        7.01861000947009660028e-03, 1.62743947309056706058e-02,
        2.53920653092620594561e-02, 3.42738629130214331033e-02,
        4.28358980222266806557e-02, 5.09980592623761761959e-02,
        5.86840934785355471448e-02, 6.58222227763618468406e-02,
        7.23457941088485062287e-02, 7.81938957870703064685e-02,
        8.33119242269467552223e-02, 8.76520930044038111450e-02,
        9.11738786957638847129e-02, 9.38443990808045656367e-02,
        9.56387200792748594185e-02, 9.65400885147278005666e-02};  

        for(i = 0;i < nradial / 2;i++) {
        y[i] =                           -x32[i];
        y[nradial - i - 1] =              x32[i];
        radial_weight[i] =                A32[i];
        radial_weight[nradial - i - 1] =  A32[i];
        }
        break;


/*
48

    3.23801709628693620343e-02,    9.70046992094626989322e-02,
    1.61222356068891718055e-01,    2.24763790394689061224e-01,
    2.87362487355455576728e-01,    3.48755886292160738148e-01,
    4.08686481990716729925e-01,    4.66902904750958404535e-01,
    5.23160974722233033658e-01,    5.77224726083972703838e-01,
    6.28867396776513624013e-01,    6.77872379632663905208e-01,
    7.24034130923814654658e-01,    7.67159032515740339276e-01,
    8.07066204029442627087e-01,    8.43588261624393530704e-01,
    8.76572020274247885885e-01,    9.05879136715569672805e-01,
    9.31386690706554333107e-01,    9.52987703160430860724e-01,
    9.70591592546247250472e-01,    9.84124583722826857765e-01,
    9.93530172266350757526e-01,    9.98771007252426118580e-01
};

static const double A[] = {
    6.47376968126839225006e-02,    6.44661644359500822082e-02,
    6.39242385846481866207e-02,    6.31141922862540256548e-02,
    6.20394231598926639029e-02,    6.07044391658938800517e-02,
    5.91148396983956357477e-02,    5.72772921004032157044e-02,
    5.51995036999841628676e-02,    5.28901894851936670964e-02,
    5.03590355538544749590e-02,    4.76166584924904748267e-02,
    4.46745608566942804201e-02,    4.15450829434647492133e-02,
    3.82413510658307063158e-02,    3.47772225647704388909e-02,
    3.11672278327980889025e-02,    2.74265097083569482001e-02,
    2.35707608393243791410e-02,    1.96161604573555278142e-02,
    1.55793157229438487279e-02,    1.14772345792345394895e-02,
    7.32755390127626210220e-03,    3.15334605230583863260e-03


64

    2.43502926634244325092e-02,    7.29931217877990394521e-02,
    1.21462819296120554468e-01,    1.69644420423992818036e-01,
    2.17423643740007084148e-01,    2.64687162208767416384e-01,
    3.11322871990210956161e-01,    3.57220158337668115941e-01,
    4.02270157963991603693e-01,    4.46366017253464087978e-01,
    4.89403145707052957474e-01,    5.31279464019894545634e-01,
    5.71895646202634034291e-01,    6.11155355172393250241e-01,
    6.48965471254657339884e-01,    6.85236313054233242559e-01,
    7.19881850171610826840e-01,    7.52819907260531896590e-01,
    7.83972358943341407619e-01,    8.13265315122797559746e-01,
    8.40629296252580362743e-01,    8.65999398154092819759e-01,
    8.89315445995114105856e-01,    9.10522137078502805780e-01,
    9.29569172131939575846e-01,    9.46411374858402816069e-01,
    9.61008799652053718944e-01,    9.73326827789910963733e-01,
    9.83336253884625956939e-01,    9.91013371476744320732e-01,
    9.96340116771955279355e-01,    9.99305041735772139465e-01
};

static const double A[] = {
    4.86909570091397203845e-02,    4.85754674415034269347e-02,
    4.83447622348029571711e-02,    4.79993885964583077283e-02,
    4.75401657148303086630e-02,    4.69681828162100173254e-02,
    4.62847965813144172952e-02,    4.54916279274181444806e-02,
    4.45905581637565630589e-02,    4.35837245293234533762e-02,
    4.24735151236535890083e-02,    4.12625632426235286114e-02,
    3.99537411327203413872e-02,    3.85501531786156291275e-02,
    3.70551285402400460401e-02,    3.54722132568823838103e-02,
    3.38051618371416093907e-02,    3.20579283548515535845e-02,
    3.02346570724024788683e-02,    2.83396726142594832269e-02,
    2.63774697150546586725e-02,    2.43527025687108733385e-02,
    2.22701738083832541592e-02,    2.01348231535302093718e-02,
    1.79517157756973430852e-02,    1.57260304760247193221e-02,
    1.34630478967186425983e-02,    1.11681394601311288184e-02,
    8.84675982636394772300e-03,    6.50445796897836285625e-03,
    4.14703326056246763525e-03,    1.78328072169643294728e-03


96

    1.62767448496029695789e-02,    4.88129851360497311115e-02,
    8.12974954644255589937e-02,    1.13695850110665920914e-01,
    1.45973714654896941992e-01,    1.78096882367618602759e-01,
    2.10031310460567203601e-01,    2.41743156163840012331e-01,
    2.73198812591049141485e-01,    3.04364944354496353015e-01,
    3.35208522892625422607e-01,    3.65696861472313635024e-01,
    3.95797649828908603298e-01,    4.25478988407300545368e-01,
    4.54709422167743008628e-01,    4.83457973920596359777e-01,
    5.11694177154667673569e-01,    5.39388108324357436218e-01,
    5.66510418561397168413e-01,    5.93032364777572080659e-01,
    6.18925840125468570377e-01,    6.44163403784967106801e-01,
    6.68718310043916153943e-01,    6.92564536642171561332e-01,
    7.15676812348967626199e-01,    7.38030643744400132876e-01,
    7.59602341176647498681e-01,    7.80369043867433217620e-01,
    8.00308744139140817216e-01,    8.19400310737931675557e-01,
    8.37623511228187121497e-01,    8.54959033434601455438e-01,
    8.71388505909296502900e-01,    8.86894517402420416068e-01,
    9.01460635315852341334e-01,    9.15071423120898074215e-01,
    9.27712456722308690977e-01,    9.39370339752755216934e-01,
    9.50032717784437635746e-01,    9.59688291448742539290e-01,
    9.68326828463264212168e-01,    9.75939174585136466459e-01,
    9.82517263563014677430e-01,    9.88054126329623799458e-01,
    9.92543900323762624555e-01,    9.95981842987209290633e-01,
    9.98364375863181677739e-01,    9.99689503883230766825e-01
};

static const double A[] = {
    3.25506144923631662418e-02,    3.25161187138688359885e-02,
    3.24471637140642693631e-02,    3.23438225685759284293e-02,
    3.22062047940302506674e-02,    3.20344562319926632176e-02,
    3.18287588944110065352e-02,    3.15893307707271685576e-02,
    3.13164255968613558124e-02,    3.10103325863138374230e-02,
    3.06713761236691490147e-02,    3.02999154208275937943e-02,
    2.98963441363283859846e-02,    2.94610899581679059697e-02,
    2.89946141505552365432e-02,    2.84974110650853856455e-02,
    2.79700076168483344400e-02,    2.74129627260292428232e-02,
    2.68268667255917621977e-02,    2.62123407356724139131e-02,
    2.55700360053493614996e-02,    2.49006332224836102884e-02,
    2.42048417923646912830e-02,    2.34833990859262198430e-02,
    2.27370696583293740018e-02,    2.19666444387443491945e-02,
    2.11729398921912989884e-02,    2.03567971543333245953e-02,
    1.95190811401450224097e-02,    1.86606796274114673859e-02,
    1.77825023160452608374e-02,    1.68854798642451724498e-02,
    1.59705629025622913804e-02,    1.50387210269949380059e-02,
    1.40909417723148609158e-02,    1.31282295669615726374e-02,
    1.21516046710883196352e-02,    1.11621020998384985916e-02,
    1.01607705350084157574e-02,    9.14867123078338663265e-03,
    8.12687692569875921698e-03,    7.09647079115386526927e-03,
    6.05854550423596168331e-03,    5.01420274292751769241e-03,
    3.96455433844468667376e-03,    2.91073181793494640833e-03,
    1.85396078894692173237e-03,    7.96792065552012429429e-04


*/

       case 64:
        static const double x64[] = {
        9.993050417357721394569e-01, 9.963401167719552793469e-01,
        9.910133714767443207394e-01, 9.833362538846259569313e-01,
        9.733268277899109637419e-01, 9.610087996520537189186e-01,
        9.464113748584028160625e-01, 9.295691721319395758215e-01,
        9.105221370785028057564e-01, 8.893154459951141058534e-01,
        8.659993981540928197608e-01, 8.406292962525803627517e-01,
        8.132653151227975597419e-01, 7.839723589433414076102e-01,
        7.528199072605318966119e-01, 7.198818501716108268489e-01,
        6.852363130542332425636e-01, 6.489654712546573398578e-01,
        6.111553551723932502489e-01, 5.718956462026340342839e-01,
        5.312794640198945456580e-01, 4.894031457070529574785e-01,
        4.463660172534640879849e-01, 4.022701579639916036958e-01,
        3.572201583376681159504e-01, 3.113228719902109561575e-01,
        2.646871622087674163740e-01, 2.174236437400070841496e-01,
        1.696444204239928180373e-01, 1.214628192961205544704e-01,
        7.299312178779903944954e-02, 2.435029266342443250896e-02};
       static const double A64[] = {
        0.0017832807216964329473,    0.00414703326056246763529,
        0.00650445796897836285612,   0.008846759826363947723,
        0.0111681394601311288186,    0.0134630478967186425981,
        0.015726030476024719322,     0.01795171577569734308505,
        0.02013482315353020937234,   0.0222701738083832541593,
        0.0243527025687108733382,    0.02637746971505465867169,
        0.0283396726142594832275,    0.03023465707240247886797,
        0.0320579283548515535855,    0.0338051618371416093916,
        0.0354722132568823838107,    0.0370551285402400460404,
        0.038550153178615629129,     0.0399537411327203413867,
        0.0412625632426235286102,    0.04247351512365358900734,
        0.04358372452932345337683,   0.04459055816375656306013,
        0.0454916279274181444798,    0.046284796581314417296,
        0.0469681828162100173253,    0.0475401657148303086623,
        0.0479993885964583077281,    0.0483447622348029571698,
        0.0485754674415034269348,    0.0486909570091397203834};

        for(i = 0;i < nradial / 2;i++) {
        y[i] = -x64[i];
        y[nradial - i - 1] = x64[i];
        radial_weight[i] =  A64[i];
        radial_weight[nradial - i - 1] =  A64[i];
        }
        break;

        default:
        if (job->taskid == 0)
        fprintf(file.out,"Number of radial intervals %d in Gauss_Legendre_weights_abscissa not allowed.\nAllowed values 4, 8, 64\n",nradial);
        MPI_Finalize();
        exit(1);

       } // close switch

}

