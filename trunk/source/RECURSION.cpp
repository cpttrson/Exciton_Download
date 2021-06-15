#include "mycomplex.h"
#include "myconstants.h"
#include "USER_DATA.h"
#include "RECURSION.h"

using namespace std;

double e(int i, int ip, int t, double p, double PAx, double PBx)

{

  double recurrence;
  int flag;
  //fprintf(file_out,"i = %d ip = %d t = %d \n",i,ip,t) ;
  flag = (t == 0 && i == 0 && ip == 0) ? 1 : 0;
  if (flag == 1) {
    return k_one;
  }

  else if

  (i > ip)

  {
    recurrence = (t >= 0 && t <= i + ip) ? (e(i - 1, ip, t - 1, p, PAx, PBx) / 2 / p + PAx * e(i - 1, ip, t, p, PAx,
        PBx) + (t + 1) * e(i - 1, ip, t + 1, p, PAx, PBx)) : k_zero;
  }
  //fprintf(file_out,"i = %d ip = %d t = %d recurrence2 = %lf \n",i,ip,t,recurrence) ; }

  else 

  {
    recurrence = (t >= 0 && t <= i + ip) ? (e(i, ip - 1, t - 1, p, PAx, PBx) / 2 / p + PBx * e(i, ip - 1, t, p, PAx,
        PBx) + (t + 1) * e(i, ip - 1, t + 1, p, PAx, PBx)) : k_zero;
  }
  //fprintf(file_out,"i = %d ip = %d t = %d recurrence3 = %lf \n",i,ip,t,recurrence) ; }

  return recurrence;

}

double ftuvn(int t, int u, int v, int n, double *f, VECTOR_DOUBLE x)

{

  double result;
  if (t < 0 || u < 0 || v < 0)
    return k_zero;
  if (t >= u && t >= v) {
    //printf("t %d %d %d %d \n",t,u,v,n) ;
    //if (t == 0 && u == 0 && v == 0){fprintf(file_out,"tuv0 %d %d %d %d \n",t,u,v,n); }
    //if (t == 0 && u == 0 && v == 0){fprintf(file_out,"tuv0 %d %d %d %d %lf \n",t,u,v,n,f[n]); }
  if (t == 0 && u == 0 && v == 0)
    return f[n];
    result = (t > 0) ? x.comp1 * ftuvn(t - 1, u, v, n + 1, f, x) + (t - 1) * ftuvn(t - 2, u, v, n + 1, f, x) : 0.0;
  } else if (u > t && u >= v) {
    //printf("u %d %d %d %d \n",t,u,v,n) ;
    result = (u > 0) ? x.comp2 * ftuvn(t, u - 1, v, n + 1, f, x) + (u - 1) * ftuvn(t, u - 2, v, n + 1, f, x) : 0.0;
  } else {
    //printf("v %d %d %d %d \n",t,u,v,n) ;
    result = (v > 0) ? x.comp3 * ftuvn(t, u, v - 1, n + 1, f, x) + (v - 1) * ftuvn(t, u, v - 2, n + 1, f, x) : 0.0;
  }
  //printf("ftuvn tuv n %d %d %d %d %lf \n",t,u,v,n,result) ;
  return result;

}

double cosfactor(int tuv, double GdotR)

{

  int cosp = tuv - 4 * (tuv / 4);
  double cosfac;

  switch (cosp) {

    case 0:
      cosfac = cos(GdotR);
      break;
    case 1:
      cosfac = -sin(GdotR);
      break;
    case 2:
      cosfac = -cos(GdotR);
      break;
    case 3:
      cosfac = sin(GdotR);
  }

  return cosfac;

}

Complex cosfactor_complex(int tuv, double GdotR)

{
  int cosp = tuv - 4 * (tuv / 4);
  Complex cosfac;

  switch (cosp) {

    case 0:
      cosfac = Complex(cos(GdotR),sin(GdotR));
      break;
    case 1:
      cosfac = Complex(-sin(GdotR),cos(GdotR));
      break;
    case 2:
      cosfac = Complex(-cos(GdotR),-sin(GdotR));
      break;
    case 3:
      cosfac = Complex(sin(GdotR),-cos(GdotR));
  }

  return cosfac;

}

void non_recursive_ftuvn(int mm, int index_R, double f[][13][13][13], double en[][55], VECTOR_DOUBLE *r_12)

{

  int n;

          for (n = 0; n <= mm; n++)
          f[0][0][0][n] = en[index_R][n];

          for (n = 0; n <= mm; n++)
          f[1][0][0][n] = r_12->comp1 * en[index_R][n + 1];
          for (n = 0; n <= mm; n++)
          f[0][1][0][n] = r_12->comp2 * en[index_R][n + 1];
          for (n = 0; n <= mm; n++)
          f[0][0][1][n] = r_12->comp3 * en[index_R][n + 1];

          for (n = 0; n <= mm - 1; n++)
          f[1][1][0][n] = r_12->comp1 * f[0][1][0][n + 1];
          for (n = 0; n <= mm - 1; n++)
          f[1][0][1][n] = r_12->comp1 * f[0][0][1][n + 1];
          for (n = 0; n <= mm - 1; n++)
          f[0][1][1][n] = r_12->comp2 * f[0][0][1][n + 1];
          for (n = 0; n <= mm - 1; n++)
          f[2][0][0][n] = r_12->comp1 * f[1][0][0][n + 1] + f[0][0][0][n + 1];
          for (n = 0; n <= mm - 1; n++)
          f[0][2][0][n] = r_12->comp2 * f[0][1][0][n + 1] + f[0][0][0][n + 1];
          for (n = 0; n <= mm - 1; n++)
          f[0][0][2][n] = r_12->comp3 * f[0][0][1][n + 1] + f[0][0][0][n + 1];

          for (n = 0; n <= mm - 2; n++)
          f[2][1][0][n] = r_12->comp1 * f[1][1][0][n + 1] + f[0][1][0][n + 1];
          for (n = 0; n <= mm - 2; n++)
          f[2][0][1][n] = r_12->comp1 * f[1][0][1][n + 1] + f[0][0][1][n + 1];
          for (n = 0; n <= mm - 2; n++)
          f[1][2][0][n] = r_12->comp2 * f[1][1][0][n + 1] + f[1][0][0][n + 1];
          for (n = 0; n <= mm - 2; n++)
          f[0][2][1][n] = r_12->comp2 * f[0][1][1][n + 1] + f[0][0][1][n + 1];
          for (n = 0; n <= mm - 2; n++)
          f[0][1][2][n] = r_12->comp3 * f[0][1][1][n + 1] + f[0][1][0][n + 1];
          for (n = 0; n <= mm - 2; n++)
          f[1][0][2][n] = r_12->comp3 * f[1][0][1][n + 1] + f[1][0][0][n + 1];
          for (n = 0; n <= mm - 2; n++)
          f[1][1][1][n] = r_12->comp1 * f[0][1][1][n + 1];
          for (n = 0; n <= mm - 2; n++)
          f[3][0][0][n] = r_12->comp1 * f[2][0][0][n + 1] + two * f[1][0][0][n + 1];
          for (n = 0; n <= mm - 2; n++)
          f[0][3][0][n] = r_12->comp2 * f[0][2][0][n + 1] + two * f[0][1][0][n + 1];
          for (n = 0; n <= mm - 2; n++)
          f[0][0][3][n] = r_12->comp3 * f[0][0][2][n + 1] + two * f[0][0][1][n + 1];

          for (n = 0; n <= mm - 3; n++)
          f[2][1][1][n] = r_12->comp1 * f[1][1][1][n + 1] + f[0][1][1][n + 1];
          for (n = 0; n <= mm - 3; n++)
          f[1][2][1][n] = r_12->comp2 * f[1][1][1][n + 1] + f[1][0][1][n + 1];
          for (n = 0; n <= mm - 3; n++)
          f[1][1][2][n] = r_12->comp3 * f[1][1][1][n + 1] + f[1][1][0][n + 1];
          for (n = 0; n <= mm - 3; n++)
          f[2][2][0][n] = r_12->comp1 * f[1][2][0][n + 1] + f[0][2][0][n + 1];
          for (n = 0; n <= mm - 3; n++)
          f[2][0][2][n] = r_12->comp1 * f[1][0][2][n + 1] + f[0][0][2][n + 1];
          for (n = 0; n <= mm - 3; n++)
          f[0][2][2][n] = r_12->comp2 * f[0][1][2][n + 1] + f[0][0][2][n + 1];
          for (n = 0; n <= mm - 3; n++)
          f[3][1][0][n] = r_12->comp1 * f[2][1][0][n + 1] + two * f[1][1][0][n + 1];
          for (n = 0; n <= mm - 3; n++)
          f[3][0][1][n] = r_12->comp1 * f[2][0][1][n + 1] + two * f[1][0][1][n + 1];
          for (n = 0; n <= mm - 3; n++)
          f[0][3][1][n] = r_12->comp2 * f[0][2][1][n + 1] + two * f[0][1][1][n + 1];
          for (n = 0; n <= mm - 3; n++)
          f[1][3][0][n] = r_12->comp2 * f[1][2][0][n + 1] + two * f[1][1][0][n + 1];
          for (n = 0; n <= mm - 3; n++)
          f[1][0][3][n] = r_12->comp3 * f[1][0][2][n + 1] + two * f[1][0][1][n + 1];
          for (n = 0; n <= mm - 3; n++)
          f[0][1][3][n] = r_12->comp3 * f[0][1][2][n + 1] + two * f[0][1][1][n + 1];
          for (n = 0; n <= mm - 3; n++)
          f[4][0][0][n] = r_12->comp1 * f[3][0][0][n + 1] + three * f[2][0][0][n + 1];
          for (n = 0; n <= mm - 3; n++)
          f[0][4][0][n] = r_12->comp2 * f[0][3][0][n + 1] + three * f[0][2][0][n + 1];
          for (n = 0; n <= mm - 3; n++)
          f[0][0][4][n] = r_12->comp3 * f[0][0][3][n + 1] + three * f[0][0][2][n + 1];

          for (n = 0; n <= mm - 4; n++)
          f[3][1][1][n] = r_12->comp1 * f[2][1][1][n + 1] + two * f[1][1][1][n + 1];
          for (n = 0; n <= mm - 4; n++)
          f[1][3][1][n] = r_12->comp2 * f[1][2][1][n + 1] + two * f[1][1][1][n + 1];
          for (n = 0; n <= mm - 4; n++)
          f[1][1][3][n] = r_12->comp3 * f[1][1][2][n + 1] + two * f[1][1][1][n + 1];
          for (n = 0; n <= mm - 4; n++)
          f[2][2][1][n] = r_12->comp1 * f[1][2][1][n + 1] + f[0][2][1][n + 1];
          for (n = 0; n <= mm - 4; n++)
          f[2][1][2][n] = r_12->comp1 * f[1][1][2][n + 1] + f[0][1][2][n + 1];
          for (n = 0; n <= mm - 4; n++)
          f[1][2][2][n] = r_12->comp2 * f[1][1][2][n + 1] + f[1][0][2][n + 1];
          for (n = 0; n <= mm - 4; n++)
          f[3][2][0][n] = r_12->comp1 * f[2][2][0][n + 1] + two * f[1][2][0][n + 1];
          for (n = 0; n <= mm - 4; n++)
          f[3][0][2][n] = r_12->comp1 * f[2][0][2][n + 1] + two * f[1][0][2][n + 1];
          for (n = 0; n <= mm - 4; n++)
          f[0][3][2][n] = r_12->comp2 * f[0][2][2][n + 1] + two * f[0][1][2][n + 1];
          for (n = 0; n <= mm - 4; n++)
          f[2][3][0][n] = r_12->comp2 * f[2][2][0][n + 1] + two * f[2][1][0][n + 1];
          for (n = 0; n <= mm - 4; n++)
          f[2][0][3][n] = r_12->comp3 * f[2][0][2][n + 1] + two * f[2][0][1][n + 1];
          for (n = 0; n <= mm - 4; n++)
          f[0][2][3][n] = r_12->comp3 * f[0][2][2][n + 1] + two * f[0][2][1][n + 1];
          for (n = 0; n <= mm - 4; n++)
          f[4][1][0][n] = r_12->comp1 * f[3][1][0][n + 1] + three * f[2][1][0][n + 1];
          for (n = 0; n <= mm - 4; n++)
          f[4][0][1][n] = r_12->comp1 * f[3][0][1][n + 1] + three * f[2][0][1][n + 1];
          for (n = 0; n <= mm - 4; n++)
          f[0][4][1][n] = r_12->comp2 * f[0][3][1][n + 1] + three * f[0][2][1][n + 1];
          for (n = 0; n <= mm - 4; n++)
          f[1][4][0][n] = r_12->comp2 * f[1][3][0][n + 1] + three * f[1][2][0][n + 1];
          for (n = 0; n <= mm - 4; n++)
          f[0][1][4][n] = r_12->comp3 * f[0][1][3][n + 1] + three * f[0][1][2][n + 1];
          for (n = 0; n <= mm - 4; n++)
          f[1][0][4][n] = r_12->comp3 * f[1][0][3][n + 1] + three * f[1][0][2][n + 1];
          for (n = 0; n <= mm - 4; n++)
          f[5][0][0][n] = r_12->comp1 * f[4][0][0][n + 1] + four * f[3][0][0][n + 1];
          for (n = 0; n <= mm - 4; n++)
          f[0][5][0][n] = r_12->comp2 * f[0][4][0][n + 1] + four * f[0][3][0][n + 1];
          for (n = 0; n <= mm - 4; n++)
          f[0][0][5][n] = r_12->comp3 * f[0][0][4][n + 1] + four * f[0][0][3][n + 1];

          for (n = 0; n <= mm - 5; n++)
          f[2][2][2][n] = r_12->comp1 * f[1][2][2][n + 1] + f[0][2][2][n + 1];
          for (n = 0; n <= mm - 5; n++)
          f[3][2][1][n] = r_12->comp1 * f[2][2][1][n + 1] + two * f[1][2][1][n + 1];
          for (n = 0; n <= mm - 5; n++)
          f[3][1][2][n] = r_12->comp1 * f[2][1][2][n + 1] + two * f[1][1][2][n + 1];
          for (n = 0; n <= mm - 5; n++)
          f[1][3][2][n] = r_12->comp2 * f[1][2][2][n + 1] + two * f[1][1][2][n + 1];
          for (n = 0; n <= mm - 5; n++)
          f[2][3][1][n] = r_12->comp2 * f[2][2][1][n + 1] + two * f[2][1][1][n + 1];
          for (n = 0; n <= mm - 5; n++)
          f[2][1][3][n] = r_12->comp3 * f[2][1][2][n + 1] + two * f[2][1][1][n + 1];
          for (n = 0; n <= mm - 5; n++)
          f[1][2][3][n] = r_12->comp3 * f[1][2][2][n + 1] + two * f[1][2][1][n + 1];
          for (n = 0; n <= mm - 5; n++)
          f[4][1][1][n] = r_12->comp1 * f[3][1][1][n + 1] + three * f[2][1][1][n + 1];
          for (n = 0; n <= mm - 5; n++)
          f[1][4][1][n] = r_12->comp2 * f[1][3][1][n + 1] + three * f[1][2][1][n + 1];
          for (n = 0; n <= mm - 5; n++)
          f[1][1][4][n] = r_12->comp3 * f[1][1][3][n + 1] + three * f[1][1][2][n + 1];
          for (n = 0; n <= mm - 5; n++)
          f[4][2][0][n] = r_12->comp1 * f[3][2][0][n + 1] + three * f[2][2][0][n + 1];
          for (n = 0; n <= mm - 5; n++)
          f[4][0][2][n] = r_12->comp1 * f[3][0][2][n + 1] + three * f[2][0][2][n + 1];
          for (n = 0; n <= mm - 5; n++)
          f[2][4][0][n] = r_12->comp2 * f[2][3][0][n + 1] + three * f[2][2][0][n + 1];
          for (n = 0; n <= mm - 5; n++)
          f[0][4][2][n] = r_12->comp2 * f[0][3][2][n + 1] + three * f[0][2][2][n + 1];
          for (n = 0; n <= mm - 5; n++)
          f[2][0][4][n] = r_12->comp3 * f[2][0][3][n + 1] + three * f[2][0][2][n + 1];
          for (n = 0; n <= mm - 5; n++)
          f[0][2][4][n] = r_12->comp3 * f[0][2][3][n + 1] + three * f[0][2][2][n + 1];
          for (n = 0; n <= mm - 5; n++)
          f[3][3][0][n] = r_12->comp1 * f[2][3][0][n + 1] + two * f[1][3][0][n + 1];
          for (n = 0; n <= mm - 5; n++)
          f[3][0][3][n] = r_12->comp1 * f[2][0][3][n + 1] + two * f[1][0][3][n + 1];
          for (n = 0; n <= mm - 5; n++)
          f[0][3][3][n] = r_12->comp3 * f[0][3][2][n + 1] + two * f[0][3][1][n + 1];
          for (n = 0; n <= mm - 5; n++)
          f[5][1][0][n] = r_12->comp1 * f[4][1][0][n + 1] + four * f[3][1][0][n + 1];
          for (n = 0; n <= mm - 5; n++)
          f[5][0][1][n] = r_12->comp1 * f[4][0][1][n + 1] + four * f[3][0][1][n + 1];
          for (n = 0; n <= mm - 5; n++)
          f[1][5][0][n] = r_12->comp2 * f[1][4][0][n + 1] + four * f[1][3][0][n + 1];
          for (n = 0; n <= mm - 5; n++)
          f[0][5][1][n] = r_12->comp2 * f[0][4][1][n + 1] + four * f[0][3][1][n + 1];
          for (n = 0; n <= mm - 5; n++)
          f[0][1][5][n] = r_12->comp3 * f[0][1][4][n + 1] + four * f[0][1][3][n + 1];
          for (n = 0; n <= mm - 5; n++)
          f[1][0][5][n] = r_12->comp3 * f[1][0][4][n + 1] + four * f[1][0][3][n + 1];
          for (n = 0; n <= mm - 5; n++)
          f[6][0][0][n] = r_12->comp1 * f[5][0][0][n + 1] + five * f[4][0][0][n + 1];
          for (n = 0; n <= mm - 5; n++)
          f[0][6][0][n] = r_12->comp2 * f[0][5][0][n + 1] + five * f[0][4][0][n + 1];
          for (n = 0; n <= mm - 5; n++)
          f[0][0][6][n] = r_12->comp3 * f[0][0][5][n + 1] + five * f[0][0][4][n + 1];

          for (n = 0; n <= mm - 6; n++)
          f[5][1][1][n] = r_12->comp1 * f[4][1][1][n + 1] + four * f[3][1][1][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[1][5][1][n] = r_12->comp2 * f[1][4][1][n + 1] + four * f[1][3][1][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[1][1][5][n] = r_12->comp3 * f[1][1][4][n + 1] + four * f[1][1][3][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[3][2][2][n] = r_12->comp1 * f[2][2][2][n + 1] + two * f[1][2][2][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[2][3][2][n] = r_12->comp2 * f[2][2][2][n + 1] + two * f[2][1][2][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[2][2][3][n] = r_12->comp3 * f[2][2][2][n + 1] + two * f[2][2][1][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[3][3][1][n] = r_12->comp1 * f[2][3][1][n + 1] + two * f[1][3][1][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[3][1][3][n] = r_12->comp1 * f[2][1][3][n + 1] + two * f[1][1][3][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[1][3][3][n] = r_12->comp2 * f[1][2][3][n + 1] + two * f[1][1][3][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[4][3][0][n] = r_12->comp1 * f[3][3][0][n + 1] + three * f[2][3][0][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[4][0][3][n] = r_12->comp1 * f[3][0][3][n + 1] + three * f[2][0][3][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[0][4][3][n] = r_12->comp2 * f[0][3][3][n + 1] + three * f[0][2][3][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[3][4][0][n] = r_12->comp2 * f[3][3][0][n + 1] + three * f[3][2][0][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[0][3][4][n] = r_12->comp3 * f[0][3][3][n + 1] + three * f[0][3][2][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[3][0][4][n] = r_12->comp3 * f[3][0][3][n + 1] + three * f[3][0][2][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[4][2][1][n] = r_12->comp1 * f[3][2][1][n + 1] + three * f[2][2][1][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[4][1][2][n] = r_12->comp1 * f[3][1][2][n + 1] + three * f[2][1][2][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[2][4][1][n] = r_12->comp2 * f[2][3][1][n + 1] + three * f[2][2][1][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[1][4][2][n] = r_12->comp2 * f[1][3][2][n + 1] + three * f[1][2][2][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[2][1][4][n] = r_12->comp3 * f[2][1][3][n + 1] + three * f[2][1][2][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[1][2][4][n] = r_12->comp3 * f[1][2][3][n + 1] + three * f[1][2][2][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[5][2][0][n] = r_12->comp1 * f[4][2][0][n + 1] + four * f[3][2][0][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[5][0][2][n] = r_12->comp1 * f[4][0][2][n + 1] + four * f[3][0][2][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[2][5][0][n] = r_12->comp2 * f[2][4][0][n + 1] + four * f[2][3][0][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[0][5][2][n] = r_12->comp2 * f[0][4][2][n + 1] + four * f[0][3][2][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[0][2][5][n] = r_12->comp3 * f[0][2][4][n + 1] + four * f[0][2][3][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[2][0][5][n] = r_12->comp3 * f[2][0][4][n + 1] + four * f[2][0][3][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[6][1][0][n] = r_12->comp1 * f[5][1][0][n + 1] + five * f[4][1][0][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[6][0][1][n] = r_12->comp1 * f[5][0][1][n + 1] + five * f[4][0][1][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[1][6][0][n] = r_12->comp2 * f[1][5][0][n + 1] + five * f[1][4][0][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[0][6][1][n] = r_12->comp2 * f[0][5][1][n + 1] + five * f[0][4][1][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[0][1][6][n] = r_12->comp3 * f[0][1][5][n + 1] + five * f[0][1][4][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[1][0][6][n] = r_12->comp3 * f[1][0][5][n + 1] + five * f[1][0][4][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[7][0][0][n] = r_12->comp1 * f[6][0][0][n + 1] + six * f[5][0][0][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[0][7][0][n] = r_12->comp2 * f[0][6][0][n + 1] + six * f[0][5][0][n + 1];
          for (n = 0; n <= mm - 6; n++)
          f[0][0][7][n] = r_12->comp3 * f[0][0][6][n + 1] + six * f[0][0][5][n + 1];

          for (n = 0; n <= mm - 7; n++)
          f[3][3][2][n] = r_12->comp1 * f[2][3][2][n + 1] + two * f[1][3][2][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[3][2][3][n] = r_12->comp1 * f[2][2][3][n + 1] + two * f[1][2][3][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[2][3][3][n] = r_12->comp2 * f[2][2][3][n + 1] + two * f[2][1][3][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[4][2][2][n] = r_12->comp1 * f[3][2][2][n + 1] + three * f[2][2][2][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[2][4][2][n] = r_12->comp2 * f[2][3][2][n + 1] + three * f[2][2][2][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[2][2][4][n] = r_12->comp3 * f[2][2][3][n + 1] + three * f[2][2][2][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[6][2][0][n] = r_12->comp1 * f[5][2][0][n + 1] + five * f[4][2][0][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[6][0][2][n] = r_12->comp1 * f[5][0][2][n + 1] + five * f[4][0][2][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[2][6][0][n] = r_12->comp2 * f[2][5][0][n + 1] + five * f[2][4][0][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[0][6][2][n] = r_12->comp2 * f[0][5][2][n + 1] + five * f[0][4][2][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[2][0][6][n] = r_12->comp3 * f[2][0][5][n + 1] + five * f[2][0][4][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[0][2][6][n] = r_12->comp3 * f[0][2][5][n + 1] + five * f[0][2][4][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[4][4][0][n] = r_12->comp1 * f[3][4][0][n + 1] + three * f[2][4][0][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[4][0][4][n] = r_12->comp1 * f[3][0][4][n + 1] + three * f[2][0][4][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[0][4][4][n] = r_12->comp2 * f[0][3][4][n + 1] + three * f[0][2][4][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[5][3][0][n] = r_12->comp1 * f[4][3][0][n + 1] + four * f[3][3][0][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[5][0][3][n] = r_12->comp1 * f[4][0][3][n + 1] + four * f[3][0][3][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[0][5][3][n] = r_12->comp2 * f[0][4][3][n + 1] + four * f[0][3][3][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[3][5][0][n] = r_12->comp2 * f[3][4][0][n + 1] + four * f[3][3][0][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[0][3][5][n] = r_12->comp3 * f[0][3][4][n + 1] + four * f[0][3][3][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[3][0][5][n] = r_12->comp3 * f[3][0][4][n + 1] + four * f[3][0][3][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[4][3][1][n] = r_12->comp1 * f[3][3][1][n + 1] + three * f[2][3][1][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[4][1][3][n] = r_12->comp1 * f[3][1][3][n + 1] + three * f[2][1][3][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[3][4][1][n] = r_12->comp2 * f[3][3][1][n + 1] + three * f[3][2][1][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[1][4][3][n] = r_12->comp2 * f[1][3][3][n + 1] + three * f[1][2][3][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[1][3][4][n] = r_12->comp3 * f[1][3][3][n + 1] + three * f[1][3][2][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[3][1][4][n] = r_12->comp3 * f[3][1][3][n + 1] + three * f[3][1][2][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[5][2][1][n] = r_12->comp1 * f[4][2][1][n + 1] + four * f[3][2][1][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[5][1][2][n] = r_12->comp1 * f[4][1][2][n + 1] + four * f[3][1][2][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[1][5][2][n] = r_12->comp2 * f[1][4][2][n + 1] + four * f[1][3][2][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[2][5][1][n] = r_12->comp2 * f[2][4][1][n + 1] + four * f[2][3][1][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[2][1][5][n] = r_12->comp3 * f[2][1][4][n + 1] + four * f[2][1][3][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[1][2][5][n] = r_12->comp3 * f[1][2][4][n + 1] + four * f[1][2][3][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[6][1][1][n] = r_12->comp1 * f[5][1][1][n + 1] + five * f[4][1][1][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[1][6][1][n] = r_12->comp2 * f[1][5][1][n + 1] + five * f[1][4][1][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[1][1][6][n] = r_12->comp3 * f[1][1][5][n + 1] + five * f[1][1][4][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[7][0][1][n] = r_12->comp1 * f[6][0][1][n + 1] + six * f[5][0][1][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[7][1][0][n] = r_12->comp1 * f[6][1][0][n + 1] + six * f[5][1][0][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[1][7][0][n] = r_12->comp2 * f[1][6][0][n + 1] + six * f[1][5][0][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[0][7][1][n] = r_12->comp2 * f[0][6][1][n + 1] + six * f[0][5][1][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[0][1][7][n] = r_12->comp3 * f[0][1][6][n + 1] + six * f[0][1][5][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[1][0][7][n] = r_12->comp3 * f[1][0][6][n + 1] + six * f[1][0][5][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[8][0][0][n] = r_12->comp1 * f[7][0][0][n + 1] + seven * f[6][0][0][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[0][8][0][n] = r_12->comp2 * f[0][7][0][n + 1] + seven * f[0][6][0][n + 1];
          for (n = 0; n <= mm - 7; n++)
          f[0][0][8][n] = r_12->comp3 * f[0][0][7][n + 1] + seven * f[0][0][6][n + 1];

}

void erf_derivative(double *derivative, int max, double a, double x, RECIPROCAL_LATTICE *G)

{

  // Evaluate derivatives of erf(a * x + b) up to nth order including zeroth order
  // Argument list: max - order of highest derivative, a in a * x + b, x - full argument of erf(a * x + b)
  // derivative array must be of length n + 1
  // Derivative of erfc(a * x + b) is -derivative of erf(a * x + b) iff n >= 1

int n, k;
double a_n, ax, ax_2, ax_k, exp_ax2, fac;
max = 5;

  ax_2 = x * x;
  exp_ax2 = exp(-ax_2);
  derivative[0] = erf(x);
  fac = two / rtpi * exp_ax2;
  a_n  = a;
  for (n = 1; n <= max; n++) {
    ax_k = k_one;
    derivative[n] = k_zero;
    for (k = 0; k <= n; k++) {
      derivative[n] += fac * G->B->a[n][k] * a_n * ax_k;
      ax_k *= x;
     }
     a_n  *= a;
    }
    //for(n = 0; n <=max;n++) {
    //printf("derivative %3d %20.10lf %20.10lf %20.10lf\n",n,derivative[n],a,x); }

}

void erf_exp_derivative(double *erf_derivative, double *exp_derivative, int max, double a, double x, RECIPROCAL_LATTICE *G)

{

  // Evaluate derivatives of exp(a * x + b) and erf(a * x + b) up to nth order including zeroth order
  // Argument list: max - order of highest derivative, a in a * x + b, x - full argument of exp/erf(a * x + b)
  // derivative array must be of length n + 1
  // Derivative of erfc(a * x + b) is -derivative of erf(a * x + b) iff n >= 1

int n, k;
double a_n, ax, exp_ax2, x_k, fac;

  // exp derivative 
  max = 5;

  a_n  = k_one;
  exp_ax2 = exp(-x * x);
  for (n = 0; n <= max; n++) {
    x_k = k_one;
    exp_derivative[n] = k_zero;
    for (k = 0; k <= n; k++) {
      exp_derivative[n] += a_n * G->B->a[n + 1][k] * x_k * exp_ax2;
      x_k *= x;
     }
     a_n  *= a;
    }

  // erf derivative 

  erf_derivative[0] = erf(x);
  fac = two / rtpi * a;
  for (n = 1; n <= max; n++) {
    erf_derivative[n] =  fac * exp_derivative[n - 1];
   }

    //for(n = 0; n <=max;n++) {
    //printf("exp erf derivative %3d %20.10lf %20.10lf %20.10lf %20.10lf\n",n,exp_derivative[n],erf_derivative[n],a,x); }

}

void erfc_derivative(double *derivative, int max, double a, double x, RECIPROCAL_LATTICE *G)

{

  // Evaluate derivatives of erfc(a * x + b) up to nth order including zeroth order
  // Argument list: max - order of highest derivative, a in a * x + b, x - full argument of erfc(a * x + b)
  // derivative array must be of length n + 1
  // Derivative of erfc(a * x + b) is -derivative of erf(a * x + b) iff n >= 1

int n, k;
double a_n, ax, ax_2, ax_k, exp_ax2, fac;

  ax_2 = x * x;
  exp_ax2 = exp(-ax_2);
  derivative[0] = erfc(x);
  fac = two / rtpi * exp_ax2;
  a_n  = a;
  for (n = 1; n <= max; n++) {
    ax_k = k_one;
    derivative[n] = k_zero;
    for (k = 0; k <= n; k++) {
      derivative[n] -= fac * G->B->a[n][k] * a_n * ax_k;
      ax_k *= x;
     }
     a_n  *= a;
    }
     //for(n = 0; n <=max;n++) {
     //printf("erfc derivative %3d %20.10lf %20.10lf %20.10lf\n",n,derivative[n],a,x); }

}

/*
double e3(int i, int j, int k, int t, double p, double PAx, double PBx, double PCx)

{

  double recurrence;
  int flag;
  //fprintf(file_out,"i = %d j = %d k = %d t = %d \n",i,j,k,t) ;
  flag = (t == 0 && i == 0 && j == 0 && k == 0) ? 1 : 0;
  if (flag == 1) {
    return k_one;
  }

  else if

  (i >= j && i >= k)

  {
  recurrence = (t >= 0 && t <= i + j + k) ? (e3(i - 1, j, k, t - 1, p, PAx, PBx, PCx) / two / p + PAx * e3(i - 1, j, k, t, p, PAx, PBx, PCx)
                                 + (t + 1) * e3(i - 1, j, k, t + 1, p, PAx, PBx, PCx)) : k_zero;
  }
  //fprintf(file_out,"i = %d j = %d k = %d t = %d recurrence2 = %lf \n",i,j,k,t,recurrence) ; }

  else if

  (j > i && j >= k)

  {
  recurrence = (t >= 0 && t <= i + j + k) ? (e3(i, j - 1, k, t - 1, p, PAx, PBx, PCx) / two / p + PBx * e3(i, j - 1, k, t, p, PAx, PBx, PCx)
                                 + (t + 1) * e3(i, j - 1, k, t + 1, p, PAx, PBx, PCx)) : k_zero;
  }
  //fprintf(file_out,"i = %d j = %d k = %d t = %d recurrence2 = %lf \n",i,j,k,t,recurrence) ; }

  else 

  {
  recurrence = (t >= 0 && t <= i + j + k) ? (e3(i, j, k - 1, t - 1, p, PAx, PBx, PCx) / two / p + PCx * e3(i, j, k - 1, t, p, PAx, PBx, PCx)
                                 + (t + 1) * e3(i, j, k - 1, t + 1, p, PAx, PBx, PCx)) : k_zero;
  }
  //fprintf(file_out,"i = %d j = %d k = %d t = %d recurrence3 = %lf \n",i,j,k,t,recurrence) ; }
  return recurrence;

}

void non_recursive_e(double e[][9][9], int mm, double p, double PAx, double PBx)

{

          // i + i' = 0   s * s

          e[0][0][0] = k_one;

          if (mm == 0) return;

          // i + i' = 1   s * p

          e[1][0][0] = PAx;
          e[0][1][0] = PBx;

          e[1][0][1] = half / p;
          e[0][1][1] = e[1][0][1];

          if (mm == 1) return;

          // i + i' = 2  p * p, s * d

          e[2][0][0] = e[1][0][0] * e[1][0][0] + e[1][0][1];
          e[2][0][1] = e[1][0][0] * e[1][0][1] * two;
          e[2][0][2] = e[1][0][1] * e[1][0][1];

          e[1][1][0] = e[1][0][0] * e[0][1][0];
          e[1][1][1] = e[1][0][0] * e[0][1][1] + e[0][1][0] * e[1][0][1];
          e[1][1][2] = e[1][0][1] * e[0][1][1];

          if (mm == 2) return;

          // i + i' = 3  p * d

          e[3][0][0] = e[2][0][0] * e[1][0][0] + e[2][0][1];
          e[3][0][1] = e[2][0][0] * e[1][0][1] + e[1][0][0] * e[2][0][1] + e[2][0][2] * two;
          e[3][0][2] = e[2][0][1] * e[1][0][1] + e[2][0][2] * e[1][0][0];
          e[3][0][3] = e[2][0][2] * e[1][0][1];

          e[2][1][0] = e[2][0][0] * e[0][1][0] + e[1][0][0] * e[1][1][0] + e[1][1][1];
          e[2][1][1] = e[1][1][0] * e[1][0][1] + e[1][1][1] * e[1][0][0] + e[1][1][2];
          e[2][1][2] = e[1][1][1] * e[1][0][1] + e[1][1][2] * e[1][0][0];
          e[2][1][3] = e[1][1][2] * e[1][0][1];

          e[1][2][0] = e[0][2][0] * e[1][0][0] + e[0][1][0] * e[1][1][0] + e[1][1][1];
          e[1][2][1] = e[1][1][0] * e[0][1][1] + e[1][1][1] * e[0][1][0] + e[1][1][2];
          e[1][2][2] = e[1][1][1] * e[0][1][1] + e[1][1][2] * e[0][1][0];
          e[1][2][3] = e[1][1][2] * e[0][1][1];

          e[0][3][0] = e[0][2][0] * e[0][1][0] + e[0][2][1];
          e[0][3][1] = e[0][2][0] * e[0][1][1] + e[0][1][0] * e[0][2][1] + e[0][2][2] * two;
          e[0][3][2] = e[0][2][1] * e[0][1][1] + e[0][2][2] * e[0][1][0];
          e[0][3][3] = e[0][2][2] * e[0][1][1];

          if (mm == 3) return;

}

void non_recursive_E(double *e, double *C_max, int imax, int jmax, int off3, double *pab_inv, double *PA, double *PB, SHELL *shells, FILES file)


{

int i, j, k;
int mm, off1, off2;

  off2    = (imax + jmax + 1) * off3;
  off1    = (jmax + 1) * off2;
  mm      =  imax + jmax;

  // i + i' = 0   s * s

  for(i = 0; i < off3; i++) {
    e[i] = k_one;
    C_max[i] = k_one;
   }

  if (mm == 0) return;

  // i + i' = 1   s * p

int offset_100 =     off1;
int offset_010 =                off2;
int offset_101 =     off1 +                off3;
int offset_011 =                off2 +     off3;

  for(i = 0; i < off3; i++) {
    e[offset_100 + i] = PA[i];
    C_max[i] = fabs(PA[i]);
 }

  for(i = 0; i < off3; i++) {
    e[offset_010 + i] = PB[i];
    C_max[i] = (C_max[i] < fabs(e[offset_010 + i]) ? C_max[i] : fabs(e[offset_010 + i]));
 }

  for(i = 0; i < off3; i++) {
    e[offset_101 + i] = half * pab_inv[i];
    C_max[i] = (C_max[i] < fabs(e[offset_101 + i]) ? C_max[i] : fabs(e[offset_101 + i]));
 }

  for(i = 0; i < off3; i++) {
    e[offset_011 + i] = e[offset_101 + i];
    C_max[i] = (C_max[i] < fabs(e[offset_011 + i]) ? C_max[i] : fabs(e[offset_011 + i]));
 }

  if (mm == 1) return;

  // i + i' = 2  p * p, s * d

int offset_200 = 2 * off1;
int offset_201 = 2 * off1 +                off3;
int offset_202 = 2 * off1 +          + 2 * off3;
int offset_110 =     off1 +     off2;
int offset_111 =     off1 +     off2 +     off3;
int offset_112 =     off1 +     off2 +     off3;

  for(i = 0; i < off3; i++) {
  e[offset_200 + i] = e[offset_100 + i] * e[offset_100 + i] + e[offset_101 + i];
 }

  for(i = 0; i < off3; i++) {
  e[offset_201 + i] = e[offset_100 + i] * e[offset_101 + i] * two;
 }

  for(i = 0; i < off3; i++) {
  e[offset_202 + i] = e[offset_101 + i] * e[offset_011 + i];
 }

  for(i = 0; i < off3; i++) {
  e[offset_110 + i] = e[offset_100 + i] * e[offset_010 + i];
 }

  for(i = 0; i < off3; i++) {
  e[offset_111 + i] = e[offset_100 + i] * e[offset_011 + i] + e[offset_010 + i] * e[offset_101 + i];
 }

  for(i = 0; i < off3; i++) {
  e[offset_112 + i] = e[offset_101 + i] * e[offset_011 + i];
 }

          //e[2][0][0] = e[1][0][0] * e[1][0][0] + e[1][0][1];
          //e[2][0][1] = e[1][0][0] * e[1][0][1] * two;
          //e[2][0][2] = e[1][0][1] * e[1][0][1];

          //e[1][1][0] = e[1][0][0] * e[0][1][0];
          //e[1][1][1] = e[1][0][0] * e[0][1][1] + e[0][1][0] * e[1][0][1];
          //e[1][1][2] = e[1][0][1] * e[0][1][1];

          if (mm == 2) return;

          // i + i' = 3  p * d
*/
/*

          e[3][0][0] = e[2][0][0] * e[1][0][0] + e[2][0][1];
          e[3][0][1] = e[2][0][0] * e[1][0][1] + e[1][0][0] * e[2][0][1] + e[2][0][2] * two;
          e[3][0][2] = e[2][0][1] * e[1][0][1] + e[2][0][2] * e[1][0][0];
          e[3][0][3] = e[2][0][2] * e[1][0][1];

          e[2][1][0] = e[2][0][0] * e[0][1][0] + e[1][0][0] * e[1][1][0] + e[1][1][1];
          e[2][1][1] = e[1][1][0] * e[1][0][1] + e[1][1][1] * e[1][0][0] + e[1][1][2];
          e[2][1][2] = e[1][1][1] * e[1][0][1] + e[1][1][2] * e[1][0][0];
          e[2][1][3] = e[1][1][2] * e[1][0][1];

          e[1][2][0] = e[0][2][0] * e[1][0][0] + e[0][1][0] * e[1][1][0] + e[1][1][1];
          e[1][2][1] = e[1][1][0] * e[0][1][1] + e[1][1][1] * e[0][1][0] + e[1][1][2];
          e[1][2][2] = e[1][1][1] * e[0][1][1] + e[1][1][2] * e[0][1][0];
          e[1][2][3] = e[1][1][2] * e[0][1][1];

          e[0][3][0] = e[0][2][0] * e[0][1][0] + e[0][2][1];
          e[0][3][1] = e[0][2][0] * e[0][1][1] + e[0][1][0] * e[0][2][1] + e[0][2][2] * two;
          e[0][3][2] = e[0][2][1] * e[0][1][1] + e[0][2][2] * e[0][1][0];
          e[0][3][3] = e[0][2][2] * e[0][1][1];

          if (mm == 3) return;

}
*/

