//  Initialize coefficients and exponents for region 1
#include <math.h>
#include "region1.h"


//---------------------------------------------------------------------------
//  August 2007 IF97: IF97-Rev.pdf: P6-9
// IAPWS -IF 97 Backware Equation for Region 1:
//            Backward (p,h)->T, (p,s)->T
//    (P,S)->T (P,H)->T
//-------------------------------------------------------------
// Backward equation T(p,h) for region 1
//--------------------------------------------------------------
double ph2T_reg1(double p, double h)
{
  // Page 11, Table6 :Initialize coefficients and exponents (P,H)->T for region 1
  IJnData IJn[] = {
      {0, 0, -238.72489924521},
      {0, 1, 404.21188637945},
      {0, 2, 113.49746881718},
      {0, 6, -5.8457616048039},
      {0, 22, -1.528548241314E-04},
      {0, 32, -1.0866707695377E-06},
      {1, 0, -13.391744872602},
      {1, 1, 43.211039183559},
      {1, 2, -54.010067170506},
      {1, 3, 30.535892203916},
      {1, 4, -6.5964749423638},
      {1, 10, 9.3965400878363E-03},
      {1, 32, 1.157364750534E-07},
      {2, 10, -2.5858641282073E-05},
      {2, 32, -4.0644363084799E-09},
      {3, 10, 6.6456186191635E-08},
      {3, 32, 8.0670734103027E-11},
      {4, 32, -9.3477771213947E-13},
      {5, 32, 5.8265442020601E-15},
      {6, 32, -1.5020185953503E-17}};

  double pi, eta;
  double theta;
  pi = p / 1.0;
  eta = h / 2500.0 + 1.0;
  theta = poly(pi, eta, 20, IJn);
  return (1.0 *theta);

  /*double T1, T2, T, f1, f2;
  T1 = (1.0 * theta);
  f1 = h - pT2h_reg1(p, T1);
  if (fabs(f1) > xacc)
  {
    if (f1 > 0)
      T2 = (1.0 + f1 / h) * T1; // TODO： 1.05 用 1+f1/h 是不是更快，没有测试
    else
      T2 = (1.0 - f1 / h) * T1;

    f2 = h - pT2h_reg1(p, T2);

    T = rtsec2(pT2h_reg1, p, h, T1, T2, f1, f2, xacc, iMAX);
  }
  else
    T = T1;

  return T;*/
}

//----------------------------------------------------------------
//  Backward equation T(p,s) for region 1
//----------------------------------------------------------------
double ps2T_reg1(double p, double s)
// Page 12, Table 8 : Initialize coefficients and exponents (P,S)->T for region 1
{
  IJnData IJn[] = {
      {0, 0, 0.17478268058307e3},
      {0, 1, 0.34806930892873e2},
      {0, 2, 0.65292584978455e1},
      {0, 3, 0.33039981775489},
      {0, 11, -0.19281382923196e-6},

      {0, 31, -0.24909197244573e-22},
      {1, 0, -0.26107636489332},
      {1, 1, 0.22592965981586},
      {1, 2, -0.64256463395226e-1},
      {1, 3, 0.78876289270526e-2},

      {1, 12, 0.35672110607366e-9},
      {1, 31, 0.17332496994895e-23},
      {2, 0, 0.56608900654837e-3},
      {2, 1, -0.32635483139717e-3},
      {2, 2, 0.44778286690632e-4},

      {2, 9, -0.51322156908507e-9},
      {2, 31, -0.42522657042207e-25},
      {3, 10, 0.26400441360689e-12},
      {3, 32, 0.78124600459723e-28},
      {4, 32, -0.30732199903668e-30}};

  double pi, sigma;
  double theta;
  pi = p / 1.0;
  sigma = s / 1.0 + 2.0;
  theta = poly(pi, sigma, 20, IJn);
   return (1.0*theta);

  /*
  // iteration: refine
  double T1, T2, T, f1, f2;
  T1 = (1.0 * theta);
  f1 = s - pT2s_reg1(p, T1);
  if (fabs(f1) > xacc)
  {
    if (f1 > 0) // pT2s_reg1(p,T1)< s ,the T1< expt T，so， T2=1.05*T1 T（T1,T2)
      T2 = (1.0 + f1 / s) * T1;
    else
      T2 = (1.0 - f1 / s) * T1;

    f2 = s - pT2s_reg1(p, T2);
    T = rtsec2(pT2s_reg1, p, s, T1, T2, f1, f2, xacc, iMAX);
  }
  else
    T = T1;

  return T;*/
}

void solo_ij_pow_reg1(double vi,double vj,double *soI_pow,double *soJ_pow)
{
    // [0, 1, 2, 3, 4, 5, 8, 21, 23, 29, 30, 31, 32];
   
    soI_pow[0] = 1.0;
    for(int k=1; k<=5;k++) {
        soI_pow[k] = soI_pow[k - 1] * vi;
    }
    double vi2 = vi * vi;
    soI_pow[6] = soI_pow[5] * vi2 * vi; //8
    soI_pow[7] = soI_pow[6] * soI_pow[6] * soI_pow[5]; //21
    soI_pow[8] = soI_pow[7] * vi2; //23
    soI_pow[9] = soI_pow[6] * soI_pow[7];
    for(int k=10; k<=12;k++) {
        soI_pow[k] = soI_pow[k - 1] * vi;
    }

    //  -2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -3, 17, -4, 6, -5, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41,
    double vj2 = vj * vj;
    soJ_pow[0] = 1.0 / vj2; //-2
    soJ_pow[1] = 1.0 / vj; //-1
    soJ_pow[2] = 1.0; // 0
    soJ_pow[3] = vj; // 1
    soJ_pow[4] = vj2; //2
    soJ_pow[5]=soJ_pow[4]*vj; //3
    soJ_pow[6]=soJ_pow[5]*vj;  //4
    soJ_pow[7]=soJ_pow[6]*vj; //5
   
    double J_pow7 = vj2 * soJ_pow[7];
    double J_pow9 = soJ_pow[6] * soJ_pow[7];
    double J_pow10 = J_pow9 *vj;
  
    soJ_pow[8] = 1.0 / J_pow9; //-9
    soJ_pow[9] = 1.0 / J_pow7; // -7
    soJ_pow[10] = 1.0 / soJ_pow[5]; // -3

    soJ_pow[11] = J_pow7 * J_pow10; // 17
    soJ_pow[12] = 1.0 / soJ_pow[6]; // -4
    soJ_pow[13] = soJ_pow[7] * vj; //6
    soJ_pow[14] = 1.0 / soJ_pow[7]; //-5
    soJ_pow[15] = J_pow10; // 10

    soJ_pow[16] = 1.0 / J_pow7/ vj; //-8
    soJ_pow[17] = 1.0 / soJ_pow[15]/ vj; //-11
    soJ_pow[18] = 1.0 / soJ_pow[13]; //-6
    soJ_pow[19] = (1.0 / J_pow10) * soJ_pow[16] * soJ_pow[17]; // -29
    soJ_pow[20] = soJ_pow[19] * soJ_pow[0]; //-31

    soJ_pow[21] = soJ_pow[20] * soJ_pow[9]; //-38
    for(int k=22; k<=24; k++) {
        //  -39, -40, -41
        soJ_pow[k] = soJ_pow[k - 1] / vj;
    }
    
    
    /*for(int k=0; k<13; k++) {
        soI_pow[k] = ipowsac(vi, soI[k]);
    }
    for(int k=0; k<25; k++) {
        soJ_pow[k] = ipowsac(vj, soJ[k]);
    }*/
}

static const int soI[13] = {0, 1, 2, 3, 4, 5, 8, 21, 23, 29, 30, 31, 32};

static int i2soI[34] = {0, 0, 0, 0, 0,
                      0, 0, 0, 1, 1,
                      1, 1, 1, 1, 2,
                      2, 2, 2, 2, 3,
                      3, 3, 4, 4, 4,
                      5, 6, 6, 7, 8,
                      9, 10, 11, 12};

static const int soJ[25] = {-2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -3, 17, -4, 6, -5, 10,
                            -8, -11, -6, -29, -31, -38, -39, -40, -41};

static int j2soJ[34] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 5, 10, 2, 3, 5, 11, 12, 2,
                      13, 14, 0, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};

static IJnData IJn[] = {
    {0, -2, 0.14632971213167E+00},
    {0, -1, -0.84548187169114E+00},
    {0, 0, -0.37563603672040E+01},
    {0, 1, 0.33855169168385E+01},
    {0, 2, -0.95791963387872E+00},

    {0, 3, 0.15772038513228E+00},
    {0, 4, -0.16616417199501E-01},
    {0, 5, 0.81214629983568E-03},
    {1, -9, 0.28319080123804E-03},
    {1, -7, -0.60706301565874E-03},

    {1, -1, -0.18990068218419E-01},
    {1, 0, -0.32529748770505E-01},
    {1, 1, -0.21841717175414E-01},
    {1, 3, -0.52838357969930E-04},
    {2, -3, -0.47184321073267E-03},

    {2, 0, -0.30001780793026E-03},
    {2, 1, 0.47661393906987E-04},
    {2, 3, -0.44141845330846E-05},
    {2, 17, -0.72694996297594E-15},
    {3, -4, -0.31679644845054E-04},

    {3, 0, -0.28270797985312E-05},
    {3, 6, -0.85205128120103E-09},
    {4, -5, -0.22425281908000E-05},
    {4, -2, -0.65171222895601E-06},
    {4, 10, -0.14341729937924E-12},

    {5, -8, -0.40516996860117E-06},
    {8, -11, -0.12734301741641E-08},
    {8, -6, -0.17424871230634E-09},
    {21, -29, -0.68762131295531E-18},
    {23, -31, 0.14478307828521E-19},

    {29, -38, 0.26335781662795E-22},
    {30, -39, -0.11947622640071E-22},
    {31, -40, 0.18228094581404E-23},
    {32, -41, -0.93537087292458E-25}};

double gamma_reg1(double pi, double tau)
{
    return poly_solo(7.1 - pi, tau - 1.222, 34, IJn, i2soI, j2soJ, solo_ij_pow_reg1);
}

double gamma_pi_reg1(double pi, double tau)
{
    return -poly_solo_i(7.1 - pi, tau - 1.222, 34, IJn, i2soI, j2soJ, solo_ij_pow_reg1);
}

double gamma_tau_reg1(double pi, double tau)
{
    return poly_solo_j(7.1 - pi, tau - 1.222, 34, IJn, i2soI, j2soJ, solo_ij_pow_reg1);
}

double gamma_tautau_reg1(double pi, double tau)
{
    return poly_solo_jj(7.1 - pi, tau - 1.222, 34, IJn, i2soI, j2soJ, solo_ij_pow_reg1);
}

double gamma_pitau_reg1(double pi, double tau)
{
    return -poly_solo_ij(7.1 - pi, tau - 1.222, 34, IJn, i2soI, j2soJ, solo_ij_pow_reg1);
}

double gamma_pipi_reg1(double pi, double tau)
{
    return -poly_solo_ii(7.1 - pi, tau - 1.222, 34, IJn, i2soI, j2soJ, solo_ij_pow_reg1);
}

//------------------------------------------------------------------------------
//
//----------------------------------------------------------------------------

void polys_solo_0_i_reg1(double pi, double tau, double *poly, double *poly_pi)
{
    polys_solo_0_i(7.1 - pi, tau - 1.222, 34, IJn, i2soI, j2soJ, solo_ij_pow_reg1, poly, poly_pi);
    *poly_pi = -(*poly_pi);
}

void polys_solo_0_j_reg1(double pi, double tau, double *poly, double *poly_tau)
{
    polys_solo_0_j(7.1 - pi, tau - 1.222, 34, IJn, i2soI, j2soJ, solo_ij_pow_reg1, poly, poly_tau);
}

void polys_solo_i_j_reg1(double pi, double tau, double *poly_pi, double *poly_tau)
{
    polys_solo_i_j(7.1 - pi, tau - 1.222, 34, IJn, i2soI, j2soJ, solo_ij_pow_reg1, poly_pi, poly_tau);
    *poly_pi = -(*poly_pi);
}

void polys_solo_i_ii_reg1(double pi, double tau, double *poly_pi, double *poly_pipi)
{
    polys_solo_i_ii(7.1 - pi, tau - 1.222, 34, IJn, i2soI, j2soJ, solo_ij_pow_reg1, poly_pi, poly_pipi);
    *poly_pi = -(*poly_pi);
}

void polys_solo_i_ij_reg1(double pi, double tau, double *poly_pi, double *poly_pitau)
{
    polys_solo_i_ij(7.1 - pi, tau - 1.222, 34, IJn, i2soI, j2soJ, solo_ij_pow_reg1, poly_pi, poly_pitau);
    *poly_pi = -(*poly_pi);
    *poly_pitau = -(*poly_pitau);
}

void polys_solo_i_ii_ij_reg1(double pi, double tau, double *poly_pi, double *poly_pipi, double *poly_pitau)
{
    polys_solo_i_ii_ij(7.1 - pi, tau - 1.222, 34, IJn, i2soI, j2soJ, solo_ij_pow_reg1, poly_pi,poly_pipi, poly_pitau);
    *poly_pi = -(*poly_pi);
    *poly_pitau = -(*poly_pitau);
}

void polys_solo_i_ij_jj_reg1(double pi, double tau, double *poly_pi, double *poly_pitau, double *poly_tautau)
{
    polys_solo_i_ij_jj(7.1 - pi, tau - 1.222, 34, IJn, i2soI, j2soJ, solo_ij_pow_reg1, poly_pi,poly_pitau, poly_tautau);
    *poly_pi = -(*poly_pi);
    *poly_pitau = -(*poly_pitau);
}

void polys_solo_i_ii_ij_jj_reg1(double pi, double tau, double *poly_pi, double *poly_pipi, double *poly_pitau, double *poly_tautau)
{
    polys_solo_i_ii_ij_jj(7.1 - pi, tau - 1.222, 34, IJn, i2soI, j2soJ, solo_ij_pow_reg1, poly_pi, poly_pipi, poly_pitau, poly_tautau);
    *poly_pi = -(*poly_pi);
    *poly_pitau = -(*poly_pitau);
}