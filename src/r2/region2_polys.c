
//---------------------------------------------------------------------------
// IAPWS -IF 97 Basic Equation for Region 2:  Release : IF97-rev  August 2007
// The dimensionless Gibbs free energy gamma and its derivatives
//     REGION 2 G(p,T) EQUATIONS ：  Eq. (15), P13
//     the ideal-gas part ： Eq. (16)
// Cheng  Maohua
//-------------------------------

//  Initialize coefficients and exponents for region 2
#include <math.h>
#include "region2.h"

void solo_i_j_power_reg2(double vi, double vj, double soI_pow[], double soJ_pow[])
{
  //  [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 16, 18, 20, 21, 22, 23, 24];
  soI_pow[0] = vi;
  for (int k = 1; k <= 9; k++)
  {
    soI_pow[k] = soI_pow[k - 1] * vi;
  }
  double vi2 = vi * vi;
  soI_pow[10] = soI_pow[9] * vi2 * vi2 * vi2;
  soI_pow[11] = soI_pow[10] * vi2;
  soI_pow[12] = soI_pow[11] * vi2;
  for (int k = 13; k <= 16; k++)
  {
    soI_pow[k] = soI_pow[k - 1] * vi;
  }

  //  0, 1, 2, 3, 6, 4, 7, 36, 35, 16, 11, 25, 8, 13, 10, 14, 29, 50, 57, 20, 48, 21, 53, 39, 26, 40,  58,

  soJ_pow[0] = 1.0;                     // 0
  soJ_pow[1] = vj;                      // 1
  soJ_pow[2] = vj * vj;                 // 2
  soJ_pow[3] = vj * soJ_pow[2];         // 3
  soJ_pow[4] = soJ_pow[3] * soJ_pow[3]; // 6

  soJ_pow[5] = soJ_pow[3] * vj; // 4
  soJ_pow[6] = soJ_pow[4] * vj; // 7

  double J_pow8 = soJ_pow[6] * vj;
  double J_pow10 = J_pow8 * soJ_pow[2];
  double J_pow11 = J_pow8 * soJ_pow[3];
  double J_pow16 = J_pow10 * soJ_pow[4];
  double J_pow20 = J_pow10 * J_pow10;

  soJ_pow[7] = J_pow20 * J_pow16; // 36
  soJ_pow[8] = soJ_pow[7] / vj;   // 35
  soJ_pow[9] = J_pow16;           // 16
  soJ_pow[10] = J_pow11;          // 11

  soJ_pow[11] = soJ_pow[8] / J_pow10; // 25
  soJ_pow[12] = J_pow8;               // 8
  soJ_pow[13] = J_pow10 * soJ_pow[3]; // 13
  soJ_pow[14] = J_pow10;              // 10
  soJ_pow[15] = J_pow10 * soJ_pow[5]; // 14

  soJ_pow[16] = J_pow16 * soJ_pow[13];    // 29
  soJ_pow[17] = soJ_pow[7] * soJ_pow[15]; // 50
  soJ_pow[18] = soJ_pow[17] * soJ_pow[6]; // 57
  soJ_pow[19] = J_pow20;                  // 20
  soJ_pow[20] = soJ_pow[8] * soJ_pow[13]; // 48

  soJ_pow[21] = J_pow20 * vj;             // 21
  soJ_pow[22] = soJ_pow[17] * soJ_pow[3]; // 53
  soJ_pow[23] = soJ_pow[7] * soJ_pow[3];  // 39
  soJ_pow[24] = soJ_pow[11] * vj;         // 26
  soJ_pow[25] = J_pow20 * J_pow20;        // 40
  soJ_pow[26] = soJ_pow[18] * vj;         // 58

  // for (int k = 0; k < 27; k++)
  //{
  //   soJ_pow[k] = ipowsac(vj, soJ[k]);
  // }
}

static const int soI[17] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 16, 18, 20, 21, 22, 23, 24};

static int i2soI[43] = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 5, 5, 5,
                        6, 6, 6, 7, 7, 8, 9, 9, 9, 10, 10, 11, 12, 12, 12,
                        13, 14, 15, 16, 16, 16};

static const int soJ[27] = {0, 1, 2, 3, 6, 4, 7,
                            36, 35, 16, 11, 25, 8, 13, 10,
                            14, 29, 50, 57, 20, 48, 21, 53, 39, 26, 40, 58};

static int j2soJ[43] = {0, 1, 2, 3, 4, 1, 2, 5, 6, 7, 0, 1, 3,
                        4, 8, 1, 2, 3, 6, 3, 9, 8, 0, 10, 11, 12, 7, 13, 5, 14, 15,
                        16, 17, 18, 19, 8, 20, 21, 22, 23, 24, 25, 26};

static int r2j0[9] = {0, 1, -5, -4, -3,
    -2, -1, 2, 3};

static double n0[9] = {-9.6927686500217,
     10.086655968018,
     -0.005608791128302,
     0.071452738081455,
     -0.40710498223928,
     1.4240819171444,
     -4.383951131945,
     -0.28408632460772,
     0.021268463753307};

// Table 11 Page 13 The residual part gamma'r for  region 2, Eq. (17)
// 43 items
static IJnData IJn[] = {
{1, 0, -0.17731742473213E-02},
{1, 1, -0.17834862292358E-01},
{1, 2, -0.45996013696365E-01},
{1, 3, -0.57581259083432E-01},
{1, 6, -0.5032527872793E-01},
{2, 1, -0.33032641670203E-04},
{2, 2, -0.18948987516315E-03},
{2, 4, -0.39392777243355E-02},
{2, 7, -0.43797295650573E-01},
{2, 36, -2.6674547914087E-05},
{3, 0, 0.20481737692309E-07},
{3, 1, 0.43870667284435E-06},
{3, 3, -0.3227767723857E-04},
{3, 6, -0.15033924542148E-02},
{3, 35, -0.040668253562649},
{4, 1, -7.8847309559367E-10},
{4, 2, 1.2790717852285E-08},
{4, 3, 4.8225372718507E-07},
{5, 7, 2.2922076337661E-06},
{6, 3, -1.6714766451061E-11},
{6, 16, -2.1171472321355E-03},
{6, 35, -23.895741934104},
{7, 0, -5.905956432427E-18},
{7, 11, -1.2621808899101E-06},
{7, 25, -0.038946842435739},
{8, 8, 1.1256211360459E-11},
{8, 36, -8.2311340897998},
{9, 13, 1.9809712802088E-08},
{10, 4, 1.0406965210174E-19},
{10, 10, -1.0234747095929E-13},
{10, 14, -1.0018179379511E-09},
{16, 29, -8.0882908646985E-11},
{16, 50, 0.10693031879409},
{18, 57, -0.33662250574171},
{20, 20, 8.9185845355421E-25},
{20, 35, 3.0629316876232E-13},
{20, 48, -4.2002467698208E-06},
{21, 21, -5.9056029685639E-26},
{22, 53, 3.7826947613457E-06},
{23, 39, -1.2768608934681E-15},
{24, 26, 7.3087610595061E-29},
{24, 40, 5.5414715350778E-17},
{24, 58, -9.436970724121E-07},
};

// Ideal-gas part of fundamental equation for region 2
// Eq16 P13
double gamma0_reg2(double pi, double tau)

{
    double gamma0 = log(pi);
    for (int i = 0; i < 9; i++)
        gamma0 += n0[i] * ipowsac(tau, r2j0[i]);
    return gamma0;
}

double gamma0_pi_reg2(double pi)
// First derivative in pi of ideal-gas part of fundamental equation for region 2
{
    return 1.0 / pi;
}

double gamma0_pipi_reg2(double pi)
// Second derivative in pi of ideal-gas part of fundamental equation for region 2
{
    return -1.0 / pi / pi;
}

double gamma0_tau_reg2(double tau)
// First derivative in tau of ideal-gas part of fundamental equation for region 2
{
    double gamma0tau = 0.0;
    for (int i = 0; i < 9; i++)
        gamma0tau += n0[i] * r2j0[i] * ipowsac(tau, r2j0[i] - 1);
    return gamma0tau;
}

double gamma0_tautau_reg2(double pi, double tau)
// Second derivative in tau of ideal-gas part of fundamental equation for region 2
{
    double gamma0tautau = 0.0;
    for (int i = 0; i < 9; i++)
        gamma0tautau += n0[i] * r2j0[i] * (r2j0[i] - 1) * ipowsac(tau, r2j0[i] - 2);
    return gamma0tautau;
}

double gamma0_pitau_reg2()
// Second derivative in pi and tau of ideal-gas part of fundamental equation for region 2
{
    return 0.0;
}

//-----------------------------------------------------------------
//  polynominal with solo i,j
//--------------------------------------------------------------------------------

double gammar_reg2(double pi, double tau)
{
    return poly_solo(pi, tau - 0.5, 43, IJn, i2soI, j2soJ, solo_i_j_power_reg2);
}

double gammar_pi_reg2(double pi, double tau)
{
    return poly_solo_i(pi, tau - 0.5, 43, IJn, i2soI, j2soJ, solo_i_j_power_reg2);
}

double gammar_pipi_reg2(double pi, double tau)
{
    return poly_solo_ii(pi, tau - 0.5, 43, IJn, i2soI, j2soJ, solo_i_j_power_reg2);
}

double gammar_pitau_reg2(double pi, double tau)
{
    return poly_solo_ij(pi, tau - 0.5, 43, IJn, i2soI, j2soJ, solo_i_j_power_reg2);
}

double gammar_tau_reg2(double pi, double tau)
{
    return poly_solo_j(pi, tau - 0.5, 43, IJn, i2soI, j2soJ, solo_i_j_power_reg2);
}

double gammar_tautau_reg2(double pi, double tau)
{
    return poly_solo_jj(pi, tau - 0.5, 43, IJn, i2soI, j2soJ, solo_i_j_power_reg2);
}

//-----------------------------------------------------------------
//  multiple polynominal with solo i,j
//--------------------------------------------------------------------------------

void polys_solo_0_i_reg2(double pi, double tau, double *poly, double *poly_pi)
{
    polys_solo_0_i(pi, tau - 0.5, 43, IJn, i2soI, j2soJ, solo_i_j_power_reg2, poly, poly_pi);
}

void polys_solo_i_j_reg2(double pi, double tau, double *poly_pi, double *poly_tau)
{
    polys_solo_i_j(pi, tau - 0.5, 43, IJn, i2soI, j2soJ, solo_i_j_power_reg2, poly_pi, poly_tau);
}

void polys_solo_0_j_reg2(double pi, double tau, double *poly, double *poly_tau)
{
    polys_solo_0_j(pi, tau - 0.5, 43, IJn, i2soI, j2soJ, solo_i_j_power_reg2, poly, poly_tau);
}

void polys_solo_i_ij_reg2(double pi, double tau, double *poly_pi, double *poly_pitau)
{
    polys_solo_i_ij(pi, tau - 0.5, 43, IJn, i2soI, j2soJ, solo_i_j_power_reg2, poly_pi, poly_pitau);
}

void polys_solo_i_ii_reg2(double pi, double tau, double *poly_pi, double *poly_pipi)
{
    polys_solo_i_ii(pi, tau - 0.5, 43, IJn, i2soI, j2soJ, solo_i_j_power_reg2, poly_pi, poly_pipi);
}

void polys_solo_i_ij_jj_reg2(double pi, double tau, double *poly_pi, double *poly_pitau, double *poly_tautau)
{
    polys_solo_i_ij_jj(pi, tau - 0.5, 43, IJn, i2soI, j2soJ, solo_i_j_power_reg2, poly_pi, poly_pitau, poly_tautau);
}

void polys_solo_i_ii_ij_reg2(double pi, double tau, double *poly_pi, double *poly_pipi, double *poly_pitau)
{
    polys_solo_i_ii_ij(pi, tau - 0.5, 43, IJn, i2soI, j2soJ, solo_i_j_power_reg2, poly_pi, poly_pipi, poly_pitau);
}

void polys_solo_i_ii_ij_jj_reg2(double pi, double tau, double *poly_pi, double *poly_pipi, double *poly_pitau, double *poly_tautau)
{
    polys_solo_i_ii_ij_jj(pi, tau - 0.5, 43, IJn, i2soI, j2soJ, solo_i_j_power_reg2, poly_pi, poly_pipi, poly_pitau, poly_tautau);
}


double pT2s_reg2(double p, double T)
//
// specific entropy in region 2
// sreg2 in kJ/(kg K)
// temperature in K
// pressure in Mpa
//
{
    double tau = r2Tstar / T;
    double pi = p;
    // return rgas_water * (tau * (gamma0_tau_reg2(tau) + gammartau_reg2(pi,tau)) - (gamma0__reg2(pi,tau) + gammar_reg2(pi,tau)));

    // gamma0_
    double gamma0_ = log(pi);
    double gamma0_tau = 0.0;
    double gamma0_item;
    for (int i = 0; i < 9; i++)
    {
        gamma0_item = n0[i] * ipowsac(tau, r2j0[i]);
        gamma0_ += gamma0_item;
        gamma0_tau += gamma0_item * r2j0[i];
    }
    gamma0_tau /= tau;

    // gammar
    double gammar = 0;
    double gammar_tau = 0.0;
    polys_solo_0_j_reg2(pi, tau, &gammar, &gammar_tau);
    return rgas_water * (tau * (gamma0_tau + gammar_tau) - (gamma0_ + gammar));
}