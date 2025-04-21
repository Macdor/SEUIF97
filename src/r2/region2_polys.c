
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

// First derivative in pi of ideal-gas part of fundamental equation for region 2
double gamma0_pi_reg2(double pi)
{
    return 1.0 / pi;
}

// Second derivative in pi of ideal-gas part of fundamental equation for region 2
double gamma0_pipi_reg2(double pi)
{
    return -1.0 / pi / pi;
}

// First derivative in tau of ideal-gas part of fundamental equation for region 2
double gamma0_tau_reg2(double tau)
{
    double gamma0tau = 0.0;
    for (int i = 0; i < 9; i++)
        gamma0tau += n0[i] * r2j0[i] * ipowsac(tau, r2j0[i] - 1);
    return gamma0tau;
}

// Second derivative in tau of ideal-gas part of fundamental equation for region 2
double gamma0_tautau_reg2(double pi, double tau)
{
    double gamma0tautau = 0.0;
    for (int i = 0; i < 9; i++)
        gamma0tautau += n0[i] * r2j0[i] * (r2j0[i] - 1) * ipowsac(tau, r2j0[i] - 2);
    return gamma0tautau;
}

// Second derivative in pi and tau of ideal-gas part of fundamental equation for region 2
double gamma0_pitau_reg2()
{
    return 0.0;
}

//-----------------------------------------------------------------
//  polynominal with solo i,j
//--------------------------------------------------------------------------------

double gammar_reg2(double pi, double tau)
{
    return poly(pi, tau - 0.5, 43, IJn);
}

double gammar_pi_reg2(double pi, double tau)
{
    return poly_i(pi, tau - 0.5, 43, IJn);
}

double gammar_pipi_reg2(double pi, double tau)
{
    return poly_ii(pi, tau - 0.5, 43, IJn);
}

double gammar_pitau_reg2(double pi, double tau)
{
    return poly_ij(pi, tau - 0.5, 43, IJn);
}

double gammar_tau_reg2(double pi, double tau)
{
    return poly_j(pi, tau - 0.5, 43, IJn);
}

double gammar_tautau_reg2(double pi, double tau)
{
    return poly_jj(pi, tau - 0.5, 43, IJn);
}

//-----------------------------------------------------------------
//  multiple polynominal with solo i,j
//--------------------------------------------------------------------------------

void polys_0_i_ii_j_jj_ij_reg2(double pi, double tau, double *gamma, double *gamma_pi, double *gamma_pipi, double *gamma_tau, double *gamma_tautau, double *gamma_pitau)
{
    polys_0_i_ii_j_jj_ij(pi, tau - 0.5, 43, IJn, gamma, gamma_pi, gamma_pipi, gamma_tau, gamma_tautau, gamma_pitau);
}
