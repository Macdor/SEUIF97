/*
 The extended pairs
   (p,v)->T
   (T,v)->p  (T,h)->p  (T,s)->p
*/
/* -------------------------------------------------------------------------
  IAPWS-IF97 Region1:
    Properties(9):
     T in K           temperature
     p in MPa         pressure
     v  in m^3/kg     volume
     u  in kJ/kg      internal energy
     s in kJ/(kg K)   entropy
     h  in kJ/kg      enthalpy
    cp  in kJ/(kg K)  specific isobaric heat capacity
    cv  in kJ/(kg K)  specific isochoric heat capacity
     w  in m/s        speed of sound
----------------------------------------------------------------*/
#include <math.h>
#include <stdint.h>
#include "region1.h"
#include "../r4/region4.h"

// Region 1  (p,v)->T using the secant method and refine adjust
// * p: pressure  MPa
// * v: specific volume m^3/kg
// * T: temperature  K
double pv2T_reg1(double p, double v)
{
    double T1, T2, f1, f, v1, v2, T;
    T1 = TMIN1;
    v1 = pT2v_reg1(p, T1);
    if ((p >= 16.5291643) && (p <= 100.0))
    {
        T2 = TMAX1;
    }
    else
    {
        T2 = TSat(p);
    };
    v2 = pT2v_reg1(p, T2);
    f = v - pT2v_reg1(p, T2);
    if ((v2 - v1) != 0.0)
    {
        T1 = T1 + (T2 - T1) * (v - v1) / (v2 - v1);
    }
    f1 = v - pT2v_reg1(p, T1);
    T = rtsec2(pT2v_reg1, p, v, T1, T2, f1, f, xacc, iMAX);
    if (T >= 300)
        if (fabs(f1) < xacc)
            return (T1);

    f = v - pT2v_reg1(p, T);
    if (fabs(f) > xacc)
    {
        int success = 0;
        int sum = 0;
        T1 = T;
        f1 = f;
        if (f < 0) // t大了
        {
            while (!success)
            {
                T1 -= 0.1;
                if (T1 < 273.15)
                {
                    T = 273.15;
                    break;
                };
                f1 = v - pT2v_reg1(p, T1);
                if (fabs(f1) < fabs(f))
                {
                    f = f1;
                    T = T1;
                }
                if (fabs(f1) < xacc)
                {
                    T = T1;
                    break;
                }
                sum += 1;
                if (sum > 10000)
                    success = 1;
            };
        }
        else
        {
            while (!success)
            {
                T1 += 0.1;
                if (T1 > 623.15)
                {
                    T = 623.15;
                    break;
                };
                f1 = v - pT2v_reg1(p, T1);
                if (fabs(f1) < fabs(f))
                {
                    f = f1;
                    T = T1;
                }
                if (100 * fabs(f1) < xacc)
                {
                    T = T1;
                    break;
                }
                sum += 1;
                if (sum > 10000)
                {
                    success = 1;
                }
            }
        }
    };
    // 1 在低温区，v变化很小，
    // 2 在低温，低压区，v变化很小，同时变化规律发生反方向变化， T加，v减
    if (T < 290)
    {
    };
    if (T < TMIN1)
    {
        T = TMIN1;
    }
    else
    {
        if (T > TMAX1)
        {
            T = TMAX1;
        };
    };
    return (T);
}

// Region 1  (T,v)->p using the secant method
//  T: temperature  K
//  v: specific volume m^3/kg
//  p: pressure  MPa
double Tv2p_reg1(double T, double v)
{
    double p1 = 0.3 * (pSat(T) + PMAX1);
    double p2 = 1.05 * p1;
    double f1 = v - pT2v_reg1(p1, T);
    double f = v - pT2v_reg1(p2, T);
    double p = rtsec1(pT2v_reg1, T, v, p1, p2, f1, f, xacc, iMAX);
    return p;
}

//----------------------------------------------
//      (T,s)->p
//----------------------------------------------
double Ts2p_reg1(double T, double s)
{
    double p, p1, p2, pmin1, s1, s2, f1, f2;
    pmin1 = pSat(T);
    p1 = pmin1;
    s1 = pT2s_reg1(p1, T);
    f1 = s - s1;
    p2 = PMAX1;
    s2 = pT2s_reg1(p2, T);
    f2 = s - s2;
    p1 = p2 - (p2 - p1) * (s - s2) / (s1 - s2);
    if (p1 < pmin1)
        p1 = pmin1;
    s1 = pT2s_reg1(p1, T);
    f1 = s - s1;
    if (fabs(f1) < xacc)
        return (p1);
    p = rtsec1(pT2s_reg1, T, s, p1, p2, f1, f2, xacc, iMAX);
    if (p > PMAX1)
        p = PMAX1;
    if (p < pmin1)
        p = pmin1;
    return (p);
}

// Region 1  (T,h)->p using the secant method
//  *  T: temperature  K
//  *  h: specific enthalpy kJ/kg
//  *  p: pressure  MPa
double Th2p_reg1(double T, double h)
{
    double pmin1 = pSat(T);
    double p1 = pmin1;
    double p2 = PMAX1; // p1 + stepa
    double h1 = pT2h_reg1(p1, T);
    if (fabs(h - h1) < xacc)
    {
        return p1;
    };
    double h2 = pT2h_reg1(p2, T);
    if (fabs(h - h2) < xacc)
    {
        return p2;
    }

    double f1 = h - pT2h_reg1(p1, T);
    double f = h - pT2h_reg1(p2, T);
    double p = rtsec1(pT2h_reg1, T, h, p1, p2, f1, f, xacc, iMAX);

    if (p > PMAX1)
    {
        p = PMAX1;
    }
    if (p < pmin1)
    {
        p = pmin1;
    }
    return p;
}

// helper for iter (h,s)->p
double ph2s_reg1(double p, double h)
{
  double T;
  T = ph2T_reg1(p, h);
  return pT2s_reg1(p, T);
}

//----------------------------------------------------------------
//  Backward equation p(h,s) for region 1
//----------------------------------------------------------------
double hs2p_reg1(double h, double s)
{
  // Page 5, Table 2 :
  // Initialize coefficients and exponents (H,S)->P for region 1
  IJnData IJn[] = {
      {0, 0, -.691997014660582},
      {0, 1, -.183612548787560e2},
      {0, 2, -.928332409297335e1},
      {0, 4, .659639569909906e2},
      {0, 5, -.162060388912024e2},

      {0, 6, .450620017338667e3},
      {0, 8, .854680678224170e3},
      {0, 14, .607523214001162e4},
      {1, 0, .326487682621856e2},
      {1, 1, -.269408844582931e2},

      {1, 4, -.319947848334300e3},
      {1, 6, -.928354307043320e3},
      {2, 0, .303634537455249e2},
      {2, 1, -.650540422444146e2},
      {2, 10, -.430991316516130e4},

      {3, 4, -.747512324096068e3},
      {4, 1, .730000345529245e3},
      {4, 4, .114284032569021e4},
      {5, 0, -.436407041874559e3}};
  double eta, sigma;
  double pi;
  eta = h / 3400.0 + 0.05;
  sigma = s / 7.6 + 0.05;
  pi= poly(eta,sigma , 19,IJn);
  return (100.0*pi);
}


// specific volume in region 1
double pT2v_reg1(double p, double T)
{
  double pi, tau;
  tau = r1Tstar / T;
  pi = p / r1pstar;
  return 0.001 * rgas_water * T * pi * gamma_pi_reg1(pi,tau) / p;
}

// specific internal energy in region 1
double pT2u_reg1(double p, double T)
{
  double pi, tau;
  tau = r1Tstar / T;
  pi = p / r1pstar;

  double gamma = 0.0;
  double gamma_pi = 0.0;
  double gamma_pipi = 0.0;
  double gamma_tau = 0.0;
  double gamma_tautau = 0.0;
  double gamma_pitau = 0.0;

  polys_0_i_ii_j_jj_ij_reg1(pi, tau, &gamma, &gamma_pi, &gamma_pipi, &gamma_tau, &gamma_tautau, &gamma_pitau);

  return rgas_water * T * (tau * gamma_tau - pi * gamma_pi);
}

// specific entropy in region 1
double pT2s_reg1(double p, double T)
{

  double pi, tau;
  tau = r1Tstar / T;
  pi = p / r1pstar;

  double gamma = 0.0;
  double gamma_pi = 0.0;
  double gamma_pipi = 0.0;
  double gamma_tau = 0.0;
  double gamma_tautau = 0.0;
  double gamma_pitau = 0.0;

  polys_0_i_ii_j_jj_ij_reg1(pi, tau, &gamma, &gamma_pi, &gamma_pipi, &gamma_tau, &gamma_tautau, &gamma_pitau); 
  return rgas_water * (tau * gamma_tau - gamma);
}

// specific enthalpy in region 1
double pT2h_reg1(double p, double T)
{
  double pi, tau;
  tau = r1Tstar / T;
  pi = p / r1pstar;
  return rgas_water * T * tau * gamma_tau_reg1(pi,tau);
}

// specific isobaric heat capacity in region 1
double pT2cp_reg1(double p, double T)
{
  double pi, tau;
  tau = r1Tstar / T;
  pi = p / r1pstar;
   return -rgas_water * tau * tau * gamma_tautau_reg1(pi, tau);
}

// specific isochoric heat capacity in region 1
// cv in kJ/(kg K), T in K, p in MPa
double pT2cv_reg1(double p, double T)
{

  double pi, tau, a, b;
  tau = r1Tstar / T;
  pi = p / r1pstar;

  double gamma = 0.0;
  double gamma_pi = 0.0;
  double gamma_pipi = 0.0;
  double gamma_tau = 0.0;
  double gamma_tautau = 0.0;
  double gamma_pitau = 0.0;

  polys_0_i_ii_j_jj_ij_reg1(pi, tau, &gamma, &gamma_pi, &gamma_pipi, &gamma_tau, &gamma_tautau, &gamma_pitau);
  a = -tau * tau * gamma_tautau;
  b = gamma_pi - tau * gamma_pitau;
  return rgas_water * (a + b * b / gamma_pipi);
}

// speed of sound in region 1
// w in m/s, T in K, p in Mpa
double pT2w_reg1(double p, double T)
{
  double pi, tau, a, b;
  tau = r1Tstar / T;
  pi = p / r1pstar;

  double gamma = 0.0;
  double gamma_pi = 0.0;
  double gamma_pipi = 0.0;
  double gamma_tau = 0.0;
  double gamma_tautau = 0.0;
  double gamma_pitau = 0.0;

  polys_0_i_ii_j_jj_ij_reg1(pi, tau, &gamma, &gamma_pi, &gamma_pipi, &gamma_tau, &gamma_tautau, &gamma_pitau);
  a = gamma_pi - tau * gamma_pitau;
  b = a * a  / (tau * tau * gamma_tautau);
  return gamma_pi * sqrt(1000.0 * rgas_water * T / (b - gamma_pipi));
}



//---------------------------------------------------------------------------
//  August 2007 IF97: IF97-Rev.pdf: P6-9
// IAPWS -IF 97 Backware Equation for Region 1:
//            Backward (p,h)->T, (p,s)->T
//    (P,S)->T (P,H)->T
//-------------------------------------------------------------
// Backward equation T(p,h) for region 1
//--------------------------------------------------------------
// Page 11, Table6 :Initialize coefficients and exponents (P,H)->T for region 1
double ph2T_reg1(double p, double h)
{
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
  eta = h / 2500.0;
  theta = poly(pi, eta + 1.0, 20, IJn);
  return (1.0 *theta);
}

//----------------------------------------------------------------
//  Backward equation T(p,s) for region 1
//----------------------------------------------------------------
// Page 12, Table 8 : Initialize coefficients and exponents (P,S)->T for region 1
double ps2T_reg1(double p, double s)
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
  sigma = s / 1.0;
  theta = poly(pi, sigma + 2.0, 20, IJn);
   return (1.0*theta);
}