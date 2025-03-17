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
    double tau, p, p1, p2, pmin1, s1, s2, f1, f2;
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

  /*// iteration: refine
  double p1, p2, p, f1, f2;
  p1 = (100.0 * pi);
  f1 = s - ph2s_reg1(p1, h);
  if (fabs(f1) > xacc)
  {
    if (f1 > 0) // pT2s_reg1(p,h)< s ,the p1< expt p，so， p2=1.05*p1 p（p1,p2)
      p2 = (1.0 + f1 / s) * p1;
    else
      p2 = (1.0 - f1 / s) * p1;

    f2 = s - ph2s_reg1(p2, h);
    p = rtsec1(ph2s_reg1, h, s, p1, p2, f1, f2, xacc, iMAX);
  }
  else
    p = p1;

  return p;*/
}


double pT2v_reg1(double p, double T)
// specific volume in region 1
{
  double pi, tau;
  tau = r1Tstar / T;
  pi = p / r1pstar;
  return 0.001 * rgas_water * T * pi * gamma_pi_reg1(pi,tau) / p;
}

double pT2u_reg1(double p, double T)
// specific internal energy in region 1
{
  double pi, tau;
  tau = r1Tstar / T;
  pi = p / r1pstar;
  // return rgas_water * T *
  //       (tau * gammatau_reg1(pi,tau) - pi * gammapi_reg1(pi,tau));

  double gammapi = 0.0;
  double gammatau = 0.0;
  //polys_solo_i_j(7.1 - pi, tau - 1.222, 34, IJn, i2soI, j2soJ, solo_ij_pow_reg1, &gammapi, &gammatau);

  polys_solo_i_j_reg1(pi, tau, &gammapi, &gammatau);
  return rgas_water * T * (tau * gammatau - pi * gammapi);
}

double pT2s_reg1(double p, double T)
// specific entropy in region 1
{

  double pi, tau;
  tau = r1Tstar / T;
  pi = p / r1pstar;
  // return rgas_water * (tau * gammatau_reg1(pi,tau) - gamma_reg1(pi,tau));

  double gamma = 0.0;
  double gammatau = 0.0;

  // polys_0_j(7.1 - pi, tau - 1.222,34,IJn,&gamma,&gammatau);
  polys_solo_0_j_reg1(pi, tau, &gamma, &gammatau);  
  return rgas_water * (tau * gammatau - gamma);
}

double pT2h_reg1(double p, double T)
// specific enthalpy in region 1
{
  double pi, tau;
  tau = r1Tstar / T;
  pi = p / r1pstar;
  return rgas_water * T * tau * gamma_tau_reg1(pi,tau);

  //double gammatau = poly_solo_j(7.1 - pi, tau - 1.222, 34, IJn, i2soI, j2soJ, solo_ij_pow_reg1);
  //return rgas_water * T * tau * gammatau;
}

double pT2cp_reg1(double p, double T)
// specific isobaric heat capacity in region 1
{
  double pi, tau;
  tau = r1Tstar / T;
  pi = p / r1pstar;
   return -rgas_water * tau * tau * gamma_tautau_reg1(pi, tau);

  //double gammatautau = poly_solo_jj(7.1 - pi, tau - 1.222, 34, IJn, i2soI, j2soJ, solo_ij_pow_reg1);
  //return -rgas_water * tau * tau * gammatautau;
}

double pT2cv_reg1(double p, double T)
// specific isochoric heat capacity in region 1
// cv in kJ/(kg K), T in K, p in MPa
{

  double pi, tau, a, b;
  tau = r1Tstar / T;
  pi = p / r1pstar;

  // a = -tau * tau * gammatautau_reg1(pi,tau);
  // b = gammapi_reg1(pi,tau) - tau * gammapitau_reg1(pi,tau);
  // b *= b;
  // return rgas_water * (a + b / gammapipi_reg1(pi,tau));

  double poly_pi = 0;
  double poly_pitau = 0;
  double poly_pipi = 0;
  double poly_tautau = 0;
  // polys_i_ii_ij_jj(7.1 - pi, tau - 1.222, 34, IJn, &poly_pi, &poly_pipi, &poly_pitau, &poly_tautau);

  //polys_solo_i_ii_ij_jj(7.1 - pi, tau - 1.222, 34, IJn, i2soI, j2soJ, solo_ij_pow_reg1, &poly_pi, &poly_pipi, &poly_pitau, &poly_tautau);
  polys_solo_i_ii_ij_jj_reg1(pi, tau, &poly_pi, &poly_pipi, &poly_pitau, &poly_tautau);
  a = -tau * tau * poly_tautau;
  b = poly_pi - tau * poly_pitau;
  b *= b;
  return rgas_water * (a + b / poly_pipi);
}

double pT2w_reg1(double p, double T)
// speed of sound in region 1
// w in m/s, T in K, p in Mpa
{
  double pi, tau, gammapi, a, b;
  tau = r1Tstar / T;
  pi = p / r1pstar;

  /*
  gammapi = gammapi_reg1(pi,tau);
  a = gammapi - tau * gammapitau_reg1(pi,tau);
  a *= a;
  b = a / (tau * tau * gammatautau_reg1(pi,tau));
  b = b - gammapipi_reg1(pi,tau);
  return gammapi * sqrt(1000.0 * rgas_water * T / b);
*/

  double poly_pi = 0;
  double poly_pitau = 0;
  double poly_pipi = 0;
  double poly_tautau = 0;
//  polys_i_ii_ij_jj(7.1 - pi, tau - 1.222, 34, IJn, &poly_pi, &poly_pipi, &poly_pitau, &poly_tautau);

  //polys_solo_i_ii_ij_jj(7.1 - pi, tau - 1.222, 34, IJn, i2soI, j2soJ, solo_ij_pow_reg1, &poly_pi, &poly_pipi, &poly_pitau, &poly_tautau);

  polys_solo_i_ii_ij_jj_reg1(pi, tau, &poly_pi, &poly_pipi, &poly_pitau, &poly_tautau);
  a = poly_pi - tau * poly_pitau;
  a *= a;
  b = a / (tau * tau * poly_tautau);
  b = b - poly_pipi;
  return poly_pi * sqrt(1000.0 * rgas_water * T / b);
}
