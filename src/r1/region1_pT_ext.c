/*  Region1: The extended Properties:
  * pT_ext_reg1(p: f64, T: f64, o_id: i32) -> f64

  Properties(12 + 5):

   *  k: Isentropic exponent OKISE 11 -
   *  ipcec: Isobaric cubic expansion coefficient  OIPCEC 17 αv  1/K
   *  kt: Isothermal compressibility, OKT kT 1/MPa

   *  e: Specific exergy    kJ/kg
   *  f: Specific Helmholtz free energy OF 12  kJ/kg
   *  g: Specific Gibbs free energy  OG 13 kJ/kg

   *  joule : Joule-Thomson coefficient  OJTC  23   K/MPa
   *  ijoule: Isothermal Joule-Thomson coefficient OIJTC 22 kJ/(kg·MPa)
   *
   *  z: Compressibility factor   -

   *  dpdtcv:  Partial derivative (dP/dT)v  MPa/K
   *  dvdpct: Partial derivative (dV/dP)T  m³/(kg·MPa)
   *  dvdtcp: Partial derivative (dV/dT)p  m³/(kg·K)
  TODO:
   *  pc:  isochoric pressure coefficient  1/K
   *  batap ：Isothermal stress coefficient, kg/m³
   *  fi: Fugacity coefficient
   *  fu: Fugacity, MPa
   * alfap: relative pressure coefficient  1/K
*/

#include <math.h>
#include <stdint.h>
#include "region1.h"

// k:Isentropic exponent OKISE 11
// * k= -(v/p)*/1000*(dp/dv)s
double pT2k_reg1(double p, double T)
{
    double v = pT2v_reg1(p, T);
    double w = pT2w_reg1(p, T);
    return 1.0E-6 * w * w / v / p;
}

// αv: Isobaric cubic expansion coefficient / Isobaric volume expansion coefficient  1/K
// * αv=(1.0/v)*(∂v/∂T)p ( The book 2019,Page 28)
double pT2ipcec_reg1(double p, double T)
{
    double tau = r1Tstar / T;
    double pi = p / r1pstar;
    
    double gamma = 0.0;
    double gammapi = 0.0;
    double gammapipi = 0.0;
    double gammatau = 0.0;
    double gammatautau = 0.0;
    double gammapitau = 0.0;
  
    polys_0_i_ii_j_jj_ij_reg1(pi, tau, &gamma, &gammapi, &gammapipi, &gammatau, &gammatautau, &gammapitau);
    return (1.0 - tau * gammapitau / gammapi) / T;
}

/// kt Isothermal compressibility 1/Mpa
// * kt=-(1.0/V)*(dv/dp)T
double pT2kt_reg1(double p, double T)
{
    double tau = r1Tstar / T;
    double pi = p / r1pstar;
    
    double gamma = 0.0;
    double gammapi = 0.0;
    double gammapipi = 0.0;
    double gammatau = 0.0;
    double gammatautau = 0.0;
    double gammapitau = 0.0;
  
    polys_0_i_ii_j_jj_ij_reg1(pi, tau, &gamma, &gammapi, &gammapipi, &gammatau, &gammatautau, &gammapitau);
    return -gammapipi / gammapi / r1pstar;
}

/// the specific Gibbs free energy, Page 6, Eq(7)
/// *  g=R*t*gamma
double pT2g_reg1(double p, double T)
{
    double tau = r1Tstar / T;
    double pi = p / r1pstar;
    return rgas_water * T * gamma_reg1(tau, pi);
}

/// the Helmholtz Specific free energy:
/// *  f=u-T*s=R*T*(gamma-gamma_pi)
double pT2f_reg1(double p, double T)
{
    double tau = r1Tstar / T;
    double pi = p / r1pstar;
    
    double gamma = 0.0;
    double gammapi = 0.0;
    double gammapipi = 0.0;
    double gammatau = 0.0;
    double gammatautau = 0.0;
    double gammapitau = 0.0;
  
    polys_0_i_ii_j_jj_ij_reg1(pi, tau, &gamma, &gammapi, &gammapipi, &gammatau, &gammatautau, &gammapitau);
    return rgas_water * T * (gamma - pi * gammapi);
}

/// joule ： Joule-Thomson coefficient    K/MPa
/// *  (dt/dp)h
double pT2joule_reg1(double p, double T)
{
    double tau = r1Tstar / T;
    double pi = p / r1pstar;

    double gamma = 0.0;
    double gammapi = 0.0;
    double gammapipi = 0.0;
    double gammatau = 0.0;
    double gammatautau = 0.0;
    double gammapitau = 0.0;
  
    polys_0_i_ii_j_jj_ij_reg1(pi, tau, &gamma, &gammapi, &gammapipi, &gammatau, &gammatautau, &gammapitau);
    double v = rgas_water * T * gammapi / r1pstar;
    double cp = -rgas_water * tau * tau * gammatautau;
    double TCex_1 = -tau * gammapitau / gammapi;
    return (v / cp) * TCex_1;
}

//  Isothermal throttling coefficient / Isothermal Joule-Thomson coefficient
///  ijoule: kJ/(kg·MPa)
/// *  (dh/dp)T
double pT2ijoule_reg1(double p, double T)
{
    double tau = r1Tstar / T;
    double pi = p / r1pstar;
    return 0.001 * rgas_water * r1Tstar * gamma_pitau_reg1(tau, pi) / r1pstar;
}

/// z: Compressibility factor
double pT2z_reg1(double p, double T)
{
    double v = pT2v_reg1(p, T);
    return (1000.0 * p * v / rgas_water / T);
}

/// e  the specific exergy   kJ/kg  in region 1
double pT2e_reg1(double p, double T)
{
    double tau = r1Tstar / T;
    double pi = p / r1pstar;
    
    double gamma = 0.0;
    double gammapi = 0.0;
    double gammapipi = 0.0;
    double gammatau = 0.0;
    double gammatautau = 0.0;
    double gammapitau = 0.0;
  
    polys_0_i_ii_j_jj_ij_reg1(pi, tau, &gamma, &gammapi, &gammapipi, &gammatau, &gammatautau, &gammapitau);

    return rgas_water * (T * gamma + (T - 273.16) * (tau * gammatau - gamma));
}

/// Region 1 - (dp/dt)v  MPa/K
double pT2dpdtcv_reg1(double p, double T)
{
    double tau = r1Tstar / T;
    double pi = p / r1pstar;
    
    double gamma = 0.0;
    double gammapi = 0.0;
    double gammapipi = 0.0;
    double gammatau = 0.0;
    double gammatautau = 0.0;
    double gammapitau = 0.0;
  
    polys_0_i_ii_j_jj_ij_reg1(pi, tau, &gamma, &gammapi, &gammapipi, &gammatau, &gammatautau, &gammapitau);
    return r1pstar * (gammapitau * r1Tstar - gammapi * T) / (T * T * gammapipi);
}

/// dvdp Partial derivative (dV/dP)T  m3/(kg·MPa)
double pT2dvdpct_reg1(double p, double T)
{
    double tau = r1Tstar / T;
    double pi = p / r1pstar;
    return 0.001 * rgas_water * T * gamma_pipi_reg1(tau, pi) / r1pstar / r1pstar;
}

// (dv/dT)p m3/(kg.K)
double pT2dvdtcp_reg1(double p, double T)
{
    double tau = r1Tstar / T;
    double pi = p / r1pstar;
    
    double gamma = 0.0;
    double gammapi = 0.0;
    double gammapipi = 0.0;
    double gammatau = 0.0;
    double gammatautau = 0.0;
    double gammapitau = 0.0;
  
    polys_0_i_ii_j_jj_ij_reg1(pi, tau, &gamma, &gammapi, &gammapipi, &gammatau, &gammatautau, &gammapitau);
    return 0.001 * rgas_water * (gammapi - tau * gammapitau) / r1pstar;
}

// (dv/dh)p m3/(kg.K)
double pT2dvdhcp_reg1(double p, double T)
{
    double tau = r1Tstar / T;
    double pi = p / r1pstar;
    
    double gamma = 0.0;
    double gammapi = 0.0;
    double gammapipi = 0.0;
    double gammatau = 0.0;
    double gammatautau = 0.0;
    double gammapitau = 0.0;
  
    polys_0_i_ii_j_jj_ij_reg1(pi, tau, &gamma, &gammapi, &gammapipi, &gammatau, &gammatautau, &gammapitau);

    double cp = -rgas_water*tau*tau*gammatautau;
    double dvdtcp = 0.001*rgas_water*(gammapi - tau*gammapitau)/r1pstar;
    double dvdhcp = dvdtcp/cp;

    return dvdhcp;
}

// (dv/dp)h m3/(kg.K)
double pT2dvdpch_reg1(double p, double T)
{
    double tau = r1Tstar / T;
    double pi = p / r1pstar;
    
    double gamma = 0.0;
    double gammapi = 0.0;
    double gammapipi = 0.0;
    double gammatau = 0.0;
    double gammatautau = 0.0;
    double gammapitau = 0.0;
  
    polys_0_i_ii_j_jj_ij_reg1(pi, tau, &gamma, &gammapi, &gammapipi, &gammatau, &gammatautau, &gammapitau);

    double cp = -rgas_water*tau*tau*gammatautau;
    double dhdpct = rgas_water*r1Tstar*gammapitau/r1pstar;
    double dvdpct = 0.001*rgas_water*T*gammapipi/r1pstar/r1pstar;
    double dvdtcp = 0.001*rgas_water*(gammapi - tau*gammapitau)/r1pstar;
    double dvdhcp = dvdtcp/cp;
    double dvdpch = dvdpct - dvdhcp*dhdpct;
    return dvdpch;
}
