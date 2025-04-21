/*  Region1: The extended Properties:
  * pT_ext_reg2(p: f64, T: f64, o_id: i32) -> f64

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
#include "region2.h"

// k:Isentropic exponent OKISE 11
// * k= -(v/p)*/1000*(dp/dv)s
double pT2k_reg2(double p, double T)
{
    double v = pT2v_reg2(p, T);
    double w = pT2w_reg2(p, T);
    return 1.0E-6 * w * w / v / p;
}

// αv: Isobaric cubic expansion coefficient / Isobaric volume expansion coefficient  1/K
// * αv=(1.0/v)*(∂v/∂T)p ( The book 2019,Page 32)
double pT2ipcec_reg2(double p, double T)
{
    double tau = r2Tstar / T;
    double pi = p / r2pstar;
    
    double gammar = 0.0;
    double gammarpi = 0.0;
    double gammarpipi = 0.0;
    double gammartau = 0.0;
    double gammartautau = 0.0;
    double gammarpitau = 0.0;
  
    polys_0_i_ii_j_jj_ij_reg1(pi, tau, &gammar, &gammarpi, &gammarpipi, &gammartau, &gammartautau, &gammarpitau);
    
    return (1.0 - pi * tau * gammarpitau / (1.0 + p * gammarpi)) / T;
}

/// kt Isothermal compressibility 1/Mpa
// * kt=-(1.0/V)*(dv/dp)T the book 2019 p32
double pT2kt_reg2(double p, double T)
{
    double tau = r2Tstar / T;
    double pi = p / r2pstar;

    double gammar = 0.0;
    double gammarpi = 0.0;
    double gammarpipi = 0.0;
    double gammartau = 0.0;
    double gammartautau = 0.0;
    double gammarpitau = 0.0;
  
    polys_0_i_ii_j_jj_ij_reg1(pi, tau, &gammar, &gammarpi, &gammarpipi, &gammartau, &gammartautau, &gammarpitau);
    return (1 - pi * pi * gammarpipi) / (1 + pi * gammarpi) / p;
}

/// the specific Gibbs free energy, Page 6, Eq(7)
/// *  g=R*t*gamma
double pT2g_reg2(double p, double T)
{
    double tau = r2Tstar / T;
    double pi = p / r2pstar;
    double gamma0 = gamma0_reg2(pi, tau);
    double gammar = gammar_reg2(pi, tau);
    return rgas_water * T * (gamma0 + gammar);
}

/// the Helmholtz Specific free energy:
/// *  f=u-T*s=R*T*(gamma-gamma_pi)
double pT2f_reg2(double p, double T)
{
    double tau = r2Tstar / T;
    double pi = p / r2pstar;

    double gamma0 = gamma0_reg2(pi, tau);
    double gamma0_pi = gamma0_pi_reg2(pi);

    double gammar = 0.0;
    double gammarpi = 0.0;
    double gammarpipi = 0.0;
    double gammartau = 0.0;
    double gammartautau = 0.0;
    double gammarpitau = 0.0;
  
    polys_0_i_ii_j_jj_ij_reg1(pi, tau, &gammar, &gammarpi, &gammarpipi, &gammartau, &gammartautau, &gammarpitau);

    return rgas_water * T * (gamma0 + gammar - pi * (gamma0_pi + gammarpi));
}

/// joule ： Joule-Thomson coefficient    K/MPa (10^-3*K/kPa)
/// *  (dt/dp)h  The book: Table 22 Page  344
double pT2joule_reg2(double p, double T)
{
    double tau = r2Tstar / T;
    double pi = p / r2pstar;

    // let g1pi: f64 = gamma_pi_reg2(tau, pi);
    // let v: f64 = 0.001 * RGAS_WATER * T * g1pi / r2pstar;
    //  let cp: f64 = -RGAS_WATER * tau * tau * gamma_tautau_reg2(tau, pi);
    //  let TCex_1: f64 = -tau * gamma_pitau_reg2(tau, pi) / g1pi;
    //  (v / cp) * TCex_1

    double gammar = 0.0;
    double gammarpi = 0.0;
    double gammarpipi = 0.0;
    double gammartau = 0.0;
    double gammartautau = 0.0;
    double gammarpitau = 0.0;
  
    polys_0_i_ii_j_jj_ij_reg1(pi, tau, &gammar, &gammarpi, &gammarpipi, &gammartau, &gammartautau, &gammarpitau);

    double gamma0pi = gamma0_pi_reg2(pi);
    double gamma0pitau = gamma0_pitau_reg2();
    double gamma0tautau = gamma0_tautau_reg2(pi, tau);

    double gammapi = gamma0pi + gammarpi;
    double gammapitau = gamma0pitau + gammarpitau;
    double gammatautau = gamma0tautau + gammartautau;

    double v = rgas_water * T * gammapi;
    double cp = -rgas_water * tau * tau * gammatautau;
    double TCex_1 = -tau * gammapitau / gammapi;
    return (v / cp) * TCex_1;
}

//  Isothermal throttling coefficient / Isothermal Joule-Thomson coefficient
///  ijoule: kJ/(kg·MPa)
/// *  (dh/dp)T
double pT2ijoule_reg2(double p, double T)
{
    double tau = r2Tstar / T;
    double pi = p / r2pstar;
    double gamma_pitau = gammar_pitau_reg2(pi, tau);

    gamma_pitau += gamma0_pitau_reg2();
    return 0.001 * rgas_water * r2Tstar * gamma_pitau / r2pstar;
}

/// z: Compressibility factor
double pT2z_reg2(double p, double T)
{
    double v = pT2v_reg2(p, T);
    return (1000.0 * p * v / rgas_water / T);
}

/// e  the specific exergy   kJ/kg  in region 1
double pT2e_reg2(double p, double T)
{
    double tau = r2Tstar / T;
    double pi = p / r2pstar;

    double gammar = 0.0;
    double gammarpi = 0.0;
    double gammarpipi = 0.0;
    double gammartau = 0.0;
    double gammartautau = 0.0;
    double gammarpitau = 0.0;
  
    polys_0_i_ii_j_jj_ij_reg1(pi, tau, &gammar, &gammarpi, &gammarpipi, &gammartau, &gammartautau, &gammarpitau);

    double gamma = gamma0_reg2(pi, tau) + gammar;
    double gammatau = gamma0_tau_reg2(tau) + gammartau;
    return rgas_water * (T * gamma + (T - 273.16) * (tau * gammatau - gamma));
}

/// Region 1 - (dp/dt)v  MPa/K
double pT2dpdtcv_reg2(double p, double T)
{
    double tau = r2Tstar / T;
    double pi = p / r2pstar;

    double gammar = 0.0;
    double gammarpi = 0.0;
    double gammarpipi = 0.0;
    double gammartau = 0.0;
    double gammartautau = 0.0;
    double gammarpitau = 0.0;
  
    polys_0_i_ii_j_jj_ij_reg1(pi, tau, &gammar, &gammarpi, &gammarpipi, &gammartau, &gammartautau, &gammarpitau);

    double gamma_pi = gamma0_pi_reg2(pi) + gammarpi;
    double gamma_pipi = gamma0_pipi_reg2(pi) + gammarpipi;
    double gamma_pitau = gamma0_pitau_reg2() + gammarpitau;
    return r2pstar * (gamma_pitau * r2Tstar - gamma_pi * T) / (T * T * gamma_pipi);
}

/// dvdp Partial derivative (dV/dP)T  m3/(kg·MPa)
double pT2dvdpct_reg2(double p, double T)
{
    double tau = r2Tstar / T;
    double pi = p / r2pstar;
    double gamma_pipi = gammar_pipi_reg2(pi, tau);
    gamma_pipi += gamma0_pipi_reg2(pi);
    return 0.001 * rgas_water * T * gamma_pipi / r2pstar / r2pstar;
}

// (dv/dT)p m3/(kg.K)
double pT2dvdtcp_reg2(double p, double T)
{
    double tau = r2Tstar / T;
    double pi = p / r2pstar;

    double gammar = 0.0;
    double gammarpi = 0.0;
    double gammarpipi = 0.0;
    double gammartau = 0.0;
    double gammartautau = 0.0;
    double gammarpitau = 0.0;
  
    polys_0_i_ii_j_jj_ij_reg1(pi, tau, &gammar, &gammarpi, &gammarpipi, &gammartau, &gammartautau, &gammarpitau);

    double gammapi = gamma0_pi_reg2(pi) + gammarpi;
    double gammapitau = gamma0_pitau_reg2() + gammarpitau;
    return 0.001 * rgas_water * (gammapi - tau * gammapitau) / r2pstar;
}

// (dv/dh)p m3/(kg.K)
double pT2dvdhcp_reg2(double p, double T)
{
    double tau = r2Tstar / T;
    double pi = p / r2pstar;
    
    double gammar = 0.0;
    double gammarpi = 0.0;
    double gammarpipi = 0.0;
    double gammartau = 0.0;
    double gammartautau = 0.0;
    double gammarpitau = 0.0;
  
    polys_0_i_ii_j_jj_ij_reg2(pi, tau, &gammar, &gammarpi, &gammarpipi, &gammartau, &gammartautau, &gammarpitau);

    double gammapi = gamma0_pi_reg2(pi) + gammarpi;
    double gammatautau = gamma0_tautau_reg2(pi, tau) + gammartautau;
    double gammapitau = gamma0_pitau_reg2() + gammarpitau;

    double cp = -rgas_water*tau*tau*gammatautau;
    double dvdtcp = 0.001*rgas_water*(gammapi - tau*gammapitau)/r2pstar;
    double dvdhcp = dvdtcp/cp;

    return dvdhcp;
}

// (dv/dp)h m3/(kg.K)
double pT2dvdpch_reg2(double p, double T)
{
    double tau = r2Tstar / T;
    double pi = p / r2pstar;
    
    double gammar = 0.0;
    double gammarpi = 0.0;
    double gammarpipi = 0.0;
    double gammartau = 0.0;
    double gammartautau = 0.0;
    double gammarpitau = 0.0;
  
    polys_0_i_ii_j_jj_ij_reg2(pi, tau, &gammar, &gammarpi, &gammarpipi, &gammartau, &gammartautau, &gammarpitau);

    double gammapi = gamma0_pi_reg2(pi) + gammarpi;
    double gammapipi = gamma0_pipi_reg2(pi) + gammarpipi;
    double gammatautau = gamma0_tautau_reg2(pi, tau) + gammartautau;
    double gammapitau = gamma0_pitau_reg2() + gammarpitau;

    double cp = -rgas_water*tau*tau*gammatautau;
    double dhdpct = rgas_water*r2Tstar*gammapitau/r2pstar;
    double dvdpct = 0.001*rgas_water*T*gammapipi/r2pstar/r2pstar;
    double dvdtcp = 0.001*rgas_water*(gammapi - tau*gammapitau)/r2pstar;
    double dvdhcp = dvdtcp/cp;
    double dvdpch = dvdpct - dvdhcp*dhdpct;
    return dvdpch;
}
