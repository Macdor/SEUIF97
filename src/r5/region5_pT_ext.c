/*
Region 5: the extended properties

*/
#include <math.h>
#include <stdint.h>
#include "region5.h"


// k: isentropic exponent  T K,	p MPa OKISE 11
//    k= -(v/p)*/1000*(dp/dv)s
double pT2k_reg5(double p, double T)
{
    double k;
    double a, b, c, d, e;
    double tau = r5Tstar / T;
    double pi = p / r5Pstar;
    double dgammar_pi = gammar_pi_reg5(pi, tau);

    a = ipowsac(1.0 + pi * dgammar_pi, 2);
    b = 1.0 - pi * pi * gammar_pipi_reg5(pi, tau);
    c = 1.0 + pi * (dgammar_pi - tau * gammar_pitau_reg5(pi, tau));
    d = tau * tau * (gamma0_tautau_reg5(tau) + gammar_tautau_reg5(pi, tau));
    e = pi * (gamma0_pi_reg5(pi) + gammar_pi_reg5(pi, tau));
    k = a / (b + (c * c) / d) / e;
    return k;
}

// kt: Isothermal compressibility, [1/MPa]
//* kt=-(1.0/V)*(dv/dp)T
double pT2kt_reg5(double p, double T)
{
    double tau = r5Tstar / T;
    double pi = p / r5Pstar;
    double gamma_pipi = gamma0_pipi_reg5(pi) + gammar_pipi_reg5(pi, tau);
    double gamma_pi = gamma0_pi_reg5(pi) + gammar_pi_reg5(pi, tau);
    return -(gamma_pipi / gamma_pi);
}

//  the specific Gibbs free energy, Page 6, Eq(7)
//   g=R*t*gamma
double pT2g_reg5(double p, double T)
{
    double tau = r5Tstar / T;
    double pi = p / r5Pstar;
    return rgas_water * T * (gamma0_reg5(pi, tau) + gammar_reg5(pi, tau));
}

// the Helmholtz Specific free energy:
//  f=u-T*s=R*T(gamma-gamma_pi)
double pT2f_reg5(double p, double T)
{
    double tau = r5Tstar / T;
    double pi = p / r5Pstar;
    return rgas_water * T * (gamma0_reg5(pi, tau) + gammar_reg5(pi, tau) - p * (gamma0_pi_reg5(pi) + gammar_pi_reg5(pi, tau)));
}

// e
double pT2e_reg5(double p, double T)
{
    double tau = r5Tstar / T;
    double pi = p / r5Pstar;
    double gamma = gamma0_reg5(pi, tau) + gammar_reg5(pi, tau);
    double gamma_tau = gamma0_tau_reg5(tau) + gammar_tau_reg5(pi, tau);
    return rgas_water * (T * gamma + (T - Tt) * (tau * gamma_tau - gamma));
}

// coefficient of thermal expansion,
// ipcec-Isobaric cubic expansion coefficient  1/K
// α=(1.0/V)*(dv/dT)p
double pT2ipcec_reg5(double p, double T)
{
    double d = 1.0 / pT2v_reg5(p, T);
    return d * pT2dvdtcp_reg5(p, T);
}

// joule ： Joule-Thomson coefficient    K/MPa
//*  (dt/dp)h
double pT2joule_reg5(double p, double T)
{
    double tau = r5Tstar / T;
    double pi = p / r5Pstar;
    double gamma_pi = gamma0_pi_reg5(pi) + gammar_pi_reg5(pi, tau);
    double gamma_tautau = gamma0_tautau_reg5(tau) + gammar_tautau_reg5(pi, tau);
    double gamma_pitau = gamma0_pitau_reg5() + gammar_pitau_reg5(pi, tau);
    double v = 0.001 * rgas_water * T * gamma_pi;
    double cp = -rgas_water * tau * tau * gamma_tautau;
    double TCex_1 = -tau * gamma_pitau / gamma_pi;
    return 1000.0 * (v / cp) * TCex_1;
}

//  Isothermal throttling coefficient
// iJTC Isothermal Joule-Thomson coefficient kJ/(kg·MPa)
//*  (dh/dp)t
double pT2ijoule_reg5(double p, double T)
{
    double tau = r5Tstar / T;
    double pi = p / r5Pstar;
    double gamma_pitau = gamma0_pitau_reg5() + gammar_pitau_reg5(pi, tau);
    return 0.001 * rgas_water * r5Tstar * gamma_pitau;
}

// z: Compressibility factor
double pT2z_reg5(double p, double T)
{
    double v = pT2v_reg5(p, T);
    return 1000.0 * p * v / rgas_water / T;
}

// (dp/dt)v MPa/K
double pT2dpdtcv_reg5(double p, double T)
{
    double tau = r5Tstar / T;
    double pi = p / r5Pstar;
    double gamma_pi = gamma0_pi_reg5(pi) + gammar_pi_reg5(pi, tau);
    double gamma_pitau = gamma0_pitau_reg5() + gammar_pitau_reg5(pi, tau);
    double gamma_pipi = gamma0_pipi_reg5(pi) + gammar_pipi_reg5(pi, tau);
    return r5Pstar * (gamma_pitau * r5Tstar - gamma_pi * T) / (T * T * gamma_pipi);
}

// (dv/dp)T  m3/MPa
double pT2dvdpct_reg5(double p, double T)
{
    double tau = r5Tstar / T;
    double pi = p / r5Pstar;
    double gamma_pipi = gamma0_pipi_reg5(pi) + gammar_pipi_reg5(pi, tau);
    return 0.001 * rgas_water * T * gamma_pipi;
}

// dvdt, cp m3/(kg.K)
double pT2dvdtcp_reg5(double p, double T)
{
    double tau = r5Tstar / T;
    double pi = p / r5Pstar;
    double gamma_pi = gamma0_pi_reg5(pi) + gammar_pi_reg5(pi, tau);
    double gamma_pitau = gamma0_pitau_reg5() + gammar_pitau_reg5(pi, tau);
    return 0.001 * rgas_water * (gamma_pi - tau * gamma_pitau);
}

// (dv/dh)p m3/(kg.K)
double pT2dvdhcp_reg5(double p, double T)
{
    double dvdtcp = pT2dvdtcp_reg5(p, T);
    double cp = pT2cp_reg5(p, T);
    double dvdhcp = dvdtcp/cp;
    return dvdhcp;
}

// (dv/dp)h m3/(kg.K)
double pT2dvdpch_reg5(double p, double T)
{
    double dvdpct = pT2dvdpct_reg5(p, T);
    double dvdtcp = pT2dvdtcp_reg5(p, T);
    double dhdpct = pT2ijoule_reg5(p, T);
    double cp = pT2cp_reg5(p, T);
    double dvdpch = dvdpct - dvdtcp*dhdpct/cp;
    return dvdpch;
}
