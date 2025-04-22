/*
 The extended pairs
   (p,v)->T
   (T,v)->p  (T,h)->p  (T,s)->p
*/
#include <math.h>
#include <stdint.h>
#include "region5.h"

// Region 5  (p,v)->T using the secant method and refine adjust
// * p: pressure  MPa
// * v: specific volume m^3/kg
// * T: temperature  K
double pv2T_reg5(double p, double v)
{
    double T1, T2, f1, f, v1, v2, T;
    T1 = TMIN5;
    v1 = pT2v_reg5(p, T1);
    T2 = TMAX5;
    v2 = pT2v_reg5(p, T2);
    f = v - pT2v_reg5(p, T2);
    if ((v2 - v1) != 0.0)
    {
        T1 = T1 + (T2 - T1) * (v - v1) / (v2 - v1);
    }
    f1 = v - pT2v_reg5(p, T1);
    T = rtsec2(pT2v_reg5, p, v, T1, T2, f1, f, xacc, iMAX);
    if (T < TMIN5)
    {
        T = TMIN5;
    }
    else
    {
        if (T > TMAX5)
        {
            T = TMAX5;
        };
    };
    return (T);
}

// Region 5 (T,v)->p using the secant method
//  T: temperature  K
//  v: specific volume m^3/kg
//  p: pressure  MPa
double Tv2p_reg5(double T, double v)
{
    double p, p1, p2, v1, v2, f1, f2;
    int vbounded = 0;
    p1 = PMIN5; //
    v1 = pT2v_reg5(p1, T);
    if (v == v1)
        return (p1);
    p2 = 10 * p1;
    v2 = pT2v_reg5(p2, T);
    if (v == v2)
        return (p2);
    if ((v > v2) && (v < v1))
        vbounded = 1;
    while (vbounded == 0)
    {
        p1 = p2;
        v1 = v2;
        p2 = 5 * p1;
        if (p2 >= PMAX5)
        {
            p2 = PMAX5;
            vbounded = 1;
        }
        v2 = pT2v_reg5(p2, T);
        if (v == v2)
            return (p2);
        if ((v > v2) && (v < v1))
            vbounded = 1;
    }
    f2 = v - v2;
    p1 = p2 - (p2 - p1) * (v - v2) / (v1 - v2);
    if (p1 < PMIN5)
        p1 = PMIN5;
    f1 = v - pT2v_reg5(p1, T);
    if (fabs(f1) < xacc)
        return (p1);
    p = rtsec1(pT2v_reg5, T, v, p1, p2, f1, f2, xacc, iMAX);
    if (p > PMAX5)
        p = PMAX5;
    if (p < PMIN5)
        p = PMIN5;
    return (p);
}

//----------------------------------------------
//      (T,s)->p
//----------------------------------------------
double Ts2p_reg5(double T, double s)
{
    double p, p1, p2, s1, s2, f1, f2;
    p1 = PMIN5; //
    s1 = pT2s_reg5(p1, T);
    f1 = s - s1;
    p2 = PMAX5;
    s2 = pT2s_reg5(p2, T);
    f2 = s - s2;
    p1 = p2 - (p2 - p1) * (s - s2) / (s1 - s2);
    if (p1 < PMIN5)
        p1 = PMIN5;
    s1 = pT2s_reg5(p1, T);
    f1 = s - s1;
    if (fabs(f1) < xacc)
        return (p1);
    p = rtsec1(pT2s_reg5, T, s, p1, p2, f1, f2, xacc, iMAX);
    if (p > PMAX5)
        p = PMAX5;
    if (p < PMIN5)
        p = PMIN5;
    return (p);
}

// Region 5 (T,h)->p using the secant method
//  *  T: temperature  K
//  *  h: specific enthalpy kJ/kg
//  *  p: pressure  MPa
double Th2p_reg5(double T, double h)
{
    double p, p1, p2, h1, h2, f1, f2;
    p1 = PMIN5;
    h1 = pT2h_reg5(p1, T);
    p2 = PMAX5; //
    h2 = pT2h_reg5(p2, T);
    f2 = h - h2;
    p1 = p2 - (p2 - p1) * fabs(h - h2) / (h1 - h2);
    f1 = h - pT2h_reg5(p1, T);
    if (fabs(f1) < xacc)
        return (p1);
    p = rtsec1(pT2h_reg5, T, h, p1, p2, f1, f2, xacc, iMAX);
    if (p > PMAX5)
        p = PMAX5;
    if (p < PMIN5)
        p = PMIN5;
    return (p);
}

double ph2T_reg5(double p, double h)
{
   double T, T1, T2;
   double f1, f2;

   T = -1000;

   T1 = 0.5 * (2273.15 + 1073.15);

   f1 = h - pT2h_reg5(p, T1);
   if (fabs(f1) > xacc)
   {
      if (f1 > 0)
         T2 = (1.0 + f1 / h) * T1;
      else
         T2 = (1.0 - f1 / h) * T1;

      f2 = h - pT2h_reg5(p, T2);

      T = rtsec2(pT2h_reg5, p, h, T1, T2, f1, f2, xacc, iMAX);
   }
   else
      T = T1;

   if (T < TMIN5)
      T = TMIN5;
   else if (T > TMAX5)
      T = TMAX5;

   return (T);
}

double ps2T_reg5(double p, double s)
{
   double T, T1, T2, f1, f;

   T = -1000;

   T1 = 0.5 * (2273.15 + 1073.15); // ,Get initial value

   f1 = s - pT2s_reg5(p, T1);

   if (fabs(f1) > xacc)
   {
      if (f1 > 0)
         T2 = (1.0 + f1 / s) * T1;
      else
         T2 = (1.0 - f1 / s) * T1;

      f = s - pT2s_reg5(p, T2);
      T = rtsec2(pT2s_reg5, p, s, T1, T2, f1, f, xacc, iMAX);
   }
   else
      T = T1;

   if (T < TMIN5)
      T = TMIN5;
   else if (T > TMAX5)
      T = TMAX5;

   return T;
}

// helper for hs2preg5
double ph2s_reg5(double p, double h)
{
   double T = ph2T_reg5(p, h);
   double s = pT2s_reg5(p, T);
   return s;
}

double hs2p_reg5(double h, double s)
{
   double p, p1, p2, f1, f2;

   // TODO:  迭代初始值，可测试那个更好?
   // 也可以计算smin,smax,2元插值得到更接近的p1

   // 测试表明：更复杂的方法计算迭代初始数值并不更好
   //double  hmin =pT2h_reg5(PMIN5,1073.15);
   //double  hmax =pT2h_reg5(PMAX5,2273.15);
   //p1=PMIN5+(PMAX5-PMIN5)*(h-hmin)/(hmax-hmin);

   p1 = 0.5 * (PMIN5 + PMAX5);

   f1 = s - ph2s_reg5(p1, h);

   if (fabs(f1) > xacc)
   {
      if (f1 > 0)
         p2 = (1.0 + f1 / s) * p1;
      else
         p2 = (1.0 - f1 / s) * p1;

      f2 = s - ph2s_reg5(p2, h);

      p = rtsec1(ph2s_reg5, h, s, p1, p2, f1, f2, xacc, iMAX);
   }
   else
      p = p1;

   return (p);
}


double pT2v_reg5(double p, double T)
{ // t[K],	p[MPa]
	double pi = p / r5Pstar;
	double tau = r5Tstar / T;
	double a, v;
	a = pi * (gamma0_pi_reg5(pi) + gammar_pi_reg5(pi, tau));
	v = rgas_water * T / p * a;
	return v * 0.001;
}

double pT2u_reg5(double p, double T)
{ // t[K],	p[MPa]
	double pi = p / r5Pstar;
	double tau = r5Tstar / T;
	double a, u;
	a = tau * (gamma0_tau_reg5(tau) + gammar_tau_reg5(pi, tau)) - pi * (gamma0_pi_reg5(pi) + gammar_pi_reg5(pi, tau));
	u = rgas_water * T * a;
	return u;
}

double pT2s_reg5(double p, double T)
{ // t[K],	p[MPa]
	double pi = p / r5Pstar;
	double tau = r5Tstar / T;
	double a, s;
	a = tau * (gamma0_tau_reg5(tau) + gammar_tau_reg5(pi, tau)) - (gamma0_reg5(pi, tau) + gammar_reg5(pi, tau));
	s = rgas_water * a;
	return s;
}

double pT2h_reg5(double p, double T)
{ // t[K],	p[MPa]
	double pi = p / r5Pstar;
	double tau = r5Tstar / T;
	double a, h;
	a = tau * (gamma0_tau_reg5(tau) + gammar_tau_reg5(pi, tau));
	h = rgas_water * T * a;
	return h;
}

double pT2cp_reg5(double p, double T)
{ // t[K],	p[MPa]
	double pi = p / r5Pstar;
	double tau = r5Tstar / T;
	double a, cp;
	a = -tau*tau * (gamma0_tautau_reg5(tau) + gammar_tautau_reg5(pi, tau));
	cp = rgas_water * a;
	return cp;
}

double pT2cv_reg5(double p, double T)
{ // t[K],	p[MPa]
	double pi = p / r5Pstar;
	double tau = r5Tstar / T;
	double a, b, c; // for temp
	double cv;
	a = (-pow(tau, 2.0)) * (gamma0_tautau_reg5(tau) + gammar_tautau_reg5(pi, tau));
	b = 1.0 + pi * gammar_pi_reg5(pi, tau) - tau * pi * gammar_pitau_reg5(pi, tau);
	c = 1.0 - pow(pi, 2.0) * gammar_pipi_reg5(pi, tau);
	cv = rgas_water * (a - pow(b, 2.0) / c);
	return cv;
}

double pT2w_reg5(double p, double T)
{
	double w;
	double a, b, c, d, w2;
	double tau = r5Tstar / T;
	double pi = p / r5Pstar;
	double dgammar_pi = gammar_pi_reg5(pi, tau);

	a = pow(1.0 + pi * dgammar_pi, 2.0);
	b = 1.0 - pi * pi * gammar_pipi_reg5(pi, tau);
	c = 1.0 + pi * dgammar_pi - tau * pi * gammar_pitau_reg5(pi, tau);
	d = tau * tau * (gamma0_tautau_reg5(tau) + gammar_tautau_reg5(pi, tau));
	w2 = rgas_water * T * a / (b + (c * c) / d) * 1000.0;
	if (w2 < 0.0)
		w2 = 0.0;
	w = sqrt(w2);
	return w;
}
