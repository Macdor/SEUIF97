/*
 The extended pairs
   (p,v)->T
   (T,v)->p  (T,h)->p  (T,s)->p
*/
#include <math.h>
#include <stdint.h>
#include "region2.h"
#include "../r4/region4.h"

/// the helper for the extended input pair
double p2Tmin_reg2(double p)
{
    double Tmin = TMIN2;
    if (p > 0.0 && p < PMIN1)
    {
        Tmin = TMIN2;
    }
    else
    {
        if (p >= PMIN1 && p <= pSat(TMAX1))
        {
            Tmin = TSat(p);
        }
        else
        {
            Tmin = B23_p2T(p);
        }
    }
    return Tmin;
}

double T2pmax_reg2(double T)
{
    double pmax = PMAX2;
    if (T >= TMIN1 && T <= TMAX1)
    {
        pmax = pSat(T);
    }
    else
    {
        if (T > TMAX1 && T <= TMAX3)
        {
            pmax = B23_T2p(T);
        }
        else
        {
            if (T > TMAX3 && T <= TMAX2)
            {
                pmax = PMAX3;
            };
        };
    };
    return pmax;
}

// Region 2  (p,v)->T using the secant method and refine adjust
// * p: pressure  MPa
// * v: specific volume m^3/kg
// * T: temperature  K
double pv2T_reg2(double p, double v)
{
    double T1, T2, f1, f, v1, v2, T;
    double Tmin2 = p2Tmin_reg2(p);
    T1 = Tmin2;
    v1 = pT2v_reg2(p, T1);
    f1 = v - v1;
    T2 = TMAX2;
    T2 = 1.1 * T1; // fast
    v2 = pT2v_reg2(p, T2);
    f = v - v2;
    T1 = T1 + (T2 - T1) * fabs(v - v1) / (v2 - v1);
    v1 = pT2v_reg2(p, T1);
    f1 = v - v1;
    if (fabs(f1) < xacc)
        return (T1);
    T = rtsec2(pT2v_reg2, p, v, T1, T2, f1, f, xacc, iMAX);
    if (T < Tmin2)
        T = Tmin2;
    if (T > TMAX2)
        T = TMAX2;

    // zoom
    double v0 = pT2v_reg2(p, T);
    int i = 0;
    int success = 1;
    if (fabs(v0 - v) > xacc)
        success = 0;
    while (success == 0)
    {
        if (v0 > v)
        {
            i += 1;
            T = T - 0.001;
            if (T < Tmin2)
                v0 = pT2v_reg2(p, T);
        }
        else
        {
            i += 1;
            T = T + 0.001;
            if (T > TMAX2)
                v0 = pT2v_reg2(p, T);
        }
        if (i < 1000)
            success = 1;
    }
    return (T);
}

// Region 2(T,v)->p using the secant method
//  T: temperature  K
//  v: specific volume m^3/kg
//  p: pressure  MPa
double Tv2p_reg2(double T, double v)
{
    double pmax2, p, p1, p2, pmid, v1, v2, f1, f2, stepa, stepm;
    int bounded = 0;
    p1 = PMIN2;
    pmax2 = T2pmax_reg2(T);
    stepa = 1.0;
    stepm = 5.0;
    v1 = pT2v_reg2(p1, T);
    if (fabs(v - v1) < xacc)
        return (p1);
    p2 = p1 * stepm;
    v2 = pT2v_reg2(p2, T);
    if (fabs(v - v2) < xacc)
        return (p2);
    if ((v > v2) && (v < v1))
        bounded = 1;
    while (bounded == 0)
    {
        p1 = p2;
        v1 = v2;
        if (p1 > 1.0)
            p2 = p1 + stepa;
        if (p1 < 1.0)
            p2 = p1 * stepm;
        if (p2 >= pmax2)
        {
            p2 = pmax2;
            bounded = 1;
        }
        v2 = pT2v_reg2(p2, T);
        if (fabs(v - v2) < xacc)
            return (p2);
        if ((v > v2) && (v < v1))
            bounded = 1;
    }
    f2 = v - v2;
    pmid = p2 - (p2 - p1) * (v - v2) / (v1 - v2);
    if (v < pT2v_reg2(pmid, T))
        p1 = pmid;
    if (p1 < PMIN2)
        p1 = PMIN2;
    f1 = v - pT2v_reg2(p1, T);
    if (fabs(f1) < xacc)
        return (p1);
    p = rtsec1(pT2v_reg2, T, v, p1, p2, f1, f2, xacc, iMAX);
    if (p < PMIN2)
        p = PMIN2;
    if (p > pmax2)
        p = pmax2;
    return (p);
}

//----------------------------------------------
//      (T,s)->p
//----------------------------------------------
double Ts2p_reg2(double T, double s)
{
    double p, p1, p2, pmax2, f1, f2, s1, s2;
    p1 = PMIN2;
    s1 = pT2s_reg2(p1, T);
    f1 = s - s1;
    pmax2 = T2pmax_reg2(T);
    p2 = pmax2;
    s2 = pT2s_reg2(p2, T);
    f2 = s - s2;
    p1 = p2 - (p2 - p1) * (s - s2) / (s1 - s2);
    if (p1 < PMIN2)
        p1 = PMIN2;
    s1 = pT2s_reg2(p1, T);
    f1 = s - s1;
    if (fabs(f1) < xacc)
        return (p1);
    p = rtsec1(pT2s_reg2, T,s, p1, p2, f1, f2, xacc, iMAX);
    if (p < PMIN2)
        p = PMIN2;
    if (p > pmax2)
        p = pmax2;
    return (p);
}

// Region 2(T,h)->p using the secant method
//  *  T: temperature  K
//  *  h: specific enthalpy kJ/kg
//  *  p: pressure  MPa
double Th2p_reg2(double T, double h)
{
    double p, p1, p2, pmax2, pmid, stepa, stepm, f1, f2, h1, h2;
    int bounded = 0;

    p1 = PMIN2; // 0.0001;
    pmax2 = T2pmax_reg2(T);
    p2 = pmax2;
    stepa = 1;
    stepm = 5.0;
    h1 = pT2h_reg2(p1, T);
    if (fabs(h - h1) < xacc)
        return (p1);
    p2 = p1 * stepm;
    h2 = pT2h_reg2(p2, T);
    if (fabs(h - h2) < xacc)
        return (p2);
    if ((h > h2) && (h < h1))
        bounded = 1;
    while (bounded == 0)
    {
        p1 = p2;
        h1 = h2;
        if (p1 > 1)
            p2 = p1 + stepa;
        if (p1 < 1)
            p2 = p1 * stepm;
        if (p2 >= pmax2)
        {
            p2 = pmax2;
            bounded = 1;
        }
        h2 = pT2h_reg2(p2, T);
        if (fabs(h - h2) < xacc)
            return (p2);
        if ((h > h2) && (h < h1))
            bounded = 1;
    }
    f2 = h - h2;
    pmid = p2 - (p2 - p1) * (h - h2) / (h1 - h2);
    if (pmid < PMIN2)
    {
        pmid = PMIN2;
    };
    h1 = pT2h_reg2(pmid, T);
    if (h < pT2h_reg2(pmid, T))
        p1 = pmid;
    h1 = pT2h_reg2(p1, T);
    f1 = h - h1;
    if (fabs(f1) < xacc)
        return (p1);
    p = rtsec1(pT2h_reg2, T, h, p1, p2, f1, f2, xacc, iMAX);
    if (p < PMIN2)
        p = PMIN2;
    if (p > pmax2)
        p = pmax2;
    return (p);
}

// specific volume in region 2
// vreg2 in m^3/kg
// temperature in K
// pressure in MPa
double pT2v_reg2(double p, double T)
{
    double tau = r2Tstar / T;
    double pi = p;
    return 0.001 * rgas_water * T * pi * (gamma0_pi_reg2(pi) + gammar_pi_reg2(pi, tau)) / p;
}

// specific enthalpy in region 2
// hreg2 in kJ/kg
// temperature in K
// pressure in Map
double pT2h_reg2(double p, double T)
{
    double tau = r2Tstar / T;
    double pi = p;
    return rgas_water * T * tau * (gamma0_tau_reg2(tau) + gammar_tau_reg2(pi, tau));
}

// specific internal energy in region 2
// ureg2 in kJ/kg
// temperature in K
// pressure in Mpa
double pT2u_reg2(double p, double T)
{
    double tau = r2Tstar / T;
    double pi = p;

    double gammar = 0.0;
    double gammarpi = 0.0;
    double gammarpipi = 0.0;
    double gammartau = 0.0;
    double gammartautau = 0.0;
    double gammarpitau = 0.0;
  
    polys_0_i_ii_j_jj_ij_reg2(pi, tau, &gammar, &gammarpi, &gammarpipi, &gammartau, &gammartautau, &gammarpitau);

    return rgas_water * T * (tau * (gamma0_tau_reg2(tau) + gammartau) - pi * (gamma0_pi_reg2(pi) + gammarpi));
}

// specific isobaric heat capacity in region 2
// cpreg2 in kJ/(kg K)
// temperature in K
// pressure in Mpa
double pT2cp_reg2(double p, double T)
{
    double tau = r2Tstar / T;
    double pi = p;
    return -rgas_water * tau * tau * (gamma0_tautau_reg2(pi, tau) + gammar_tautau_reg2(pi, tau));
}

// specific isochoric heat capacity in region 2
// cvreg2 in kJ/(kg K)
// temperature in K
// pressure in bar
double pT2cv_reg2(double p, double T)
{
    double tau = r2Tstar / T;
    double pi = p;
    
    double gammar = 0.0;
    double gammarpi = 0.0;
    double gammarpipi = 0.0;
    double gammartau = 0.0;
    double gammartautau = 0.0;
    double gammarpitau = 0.0;
  
    polys_0_i_ii_j_jj_ij_reg2(pi, tau, &gammar, &gammarpi, &gammarpipi, &gammartau, &gammartautau, &gammarpitau);

    double a = 1 + pi * (gammarpi - tau * gammarpitau);
    return rgas_water * (-tau * tau * (gamma0_tautau_reg2(pi, tau) + gammartautau) - a * a / (1 - pi * pi * gammarpipi));
}

// speed of sound in region 2
// wreg2 in m/s
// temperature in K
// pressure in Mpa
double pT2w_reg2(double p, double T)
{
    double tau = r2Tstar / T;
    double pi = p;

    double gammar = 0.0;
    double gammarpi = 0.0;
    double gammarpipi = 0.0;
    double gammartau = 0.0;
    double gammartautau = 0.0;
    double gammarpitau = 0.0;
  
    polys_0_i_ii_j_jj_ij_reg2(pi, tau, &gammar, &gammarpi, &gammarpipi, &gammartau, &gammartautau, &gammarpitau);

    double a = pi * gammarpi;
    double b = 1 + pi * (gammarpi - tau * gammarpitau);
    return sqrt(1000.0 * rgas_water * T * (1 + 2 * pi * gammarpi + a * a) /
                ((1 - pi * pi * gammarpipi) + b * b / (tau * tau * (gamma0_tautau_reg2(pi, tau) + gammartautau))));
}


// for iter (h,s)->p
double ph2s_reg2(double p, double h)
{
  double T;
  T = ph2T_reg2(p, h);
  return pT2s_reg2(p, T);
}

double pi2aHS(double eta, double sigma)
{
  // Initialize coefficients and exponents (H,S)->P for region 2a
  IJnData IJn[] = {{0, 1, -0.182575361923032E-01},
                   {0, 3, -0.125229548799536},
                   {0, 6, 0.592290437320145},
                   {0, 16, 0.604769706185122E+01},
                   {0, 20, 0.238624965444474E+03},

                   {0, 22, -0.298639090222922E+03},
                   {1, 0, 0.512250813040750E-01},
                   {1, 1, -0.437266515606486},
                   {1, 2, 0.413336902999504},
                   {1, 3, -0.516468254574773E+01},

                   {1, 5, -0.557014838445711E+01},
                   {1, 6, 0.128555037824478E+02},
                   {1, 10, 0.114144108953290E+02},
                   {1, 16, -0.119504225652714E+03},
                   {1, 20, -0.284777985961560E+04},

                   {1, 22, 0.431757846408006E+04},
                   {2, 3, 0.112894040802650E+01},
                   {2, 16, 0.197409186206319E+04},
                   {2, 20, 0.151612444706087E+04},
                   {3, 0, 0.141324451421235E-01},

                   {3, 2, 0.585501282219601},
                   {3, 3, -0.297258075863012E+01},
                   {3, 6, 0.594567314847319E+01},
                   {3, 16, -0.623656565798905E+04},
                   {4, 16, 0.965986235133332E+04},

                   {5, 3, 0.681500934948134E+01},
                   {5, 16, -0.633207286824489E+04},
                   {6, 3, -0.558919224465760E+01},
                   {7, 1, 0.400645798472063E-01}};
  double pi, pi2;

  pi = poly(eta - 0.5, sigma - 1.2, 29, IJn);
  pi2 = pi * pi;
  return (pi2 * pi2);
}

double hs2p_reg2a(double h, double s)
{
  double eta, sigma;
  eta = h / 4200.0;
  sigma = s / 12.0;
  return (4.0 * pi2aHS(eta, sigma));
}
//=========================================================
// Initialize coefficients and exponents (H,S)->P for region 2b
//  Table 7
//=========================================================
double pi2bHS(double eta, double sigma)
{
  IJnData IJn[] = {{0, 0, 0.801496989929495E-01},
                   {0, 1, -0.543862807146111},
                   {0, 2, 0.337455597421283},
                   {0, 4, 0.890555451157450E+01},
                   {0, 8, 0.313840736431485E+03},
                   {1, 0, 0.797367065977789},
                   {1, 1, -0.121616973556240E+01},
                   {1, 2, 0.872803386937477E+01},
                   {1, 3, -0.169769781757602E+02},
                   {1, 5, -0.186552827328416E+03},
                   {1, 12, 0.951159274344237E+05},
                   {2, 1, -0.189168510120494E+02},
                   {2, 6, -0.433407037194840E+04},
                   {2, 18, 0.543212633012715E+09},
                   {3, 0, 0.144793408386013},
                   {3, 1, 0.128024559637516E+03},
                   {3, 7, -0.672309534071268E+05},
                   {3, 12, 0.336972380095287E+08},
                   {4, 1, -0.586634196762720E+03},
                   {4, 16, -0.221403224769889E+11},
                   {5, 1, 0.171606668708389E+04},
                   {5, 12, -0.570817595806302E+09},
                   {6, 1, -0.312109693178482E+04},
                   {6, 8, -0.207841384633010E+07},
                   {6, 18, 0.305605946157786E+13},
                   {7, 1, 0.322157004314333E+04},
                   {7, 16, 0.326810259797295E+12},
                   {8, 1, -0.144104158934487E+04},
                   {8, 3, 0.410694867802691E+03},
                   {8, 14, 0.109077066873024E+12},
                   {8, 18, -0.247964654258893E+14},
                   {12, 10, 0.188801906865134E+10},
                   {14, 16, -0.123651009018773E+15}};
  double pi, pi2;

  pi = poly(eta - 0.6, sigma - 1.01, 33, IJn);
  pi2 = pi * pi;
  return (pi2 * pi2);
}

double hs2p_reg2b(double h, double s)
{
  double eta, sigma;
  eta = h / 4100.0;
  sigma = s / 7.9;
  return (100.0 * pi2bHS(eta, sigma));
}

//=========================================================
// Initialize coefficients and exponents (H,S)->P for region 2c
//  Table 8
//=========================================================
double pi2cHS(double eta, double sigma)
{
  IJnData IJn[] = {{0, 0, 0.112225607199012E+00},
                   {0, 1, -0.339005953606712E+01},
                   {0, 2, -0.320503911730094E+02},
                   {0, 3, -0.197597305104900E+03},
                   {0, 4, -0.407693861553446E+03},
                   {0, 8, 0.132943775222331E+05},
                   {1, 0, 0.170846839774007E+01},
                   {1, 2, 0.373694198142245E+02},
                   {1, 5, 0.358144365815434E+04},
                   {1, 8, 0.423014446424664E+06},
                   {1, 14, -0.751071025760063E+09},
                   {2, 2, 0.523446127607898E+02},
                   {2, 3, -0.228351290812417E+03},
                   {2, 7, -0.960652417056937E+06},
                   {2, 10, -0.807059292526074E+08},
                   {2, 18, 0.162698017225669E+13},
                   {3, 0, 0.772465073604171},
                   {3, 5, 0.463929973837746E+05},
                   {3, 8, -0.137317885134128E+08},
                   {3, 16, 0.170470392630512E+13},
                   {3, 18, -0.251104628187308E+14},
                   {4, 18, 0.317748830835520E+14},
                   {5, 1, 0.538685623675312E+02},
                   {5, 4, -0.553089094625169E+05},
                   {5, 6, -0.102861522421405E+07},
                   {5, 14, 0.204249418756234E+13},
                   {6, 8, 0.273918446626977E+09},
                   {6, 18, -0.263963146312685E+16},
                   {10, 7, -0.107890854108088E+10},
                   {12, 7, -0.296492620980124E+11},
                   {16, 10, -0.111754907323424E+16}};
  double pi, pi2;

  pi = poly(eta - 0.7, sigma - 1.1, 31, IJn);
  pi2 = pi * pi;
  return (pi2 * pi2);
}

double hs2p_reg2c(double h, double s)
{
  double eta, sigma;
  eta = h / 3500.0;
  sigma = s / 5.9;
  return (100.0 * pi2cHS(eta, sigma));
}
//------------------------------------------------------------
// Region 2  HS
//------------------------------------------------------------

double s2hreg2ab(double s)
/*  Define the boundary between Region 2a and 2b, h=f(s)
   Water and Steam, http://www.iapws.org/relguide/Supp-PHS12-2014.pdf, Eq 2
     >>> _hab_s(7)   3376.437884
   """
*/
{
  static double n[4] = {-0.349898083432139E+04,
                        0.257560716905876E+04,
                        -0.421073558227969E+03,
                        0.276349063799944E+02};
  double eta, sigma;
  sigma = s / 1.0;
  eta = n[0] + sigma * (n[1] + sigma * (n[2] + n[3] * sigma));
  return (eta * 1.0);
}

double hs2p_reg2(double h, double s)
{
  double h2ab, p;
  h2ab = s2hreg2ab(s);
  if (h > h2ab)
  {
    if (s >= 5.85)
    {
      p = hs2p_reg2b(h, s);
    }
    else
    {
      p = hs2p_reg2c(h, s);
    };
  }
  else
  {
    p = hs2p_reg2a(h, s);
  };
  return (p);

  /*  double p1,p2,f1,f2;
    p1 =p;
    f1=s-ph2s_reg2(p1,h);
    if (fabs(f1)>xacc)
    {
      if (f1>0) // pT2sreg1(p,h)< s ,the p1< expt p，so， p2=1.05*p1 p（p1,p2)
        p2=(1.0+f1/s)*p1;
      else
        p2=(1.0-f1/s)*p1;

      f2=s-ph2s_reg2(p2,h);

      p=rtsec1(ph2s_reg2,h,s,p1,p2,f1,f2,xacc,iMAX);
    }
    else
      p=p1;

   return p;   */
}

double ph2T_reg2a(double p, double h)
{
  // Table 20.
  // Numerical values of the coefficients and exponents of
  // the backward equation T ( p,h ) for
  // subregion 2a, Eq. (22)
  IJnData IJn[34] = {
      {0, 0, 0.10898952318288E+04},
      {0, 1, 0.84951654495535E+03},
      {0, 2, -0.10781748091826E+03},
      {0, 3, 0.33153654801263E+02},
      {0, 7, -0.74232016790248E+01},

      {0, 20, 0.11765048724356E+02},
      {1, 0, 0.18445749355790E+01},
      {1, 1, -0.41792700549624E+01},
      {1, 2, 0.62478196935812E+01},
      {1, 3, -0.17344563108114E+02},

      {1, 7, -0.20058176862096E+03},
      {1, 9, 0.27196065473796E+03},
      {1, 11, -0.45511318285818E+03},
      {1, 18, 0.30919688604755E+04},
      {1, 44, 0.25226640357872E+06},

      {2, 0, -0.61707422868339E-02},
      {2, 2, -0.31078046629583E+00},
      {2, 7, 0.11670873077107E+02},
      {2, 36, 0.12812798404046E+09},
      {2, 38, -0.98554909623276E+09}

      ,
      {2, 40, 0.28224546973002E+10},
      {2, 42, -0.35948971410703E+10},
      {2, 44, 0.17227349913197E+10},
      {3, 24, -0.13551334240775E+05},
      {3, 44, 0.12848734664650E+08},

      {4, 12, 0.13865724283226E+01},
      {4, 32, 0.23598832556514E+06},
      {4, 44, -0.13105236545054E+08},
      {5, 32, 0.73999835474766E+04},
      {5, 36, -0.55196697030060E+06},

      {5, 42, 0.37154085996233E+07},
      {6, 34, 0.19127729239660E+05},
      {6, 44, -0.41535164835634E+06},
      {7, 28, -0.62459855192507E+02}};

  double pi, eta, theta;
  pi = p / 1.0;
  eta = h / 2000.0 - 2.1;

  return poly(pi, eta, 34, IJn);
}

double ph2T_reg2b(double p, double h)
{
  // Table 21. Numerical values of the coefficients and exponents
  // of the backward equation T ( p,h ) for
  IJnData IJn[38] = {
      {0, 0, 1.4895041079516e3},
      {0, 1, 7.4307798314034e2},
      {0, 2, -9.7708318797837e1},
      {0, 12, 2.4742464705674},
      {0, 18, -6.3281320016026e-1},

      {0, 24, 1.1385952129658},
      {0, 28, -4.7811863648625e-1},
      {0, 40, 8.5208123431544e-3},
      {1, 0, 9.3747147377932e-1},
      {1, 2, 3.3593118604916},

      {1, 6, 3.3809355601454},
      {1, 12, 1.6844539671904e-1},
      {1, 18, 7.3875745236695e-1},
      {1, 24, -4.7128737436186e-1},
      {1, 28, 1.5020273139707e-1},

      {1, 40, -2.1764114219750e-3},
      {2, 2, -2.1810755324761e-2},
      {2, 8, -1.0829784403677e-1},
      {2, 18, -4.6333324635812e-2},
      {2, 40, 7.1280351959551e-5},

      {3, 1, 1.1032831789999e-4},
      {3, 2, 1.8955248387902e-4},
      {3, 12, 3.0891541160537e-3},
      {3, 24, 1.3555504554949e-3},
      {4, 2, 2.8640237477456e-7},

      {4, 12, -1.0779857357512e-5},
      {4, 18, -7.6462712454814e-5},
      {4, 24, 1.4052392818316e-5},
      {4, 28, -3.1083814331434e-5},
      {4, 40, -1.0302738212103e-6},

      {5, 18, 2.8217281635040e-7},
      {5, 24, 1.2704902271945e-6},
      {5, 40, 7.3803353468292e-8},
      {6, 28, -1.1030139238909e-8},
      {7, 2, -8.1456365207833e-14},

      {7, 28, -2.5180545682962e-11},
      {9, 1, -1.7565233969407e-18},
      {9, 40, 8.6934156344163e-15}};
  double pi, eta, theta;
  pi = p / 1.0 - 2;
  eta = h / 2000.0 - 2.6;

  return poly(pi, eta, 38, IJn);
};

// Table 22. Numerical values of the coefficients and exponents of
// the backward   #equation T ( p,h ) for subregion 2c, Eq. (24)
double ph2T_reg2c(double p, double h)
{
  IJnData IJn[23] = {
      {-7, 0, -3236839855524.2},
      {-7, 4, 7326335090218.1},
      {-6, 0, 358250899454.47},
      {-6, 2, -583401318515.9},
      {-5, 0, -10783068217.47},

      {-5, 2, 20825544563.171},
      {-2, 0, 610747.83564516},
      {-2, 1, 859777.2253558},
      {-1, 0, -25745.72360417},
      {-1, 2, 31081.088422714},

      {0, 0, 1208.2315865936},
      {0, 1, 482.19755109255},
      {1, 4, 3.7966001272486},
      {1, 8, -10.842984880077},
      {2, 4, -0.04536417267666},

      {6, 0, 1.4559115658698E-13},
      {6, 1, 1.126159740723E-12},
      {6, 4, -1.7804982240686E-11},
      {6, 10, 1.2324579690832E-07},
      {6, 12, -1.1606921130984E-06},

      {6, 16, 2.7846367088554E-05},
      {6, 20, -5.9270038474176E-04},
      {6, 22, 1.2918582991878E-03}};

  double pi, eta, theta;
  pi = p / 1.0 + 25;
  eta = h / 2000.0 - 1.8;

  return poly(pi, eta, 23, IJn);
}

// http://www.iapws.org/relguide/IF97-Rev.html, Eq 21
double enthalpy2bc(double p)
{
  const double n[] = {0.12809002730136e-3,
                      0.26526571908428e4,
                      0.45257578905948e1};
  return n[1] + sqrt((p - n[2]) / n[0]);
}

double ph2T_reg2(double p, double h)
{
  double T;
  if (p > 4)
    if (h < enthalpy2bc(p))
      T = ph2T_reg2c(p, h);
    else
      T = ph2T_reg2b(p, h);
  else
    T = ph2T_reg2a(p, h);

  /*double T1, T2, f1, f2;
  T1 = T;
  f1 = h - pT2h_reg2(p, T1);
  if (fabs(f1) > xacc)
  {
    if (f1 > 0)
      T2 = (1.0 + f1 / h) * T1; // TODO： 1.05 用 1+f1/h 是不是更快，没有测试
    else
      T2 = (1.0 - f1 / h) * T1;

    f2 = h - pT2h_reg2(p, T2);

    T = rtsec2(pT2h_reg2, p, h, T1, T2, f1, f2, xacc, iMAX);
  }
  else
    T = T1;*/

  return T;
}

typedef struct
{
    double I;
    int J;
    double n;
} dIJnData;

//-----  the backward equation T( p,s ) ---------------------

// Table 25. Numerical values of the coefficients and exponents of
// the backward equation T( p,s ) for
// subregion 2a, Eq. (25)
double ps2T_reg2a(double p, double s)
{
    dIJnData IJn[46] = {
        {-1.5, -24, -0.39235983861984E+06},
        {-1.5, -23, 0.51526573827270E+06},
        {-1.5, -19, 0.40482443161048E+05},
        {-1.5, -13, -0.32193790923902E+03},
        {-1.5, -11, 0.96961424218694E+02},

        {-1.5, -10, -0.22867846371773E+02},
        {-1.25, -19, -0.44942914124357E+06},
        {-1.25, -15, -0.50118336020166E+04},
        {-1.25, -6, 0.35684463560015E+00},
        {-1.0, -26, 0.44235335848190E+05},

        {-1.0, -21, -0.13673388811708E+05},
        {-1.0, -17, 0.42163260207864E+06},
        {-1.0, -16, 0.22516925837475E+05},
        {-1.0, -9, 0.47442144865646E+03},
        {-1.0, -8, -0.14931130797647E+03},

        {-0.75, -15, -0.19781126320452E+06},
        {-0.75, -14, -0.23554399470760E+05},
        {-0.5, -26, -0.19070616302076E+05},
        {-0.5, -13, 0.55375669883164E+05},
        {-0.5, -9, 0.38293691437363E+04},

        {-0.5, -7, -0.60391860580567E+03},
        {-0.25, -27, 0.19363102620331E+04},
        {-0.25, -25, 0.42660643698610E+04},
        {-0.25, -11, -0.59780638872718E+04},
        {-0.25, -6, -0.70401463926862E+03},

        {0.25, 1, 0.33836784107553E+03},
        {0.25, 4, 0.20862786635187E+02},
        {0.25, 8, 0.33834172656196E-01},
        {0.25, 11, -0.43124428414893E-04},
        {0.5, 0, 0.16653791356412E+03},

        {0.5, 1, -0.13986292055898E+03},
        {0.5, 5, -0.78849547999872E+00},
        {0.5, 6, 0.72132411753872E-01},
        {0.5, 10, -0.59754839398283E-02},
        {0.5, 14, -0.12141358953904E-04},

        {0.5, 16, 0.23227096733871E-06},
        {0.75, 0, -0.10538463566194E+02},
        {0.75, 4, 0.20718925496502E+01},
        {0.75, 9, -0.72193155260427E-01},
        {0.75, 17, 0.20749887081120E-06},

        {1.0, 7, -0.18340657911379E-01},
        {1.0, 18, 0.29036272348696E-06},
        {1.25, 3, 0.21037527893619E+00},
        {1.25, 15, 0.2568123972999E-03},
        {1.5, 5, -0.1279900293381E-01},

        {1.5, 18, -0.82198102652018E-05}};

    double pi, sigma, theta;
    pi = p / 1.0;
    sigma = s / 2.0 - 2.0;

    theta = 0.0;
    for (int k = 0; k < 46; k++)
    {
        theta += IJn[k].n * pow(pi, IJn[k].I) * pow(sigma, IJn[k].J);
    }
    return 1.0 * theta;
}

// Table 26. Numerical values of the coefficients and exponents of
// the backward equation T( p,s ) for
// subregion 2b, Eq. (26)
double ps2T_reg2b(double p, double s)
{
    IJnData IJn[44] = {
        {-6, 0, 0.31687665083497e6},
        {-6, 11, 0.20864175881858e2},
        {-5, 0, -0.39859399803599e6},
        {-5, 11, -0.21816058518877e2},
        {-4, 0, 0.22369785194242e6},

        {-4, 1, -0.27841703445817e4},
        {-4, 11, 0.99207436071480e1},
        {-3, 0, -0.75197512299157e5},
        {-3, 1, 0.29708605951158e4},
        {-3, 11, -0.34406878548526e1},

        {-3, 12, 0.38815564249115},
        {-2, 0, 0.17511295085750e5},
        {-2, 1, -0.14237112854449e4},
        {-2, 6, 0.10943803364167e1},
        {-2, 10, 0.89971619308495},

        {-1, 0, -0.33759740098958e4},
        {-1, 1, 0.47162885818355e3},
        {-1, 5, -0.19188241993679e1},
        {-1, 8, 0.41078580492196},
        {-1, 9, -0.33465378172097},

        {0, 0, 0.13870034777505e4},
        {0, 1, -0.40663326195838e3},
        {0, 2, 0.41727347159610e2},
        {0, 4, 0.21932549434532e1},
        {0, 5, -0.10320050009077e1},

        {0, 6, 0.35882943516703},
        {0, 9, 0.52511453726066e-2},
        {1, 0, 0.12838916450705e2},
        {1, 1, -0.28642437219381e1},
        {1, 2, 0.56912683664855},

        {1, 3, -0.99962954584931e-1},
        {1, 7, -0.32632037778459e-2},
        {1, 8, 0.23320922576723e-3},
        {2, 0, -0.15334809857450},
        {2, 1, 0.29072288239902e-1},

        {2, 5, 0.37534702741167e-3},
        {3, 0, 0.17296691702411e-2},
        {3, 1, -0.38556050844504e-3},
        {3, 3, -0.35017712292608e-4},
        {4, 0, -0.14566393631492e-4},

        {4, 1, 0.56420857267269e-5},
        {5, 0, 0.41286150074605e-7},
        {5, 1, -0.20684671118824e-7},
        {5, 2, 0.16409393674725e-8}};

    double pi, sigma;
    pi = p / 1.0;
    sigma = 10 - s / 0.7853;
    
    double theta = poly(pi, sigma, 44, IJn);
    return 1.0 * theta;
}

// Table 27. Numerical values of the coefficient s and exponents of
//  the backward equation T( p,s ) for
// subregion 2c, Eq. (27)
double ps2T_reg2c(double p, double s)
{
    IJnData IJn[30] = {
        {-2, 0, 0.90968501005365e3},
        {-2, 1, 0.24045667088420e4},
        {-1, 0, -0.59162326387130e3},
        {0, 0, 0.54145404128074e3},
        {0, 1, -0.27098308411192e3},

        {0, 2, 0.97976525097926e3},
        {0, 3, -0.46966772959435e3},
        {1, 0, 0.14399274604723e2},
        {1, 1, -0.19104204230429e2},
        {1, 3, 0.53299167111971e1},

        {1, 4, -0.21252975375934e2},
        {2, 0, -0.31147334413760},
        {2, 1, 0.60334840894623},
        {2, 2, -0.42764839702509e-1},
        {3, 0, 0.58185597255259e-2},

        {3, 1, -0.14597008284753e-1},
        {3, 5, 0.56631175631027e-2},
        {4, 0, -0.76155864584577e-4},
        {4, 1, 0.22440342919332e-3},
        {4, 4, -0.12561095013413e-4},

        {5, 0, 0.63323132660934e-6},
        {5, 1, -0.20541989675375e-5},
        {5, 2, 0.36405370390082e-7},
        {6, 0, -0.29759897789215e-8},
        {6, 1, 0.10136618529763e-7},

        {7, 0, 0.59925719692351e-11},
        {7, 1, -0.20677870105164e-10},
        {7, 3, -0.20874278181886e-10},
        {7, 4, 0.10162166825089e-9},
        {7, 5, -0.16429828281347e-9}};

    double pi, sigma;
    pi = p / 1.0;
    sigma = 2 - s / 2.9251;
    
    double theta = poly(pi, sigma, 30, IJn);
    return 1.0 * theta;
}

double ps2T_reg2(double p, double s)
{
    double T;
    if (p > 4)
    {
        if (s < 5.85)
        {
            T = ps2T_reg2c(p, s);
        }
        else
        {
            T = ps2T_reg2b(p, s);
        }
    }
    else
    {
        T = ps2T_reg2a(p, s);
    }
    return T;
    /*
        double T1, T2, f1, f2;
        T1 = T;
        f1 = s - pT2s_reg2(p, T1);
        if (fabs(f1) > xacc)
        {
            if (f1 > 0) // pT2sreg1(p,T1)< s ,the T1< expt T，so， T2=1.05*T1 T（T1,T2)
                T2 = (1.0 + f1 / s) * T1;
            else
                T2 = (1.0 - f1 / s) * T1;

            f2 = s - pT2s_reg2(p, T2);

            T = rtsec2(pT2s_reg2, p, s, T1, T2, f1, f2, xacc, iMAX);
        }
        else
            T = T1;

        return T;*/
}

// specific entropy in region 2
// sreg2 in kJ/(kg K)
// temperature in K
// pressure in Mpa
double pT2s_reg2(double p, double T)
{
    double tau = r2Tstar / T;
    double pi = p;

    // gamma0_
    double gamma0 = gamma0_reg2(pi, tau);
    double gamma0tau = gamma0_tau_reg2(tau);

    // gammar
    double gammar = 0.0;
    double gammarpi = 0.0;
    double gammarpipi = 0.0;
    double gammartau = 0.0;
    double gammartautau = 0.0;
    double gammarpitau = 0.0;
  
    polys_0_i_ii_j_jj_ij_reg2(pi, tau, &gammar, &gammarpi, &gammarpipi, &gammartau, &gammartautau, &gammarpitau);

    return rgas_water * (tau * (gamma0tau + gammartau) - (gamma0 + gammar));
}
