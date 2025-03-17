/*
 Region 4 - The extended input pairs: (p,v)->x, (t,v),(t,h),(t,s) ->x
    x: Steam quality
*/

/*
 Region4 : http://www.iapws.org/relguide/Supp-phs3-2014.pdf. Eq 9
       Page30  Supp-phs3-2014.pdf Page25
           5.3 Backward Equation Tsat(h,s)
               s> 5.210 887 825
         temperature range is T （273.15 ， 623.15 ）
*/
#include <math.h>
#include <stdio.h>
#include "../r1/region1.h"
#include "../r2/region2.h"
#include "../r3/region3.h"
#include "region4.h"

/// Region 4: (p,v)-> x
double pv2x_reg4(double p, double v)
{
    double v_sat_w = p2SatWater(p, OV);
    double v_sat_s = p2SatSteam(p, OV);
    return (v - v_sat_w) / (v_sat_s - v_sat_w);
}

/// Region 4: (T,v)-> x
double Tv2x_reg4(double T, double v)
{
    double v_sat_w = T2SatWater(T, OV);
    double v_sat_s = T2SatSteam(T, OV);
    return (v - v_sat_w) / (v_sat_s - v_sat_w);
}

/// Region 4: (T,h)-> x
double Th2x_reg4(double T, double h)
{
    double h_sat_w = T2SatWater(T, OH);
    double h_sat_s = T2SatSteam(T, OH);
    return (h - h_sat_w) / (h_sat_s - h_sat_w);
}

/// Region 4: (T,s)-> x
double Ts2x_reg4(double T, double s)
{
    double s_sat_w = T2SatWater(T, OS);
    double s_sat_s = T2SatSteam(T, OS);
    return (s - s_sat_w) / (s_sat_s - s_sat_w);
}

//
//  Initialize coefficients for region 4
//
static double n[11] = {0, 0.11670521452767E+04, -0.72421316703206E+06, -0.17073846940092E+02,
    0.12020824702470E+05, -0.32325550322333E+07, 0.14915108613530E+02,
    -0.48232657361591E+04, 0.40511340542057E+06, -0.23855557567849E+00,
    0.65017534844798E+03};

double pSat(double T)
// saturation pressure of water
// pSatW in bar
// T :temperaturein K
//
// pSat = -1: temperature outside range
//
{
double pS;
if (T < TMIN4 || T > tc_water) // tc_water=647.096
pS = -1.0;
else
{
double del = T + n[9] / (T - n[10]);
double aco = del * (del + n[1]) + n[2];
double bco = del * (n[3] * del + n[4]) + n[5];
double cco = del * (n[6] * del + n[7]) + n[8];
pS = ipowsac(2 * cco / (-bco + sqrt(bco * bco - 4 * aco * cco)), 4);
}
return pS;
}

double TSat(double p)
//
// saturation temperature of water
// tSatW in K
// p :pressure in bar
//
// tSatW=-1: pressure outside range
//
{
double TS;
if (p < Pmin || p > pc_water)
TS = -1.0;
else
{
double bet = pow(p, 0.25);
double eco = bet * (bet + n[3]) + n[6];
double fco = bet * (n[1] * bet + n[4]) + n[7];
double gco = bet * (n[2] * bet + n[5]) + n[8];
double dco = 2.0 * gco / (-fco - sqrt(fco * fco - 4.0 * eco * gco));
TS = 0.5 * (n[10] + dco - sqrt((n[10] + dco) * (n[10] + dco) - 4.0 * (n[9] + n[10] * dco)));
}
return TS;
}


static double s4L_623 = 3.778281340;       // T=623.16 ，Sature Liquid s=3.778281340
static double s4V_623 = 5.210887825;       // T=623.16 ，Sature steam s=5.210887825
static double h4V_623 = 2.5635920043 + 03; // Page 25  Sature steam h=2.5635920043 + 03
static double s4L_273 = -1.545495919e-04;  // T=273.15 ，Sature Liquid s=-1.545495919e-04
static double s4V_273 = 9.155759395;       // T=273.15 ，Sature steam s=9.155759395

double hs2T_reg43(double h, double s)
/* Backward equation for region 4, T=f(h,s)
    h :        Specific enthalpy [kJ/kg]
    s :        Specific entropy [kJ/kgK]
 */
{
    IJnData IJn[] = {
        {0, 0, .179882673606601},
        {0, 3, -.267507455199603},
        {0, 12, .116276722612600e1},
        {1, 0, .147545428713616},
        {1, 1, -.512871635973248},
        {1, 2, .421333567697984},
        {1, 5, .563749522189870},
        {2, 0, .429274443819153},
        {2, 5, -.335704552142140e1},
        {2, 8, .108890916499278e2},
        {3, 0, -.248483390456012},
        {3, 2, .304153221906390},
        {3, 3, -.494819763939905},
        {3, 4, .107551674933261e1},
        {4, 0, .733888415457688e-1},
        {4, 1, .140170545411085e-1},
        {5, 1, -.106110975998808},
        {5, 2, .168324361811875e-1},
        {5, 4, .125028363714877e1},
        {5, 16, .101316840309509e4},
        {6, 6, -.151791558000712e1},
        {6, 8, .524277865990866e2},
        {6, 22, .230495545563912e5},
        {8, 1, .249459806365456e-1},
        {10, 20, .210796467412137e7},
        {10, 36, .366836848613065e9},
        {12, 24, -.144814105365163e9},
        {14, 1, -.179276373003590e-2},
        {14, 28, .489955602100459e10},
        {16, 12, .471262212070518e3},
        {16, 32, -.829294390198652e11},
        {18, 14, -.171545662263191e4},
        {18, 22, .355777682973575e7},
        {18, 36, .586062760258436e12},
        {20, 24, -.129887635078195e8},
        {28, 36, .317247449371057e11}};

    if (s < s4V_623)
        return INVALID_S;

    double nu = h / 2800;
    double sigma = s / 9.2;
    // double suma = 0;
    // for (int i = 0; i < 36; i++)
    //    suma += IJn[i].n * ipowsac(nu - 0.119, IJn[i].I) * ipowsac(sigma - 1.07, IJn[i].J);
    return 550.0 * poly(nu - 0.119, sigma - 1.07, 36, IJn);
}

double hs2T_reg4(double h, double s)
{

    double T;
    if (s > s4V_623 && s < s4V_273)
    {
        T = hs2T_reg43(h, s);
        return T;
    };

    // The if97 function hs2Treg43 is only valid for part of region4. Use iteration outsida.
    double Low_Bound;
    double High_Bound;
    double PL, Ts;

    if (s > s4L_273 && s <= s4L_623)
    {
        Low_Bound = Pmin;
        High_Bound = Ps_623;

        double hL = INVALID_OUTID;
        while (fabs(hL - h) > 1.0e-04 && fabs(High_Bound - Low_Bound) > 1.0e-4)
        {
            PL = (Low_Bound + High_Bound) / 2;
            Ts = TSat(PL);
            hL = pT2h_reg1(PL, Ts);
            if (hL > h)
                High_Bound = PL;
            else
                Low_Bound = PL;
        }
    };

    if (s > s4L_623 && s <= sc_water)
    {
        PL = h2pSat_reg3(h); // liquid
        Low_Bound = Pmin;
        High_Bound = PL;
    }
    if (s > sc_water && s <= s4V_623)
    {
        PL = h2pSat_reg3(h); // steam
        Low_Bound = Pmin;
        High_Bound = PL;
    }

    double sss = INVALID_OUTID;
    double p, xs, s4v, s4L, v4v, v4L;

    while (fabs(s - sss) > 1.0e-6)
    {
        p = 0.5 * (Low_Bound + High_Bound);

        Ts = TSat(p);
        xs = ph_reg4(p, h, OX);

        if (p < Ps_623)
        {
            s4v = pT2s_reg2(p, Ts);
            s4L = pT2s_reg1(p, Ts);
        }
        else
        {
            v4v = ph_reg3(p, p2SatSteam(p, OH), OV);
            s4v = Td2s_reg3(Ts, 1 / v4v);
            v4L = ph_reg3(p, p2SatWater(p, OH), OV);
            s4L = Td2s_reg3(Ts, 1 / v4L);
        };

        sss = (xs * s4v + (1 - xs) * s4L);

        if (sss < s)
        {
            High_Bound = p;
            Low_Bound = (1 + (sss - s) / s) * p;
        }
        else
        {
            Low_Bound = p;
            High_Bound = (1 + (sss - s) / s) * p;
        }
    } // end of while  (fabs(s - sss) > 1.0e-6)

    T = TSat(p);
    return T;
}