/*
   The API of region 4

*/
#include <math.h>
#include "../r1/region1.h"
#include "../r2/region2.h"
#include "../r3/region3.h"
#include "region4.h"

// TODO(cmh@seu.edu.cn): testing
double p2SatWater(double p, int o_id)
{
   // include region 3
   double T, v, r;

   if (o_id == OP)
      return p;

   T = TSat(p);
   if (o_id == OT)
      return T;

   if (p >= Pmin && p <= Ps_623)
   {
      r = pT_reg1(p, T, o_id);
   }
   else
   {
      if (p == pc_water)
         v = 1.0 / 322.0;
      else
         v = pT2vSat_reg3(p, T, 0);

      if (o_id == OV)
         return v;
      if (o_id == OD)
         return (1.0 / v);

      r = Td_reg3(T, 1 / v, o_id);
   }
   return r;
};

double p2SatSteam(double p, int o_id)
{
   double T, v, r;
   if (o_id == OP)
      return p;

   T = TSat(p);
   if (o_id == OT)
      return T;

   if (p >= Pmin && p <= Ps_623)
   {
      r = pT_reg2(p, T, o_id);
   }
   else
   { // reg3d =ss
      if (p == pc_water)
         v = 1.0 / 322.0;
      else
         v = pT2vSat_reg3(p, T, 1);

      if (o_id == OV)
         return v;

      if (o_id == OD)
         return (1.0 / v);

      r = Td_reg3(T, 1.0 / v, o_id);
   }

   return r;
}

//--------------------------------------------------
double T2SatWater(double T, int o_id)
{
   double p, v, r;
   if (o_id == OT)
      return T;

   if (T == tc_water)
      p = pc_water;
   else
      p = pSat(T);
   if (o_id == OP)
      return p;

   if (T <= 623.15)
      r = pT_reg1(p, T, o_id);
   else
   { // region3
      if (T == tc_water)
      {
         if (o_id == OS)
            return sc_water;
         if (o_id == OH)
            return hc_water;
         v = 1.0 / 322.0;
      }
      else
         v = pT2vSat_reg3(p, T, 0);

      if (o_id == OV)
         return v;
      if (o_id == OD)
         return (1.0 / v);

      r = Td_reg3(T, 1 / v, o_id);
   }
   return (r);
}

double T2SatSteam(double T, int o_id)
{
   double p, v, r;
   if (o_id == OT)
      return T;

   if (T == tc_water)
      p = pc_water;
   else
      p = pSat(T);
   if (o_id == OP)
      return p;

   if (T <= 623.15)
      r = pT_reg2(p, T, o_id);
   else
   { // region3
      if (T == tc_water)
      {
         if (o_id == OS)
            return sc_water;
         if (o_id == OH)
            return hc_water;
         v = 1.0 / 322.0;
         if (o_id == OV)
            return v;
      }
      else
         v = pT2vSat_reg3(p, T, 1);

      r = Td_reg3(T, 1 / v, o_id);
   }
   return (r);
}

double px_reg4(double p, double x, int o_id)
/*
* T: Saturated temperature [K]
 * P: Saturated pressure [MPa]
 * x: Vapor quality [-]
 * v: Specific volume [m³/kg]
 * h: Specific enthalpy [kJ/kg]
 * s: Specific entropy [kJ/kgK]
 */
{
   // TODO: testing
   double r = INVALID_OUTID;

   if (o_id == OP)
   return p;
   if (o_id == OT)
   return TSat(p);
   if (o_id == OR)
   return 4.0;

   if (x == 0)
   {
      return p2SatWater(p, o_id);
   }
   else if (x == 1)
   {
      return p2SatSteam(p, o_id);
   }

   if (x > 0 && x < 1)
   {
      int cond = o_id == OV || o_id == OD || o_id == OH || o_id == OS || o_id == OU || o_id == OF || o_id == OG;
      if(cond)
      {
         double rL = p2SatWater(p, o_id);
         double rV = p2SatSteam(p, o_id);
         r = rL + x * (rV - rL);
      }
      cond = o_id == ODVDPCH;
      if(o_id == cond)
      {
         double hL = p2SatWater(p, OH);
         double hV = p2SatSteam(p, OH);
         double rL = p2SatWater(p, o_id);
         double rV = p2SatSteam(p, o_id);
         double dx = (x*rV + (1-x)*rL)/(hL-hV);
         r = rL + dx * (hV - hL) + x * (rV - rL);
      }
      cond = o_id == ODVDHCP;
      if(o_id == cond)
      {
         double hL = p2SatWater(p, OH);
         double hV = p2SatSteam(p, OH);
         double vL = p2SatWater(p, OV);
         double vV = p2SatSteam(p, OV);
         r = (vV - vL)/(hV - hL);
      }
   }
   return r;
}

double Tx_reg4(double T, double x, int o_id)
{
   double r = INVALID_OUTID;

   if (o_id == OT)
   return T;
   if (o_id == OT)
   return pSat(T);
   if (o_id == OR)
   return 4.0;

   if (x == 0)
   {
      return T2SatWater(T, o_id);
   }
   else if (x == 1)
   {
      return T2SatSteam(T, o_id);
   }

   if (x > 0 && x < 1)
   {
      int cond = o_id == OV || o_id == OD || o_id == OH || o_id == OS || o_id == OU || o_id == OF || o_id == OG;
      if(cond)
      {
         double r1 = T2SatWater(T, o_id);
         double r2 = T2SatSteam(T, o_id);
         r = r1 + x * (r2 - r1);
      }
      cond = o_id == ODVDPCH;
      if(o_id == cond)
      {
         double hL = T2SatWater(T, OH);
         double hV = T2SatSteam(T, OH);
         double rL = T2SatWater(T, o_id);
         double rV = T2SatSteam(T, o_id);
         double dx = (x*rV + (1-x)*rL)/(hL-hV);
         r = rL + dx * (hV - hL) + x * (rV - rL);
      }
      cond = o_id == ODVDHCP;
      if(o_id == cond)
      {
         double hL = T2SatWater(T, OH);
         double hV = T2SatSteam(T, OH);
         double vL = T2SatWater(T, OV);
         double vV = T2SatSteam(T, OV);
         r = (vV - vL)/(hV - hL);
      }
   }
   return r;
}

double ph_reg4(double p, double h, int o_id)
{
   if (o_id == OP)
      return p;
   if (o_id == OH)
      return h;
   if (o_id == OR)
      return 4.0;
   double h1, h2, x;
   h1 = p2SatWater(p, OH);
   h2 = p2SatSteam(p, OH);
   x = (h - h1) / (h2 - h1);
   if (o_id == OX)
      return x;
   return px_reg4(p, x, o_id);
}

double ps_reg4(double p, double s, int o_id)
{
   if (o_id == OP)
      return p;
   if (o_id == OS)
      return s;
   if (o_id == OR)
      return 4.0;
   double s1, s2, x;
   s1 = p2SatWater(p, OS);
   s2 = p2SatSteam(p, OS);
   x = (s - s1) / (s2 - s1);
   if (o_id == OX)
      return x;
   return px_reg4(p, x, o_id);
}

double hs_reg4(double h, double s, int o_id)
{
   if (o_id == OH)
      return h;
   if (o_id == OS)
      return s;
   if (o_id == OR)
      return 4.0;
   double T, p;
   T = hs2T_reg4(h, s);
   if (o_id == OT)
      return T;
   p = pSat(T);
   if (o_id == OP)
      return p;

   double h1 = p2SatWater(p, OH);
   double h2 = p2SatSteam(p, OH);

   double x = (h - h1) / (h2 - h1);
   if (o_id == OX)
      return x;

   return px_reg4(p, x, o_id);
}

double pT_reg4(double p, double T, int o_id)
{
   if (o_id == OP)
      return p;
   if (o_id == OT)
      return T;
   if (o_id == OR)
      return 4.0;
   double value;
   if (T == tc_water && p == pc_water)
      value = Td_reg3(T, dc_water, o_id);
   else
      value = INVALID_VALUE;
}

// the extend input pairs
//  (p,v) (T,v)p (T,s) (T,h)
double pv_reg4(double p, double v, int o_id)
{
   if (o_id == OP)
      return p;
   if (o_id == OV)
      return v;
   if (o_id == OD)
      return 1/v;
   if (o_id == OR)
      return 4.0;
   double x = pv2x_reg4(p, v);
   if (o_id == OX)
   {
      return x;
   }
   return px_reg4(p, x, o_id);
};
double Tv_reg4(double T, double v, int o_id)
{
   if (o_id == OT)
      return T;
   if (o_id == OV)
      return v;
   if (o_id == OD)
      return 1/v;
   if (o_id == OR)
      return 4.0;
   double x = Tv2x_reg4(T, v);
   if (o_id == OX)
   {
      return x;
   }
   return Tx_reg4(T, x, o_id);
};
double Ts_reg4(double T, double s, int o_id)
{
   if (o_id == OT)
      return T;
   if (o_id == OS)
      return s;
   if (o_id == OR)
      return 4.0;
   double x = Ts2x_reg4(T, s);
   if (o_id == OX)
   {
      return x;
   }
   return Tx_reg4(T, x, o_id);
};
double Th_reg4(double T, double h, int o_id)
{
   if (o_id == OT)
      return T;
   if (o_id == OH)
      return h;
   if (o_id == OR)
      return 4.0;
   double x = Th2x_reg4(T, h);
   if (o_id == OX)
   {
      return x;
   }
   return Tx_reg4(T, x, o_id);
};

// (h,x,o_id) T
double hx_reg4(double h, double x, int o_id)
{
    if (o_id == OH)
        return h;
    if (o_id == OX)
        return x;
    if (o_id == OR)
        return 4.0;
    double Tl = TMIN4;
    double Tr = TMAX4;
    double T = bisection_reg4(x, h, OH, Tl, Tr, 0.00001, 1000);
    if (o_id == OT)
        return T;
    return Tx_reg4(T, x, o_id);
};

/// (s,x,o_id)
double sx_reg4(double s, double x, int o_id)
{
    if (o_id == OS)
        return s;
    if (o_id == OX)
        return x;
    if (o_id == OR)
        return 4.0;
    double Tl = TMIN4;
    double Tr = TMAX4;
    double T = bisection_reg4(x, s, OS, Tl, Tr, 0.00001, 1000);
    if (o_id == OT)
        return T;
    return Tx_reg4(T, x, o_id);
}
