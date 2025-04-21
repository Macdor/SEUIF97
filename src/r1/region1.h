/* ----------------------------------------------------
  IAPWS-IF97 Region1:
 1: IAPWS, R7-97(2012)
      IF97-Rev.pdf: P6-9
               1)fundamental: (p,t)->v,u,s,h,cp,cv,w
               2)backward: (p,h)->T, (p,s)->T

2: IAPWS, SR2-01(2014)
             Supp-PHS12-2014.pdf  (h,s)->p

  Author： Cheng Maohua
  Email: cmh@seu.edu.cn
*/
#pragma once
#include "../common/common.h"

static const double r1pstar = 16.53;  // MPa
static const double r1Tstar = 1386.0; // K

double gamma_reg1(double pi, double tau);
double gamma_pi_reg1(double pi, double tau);
double gamma_pipi_reg1(double pi, double tau);
double gamma_pitau_reg1(double pi, double tau);
double gamma_tau_reg1(double pi, double tau);
double gamma_tautau_reg1(double pi, double tau);

void polys_0_i_ii_j_jj_ij_reg1(double pi, double tau, double *gamma, double *gamma_pi, double *gamma_pipi, double *gamma_tau, double *gamma_tautau, double *gamma_pitau);

// IF97 fundamental: (p,t)
double pT2v_reg1(double p, double T);
double pT2u_reg1(double p, double T);
double pT2s_reg1(double p, double T);
double pT2h_reg1(double p, double T);
double pT2cv_reg1(double p, double T);
double pT2cp_reg1(double p, double T);
double pT2w_reg1(double p, double T);

// IF97 extended: (p,t)
double pT2k_reg1(double p, double T);
double pT2ipcec_reg1(double p, double T);
double pT2kt_reg1(double p, double T);

double pT2g_reg1(double p, double T);
double pT2f_reg1(double p, double T);

double pT2joule_reg1(double p, double T);
double pT2ijoule_reg1(double p, double T);

double pT2z_reg1(double p, double T);
double pT2e_reg1(double p, double T);

double pT2dpdtcv_reg1(double p, double T);
double pT2dvdpct_reg1(double p, double T);
double pT2dvdtcp_reg1(double p, double T);
double pT2dvdhcp_reg1(double p, double T);
double pT2dvdpch_reg1(double p, double T);

// IF97 backward: (p,h)->T, (p,s)->T,（h,s)->p
double ph2T_reg1(double p, double h);
double ps2T_reg1(double p, double s);
double hs2p_reg1(double h, double s);

// the extend input pairs
//  (p,v)->T (T,v)->p (T,s)->p
double pv2T_reg1(double p, double v);
double Tv2p_reg1(double T,double v);
double Ts2p_reg1(double T, double s);
double Th2p_reg1(double T, double h);

// (p,T)-> properties
double pT_ext_reg1(double p, double T, int o_id);
double pT_thermal_reg1(double p, double T, int o_id);

// pairs
double pT_reg1(double p, double T, int o_id);
double ph_reg1(double p, double h, int o_id);
double ps_reg1(double p, double h, int o_id);
double hs_reg1(double p, double h, int o_id);
// the extend input pairs
//  (p,v) (T,v)p (T,s)
double pv_reg1(double p, double v, int o_id);
double Tv_reg1(double T, double v, int o_id);
double Ts_reg1(double T, double s, int o_id);
double Th_reg1(double T, double h, int o_id);
