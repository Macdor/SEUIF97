/*
  The algorithms
*/
#pragma once

typedef struct
{
  int I, J;
  double n;
} IJnData;

double possac(double x, int n);
double ipowsac(double x, int n);
double posipowrqm(double x, unsigned int n);
double ipowrqm(double x, int n);

double poly(double vi, double vj, int size, IJnData *IJn);
double poly_i(double vi, double vj, int size, IJnData *IJn);
double poly_ii(double vi, double vj, int size, IJnData *IJn);
double poly_j(double vi, double vj, int size, IJnData *IJn);
double poly_jj(double vi, double vj, int size, IJnData *IJn);
double poly_ij(double vi, double vj, int size, IJnData *IJn);
// multiple
void polys_0_i_ii_j_jj_ij(double vi, double vj, int size, IJnData *IJn, double *poly_0,
  double *poly_i, double *poly_ii, double *poly_j, double *poly_jj, double *poly_ij);

static double xacc = 1.0E-08;
static int iMAX = 100;
typedef double (*callfunc)(double, double);

double rtsec1(callfunc func, double cVar2, double fr, double x1,
              double x2, double ft, double f, double xacc, int iMAX);
double rtsec2(callfunc func, double cVar1, double fr, double x1,
              double x2, double fl, double f, double xacc, int iMAX);
