#include <math.h>
#include "algorithm.h"

double poly(double vi, double vj, int size, IJnData *IJn)
{
  double v = 0.0;
  for (int k = 0; k < size; k++)
    v += IJn[k].n * ipowsac(vi, IJn[k].I) * ipowsac(vj, IJn[k].J);
  return v;
}

double poly_i(double vi, double vj, int size, IJnData *IJn)
{
  double v = 0.0;
  for (int k = 0; k < size; k++)
    v += IJn[k].n * IJn[k].I * ipowsac(vi, IJn[k].I - 1) * ipowsac(vj, IJn[k].J);
  return v;
}

double poly_ii(double vi, double vj, int size, IJnData *IJn)
{
  double v = 0.0;
  for (int k = 0; k < size; k++)
    v += IJn[k].n * IJn[k].I * (IJn[k].I - 1) * ipowsac(vi, IJn[k].I - 2) * ipowsac(vj, IJn[k].J);
  return v;
}

double poly_j(double vi, double vj, int size, IJnData *IJn)
{
  double v = 0.0;
  for (int k = 0; k < size; k++)
    v += IJn[k].n * ipowsac(vi, IJn[k].I) * IJn[k].J * ipowsac(vj, IJn[k].J - 1);
  return v;
}

double poly_jj(double vi, double vj, int size, IJnData *IJn)
{
  double v = 0.0;
  for (int k = 0; k < size; k++)
    v += IJn[k].n * ipowsac(vi, IJn[k].I) * IJn[k].J * (IJn[k].J - 1) * ipowsac(vj, IJn[k].J - 2);
  return v;
}

double poly_ij(double vi, double vj, int size, IJnData *IJn)
{
  double v = 0.0;
  for (int k = 0; k < size; k++)
    v += IJn[k].n * IJn[k].I * ipowsac(vi, IJn[k].I - 1) * IJn[k].J * ipowsac(vj, IJn[k].J - 1);
  return v;
}

void polys_0_i_ii_j_jj_ij(double vi, double vj, int size, IJnData *IJn, double *poly_0,
    double *poly_i, double *poly_ii, double *poly_j, double *poly_jj, double *poly_ij)
{
  double item = 0.0;
 
  for (int k = 0; k < size; k++)
  {
    item = IJn[k].n * ipowsac(vi, IJn[k].I) * ipowsac(vj, IJn[k].J);
    *poly_0 += item;
    *poly_i += item*IJn[k].I;
    *poly_ii += item*IJn[k].I*(IJn[k].I-1);
    *poly_j += item*IJn[k].J;
    *poly_jj += item*IJn[k].J*(IJn[k].J-1);
    *poly_ij += item*IJn[k].I * IJn[k].J;
  }
}