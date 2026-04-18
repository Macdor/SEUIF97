# The code snippets of the acceleration methods

## Key Acceleration Methods

* Shortest Addition Chain Algorithm: Used for the rapid computation of integer powers via optimal multiplication sequences.
* Power Caching for Multi-Polynomial Evaluation: Precomputes powers to avoid redundant calculations, directly retrieving them to boost performance across multiple polynomials.
* Recurrence Method for Multi-Polynomial Evaluation：Utilizes the derivative relationship to compute a base polynomial, then derives the rest via base scaling  to maximize efficiency

## Shortest Addition Chain Algorithm

The code snippets demonstrate the acceleration methods to calculate the $T(p,h)$ in region 1, illustrating the flow from the optimized kernel to the final physical property calculation:

###  The IAPWS-IF97 Backward Equations

The backward equation $T(p,h)$ for region 1 has the following dimensionless:

$$\frac{T(p,h)}{T^*} = \theta(\pi, \eta) = \sum_{i=1}^{20} n_i \pi^{I_i} (\eta + 1)^{J_i}$$

where

$$\theta = T/T^*, \ \pi = p/p^*, \text{ and } \eta = h/h^* \text{ with } T^* = 1 \text{ K}, \ p^* = 1 \text{ MPa}, \text{ and } h^* = 2500 \text{ kJ kg}^{-1}.$$

### C Code

```c
// algo/algotithm.h
#define IPOW ipowsac

// algl/asc.c
double possac(double x, int n)
{
...
}    

double ipowsac(double x, int n)
{
...
}

// algo/polynomial_solo.c
double poly(double vi, double vj, int size, IJnData *IJn)
{
  double v = 0.0;
  for (int k = 0; k < size; k++)
    v += IJn[k].n * IPOW(vi, IJn[k].I) * IPOW(vj, IJn[k].J);
  return v;
}

// r1/region1_T_phps.c
double ph2T_reg1(double p, double h)
{
  // Page 11, Table6 :Initialize coefficients and exponents (P,H)->T for region 1
  IJnData IJn[] = {
    ...
    };

  double pi, eta;
  double theta;
  pi = p / 1.0;
  eta = h / 2500.0 + 1.0;
  theta = poly(pi, eta, 20, IJn);
  return (1.0 *theta);
}
```

## Power Caching and Recurrence Method for Multi-Polynomial Evaluation

The code snippets demonstrate the acceleration methods to calculate the specific internal energy in region 1, illustrating the flow from the optimized kernel to the final physical property calculation:

### The IAPWS-IF97 Equations

The basic equation for this region 1 is a fundamental equation for the specific **Gibbs free energy** $g$. 

$$\frac{g(p,T)}{RT} = \gamma(\pi,\tau) = \sum_{i=1}^{34} n_i (7.1-\pi)^{I_i} (\tau-1.222)^{J_i}$$

where 

$\pi = p/p^{*}$

$\tau = T^{*}/T$

To derive the **specific internal energy**

$$u = g - T \left( \frac{\partial g}{\partial T} \right)_p - p \left( \frac{\partial g}{\partial p} \right)_T$$

$$\frac{u(\pi, \tau)}{RT} = \tau \gamma_{\tau} - \pi \gamma_{\pi}$$

### C code

```c
// algo/polynomial_solo.c
void polys_solo_i_j(double vi, double vj, int size, IJnData *IJn, int *i2soI, int *j2soJ, solo_power_fn solo_i_j_power, double *poly_i, double *poly_j)
{
  (*solo_i_j_power)(vi, vj, soI_pow, soJ_pow);
  double item = 0.0;
  for (int k = 0; k < size; k++)
  {
    item = IJn[k].n * soI_pow[i2soI[k]] * soJ_pow[j2soJ[k]];
    *poly_i += item * IJn[k].I;
    *poly_j += item * IJn[k].J;
  }
  *poly_i /= vi;
  *poly_j /= vj;
}

// r1/region1_solo_ij.h
static const int soI[13] = {0, 1, 2, 3, 4, 5, 8, 21, 23, 29, 30, 31, 32};
static int i2soI[34] = {0, 0, 0, 0, 0,
                      0, 0, 0, 1, 1,
                      1, 1, 1, 1, 2,
                      2, 2, 2, 2, 3,
                      3, 3, 4, 4, 4,
                      5, 6, 6, 7, 8,
                      9, 10, 11, 12};

static const int soJ[25] = {-2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -3, 17, -4, 6, -5, 10,
                            -8, -11, -6, -29, -31, -38, -39, -40, -41};

static int j2soJ[34] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 5, 10, 2, 3, 5, 11, 12, 2,
                      13, 14, 0, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};

void solo_ij_pow_reg1(double vi, double vj, double *soI_pow, double *soJ_pow);

// r1/region1_gef.c
void polys_solo_i_j_reg1(double pi, double tau, double *poly_pi, double *poly_tau)
{
    polys_solo_i_j(7.1 - pi, tau - 1.222, 34, IJn, i2soI, j2soJ, solo_ij_pow_reg1, poly_pi, poly_tau);
    *poly_pi = -(*poly_pi);
}

// r1/region1_pT.c
double pT2u_reg1(double p, double T)
// specific internal energy in region 1
{
  double pi, tau;
  tau = r1Tstar / T;
  pi = p / r1pstar;

  double gammapi = 0.0;
  double gammatau = 0.0;

  polys_solo_i_j_reg1(pi, tau, &gammapi, &gammatau);
  
  return rgas_water * T * (tau * gammatau - pi * gammapi);
}
```
