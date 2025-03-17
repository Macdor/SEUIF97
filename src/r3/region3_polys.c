/*
 Speciphic Helmholtz free energy and derivatives

*/
#include <math.h>
#include "region3.h"

/*
  The solo power in region 3
   -   use Horner’s rule to get the soI_pow,soJ_pow quickly
 */ 
void solo_i_j_power_reg3(double vi, double vj, double soI_pow[], double soJ_pow[])
{
    // [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11];
    soI_pow[0] = 1.0;
    for (int k = 1; k <= 12; k++)
    {
        soI_pow[k] = soI_pow[k - 1] * vi;
    }

    // j  [0, 1, 2, 7, 10, 12, 23, 6, 15, 17, 22, 26, 4, 16, 3];
    soJ_pow[0] = 1.0;
    soJ_pow[1] = vj;
    soJ_pow[2] = vj * vj;

    double J_pow3 = soJ_pow[1] * soJ_pow[2];
    double J_pow4 = soJ_pow[2] * soJ_pow[2];

    soJ_pow[3] = J_pow3 * J_pow4;         // 7
    soJ_pow[4] = J_pow3 * soJ_pow[3];     // 10
    soJ_pow[5] = soJ_pow[4] * soJ_pow[2]; // 12

    double J_pow22 = soJ_pow[5] * soJ_pow[4];

    soJ_pow[6] = J_pow22 * vj;            // 23
    soJ_pow[7] = J_pow3 * J_pow3;         // 6
    soJ_pow[8] = J_pow3 * soJ_pow[5];     // 15
    soJ_pow[9] = soJ_pow[8] * soJ_pow[2]; // 17
    soJ_pow[10] = J_pow22;                // 22

    soJ_pow[11] = J_pow22 * J_pow4; // 26
    soJ_pow[12] = J_pow4;           // 4
    soJ_pow[13] = soJ_pow[8] * vj;  // 16
    soJ_pow[14] = J_pow3;           // 3
}

static const int soI[12] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
static int i2soI[39] = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
    2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4,
    5, 5, 5, 6, 6, 6, 7, 8, 9, 9, 10, 10, 11};

static const int soJ[15] = {0, 1, 2, 7, 10, 12, 23, 6, 15, 17, 22, 26, 4, 16, 3};
static int j2soJ[39] = {0, 1, 2, 3, 4, 5, 6, 2, 7, 8, 9, 0, 2,
    7, 3, 10, 11, 0, 2, 12, 13, 11, 0, 2, 12, 11, 1, 14, 11,
    0, 2, 11, 2, 11, 2, 11, 0, 1, 11};

	static const double n1 = 0.10658070028513e1;
	// 39
	static IJnData IJn[] = {
		{0, 0, -0.15732845290239E+02},
		{0, 1, 0.20944396974307E+02},
		{0, 2, -0.76867707878716E+01},
		{0, 7, 0.26185947787954E+01}
	
		,
		{0, 10, -0.28080781148620E+01},
		{0, 12, 0.12053369696517E+01},
		{0, 23, -0.84566812812502E-02},
		{1, 2, -0.12654315477714E+01},
		{1, 6, -0.11524407806681E+01}
	
		,
		{1, 15, 0.88521043984318E+00},
		{1, 17, -0.64207765181607E+00},
		{2, 0, 0.38493460186671E+00},
		{2, 2, -0.85214708824206E+00},
		{2, 6, 0.48972281541877E+01}
	
		,
		{2, 7, -0.30502617256965E+01},
		{2, 22, 0.39420536879154E-01},
		{2, 26, 0.12558408424308E+00},
		{3, 0, -0.27999329698710E+00},
		{3, 2, 0.13899799569460E+01}
	
		,
		{3, 4, -0.20189915023570E+01},
		{3, 16, -0.82147637173963E-02},
		{3, 26, -0.47596035734923E+00},
		{4, 0, 0.43984074473500E-01},
		{4, 2, -0.44476435428739E+00}
	
		,
		{4, 4, 0.90572070719733E+00},
		{4, 26, 0.70522450087967E+00},
		{5, 1, 0.10770512626332E+00},
		{5, 3, -0.32913623258954E+00},
		{5, 26, -0.50871062041158E+00}
	
		,
		{6, 0, -0.22175400873096E-01},
		{6, 2, 0.94260751665092E-01},
		{6, 26, 0.16436278447961E+00},
		{7, 2, -0.13503372241348E-01},
		{8, 26, -0.14834345352472E-01}
	
		,
		{9, 2, 0.57922953628084E-03},
		{9, 26, 0.32308904703711E-02},
		{10, 0, 0.80964802996215E-04},
		{10, 1, -0.16557679795037E-03},
		{11, 26, -0.44923899061815E-04}};

/*
	 * Speciphic Helmholtz free energy.
	 * tau :dimensionless temperature [K]
	 * delta: dimensionless density [kg/m3]

*/
double phi_reg3(double delta, double tau)
// Fundamental equation for region 3
{
	// double phi = n1 * log(delta);
	// for (int k = 0; k < 39; k++)
	//		phi += IJn[k].n * pow(delta, IJn[k].I) * pow(tau, IJn[k].J);
	// return phi;
	return n1 * log(delta) + poly_solo(delta, tau, 39, IJn, i2soI, j2soJ, solo_i_j_power_reg3);
}

double phi_delta_reg3(double delta, double tau)
// first derivative in delta of fundamental equation for region 3
{
	// double phi_delta = n1 / delta;
	// for (int k = 0; k < 39; k++)
	//		phi_delta += IJn[k].n * IJn[k].I * pow(delta, IJn[k].I - 1) * pow(tau, IJn[k].J);
	//	return phi_delta;
	return (n1 / delta) + poly_solo_i(delta, tau, 39, IJn, i2soI, j2soJ, solo_i_j_power_reg3);
}

double phi_deltadelta_reg3(double delta, double tau)
// Second derivative in delta of fundamental equation for region 3
{
	// double phi_deltadelta = -n1 / delta / delta;
	// for (int k = 0; k < 39; k++)
	//		phi_deltadelta += IJn[k].n * IJn[k].I * (IJn[k].I - 1) * pow(delta, IJn[k].I - 2) * pow(tau, IJn[k].J);
	//	return phi_deltadelta;

	return (-n1 / delta / delta) + poly_solo_ii(delta, tau, 39, IJn, i2soI, j2soJ, solo_i_j_power_reg3);
}

double phi_tau_reg3(double delta, double tau)
// phirst derivative in tau of fundamental equation for region 3
{
	// double phi_tau = 0.0;
	// for (int k = 0; k < 39; k++)
	//	phi_tau += IJn[k].n * pow(delta, IJn[k].I) * IJn[k].J * pow(tau, IJn[k].J - 1);
	// return phi_tau;

	return poly_solo_j(delta, tau, 39, IJn, i2soI, j2soJ, solo_i_j_power_reg3);
}

double phi_tautau_reg3(double delta, double tau)
// Second derivative in tau of fundamental equation for region 3
{
	// double phi_tautau = 0.0;
	// for (int k = 0; k < 39; k++)
	//	phi_tautau += IJn[k].n * pow(delta, IJn[k].I) * IJn[k].J * (IJn[k].J - 1) * pow(tau, IJn[k].J - 2);
	//	return phi_tautau;

	return poly_solo_jj(delta, tau, 39, IJn, i2soI, j2soJ, solo_i_j_power_reg3);
}

double phi_deltatau_reg3(double delta, double tau)
// Second derivative in delta and tau of fundamental equation for region 3
{
	// double phi_deltatau = 0.0;
	// for (int k = 0; k < 39; k++)
	//	phi_deltatau += IJn[k].n * IJn[k].I * pow(delta, IJn[k].I - 1) * IJn[k].J * pow(tau, IJn[k].J - 1);
	// return phi_deltatau;

	return poly_solo_ij(delta, tau, 39, IJn, i2soI, j2soJ, solo_i_j_power_reg3);
}

//------------------------------------------------------------------------------------
// mutiple polys
//-----------------------------------------------------------------------------------------
void polys_solo_0_j_reg3(double delta, double tau,
						 double *poly,
						 double *poly_tau)
{
	polys_solo_0_j(delta, tau, 39, IJn, i2soI, j2soJ, solo_i_j_power_reg3, poly, poly_tau);
	*poly += n1 * log(delta);
}

void polys_solo_i_j_reg3(double delta, double tau,
						 double *poly_delta,
						 double *poly_tau)
{
	polys_solo_i_j(delta, tau, 39, IJn, i2soI, j2soJ, solo_i_j_power_reg3, poly_delta, poly_tau);

	*poly_delta += n1 / delta;
}

void polys_solo_i_ii_ij_jj_reg3(double delta, double tau,
								double *poly_delta, double *poly_deltadelta, double *poly_deltatau,
								double *poly_tautau)
{
	polys_solo_i_ii_ij_jj(delta, tau, 39, IJn, i2soI, j2soJ, solo_i_j_power_reg3, poly_delta, poly_deltadelta,
						  poly_deltatau, poly_tautau);

	*poly_delta += n1 / delta;
	*poly_deltadelta += (-n1 / delta / delta);
}