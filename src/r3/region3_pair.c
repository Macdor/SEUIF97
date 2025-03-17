/*
  The extended pairs:
     (p,v) ->t
    (T,h) ->d  (T,s)->d
*/
#include <math.h>
#include <stdint.h>
#include <string.h>
#include "region3.h"

// Region 3  (p,v)->T using the secant method and refine adjust
// * p: pressure  MPa
// * v: specific volume m^3/kg
// * T: temperature  K
double pv2T_reg3(double p, double v)
{
  double T1, T2, d, f1, f2;
  T1 = TMIN3;
  T2 = B23_p2T(p);
  d = 1.0 / v;
  f1 = p - Td2p_reg3(T1, d);
  f2 = p - Td2p_reg3(T2, d);
  return rtsec1(Td2p_reg3, d, p, T1, T2, f1, f2, xacc, iMAX);
}

//     (t,h) ->d
double Th2d_reg3(double T, double h)
{
  double p1 = B23_T2p(T);
  double d1 = 1.0 / pT2v_reg3(p1, T);
  double p2 = PMAX3;
  double d2 = 1.0 / pT2v_reg3(p2, T);
  double f1 = h - Td2h_reg3(T, d1);
  double f2 = h - Td2h_reg3(T, d2);
  return rtsec2(Td2h_reg3, T, h, d1, d2, f1, f2, xacc, iMAX);
}

// Region 3  (T,s)->d using the secant method
//  * T: temperature  K
//  * s: specific entropy  kJ/(kg K)
//  * d: density    kg/m^3
double Ts2d_reg3(double T, double s)
{
  double d1 = 100.0;
  double d2 = 1.1 * d1;
  double ft = s - Td2s_reg3(T, d1);
  double f = s - Td2s_reg3(T, d2);
  return rtsec2(Td2s_reg3, T, s, d1, d2, ft, f, xacc, iMAX);
}

double hs2p3a_reg3(double h, double s)
/*  Backward equation for region 3a, P=f(h,s)
     h : Specific enthalpy [kJ/kg]
     s : Specific entropy [kJ/kgK]
      P : Pressure [MPa] */
{
    IJnData IJn[33] = {{0, 0, 0.770889828326934e1},
                       {0, 1, -0.260835009128688e2},
                       {0, 5, 0.267416218930389e3},
                       {1, 0, 0.172221089496844e2},
                       {1, 3, -0.293542332145970e3},
                       {1, 4, 0.614135601882478e3},
                       {1, 8, -0.610562757725674e5},
                       {1, 14, -0.651272251118219e8},
                       {2, 6, 0.735919313521937e5},
                       {2, 16, -0.116646505914191e11},
                       {3, 0, 0.355267086434461e2},
                       {3, 2, -0.596144543825955e3},
                       {3, 3, -0.475842430145708e3},
                       {4, 0, 0.696781965359503e2},
                       {4, 1, 0.335674250377312e3},
                       {4, 4, 0.250526809130882e5},
                       {4, 5, 0.146997380630766e6},
                       {5, 28, 0.538069315091534e20},
                       {6, 28, 0.143619827291346e22},
                       {7, 24, 0.364985866165994e20},
                       {8, 1, -0.254741561156775e4},
                       {10, 32, 0.240120197096563e28},
                       {10, 36, -0.393847464679496e30},
                       {14, 22, 0.147073407024852e25},
                       {18, 28, -0.426391250432059e32},
                       {20, 36, 0.194509340621077e39},
                       {22, 16, 0.666212132114896e24},
                       {22, 28, 0.706777016552858e34},
                       {24, 36, 0.175563621975576e42},
                       {28, 16, 0.108408607429124e29},
                       {28, 36, 0.730872705175151e44},
                       {32, 10, 0.159145847398870e25},
                       {32, 28, 0.377121605943324e41}};

    double nu = h / 2300 - 1.01;
    double sigma = s / 4.4 - 0.75;
    // double suma = 0;
    // for (int i = 0; i < 33; i++)
    // {
    //    suma += n[i] * pow(nu, I[i]) * pow(sigma, J[i]);
    // }
    double suma = poly(nu, sigma, 33, IJn);
    return (99 * suma);
}

double hs2p3b_reg3(double h, double s)
/*Backward equation for region 3b, P=f(h,s)
    h :  Specific enthalpy [kJ/kg]
    s : Specific entropy [kJ/kgK]
    P : Pressure [MPa]*/
{
    IJnData IJn[35] = {
        {-12, 2, 0.125244360717979e-12},
        {-12, 10, -0.126599322553713e-1},
        {-12, 12, 0.506878030140626e1},
        {-12, 14, 0.317847171154202e2},
        {-12, 20, -0.391041161399932e6},
        {-10, 2, -0.975733406392044e-10},
        {-10, 10, -0.186312419488279e2},
        {-10, 14, 0.510973543414101e3},
        {-10, 18, 0.373847005822362e6},
        {-8, 2, 0.299804024666572e-7},
        {-8, 8, 0.200544393820342e2},
        {-6, 2, -0.498030487662829e-5},
        {-6, 6, -0.102301806360030e2},
        {-6, 7, 0.552819126990325e2},
        {-6, 8, -0.206211367510878e3},
        {-5, 10, -0.794012232324823e4},
        {-4, 4, 0.782248472028153e1},
        {-4, 5, -0.586544326902468e2},
        {-4, 8, 0.355073647696481e4},
        {-3, 1, -0.115303107290162e-3},
        {-3, 3, -0.175092403171802e1},
        {-3, 5, 0.257981687748160e3},
        {-3, 6, -0.727048374179467e3},
        {-2, 0, 0.121644822609198e-3},
        {-2, 1, 0.393137871762692e-1},
        {-1, 0, 0.704181005909296e-2},
        {0, 3, -0.829108200698110e2},
        {2, 0, -0.265178818131250},
        {2, 1, 0.137531682453991e2},
        {5, 0, -0.522394090753046e2},
        {6, 1, 0.240556298941048e4},
        {8, 1, -0.227361631268929e5},
        {10, 1, 0.890746343932567e5},
        {14, 3, -0.239234565822486e8},
        {14, 7, 0.568795808129714e10},
    };

    double nu = h / 2800 - 0.681;
    double sigma = s / 5.3 - 0.792;
    // double suma = 0;
    // for (int i = 0; i < 35; i++)
    //{
    //     suma += n[i] * pow(nu, I[i]) * pow(sigma, J[i]);
    // }
    double suma = poly(nu, sigma, 35, IJn);
    return (16.6 / suma);
}

double hs2p_reg3(double h, double s)
/*  Backward equation for region 3, P=f(h,s)
        h : Specific enthalpy [kJ/kg]
        s : Specific entropy [kJ/kgK]
        P : Pressure [MPa]
  */
{

    double p;
    if (s <= sc_water)
        p = hs2p3a_reg3(h, s);
    else
        p = hs2p3b_reg3(h, s);

    return (p);
}

double Td2p_reg3(double T, double d)
//
// pressure in region 3
// preg3 in bar
// T: temperature in K
// density in kg/m^3
//
{
    double tau = tc_water / T;
    double delta = d / dc_water;
    // return 0.001 * d * rgas_water * T * delta * phidelta_reg3(delta,tau);
    double phidelta = phi_delta_reg3(delta, tau);
    return 0.001 * d * rgas_water * T * delta * phidelta;
}

double Td2u_reg3(double T, double d)
// speciphic internal energy in region 3
// energyreg3 in kJ/kg
{
    double tau = tc_water / T;
    double delta = d / dc_water;
    // return rgas_water * T * tau * phi_tau_reg3(delta,tau);
    double phi_tau = phi_tau_reg3(delta, tau);
    return rgas_water * T * tau * phi_tau;
}

double Td2s_reg3(double T, double d)
// speciphic entropy in region 3
// entropyreg3 in kJ/(kg K)
{
    double tau = tc_water / T;
    double delta = d / dc_water;
    // return rgas_water * (tau * phi_tau_reg3(delta,tau) - phi_reg3(delta,tau));

    double phi = 0.0;
    double phi_tau = 0.0;
    polys_solo_0_j_reg3(delta, tau, &phi, &phi_tau);
    return rgas_water * (tau * phi_tau - phi);
}

double Td2h_reg3(double T, double d)
// speciphic enthalpy in region 3
// enthalpyreg3 in kJ/kg
{
    double tau = tc_water / T;
    double delta = d / dc_water;
    // return rgas_water * T * (tau * phi_tau_reg3(delta,tau) + delta * phidelta_reg3(delta,tau));
    double phidelta = 0.0;
    double phi_tau = 0.0;
    polys_solo_i_j_reg3(delta, tau, &phidelta, &phi_tau);
    return rgas_water * T * (tau * phi_tau + delta * phidelta);
}

double Td2cp_reg3(double T, double d)
// speciphic isobaric heat capacity in region 3
// cpreg3 in kJ/(kg K)
{
    double tau = tc_water / T;
    double delta = d / dc_water;

    // double a = delta * (phidelta_reg3(delta,tau) - tau * phideltatau_reg3(delta,tau));
    // a *= a;
    // double b = delta * (2.0 * phidelta_reg3(delta,tau) + delta * phideltadelta_reg3(delta,tau));
    // return rgas_water * (-tau * tau * phi_tautau_reg3(delta,tau) + a / b);

    double poly_delta = 0;
    double poly_deltatau = 0;
    double poly_deltadelta = 0;
    double poly_tautau = 0;
    polys_solo_i_ii_ij_jj_reg3(delta, tau, &poly_delta, &poly_deltadelta, &poly_deltatau, &poly_tautau);

    double a = delta * (poly_delta - tau * poly_deltatau);
    a *= a;
    double b = delta * (2.0 * poly_delta + delta * poly_deltadelta);
    return rgas_water * (-tau * tau * poly_tautau + a / b);
}

double Td2cv_reg3(double T, double d)
// speciphic isochoric heat capacity in region 3
// cvreg3 in kJ/(kg K)
{
    double tau = tc_water / T;
    double delta = d / dc_water;
    // return rgas_water * (-tau * tau * phi_tautau_reg3(delta,tau));

    double phi_tautau = phi_tautau_reg3(delta, tau);
    return rgas_water * (-tau * tau * phi_tautau);
}

double Td2w_reg3(double T, double d)
// speed of sound in region 3 in m/s
{
    double tau = tc_water / T;
    double delta = d / dc_water;
    /*

        double a = delta * phidelta_reg3(delta,tau) - delta * tau * phideltatau_reg3(delta,tau);
        a *= a;
        double r=1000.0 * rgas_water * T * (2 * delta * phidelta_reg3(delta,tau) + delta * delta * phideltadelta_reg3(delta,tau) - a / (tau * tau * phi_tautau_reg3(delta,tau)));
        return sqrt(1000.0 * rgas_water * T * (2 * delta * phidelta_reg3(delta,tau) + delta * delta * phideltadelta_reg3(delta,tau) - a / (tau * tau * phi_tautau_reg3(delta,tau))));
    */
    double poly_delta = 0;
    double poly_deltatau = 0;
    double poly_deltadelta = 0;
    double poly_tautau = 0;

    polys_solo_i_ii_ij_jj_reg3(delta, tau, &poly_delta, &poly_deltadelta, &poly_deltatau, &poly_tautau);
    double a = delta * (poly_delta - tau * poly_deltatau);
    a *= a;
    double r = 1000.0 * rgas_water * T * (delta * (2.0 * poly_delta + delta * poly_deltadelta) - a / (tau * tau * poly_tautau));
    return sqrt(r);
}


double T_atRegionBoundary(double p, char *boundary)
{
    // p is pressure in MPa
    // boundary is one of 3ab, 3cd, ...
    // returns the temperature at the boundary in K

    double t = 0.0; // temperature at the boundary in K

    int I3ab[] = {0, 0, 1, 2, -1, -2};
    int I3cd[] = {0, 0, 1, 2, 3};
    int I3gh[] = {0, 0, 1, 2, 3, 4};
    int I3ij[] = {0, 0, 1, 2, 3, 4};
    int I3jk[] = {0, 0, 1, 2, 3, 4};
    int I3mn[] = {0, 0, 1, 2, 3};
    int I3op[] = {0, 0, 1, 2, -1, -2};
    int I3qu[] = {0, 0, 1, 2, 3};
    int I3rx[] = {0, 0, 1, 2, 3};
    int I3uv[] = {0, 0, 1, 2, 3};
    int I3wx[] = {0, 0, 1, 2, -1, -2};

    double n3ab[] = {0.0, 0.154793642129415e4, -0.187661219490113e3, 0.213144632222113e2,
                     -0.191887498864292e4, 0.918419702359447e3};
    double n3cd[] = {0.0, 0.585276966696349e3, 0.278233532206915e1, -0.127283549295878e-1,
                     0.159090746562729e-3};
    double n3gh[] = {0.0, -0.249284240900418e5, 0.428143584791546e4, -0.269029173140130e3,
                     0.751608051114157e1, -0.787105249910383e-1};
    double n3ij[] = {0.0, 0.584814781649163e3, -0.616179320924617, 0.260763050899562,
                     -0.587071076864459e-2, 0.515308185433082e-4};
    double n3jk[] = {0.0, 0.617229772068439e3, -0.770600270141675e1, 0.697072596851896,
                     -0.157391839848015e-1, 0.137897492684194e-3};
    double n3mn[] = {0.0, 0.535339483742384e3, 0.761978122720128e1, -0.158365725441648,
                     0.192871054508108e-2};
    double n3op[] = {0.0, 0.969461372400213e3, -0.332500170441278e3, 0.642859598466067e2,
                     0.773845935768222e3, -0.152313732937084e4};
    double n3qu[] = {0.0, 0.565603648239126e3, 0.529062258221222e1, -0.102020639611016,
                     0.122240301070145e-2};
    double n3rx[] = {0.0, 0.584561202520006e3, -0.102961025163669e1, 0.243293362700452,
                     -0.294905044740799e-2};
    double n3uv[] = {0.0, 0.528199646263062e3, 0.890579602135307e1, -0.222814134903755,
                     0.286791682263697e-2};
    double n3wx[] = {0.0, 0.728052609145380e1, 0.973505869861952e2, 0.147370491183191e2,
                     0.329196213998375e3, 0.873371668682417e3};

    if (strcmp(boundary, "3ab") == 0)
    {
        for (int i = 1; i <= 5; i++)
            t += n3ab[i] * ipowsac(log(p), I3ab[i]);
    }
    else if (strcmp(boundary, "3op") == 0)
    {
        for (int i = 1; i <= 5; i++)
            t += n3op[i] * ipowsac(log(p), I3op[i]);
    }
    else if (strcmp(boundary, "3ef") == 0)
    {
        t = 3.727888004 * (p - 22.064) + 647.096;
    }
    else if (strcmp(boundary, "3cd") == 0)
    {
        for (int i = 1; i <= 4; i++)
            t += n3cd[i] * ipowsac(p, I3cd[i]);
    }
    else if (strcmp(boundary, "3gh") == 0)
    {
        for (int i = 1; i <= 5; i++)
            t += n3gh[i] * ipowsac(p, I3gh[i]);
    }
    else if (strcmp(boundary, "3ij") == 0)
    {
        for (int i = 1; i <= 5; i++)
            t += n3ij[i] * ipowsac(p, I3ij[i]);
    }
    else if (strcmp(boundary, "3jk") == 0)
    {
        for (int i = 1; i <= 5; i++)
            t += n3jk[i] * ipowsac(p, I3jk[i]);
    }
    else if (strcmp(boundary, "3mn") == 0)
    {
        for (int i = 1; i <= 4; i++)
            t += n3mn[i] * ipowsac(p, I3mn[i]);
    }
    else if (strcmp(boundary, "3qu") == 0)
    {
        for (int i = 1; i <= 4; i++)
            t += n3qu[i] * ipowsac(p, I3qu[i]);
    }
    else if (strcmp(boundary, "3rx") == 0)
    {
        for (int i = 1; i <= 4; i++)
            t += n3rx[i] * ipowsac(p, I3rx[i]);
    }
    else if (strcmp(boundary, "3uv") == 0)
    {
        for (int i = 1; i <= 4; i++)
            t += n3uv[i] * ipowsac(p, I3uv[i]);
    }
    else if (strcmp(boundary, "3wx") == 0)
    {
        for (int i = 1; i <= 5; i++)
            t += n3wx[i] * ipowsac(log(p), I3wx[i]);
    }

    return t;
}

double pT2vSat_reg3(double p, double T, double x)
{
    // p is pressure in MPa
    // t is temperature in K
    // x = 0 is liquid, =1 is steam]  x : integer    Vapor quality [-]
    // returns density in kg/m3

    double v;
    char subRegion;
    // set the region
    // force the region if density of saturated water is asked for
    if (x == 0)
    {
        if (p < 19.0088)
        {
            subRegion = 'c';
        }
        else
        {
            if (p < 21.0434)
            {
                subRegion = 's';
            }
            else
            {
                if (p < 21.9316)
                    subRegion = 'u';
                else
                    subRegion = 'y';
            }
        }
    }
    else
    {
        if (p < 20.5)
        {
            subRegion = 't';
        }
        else
        {
            if (p < 21.0434)
                subRegion = 'r';
            else if (p < 21.9009)
                subRegion = 'x';
            else
                subRegion = 'z';
        }
    }

    v = Vpt_3subreg(p, T, subRegion);
    return v;
}

//-------------------------------------------------------------------------

char SubRegion3(double p, double t)
{
    // p is pressure in MPa
    // T is temperature in K
    // sets the subregion

    char subRegion;

    if (p > 40.0 && p <= 100.0)
    {
        double tB = T_atRegionBoundary(p, "3ab");
        if (t > tB)
            subRegion = 'b';
        else
            subRegion = 'a';
    }
    else if (p > 25.0)
    {
        double tBab = T_atRegionBoundary(p, "3ab");
        double tBcd = T_atRegionBoundary(p, "3cd");
        double tBef = T_atRegionBoundary(p, "3ef");

        if (t <= tBcd)
            subRegion = 'c';
        else if (t <= tBab)
            subRegion = 'd';
        else if (t <= tBef)
            subRegion = 'e';
        else
            subRegion = 'f';
    }
    else if (p > 23.5)
    {
        double tBcd = T_atRegionBoundary(p, "3cd");
        double tBef = T_atRegionBoundary(p, "3ef");
        double tBgh = T_atRegionBoundary(p, "3gh");
        double tBij = T_atRegionBoundary(p, "3ij");
        double tBjk = T_atRegionBoundary(p, "3jk");

        if (t <= tBcd)
            subRegion = 'c';
        else if (t <= tBgh)
            subRegion = 'g';
        else if (t <= tBef)
            subRegion = 'h';
        else if (t <= tBij)
            subRegion = 'i';
        else if (t <= tBjk)
            subRegion = 'j';
        else
            subRegion = 'k';
    }
    else if (p > 23.0)
    {
        double tBcd = T_atRegionBoundary(p, "3cd");
        double tBef = T_atRegionBoundary(p, "3ef");
        double tBgh = T_atRegionBoundary(p, "3gh");
        double tBij = T_atRegionBoundary(p, "3ij");
        double tBjk = T_atRegionBoundary(p, "3jk");

        if (t <= tBcd)
            subRegion = 'c';
        else if (t <= tBgh)
            subRegion = 'l';
        else if (t <= tBef)
            subRegion = 'h';
        else if (t <= tBij)
            subRegion = 'i';
        else if (t <= tBjk)
            subRegion = 'j';
        else
            subRegion = 'k';
    }
    else if (p > 22.5)
    {
        double tBcd = T_atRegionBoundary(p, "3cd");
        double tBef = T_atRegionBoundary(p, "3ef");
        double tBgh = T_atRegionBoundary(p, "3gh");
        double tBij = T_atRegionBoundary(p, "3ij");
        double tBjk = T_atRegionBoundary(p, "3jk");
        double tBmn = T_atRegionBoundary(p, "3mn");
        double tBop = T_atRegionBoundary(p, "3op");

        if (t <= tBcd)
            subRegion = 'c';
        else if (t <= tBgh)
            subRegion = 'l';
        else if (t <= tBmn)
            subRegion = 'm';
        else if (t <= tBef)
            subRegion = 'n';
        else if (t <= tBop)
            subRegion = 'o';
        else if (t <= tBij)
            subRegion = 'p';
        else if (t <= tBjk)
            subRegion = 'j';
        else
            subRegion = 'k';
    }
    else
    {
        double pSat97 = pSat(643.15);
        if (p > pSat97)
        {
            double tBcd = T_atRegionBoundary(p, "3cd");
            double tBqu = T_atRegionBoundary(p, "3qu");
            double tBrx = T_atRegionBoundary(p, "3rx");
            double tBjk = T_atRegionBoundary(p, "3jk");

            if (t <= tBcd)
                subRegion = 'c';
            else if (t <= tBqu)
            {
                subRegion = 'q';
                // 这里的判断算法有点问题，下面这个 subRegion会被重新设定u
                //  简便处理，这个判断后，立即返回
                return subRegion;
            }

            if (t > tBrx && t <= tBjk)
            {
                subRegion = 'r';
                // 这里的判断算法有点问题，下面这个 subRegion会被重新设定u
                // 简便处理，这个判断后，立即返回
                return subRegion;
            }

            if (t > tBjk)
                subRegion = 'k';

            // tBqu < T <= tBrx
            // small regions right around critical point
            // 3u, 3x, 3y, 3z, 3v and 3w

            double tBuv = T_atRegionBoundary(p, "3uv");
            double tBwx = T_atRegionBoundary(p, "3wx");
            double tBef = T_atRegionBoundary(p, "3ef");

            if (p > 22.11)
            {
                if (t <= tBuv)
                    subRegion = 'u';
                else if (t < tBef)
                    subRegion = 'v';
                else if (t < tBwx)
                    subRegion = 'w';
                else
                    subRegion = 'x';
            }
            else if (p > 22.064)
            {
                double tBwx = T_atRegionBoundary(p, "3wx");
                double tBuv = T_atRegionBoundary(p, "3uv");
                if (t <= tBuv)
                    subRegion = 'u';
                else if (t <= tBef)
                    subRegion = 'y';
                else if (t <= tBwx)
                    subRegion = 'z';
                else
                    subRegion = 'x';
            }
            else
            {
                double tSat97 = TSat(p);

                if (t <= tSat97)
                {
                    if (p > 21.93161551)
                    {
                        if (t < tBuv)
                            subRegion = 'u';
                        else
                            subRegion = 'y';
                    }
                    else
                        subRegion = 'u';
                }
                else
                {
                    if (p > 21.90096265)
                    {
                        if (t < tBwx)
                            subRegion = 'z';
                        else
                            subRegion = 'x';
                    }
                    else
                        subRegion = 'x';
                }
            }
        }
        else if (p > 20.5)
        {
            double tBcd = T_atRegionBoundary(p, "3cd");
            double tBjk = T_atRegionBoundary(p, "3jk");
            double tSat97 = TSat(p);

            if (t <= tBcd)
                subRegion = 'c';
            else if (t <= tSat97)
                subRegion = 's';
            else if (t <= tBjk)
                subRegion = 'r';
            else
                subRegion = 'k';
        }
        else if (p > 19.00881189173929)
        {
            double tBcd = T_atRegionBoundary(p, "3cd");
            double tSat97 = TSat(p);

            if (t <= tBcd)
                subRegion = 'c';
            else if (t <= tSat97)
                subRegion = 's';
            else
                subRegion = 't';
        }
        else
        {
            double tSat97 = TSat(p);
            if (t <= tSat97)
                subRegion = 'c';
            else
                subRegion = 't';
        }
    }
    return subRegion;
}

double pT2v_reg3(double p, double T)
{
    // p is pressure in MPa
    // t is temperature in K
    // returns density in kg/m3
    char subRegion;
    // set the region
    subRegion = SubRegion3(p, T);
    double v = Vpt_3subreg(p, T, subRegion);
    return v;
}


//  Initialize coefphicients and exponents for region 3

double p2h3ab(double p)
{
  // Page 6 Table 2. Numerical values of the coefficients of the equation h3ab(p) in
  // its dimensionless form, Eq. (1), for defining the boundary
  // between subregions 3a and 3b
  static double n[4] = {0.201464004206875E+04,
                        0.374696550136983E+01,
                        -0.219921901054187E-01,
                        0.875131686009950E-04};
  double eta, pi;
  pi = p / 1.0;
  eta = n[0] + pi * (n[1] + pi * (n[2] + n[3] * pi));
  return (eta * 1.0);
}

double theta3aph(double pi, double eta)
{
  // (p,h)->T for region 3a
  // Page7 Table 3. Coefficients and exponents of the backward equation T3a(p,h) for subregion 3a in its
  // dimensionless form, Eq. (2)
  IJnData IJn[31] = {
      {-12, 0, -1.33645667811215e-7},
      {-12, 1, 4.55912656802978e-6},
      {-12, 2, -1.46294640700979e-5},
      {-12, 6, 6.39341312970080e-3},
      {-12, 14, 3.72783927268847e2},

      {-12, 16, -7.18654377460447e3},
      {-12, 20, 5.73494752103400e5},
      {-12, 22, -2.67569329111439e6},
      {-10, 1, -3.34066283302614e-5},
      {-10, 5, -2.45479214069597e-2},

      {-10, 12, 4.78087847764996e1},
      {-8, 0, 7.64664131818904e-6},
      {-8, 2, 1.28350627676972e-3},
      {-8, 4, 1.71219081377331e-2},
      {-8, 10, -8.51007304583213},

      {-5, 2, -1.36513461629781e-2},
      {-3, 0, -3.84460997596657e-6},
      {-2, 1, 3.37423807911655e-3},
      {-2, 3, -5.51624873066791e-1},
      {-2, 4, 7.29202277107470e-1},

      {-1, 0, -9.92522757376041e-3},
      {-1, 2, -1.19308831407288e-1},
      {0, 0, 7.93929190615421e-1},
      {0, 1, 4.54270731799386e-1},
      {1, 1, 2.09998591259910e-1},

      {3, 0, -6.42109823904738e-3},
      {3, 1, -2.35155868604540e-2},
      {4, 0, 2.52233108341612e-3},
      {4, 3, -7.64885133368119e-3},
      {10, 4, 1.36176427574291e-2},

      {12, 5, -1.33027883575669e-2}};

  pi = pi + 0.240;
  eta = eta - 0.615;
  // theta = 0.0;
  // for (int k = 0; k < 31; k++)
  //   theta += IJn[k].n * pow(pi, IJn[k].I) * pow(eta, IJn[k].J);

  double theta = poly(pi, eta, 31, IJn);
  return (theta);
}

double ph2T3a_reg3(double p, double h)
{
  double eta, pi;
  pi = p / 100.0;
  eta = h / 2300.0;
  return (760.0 * theta3aph(pi, eta));
}
//-------------------------------------------------
// (p,h)   3b
//--------------------------------------------------
double theta3bph(double pi, double eta)
{
  // Table 4. Coefficients and exponents of the backward equation T3b=(p,h) for subregion 3b in its
  // dimensionless form, Eq. (3)
  IJnData IJn[33] = {
      {-12, 0, 3.23254573644920e-5},
      {-12, 1, -1.27575556587181e-4},
      {-10, 0, -4.75851877356068e-4},
      {-10, 1, 1.56183014181602e-3},
      {-10, 5, 1.05724860113781e-1},

      {-10, 10, -8.58514221132534e1},
      {-10, 12, 7.24140095480911e2},
      {-8, 0, 2.96475810273257e-3},
      {-8, 1, -5.92721983365988e-3},
      {-8, 2, -1.26305422818666e-2},

      {-8, 4, -1.15716196364853e-1},
      {-8, 10, 8.49000969739595e1},
      {-6, 0, -1.08602260086615e-2},
      {-6, 1, 1.54304475328851e-2},
      {-6, 2, 7.50455441524466e-2},

      {-4, 0, 2.52520973612982e-2},
      {-4, 1, -6.02507901232996e-2},
      {-3, 5, -3.07622221350501},
      {-2, 0, -5.74011959864879e-2},
      {-2, 4, 5.03471360939849},

      {-1, 2, -9.25081888584834e-1},
      {-1, 4, 3.91733882917546},
      {-1, 6, -7.73146007130190e1},
      {-1, 10, 9.49308762098587e3},
      {-1, 14, -1.41043719679409e6},

      {-1, 16, 8.49166230819026e6},
      {0, 0, 8.61095729446704e-1},
      {0, 2, 3.23346442811720e-1},
      {1, 1, 8.73281936020439e-1},
      {3, 1, -4.36653048526683e-1},

      {5, 1, 2.86596714529479e-1},
      {6, 1, -1.31778331276228e-1},
      {8, 1, 6.76682064330275e-3}};

  pi = pi + 0.298;
  eta = eta - 0.720;
  // theta = 0.0;
  // for (int k = 0; k < 33; k++)
  // theta += IJn[k].n * pow(pi, IJn[k].I) * pow(eta, IJn[k].J);
  double theta = poly(pi, eta, 33, IJn);
  return (theta);
}

double ph2T3b_reg3(double p, double h)
{
  double eta, pi;
  pi = p / 100.0;
  eta = h / 2800.0;
  return (860.0 * theta3bph(pi, eta));
}
//----------------------------------------------
// Region 3a (p,h)->v
//----------------------------------------------
double omega3aph(double pi, double eta)
{
  // Page9 Table 6. Coefficients and exponents of the backward equation v3a(p,h) for subregion 3a in its
  // dimensionless form, Eq. (4)
  IJnData IJn[] = {
      {-12, 6, 5.29944062966028e-3},
      {-12, 8, -1.70099690234461e-1},
      {-12, 12, 1.11323814312927e1},
      {-12, 18, -2.17898123145125e3},
      {-10, 4, -5.06061827980875e-4},
      {-10, 7, 5.56495239685324e-1},
      {-10, 10, -9.43672726094016},
      {-8, 5, -2.97856807561527e-1},
      {-8, 12, 9.39353943717186e1},
      {-6, 3, 1.92944939465981e-2},
      {-6, 4, 4.21740664704763e-1},
      {-6, 22, -3.68914126282330e6},
      {-4, 2, -7.37566847600639e-3},
      {-4, 3, -3.54753242424366e-1},
      {-3, 7, -1.99768169338727},
      {-2, 3, 1.15456297059049},
      {-2, 16, 5.68366875815960e3},
      {-1, 0, 8.08169540124668e-3},
      {-1, 1, 1.72416341519307e-1},
      {-1, 2, 1.04270175292927},
      {-1, 3, -2.97691372792847e-1},
      {0, 0, 5.60394465163593e-1},
      {0, 1, 2.75234661176914e-1},
      {1, 0, -1.48347894866012e-1},
      {1, 1, -6.51142513478515e-2},
      {1, 2, -2.92468715386302},
      {2, 0, 6.64876096952665e-2},
      {2, 2, 3.52335014263844},
      {3, 0, -1.46340792313332e-2},
      {4, 2, -2.24503486668184},
      {5, 2, 1.10533464706142},
      {8, 2, -4.08757344495612e-2}};
  pi = pi + 0.128;
  eta = eta - 0.727;
  // omega = 0.0;
  // for (int k = 0; k < 32; k++)
  //   omega += IJn[k].n * pow(pi, IJn[k].I) * pow(eta, IJn[k].J);
  double omega = poly(pi, eta, 32, IJn);
  return (omega);
}

double ph2v3a_reg3(double p, double h)
{
  double eta, pi;
  pi = p / 100.0;
  eta = h / 2100.0;
  return (0.0028 * omega3aph(pi, eta));
}
//-------------------------------------------------
// (p,h)->v 3b
//--------------------------------------------------
double omega3bph(double pi, double eta)
{
  // Page 9
  // Table 7. Coefficients and exponents of the backward equation v3b (p,h) for subregion 3b in its
  // dimensionless form, Eq. (5)
  IJnData IJn[30] = {
      {-12, 0, -2.25196934336318e-9},
      {-12, 1, 1.40674363313486e-8},
      {-8, 0, 2.33784085280560e-6},
      {-8, 1, -3.31833715229001e-5},
      {-8, 3, 1.07956778514318e-3},
      {-8, 6, -2.71382067378863e-1},
      {-8, 7, 1.07202262490333},
      {-8, 8, -8.53821329075382e-1},
      {-6, 0, -2.15214194340526e-5},
      {-6, 1, 7.69656088222730e-4},
      {-6, 2, -4.31136580433864e-3},
      {-6, 5, 4.53342167309331e-1},
      {-6, 6, -5.07749535873652e-1},
      {-6, 10, -1.00475154528389e2},
      {-4, 3, -2.19201924648793e-1},
      {-4, 6, -3.21087965668917},
      {-4, 10, 6.07567815637771e2},
      {-3, 0, 5.57686450685932e-4},
      {-3, 2, 1.87499040029550e-1},
      {-2, 1, 9.05368030448107e-3},
      {-2, 2, 2.85417173048685e-1},
      {-1, 0, 3.29924030996098e-2},
      {-1, 1, 2.39897419685483e-1},
      {-1, 4, 4.82754995951394},
      {-1, 5, -1.18035753702231e1},
      {0, 0, 1.69490044091791e-1},
      {1, 0, -1.79967222507787e-2},
      {1, 1, 3.71810116332674e-2},
      {2, 2, -5.36288335065096e-2},
      {2, 6, 1.60697101092520}};
  pi = pi + 0.0661;
  eta = eta - 0.720;
  // omega = 0.0;
  // for (int k = 0; k < 30; k++)
  //  omega += IJn[k].n * pow(pi, IJn[k].I) * pow(eta, IJn[k].J);
  double omega = poly(pi, eta, 30, IJn);
  return (omega);
}

double ph2v3b_reg3(double p, double h)
{
  double eta, pi;
  pi = p / 100.0;
  eta = h / 2800.0;
  return (0.0088 * omega3bph(pi, eta));
}

//------------------------------------------------------------
// Region 3(3a,3b)  (p,h)->T
//------------------------------------------------------------

double ph2T_reg3(double p, double h)
{
  double T;
  double h3ab = p2h3ab(p);
  if (h <= h3ab)
  {
    T = ph2T3a_reg3(p, h);
  }
  else
  {
    T = ph2T3b_reg3(p, h);
  }
  return (T);
}

double ph2v_reg3(double p, double h)
{
  double v;
  double h3ab = p2h3ab(p);
  if (h <= h3ab)
  {
    v = ph2v3a_reg3(p, h);
  }
  else
  {
    v = ph2v3b_reg3(p, h);
  }
  return (v);
}

//------------------------------------------------------------
// (p,s)->T   3a Page 13, Table 11  Page 12 Eq(7)
//---------------------------------------------------------------
double theta3aps(double pi, double sigma)
{
  IJnData IJn[] = {{-12, 28, 0.150042008263875E+10},
                   {-12, 32, -0.159397258480424E+12},
                   {-10, 4, 0.502181140217975E-03},
                   {-10, 10, -0.672057767855466E+02},
                   {-10, 12, 0.145058545404456E+04},
                   {-10, 14, -0.823889534888890E+04},
                   {-8, 5, -0.154852214233853E+00},
                   {-8, 7, 0.112305046746695E+02},
                   {-8, 8, -0.297000213482822E+02},
                   {-8, 28, 0.438565132635495E+11},
                   {-6, 2, 0.137837838635464E-02},
                   {-6, 6, -0.297478527157462E+01},
                   {-6, 32, 0.971777947349413E+13},
                   {-5, 0, -0.571527767052398E-04},
                   {-5, 14, 0.288307949778420E+05},
                   {-5, 32, -0.744428289262703E+14},
                   {-4, 6, 0.128017324848921E+02},
                   {-4, 10, -0.368275545889071E+03},
                   {-4, 36, 0.664768904779177E+16},
                   {-2, 1, 0.449359251958880E-01},
                   {-2, 4, -0.422897836099655E+01},
                   {-1, 1, -0.240614376434179E+00},
                   {-1, 6, -0.474341365254924E+01},
                   {0, 0, 0.724093999126110E+00},
                   {0, 1, 0.923874349695897E+00},
                   {0, 4, 0.399043655281015E+01},
                   {1, 0, 0.384066651868009E-01},
                   {2, 0, -0.359344365571848E-02},
                   {2, 3, -0.735196448821653E+00},
                   {3, 2, 0.188367048396131E+00},
                   {8, 0, 0.141064266818704E-03},
                   {8, 1, -0.257418501496337E-02},
                   {10, 2, 0.123220024851555E-02}};
  pi = pi + 0.240;
  sigma = sigma - 0.703;
  // theta = 0.0;
  // for (int k = 0; k < 33; k++)
  //   theta += n[k] * pow(pi, i[k]) * pow(sigma, j[k]);
  double theta = poly(pi, sigma, 33, IJn);
  return (theta);
}

double ps2T3a_reg3(double p, double s)
{
  double pi, sigma;
  pi = p / 100.0;
  sigma = s / 4.4;
  return (760.0 * theta3aps(pi, sigma));
}
//------------------------------------------------------------
// (p,s)->T   3b  Page 13, Table 12, Page 12 Eq(7)
//---------------------------------------------------------------
double theta3bps(double pi, double sigma)
{
  IJnData IJn[28] = {{-12, 1, 0.527111701601660E+00},
                     {-12, 3, -0.401317830052742E+02},
                     {-12, 4, 0.153020073134484E+03},
                     {-12, 7, -0.224799398218827E+04},
                     {-8, 0, -0.193993484669048E+00},
                     {-8, 1, -0.140467557893768E+01},
                     {-8, 3, 0.426799878114024E+02},
                     {-6, 0, 0.752810643416743E+00},
                     {-6, 2, 0.226657238616417E+02},
                     {-6, 4, -0.622873556909932E+03},
                     {-5, 0, -0.660823667935396E+00},
                     {-5, 1, 0.841267087271658E+00},
                     {-5, 2, -0.253717501764397E+02},
                     {-5, 4, 0.485708963532948E+03},
                     {-5, 6, 0.880531517490555E+03},
                     {-4, 12, 0.265015592794626E+07},
                     {-3, 1, -0.359287150025783E+00},
                     {-3, 6, -0.656991567673753E+03},
                     {-2, 2, 0.241768149185367E+01},
                     {0, 0, 0.856873461222588E+00},
                     {2, 1, 0.655143675313458E+00},
                     {3, 1, -0.213535213206406E+00},
                     {4, 0, 0.562974957606348E-02},
                     {5, 24, -0.316955725450471E+15},
                     {6, 0, -0.699997000152457E-03},
                     {8, 3, 0.119845803210767E-01},
                     {12, 1, 0.193848122022095E-04},
                     {14, 2, -0.215095749182309E-04}};
  pi = pi + 0.760;
  sigma = sigma - 0.818;
  // theta = 0.0;
  // for (int k = 0; k < 28; k++)
  //   theta += n[k] * pow(pi, i[k]) * pow(sigma, j[k]);
  double theta = poly(pi, sigma, 28, IJn);
  return (theta);
}

double ps2T3b_reg3(double p, double s)
{
  double pi, sigma;
  pi = p / 100.0;
  sigma = s / 5.3;
  return (860.0 * theta3bps(pi, sigma));
}
//------------------------------------------------------------
// (p,s)->v  3a Page 15, Table 14  Page 14, Eq(8), Ok! 2003.12.18
//---------------------------------------------------------------
double omega3aps(double pi, double sigma)
{
  IJnData IJn[28] = {{-12, 10, 0.795544074093975E+02},
                     {-12, 12, -0.238261242984590E+04},
                     {-12, 14, 0.176813100617787E+05},
                     {-10, 4, -0.110524727080379E-02},
                     {-10, 8, -0.153213833655326E+02},
                     {-10, 10, 0.297544599376982E+03},
                     {-10, 20, -0.350315206871242E+08},
                     {-8, 5, 0.277513761062119},
                     {-8, 6, -0.523964271036888},
                     {-8, 14, -0.148011182995403E+06},
                     {-8, 16, 0.160014899374266E+07},
                     {-6, 28, 0.170802322663427E+13},
                     {-5, 1, 0.246866996006494E-03},
                     {-4, 5, 0.165326084797980E+01},
                     {-3, 2, -0.118008384666987},
                     {-3, 4, 0.253798642355900E+01},
                     {-2, 3, 0.965127704669424},
                     {-2, 8, -0.282172420532826E+02},
                     {-1, 1, 0.203224612353823},
                     {-1, 2, 0.110648186063513E+01},
                     {0, 0, 0.526127948451280},
                     {0, 1, 0.277000018736321},
                     {0, 3, 0.108153340501132E+01},
                     {1, 0, -0.744127885357893E-01},
                     {2, 0, 0.164094443541384E-01},
                     {4, 2, -0.680468275301065E-01},
                     {5, 2, 0.257988576101640E-01},
                     {6, 0, -0.145749861944416E-03}};

  pi = pi + 0.187;
  sigma = sigma - 0.755;
  // omega = 0.0;
  // for (int k = 0; k < 28; k++)
  //   omega += n[k] * pow(pi, i[k]) * pow(sigma, j[k]);
  double omega = poly(pi, sigma, 28, IJn);
  return (omega);
}

double ps2v3a_reg3(double p, double s)
{
  double pi, sigma;
  pi = p / 100.0;
  sigma = s / 4.4;
  return (0.0028 * omega3aps(pi, sigma));
}
//-------------------------------------------------------------------
// (p,s)->v  3b Page 15, Table 15  Page 14, Eq(9)   OK! 2003.12.18
//-------------------------------------------------------------------
double omega3bps(double pi, double sigma)
{
  IJnData IJn[] = {{-12, 0, 0.591599780322238E-04},
                   {-12, 1, -0.185465997137856E-02},
                   {-12, 2, 0.104190510480013E-01},
                   {-12, 3, 0.598647302038590E-02},
                   {-12, 5, -0.771391189901699},
                   {-12, 6, 0.172549765557036E+01},
                   {-10, 0, -0.467076079846526E-03},
                   {-10, 1, 0.134533823384439E-01},
                   {-10, 2, -0.808094336805495E-01},
                   {-10, 4, 0.508139374365767},
                   {-8, 0, 0.128584643361683E-02},
                   {-5, 1, -0.163899353915435E+01},
                   {-5, 2, 0.586938199318063E+01},
                   {-5, 3, -0.292466667918613E+01},
                   {-4, 0, -0.614076301499537E-02},
                   {-4, 1, 0.576199014049172E+01},
                   {-4, 2, -0.121613320606788E+02},
                   {-4, 3, 0.167637540957944E+01},
                   {-3, 1, -0.744135838773463E+01},
                   {-2, 0, 0.378168091437659E-01},
                   {-2, 1, 0.401432203027688E+01},
                   {-2, 2, 0.160279837479185E+02},
                   {-2, 3, 0.317848779347728E+01},
                   {-2, 4, -0.358362310304853E+01},
                   {-2, 12, -0.115995260446827E+07},
                   {0, 0, 0.199256573577909},
                   {0, 1, -0.122270624794624},
                   {0, 2, -0.191449143716586E+02},
                   {1, 0, -0.150448002905284E-01},
                   {1, 2, 0.146407900162154E+02},
                   {2, 2, -0.327477787188230E+01}};

  pi = pi + 0.298;
  sigma = sigma - 0.816;
  // omega = 0.0;
  // for (int k = 0; k < 31; k++)
  //   omega += n[k] * pow(pi, i[k]) * pow(sigma, j[k]);
  double omega = poly(pi, sigma, 31, IJn);
  return (omega);
}

double ps2v3b_reg3(double p, double s)
{
  double pi, sigma;
  pi = p / 100.0;
  sigma = s / 5.3;
  return (0.0088 * omega3bps(pi, sigma));
}

double ps2T_reg3(double p, double s)
{
  double T;
  if (s <= sc_water)
  {
    T = ps2T3a_reg3(p, s);
  }
  else
  {
    T = ps2T3b_reg3(p, s);
  }
  return (T);
}

double ps2v_reg3(double p, double s)
{
  double v;

  if (s <= sc_water)
  {
    v = ps2v3a_reg3(p, s);
  }
  else
  {
    v = ps2v3b_reg3(p, s);
  }
  return (v);
}


// The helper function to sum the sub-region (p,T)->v
// p is pressure in MPa
// t is temperature in K
// returns v is specific volume in Region 3a in m3/kg

typedef struct
{
    double PS, TS, VS;
    double a, b, c, d, e;
} para;

double pT2v_sum(double p, double T, int size, IJnData *IJn, para d)
{
    double p1 = pow((p / d.PS - d.a), d.c);
    double t1 = pow((T / d.TS - d.b), d.d);
    double v = poly(p1, t1, size, IJn);
    return powf(v, d.e) * d.VS;
}

double Vpt_3a(double p, double T)
{
    IJnData IJn[] = {
        {-12, 5, 0.110879558823853e-2},
        {-12, 10, 0.572616740810616e3},
        {-12, 12, -0.767051948380852e5},
        {-10, 5, -0.253321069529674e-1},
        {-10, 10, 0.628008049345689e4},
        {-10, 12, 0.234105654131876e6},
        {-8, 5, 0.216867826045856},
        {-8, 8, -0.156237904341963e3},
        {-8, 10, -0.269893956176613e5},
        {-6, 1, -0.180407100085505e-3},
        {-5, 1, 0.116732227668261e-2},
        {-5, 5, 0.266987040856040e2},
        {-5, 10, 0.282776617243286e5},
        {-4, 8, -0.242431520029523e4},
        {-3, 0, 0.435217323022733e-3},
        {-3, 1, -0.122494831387441e-1},
        {-3, 3, 0.179357604019989e1},
        {-3, 6, 0.442729521058314e2},
        {-2, 0, -0.593223489018342e-2},
        {-2, 2, 0.453186261685774},
        {-2, 3, 0.135825703129140e1},
        {-1, 0, 0.408748415856745e-1},
        {-1, 1, 0.474686397863312},
        {-1, 2, 0.118646814997915e1},
        {0, 0, 0.546987265727549},
        {0, 1, 0.195266770452643},
        {1, 0, -0.502268790869663e-1},
        {1, 2, -0.369645308193377},
        {2, 0, 0.633828037528420e-2},
        {2, 2, 0.797441793901017e-1},
    };

    para d =
        {
            .PS = 100.0,
            .TS = 760.0,
            .VS = 0.0024,
            .a = 0.085,
            .b = 0.817,
            .c = 1.0,
            .d = 1.0,
            .e = 1.0,
        };
    return pT2v_sum(p, T, 30, IJn, d);
}

double Vpt_3b(double p, double T)
{
    IJnData IJn[] = {
        {-12, 10, -0.827670470003621e-1},
        {-12, 12, 0.416887126010565e2},
        {-10, 8, 0.483651982197059e-1},
        {-10, 14, -0.291032084950276e5},
        {-8, 8, -0.111422582236948e3},
        {-6, 5, -0.202300083904014e-1},
        {-6, 6, 0.294002509338515e3},
        {-6, 8, 0.140244997609658e3},
        {-5, 5, -0.344384158811459e3},
        {-5, 8, 0.361182452612149e3},
        {-5, 10, -0.140699677420738e4},
        {-4, 2, -0.202023902676481e-2},
        {-4, 4, 0.171346792457471e3},
        {-4, 5, -0.425597804058632e1},
        {-3, 0, 0.691346085000334e-5},
        {-3, 1, 0.151140509678925e-2},
        {-3, 2, -0.416375290166236e-1},
        {-3, 3, -0.413754957011042e2},
        {-3, 5, -0.506673295721637e2},
        {-2, 0, -0.572212965569023e-3},
        {-2, 2, 0.608817368401785e1},
        {-2, 5, 0.239600660256161e2},
        {-1, 0, 0.122261479925384e-1},
        {-1, 2, 0.216356057692938e1},
        {0, 0, 0.398198903368642},
        {0, 1, -0.116892827834085},
        {1, 0, -0.102845919373532},
        {1, 2, -0.492676637589284},
        {2, 0, 0.655540456406790e-1},
        {3, 2, -0.240462535078530},
        {4, 0, -0.269798180310075e-1},
        {4, 1, 0.128369435967012}};

    para d =
        {
            .PS = 100.0,
            .TS = 860.0,
            .VS = 0.0041,
            .a = 0.280,
            .b = 0.779,
            .c = 1.0,
            .d = 1.0,
            .e = 1.0,
        };
    return pT2v_sum(p, T, 32, IJn, d);
}

double Vpt_3c(double p, double T)
{
    IJnData IJn[] = {
        {-12, 6, 3.11967788763030},
        {-12, 8, 2.76713458847564e+04},
        {-12, 10, 3.22583103403269e+07},
        {-10, 6, -3.42416065095363e+02},
        {-10, 8, -8.99732529907377e+05},
        {-10, 10, -7.93892049821251e+07},
        {-8, 5, 9.53193003217388e+01},
        {-8, 6, 2.29784742345072e+03},
        {-8, 7, 1.75336675322499e+05},
        {-6, 8, 7.91214365222792e+06},
        {-5, 1, 3.19933345844209e-05},
        {-5, 4, -6.59508863555767e+01},
        {-5, 7, -8.33426563212851e+05},
        {-4, 2, 6.45734680583292e-02},
        {-4, 8, -3.82031020570813e+06},
        {-3, 0, 4.06398848470079e-05},
        {-3, 3, 3.10327498492008e+01},
        {-2, 0, -8.92996718483724e-04},
        {-2, 4, 2.34604891591616e+02},
        {-2, 5, 3.77515668966951e+03},
        {-1, 0, 1.58646812591361e-02},
        {-1, 1, 7.07906336241843e-01},
        {-1, 2, 1.26016225146570e+01},
        {0, 0, 7.36143655772152e-01},
        {0, 1, 6.76544268999101e-01},
        {0, 2, -1.78100588189137e+01},
        {1, 0, -1.56531975531713e-01},
        {1, 2, 1.17707430048158e+01},
        {2, 0, 8.40143653860447e-02},
        {2, 1, -1.86442467471949e-01},
        {2, 3, -4.40170203949645e+01},
        {2, 7, 1.23290423502494e+06},
        {3, 0, -2.40650039730845e-02},
        {3, 7, -1.07077716660869e+06},
        {8, 1, 4.38319858566475e-02}};
    para d =
        {
            .PS = 40.0,
            .TS = 690.0,
            .VS = 0.0022,
            .a = 0.259,
            .b = 0.903,
            .c = 1.0,
            .d = 1.0,
            .e = 1.0,
        };
    return pT2v_sum(p, T, 35, IJn, d);
}

double Vpt_3d(double p, double T)
{
    IJnData IJn[] = {
        {-12, 4, -4.52484847171645e-10},
        {-12, 6, 3.15210389538801e-05},
        {-12, 7, -2.14991352047545e-03},
        {-12, 10, 5.08058874808345e+02},
        {-12, 12, -1.27123036845932e+07},
        {-12, 16, 1.15371133120497e+12},
        {-10, 0, -1.97805728776273e-16},
        {-10, 2, 2.41554806033972e-11},
        {-10, 4, -1.56481703640525e-06},
        {-10, 6, 2.77211346836625e-03},
        {-10, 8, -2.03578994462286e+01},
        {-10, 10, 1.44369489909053e+06},
        {-10, 14, -4.11254217946539e+10},
        {-8, 3, 6.23449786243773e-06},
        {-8, 7, -2.21774281146038e+01},
        {-8, 8, -6.89315087933158e+04},
        {-8, 10, -1.95419525060713e+07},
        {-6, 6, 3.16373510564015e+03},
        {-6, 8, 2.24040754426988e+06},
        {-5, 1, -4.36701347922356e-06},
        {-5, 2, -4.04213852833996e-04},
        {-5, 5, -3.48153203414663e+02},
        {-5, 7, -3.85294213555289e+05},
        {-4, 0, 1.35203700099403e-07},
        {-4, 1, 1.34648383271089e-04},
        {-4, 7, 1.25031835351736e+05},
        {-3, 2, 9.68123678455841e-02},
        {-3, 4, 2.25660517512438e+02},
        {-2, 0, -1.90102435341872e-04},
        {-2, 1, -2.99628410819229e-02},
        {-1, 0, 5.00833915372121e-03},
        {-1, 1, 3.87842482998411e-01},
        {-1, 5, -1.38535367777182e+03},
        {0, 0, 8.70745245971773e-01},
        {0, 2, 1.71946252068742},
        {1, 0, -3.26650121426383e-02},
        {1, 6, 4.98044171727877e+03},
        {3, 0, 5.51478022765087e-03},
    };
    para d =
        {
            .PS = 40.0,
            .TS = 690.0,
            .VS = 0.0029,
            .a = 0.559,
            .b = 0.939,
            .c = 1.0,
            .d = 1.0,
            .e = 4.0,
        };
    return pT2v_sum(p, T, 38, IJn, d);
}

double Vpt_3e(double p, double T)
{
    IJnData IJn[] = {
        {-12, 14, 7.15815808404721e+08},
        {-12, 16, -1.14328360753449e+11},
        {-10, 3, 3.76531002015720e-12},
        {-10, 6, -9.03983668691157e-05},
        {-10, 10, 6.65695908836252e+05},
        {-10, 14, 5.35364174960127e+09},
        {-10, 16, 7.94977402335603e+10},
        {-8, 7, 9.22230563421437e+01},
        {-8, 8, -1.42586073991215e+05},
        {-8, 10, -1.11796381424162e+06},
        {-6, 6, 8.96121629640760e+03},
        {-5, 6, -6.69989239070491e+03},
        {-4, 2, 4.51242538486834e-03},
        {-4, 4, -3.39731325977713e+01},
        {-3, 2, -1.20523111552278},
        {-3, 6, 4.75992667717124e+04},
        {-3, 7, -2.66627750390341e+05},
        {-2, 0, -1.53314954386524e-04},
        {-2, 1, 3.05638404828265e-01},
        {-2, 3, 1.23654999499486e+02},
        {-2, 4, -1.04390794213011e+03},
        {-1, 0, -1.57496516174308e-02},
        {0, 0, 6.85331118940253e-01},
        {0, 1, 1.78373462873903},
        {1, 0, -5.44674124878910e-01},
        {1, 4, 2.04529931318843e+03},
        {1, 6, -2.28342359328752e+04},
        {2, 0, 4.13197481515899e-01},
        {2, 2, -3.41931835910405e+01}};
    para d =
        {
            .PS = 40.0,
            .TS = 710.0,
            .VS = 0.0032,
            .a = 0.587,
            .b = 0.918,
            .c = 1.0,
            .d = 1.0,
            .e = 1.0,
        };
    return pT2v_sum(p, T, 29, IJn, d);
}

double Vpt_3f(double p, double T)
{
    IJnData IJn[] = {
        {0, -3, -2.51756547792325e-08},
        {0, -2, 6.01307193668763e-06},
        {0, -1, -1.00615977450049e-03},
        {0, 0, 9.99969140252192e-01},
        {0, 1, 2.14107759236486},
        {0, 2, -1.65175571959086e+01},
        {1, -1, -1.41987303638727e-03},
        {1, 1, 2.69251915156554},
        {1, 2, 3.49741815858722e+01},
        {1, 3, -3.00208695771783e+01},
        {2, 0, -1.31546288252539},
        {2, 1, -8.39091277286169},
        {3, -5, 1.81545608337015e-10},
        {3, -2, -5.91099206478909e-04},
        {3, 0, 1.52115067087106},
        {4, -3, 2.52956470663225e-05},
        {5, -8, 1.00726265203786e-15},
        {5, 1, -1.49774533860650},
        {6, -6, -7.93940970562969e-10},
        {7, -4, -1.50290891264717e-04},
        {7, 1, 1.51205531275133},
        {10, -6, 4.70942606221652e-06},
        {12, -10, 1.95049710391712e-13},
        {12, -8, -9.11627886266077e-09},
        {12, -4, 6.04374640201265e-04},
        {14, -12, -2.25132933900136e-16},
        {14, -10, 6.10916973582981e-12},
        {14, -8, -3.03063908043404e-07},
        {14, -6, -1.37796070798409e-05},
        {14, -4, -9.19296736666106e-04},
        {16, -10, 6.39288223132545e-10},
        {16, -8, 7.53259479898699e-07},
        {18, -12, -4.00321478682929e-13},
        {18, -10, 7.56140294351614e-09},
        {20, -12, -9.12082054034891e-12},
        {20, -10, -2.37612381140539e-08},
        {20, -6, 2.69586010591874e-05},
        {22, -12, -7.32828135157839e-11},
        {24, -12, 2.41995578306660e-10},
        {24, -4, -4.05735532730322e-04},
        {28, -12, 1.89424143498011e-10},
        {32, -12, -4.86632965074563e-10},
    };

    para d =
        {
            .PS = 40.0,
            .TS = 730.0,
            .VS = 0.0064,
            .a = 0.587,
            .b = 0.891,
            .c = 0.5,
            .d = 1.0,
            .e = 4.0,
        };
    return pT2v_sum(p, T, 42, IJn, d);
}

double Vpt_3g(double p, double T)
{
    IJnData IJn[] = {
        {-12, 7, 4.12209020652996e-05},
        {-12, 12, -1.14987238280587e+06},
        {-12, 14, 9.48180885032080e+09},
        {-12, 18, -1.95788865718971e+17},
        {-12, 22, 4.96250704871300e+24},
        {-12, 24, -1.05549884548496e+28},
        {-10, 14, -7.58642165988278e+11},
        {-10, 20, -9.22172769596101e+22},
        {-10, 24, 7.25379072059348e+29},
        {-8, 7, -6.17718249205859e+01},
        {-8, 8, 1.07555033344858e+04},
        {-8, 10, -3.79545802336487e+07},
        {-8, 12, 2.28646846221831e+11},
        {-6, 8, -4.99741093010619e+06},
        {-6, 22, -2.80214310054101e+30},
        {-5, 7, 1.04915406769586e+06},
        {-5, 20, 6.13754229168619e+27},
        {-4, 22, 8.02056715528378e+31},
        {-3, 7, -2.98617819828065e+07},
        {-2, 3, -9.10782540134681e+01},
        {-2, 5, 1.35033227281565e+05},
        {-2, 14, -7.12949383408211e+18},
        {-2, 24, -1.04578785289542e+36},
        {-1, 2, 3.04331584444093e+01},
        {-1, 8, 5.93250797959445e+09},
        {-1, 18, -3.64174062110798e+27},
        {0, 0, 9.21791403532461e-01},
        {0, 1, -3.37693609657471e-01},
        {0, 2, -7.24644143758508e+01},
        {1, 0, -1.10480239272601e-01},
        {1, 1, 5.36516031875059},
        {1, 3, -2.91441872156205e+03},
        {3, 24, 6.16338176535305e+39},
        {5, 22, -1.20889175861180e+38},
        {6, 12, 8.18396024524612e+22},
        {8, 3, 9.40781944835829e+08},
        {10, 0, -3.67279669545448e+04},
        {10, 6, -8.37513931798655e+15},
    };
    para d =
        {
            .PS = 25.0,
            .TS = 660.0,
            .VS = 0.0027,
            .a = 0.872,
            .b = 0.971,
            .c = 1.0,
            .d = 1.0,
            .e = 4.0,
        };
    return pT2v_sum(p, T, 38, IJn, d);
}

double Vpt_3h(double p, double T)
{
    IJnData IJn[] = {{-12, 8, 5.61379678887577e-02},
                     {-12, 12, 7.74135421587083e+09},
                     {-10, 4, 1.11482975877938e-09},
                     {-10, 6, -1.43987128208183e-03},
                     {-10, 8, 1.93696558764920e+03},
                     {-10, 10, -6.05971823585005e+08},
                     {-10, 14, 1.71951568124337e+13},
                     {-10, 16, -1.85461154985145e+16},
                     {-8, 0, 3.87851168078010e-17},
                     {-8, 1, -3.95464327846105e-14},
                     {-8, 6, -1.70875935679023e+02},
                     {-8, 7, -2.12010620701220e+03},
                     {-8, 8, 1.77683337348191e+07},
                     {-6, 4, 1.10177443629575e+01},
                     {-6, 6, -2.34396091693313e+05},
                     {-6, 8, -6.56174421999594e+06},
                     {-5, 2, 1.56362212977396e-05},
                     {-5, 3, -2.12946257021400},
                     {-5, 4, 1.35249306374858e+01},
                     {-4, 2, 1.77189164145813e-01},
                     {-4, 4, 1.39499167345464e+03},
                     {-3, 1, -7.03670932036388e-03},
                     {-3, 2, -1.52011044389648e-01},
                     {-2, 0, 9.81916922991113e-05},
                     {-1, 0, 1.47199658618076e-03},
                     {-1, 2, 2.02618487025578e+01},
                     {0, 0, 8.99345518944240e-01},
                     {1, 0, -2.11346402240858e-01},
                     {1, 2, 2.49971752957491e+01}};
    para d =
        {
            .PS = 25.0,
            .TS = 660.0,
            .VS = 0.0032,
            .a = 0.898,
            .b = 0.983,
            .c = 1.0,
            .d = 1.0,
            .e = 4.0,
        };
    return pT2v_sum(p, T, 29, IJn, d);
}

double Vpt_3i(double p, double T)
{
    IJnData IJn[] = {{0, 0, 1.06905684359136},
                     {0, 1, -1.48620857922333},
                     {0, 10, 2.59862256980408e+14},
                     {1, -4, -4.46352055678749e-12},
                     {1, -2, -5.66620757170032e-07},
                     {1, -1, -2.35302885736849e-03},
                     {1, 0, -2.69226321968839e-01},
                     {2, 0, 9.22024992944392},
                     {3, -5, 3.57633505503772e-12},
                     {3, 0, -1.73942565562222e+01},
                     {4, -3, 7.00681785556229e-06},
                     {4, -2, -2.67050351075768e-04},
                     {4, -1, -2.31779669675624},
                     {5, -6, -7.53533046979752e-13},
                     {5, -1, 4.81337131452891},
                     {5, 12, -2.23286270422356e+21},
                     {7, -4, -1.18746004987383e-05},
                     {7, -3, 6.46412934136496e-03},
                     {8, -6, -4.10588536330937e-10},
                     {8, 10, 4.22739537057241e+19},
                     {10, -8, 3.13698180473812e-13},
                     {12, -12, 1.64395334345040e-24},
                     {12, -6, -3.39823323754373e-06},
                     {12, -4, -1.35268639905021e-02},
                     {14, -10, -7.23252514211625e-15},
                     {14, -8, 1.84386437538366e-09},
                     {14, -4, -4.63959533752385e-02},
                     {14, 5, -9.92263100376750e+13},
                     {18, -12, 6.88169154439335e-17},
                     {18, -10, -2.22620998452197e-11},
                     {18, -8, -5.40843018624083e-08},
                     {18, -6, 3.45570606200257e-03},
                     {18, 2, 4.22275800304086e+10},
                     {20, -12, -1.26974478770487e-15},
                     {20, -10, 9.27237985153679e-10},
                     {22, -12, 6.12670812016489e-14},
                     {24, -12, -7.22693924063497e-12},
                     {24, -8, -3.83669502636822e-04},
                     {32, -10, 3.74684572410204e-04},
                     {32, -5, -9.31976897511086e+04},
                     {36, -10, -2.47690616026922e-02},
                     {36, -8, 6.58110546759474e+01}};
    para d =
        {
            .PS = 25.0,
            .TS = 660.0,
            .VS = 0.0041,
            .a = 0.910,
            .b = 0.984,
            .c = 0.5,
            .d = 1.0,
            .e = 4.0,
        };
    return pT2v_sum(p, T, 42, IJn, d);
}

double Vpt_3j(double p, double T)
{
    IJnData IJn[] = {{0, -1, -1.11371317395540e-04},
                     {0, 0, 1.00342892423685},
                     {0, 1, 5.30615581928979},
                     {1, -2, 1.79058760078792e-06},
                     {1, -1, -7.28541958464774e-04},
                     {1, 1, -1.87576133371704e+01},
                     {2, -1, 1.99060874071849e-03},
                     {2, 1, 2.43574755377290e+01},
                     {3, -2, -1.77040785499444e-04},
                     {4, -2, -2.59680385227130e-03},
                     {4, 2, -1.98704578406823e+02},
                     {5, -3, 7.38627790224287e-05},
                     {5, -2, -2.36264692844138e-03},
                     {5, 0, -1.61023121314333},
                     {6, 3, 6.22322971786473e+03},
                     {10, -6, -9.60754116701669e-09},
                     {12, -8, -5.10572269720488e-11},
                     {12, -3, 7.67373781404211e-03},
                     {14, -10, 6.63855469485254e-15},
                     {14, -8, -7.17590735526745e-10},
                     {14, -5, 1.46564542926508e-05},
                     {16, -10, 3.09029474277013e-12},
                     {18, -12, -4.64216300971708e-16},
                     {20, -12, -3.90499637961161e-14},
                     {20, -10, -2.36716126781431e-10},
                     {24, -12, 4.54652854268717e-12},
                     {24, -6, -4.22271787482497e-03},
                     {28, -12, 2.83911742354706e-11},
                     {28, -5, 2.70929002720228}};
    para d =
        {
            .PS = 25.0,
            .TS = 670.0,
            .VS = 0.0054,
            .a = 0.875,
            .b = 0.964,
            .c = 0.5,
            .d = 1.0,
            .e = 4.0,
        };
    return pT2v_sum(p, T, 29, IJn, d);
}

double Vpt_3k(double p, double T)
{
    IJnData IJn[] = {{-2, 10, -4.01215699576099e+08},
                     {-2, 12, 4.84501478318406e+10},
                     {-1, -5, 3.94721471363678e-15},
                     {-1, 6, 3.72629967374147e+04},
                     {0, -12, -3.69794374168666e-30},
                     {0, -6, -3.80436407012452e-15},
                     {0, -2, 4.75361629970233e-07},
                     {0, -1, -8.79148916140706e-04},
                     {0, 0, 8.44317863844331e-01},
                     {0, 1, 1.22433162656600e+01},
                     {0, 2, -1.04529634830279e+02},
                     {0, 3, 5.89702771277429e+02},
                     {0, 14, -2.91026851164444e+13},
                     {1, -3, 1.70343072841850e-06},
                     {1, -2, -2.77617606975748e-04},
                     {1, 0, -3.44709605486686},
                     {1, 1, 2.21333862447095e+01},
                     {1, 2, -1.94646110037079e+02},
                     {2, -8, 8.08354639772825e-16},
                     {2, -6, -1.80845209145470e-11},
                     {2, -3, -6.96664158132412e-06},
                     {2, -2, -1.81057560300994e-03},
                     {2, 0, 2.55830298579027},
                     {2, 4, 3.28913873658481e+03},
                     {5, -12, -1.73270241249904e-19},
                     {5, -6, -6.61876792558034e-07},
                     {5, -3, -3.95688923421250e-03},
                     {6, -12, 6.04203299819132e-18},
                     {6, -10, -4.00879935920517e-14},
                     {6, -8, 1.60751107464958e-09},
                     {6, -5, 3.83719409025556e-05},
                     {8, -12, -6.49565446702457e-15},
                     {10, -12, -1.49095328506000e-12},
                     {12, -10, 5.41449377329581e-09}};
    para d =
        {
            .PS = 25.0,
            .TS = 680.0,
            .VS = 0.0077,
            .a = 0.802,
            .b = 0.935,
            .c = 1.0,
            .d = 1.0,
            .e = 1.0,
        };
    return pT2v_sum(p, T, 34, IJn, d);
}

double Vpt_3l(double p, double T)
{
    IJnData IJn[] = {
        {-12, 14, 2.60702058647537e+09},
        {-12, 16, -1.88277213604704e+14},
        {-12, 18, 5.54923870289667e+18},
        {-12, 20, -7.58966946387758e+22},
        {-12, 22, 4.13865186848908e+26},
        {-10, 14, -8.15038000738060e+11},
        {-10, 24, -3.81458260489955e+32},
        {-8, 6, -1.23239564600519e-02},
        {-8, 10, 2.26095631437174e+07},
        {-8, 12, -4.95017809506720e+11},
        {-8, 14, 5.29482996422863e+15},
        {-8, 18, -4.44359478746295e+22},
        {-8, 24, 5.21635864527315e+34},
        {-8, 36, -4.87095672740742e+54},
        {-6, 8, -7.14430209937547e+05},
        {-5, 4, 1.27868634615495e-01},
        {-5, 5, -1.00752127917598e+01},
        {-4, 7, 7.77451437960990e+06},
        {-4, 16, -1.08105480796471e+24},
        {-3, 1, -3.57578581169659e-06},
        {-3, 3, -2.12857169423484},
        {-3, 18, 2.70706111085238e+29},
        {-3, 20, -6.95953622348829e+32},
        {-2, 2, 1.10609027472280e-01},
        {-2, 3, 7.21559163361354e+01},
        {-2, 10, -3.06367307532219e+14},
        {-1, 0, 2.65839618885530e-05},
        {-1, 1, 2.53392392889754e-02},
        {-1, 3, -2.14443041836579e+02},
        {0, 0, 9.37846601489667e-01},
        {0, 1, 2.23184043101700},
        {0, 2, 3.38401222509191e+01},
        {0, 12, 4.94237237179718e+20},
        {1, 0, -1.98068404154428e-01},
        {1, 16, -1.41415349881140e+30},
        {2, 1, -9.93862421613651e+01},
        {4, 0, 1.25070534142731e+02},
        {5, 0, -9.96473529004439e+02},
        {5, 1, 4.73137909872765e+04},
        {6, 14, 1.16662121219322e+32},
        {10, 4, -3.15874976271533e+15},
        {10, 12, -4.45703369196945e+32},
        {14, 10, 6.42794932373694e+32},
    };
    para d =
        {
            .PS = 24.0,
            .TS = 650.0,
            .VS = 0.0026,
            .a = 0.908,
            .b = 0.989,
            .c = 1.0,
            .d = 1.0,
            .e = 4.0,
        };
    return pT2v_sum(p, T, 43, IJn, d);
}

double Vpt_3m(double p, double T)
{
    IJnData IJn[] = {
        {0, 0, 8.11384363481847e-01},
        {3, 0, -5.68199310990094e+03},
        {8, 0, -1.78657198172556e+10},
        {20, 2, 7.95537657613427e+31},
        {1, 5, -8.14568209346872e+04},
        {3, 5, -6.59774567602874e+07},
        {4, 5, -1.52861148659302e+10},
        {5, 5, -5.60165667510446e+11},
        {1, 6, 4.58384828593949e+05},
        {6, 6, -3.85754000383848e+13},
        {2, 7, 4.53735800004273e+07},
        {4, 8, 9.39454935735563e+11},
        {14, 8, 2.66572856432938e+27},
        {2, 10, -5.47578313899097e+09},
        {5, 10, 2.00725701112386e+14},
        {3, 12, 1.85007245563239e+12},
        {0, 14, 1.85135446828337e+08},
        {1, 14, -1.70451090076385e+11},
        {1, 18, 1.57890366037614e+14},
        {1, 20, -2.02530509748774e+15},
        {28, 20, 3.68193926183570e+59},
        {2, 22, 1.70215539458936e+17},
        {16, 22, 6.39234909918741e+41},
        {0, 24, -8.21698160721956e+14},
        {5, 24, -7.95260241872306e+23},
        {0, 28, 2.33415869478510e+17},
        {3, 28, -6.00079934586803e+22},
        {4, 28, 5.94584382273384e+24},
        {12, 28, 1.89461279349492e+39},
        {16, 28, -8.10093428842645e+45},
        {1, 32, 1.88813911076809e+21},
        {8, 32, 1.11052244098768e+35},
        {14, 32, 2.91133958602503e+45},
        {0, 36, -3.29421923951460e+21},
        {2, 36, -1.37570282536696e+25},
        {3, 36, 1.81508996303902e+27},
        {4, 36, -3.46865122768353e+29},
        {8, 36, -2.11961148774260e+37},
        {14, 36, -1.28617899887675e+48},
        {24, 36, 4.79817895699239e+64}};
    para d =
        {
            .PS = 23.0,
            .TS = 650.0,
            .VS = 0.0028,
            .a = 1.0,
            .b = 0.997,
            .c = 1.0,
            .d = 0.25,
            .e = 1.0,
        };
    return pT2v_sum(p, T, 40, IJn, d);
}

double Vpt_3n(double p, double T)
{
    IJnData IJn[] = {{0, -12, 2.80967799943151e-39},
                     {3, -12, 6.14869006573609e-31},
                     {4, -12, 5.82238667048942e-28},
                     {6, -12, 3.90628369238462e-23},
                     {7, -12, 8.21445758255119e-21},
                     {10, -12, 4.02137961842776e-15},
                     {12, -12, 6.51718171878301e-13},
                     {14, -12, -2.11773355803058e-08},
                     {18, -12, 2.64953354380072e-03},
                     {0, -10, -1.35031446451331e-32},
                     {3, -10, -6.07246643970893e-24},
                     {5, -10, -4.02352115234494e-19},
                     {6, -10, -7.44938506925544e-17},
                     {8, -10, 1.89917206526237e-13},
                     {12, -10, 3.64975183508473e-06},
                     {0, -8, 1.77274872361946e-26},
                     {3, -8, -3.34952758812999e-19},
                     {7, -8, -4.21537726098389e-09},
                     {12, -8, -3.91048167929649e-02},
                     {2, -6, 5.41276911564176e-14},
                     {3, -6, 7.05412100773699e-12},
                     {4, -6, 2.58585887897486e-09},
                     {2, -5, -4.93111362030162e-11},
                     {4, -5, -1.58649699894543e-06},
                     {7, -5, -5.25037427886100e-01},
                     {4, -4, 2.20019901729615e-03},
                     {3, -3, -6.43064132636925e-03},
                     {5, -3, 6.29154149015048e+01},
                     {6, -3, 1.35147318617061e+02},
                     {0, -2, 2.40560808321713e-07},
                     {0, -1, -8.90763306701305e-04},
                     {3, -1, -4.40209599407714e+03},
                     {1, 0, -3.02807107747776e+02},
                     {0, 1, 1.59158748314599e+03},
                     {1, 1, 2.32534272709876e+05},
                     {0, 2, -7.92681207132600e+05},
                     {1, 4, -8.69871364662769e+10},
                     {0, 5, 3.54542769185671e+11},
                     {1, 6, 4.00849240129329e+14}};

    const double vStar = 0.0031;
    const double a = 0.976;
    const double b = 0.997;
    p = p / 23.0 - a;
    T = T / 650.0 - b;
    double v = poly(p, T, 39, IJn);
    return exp(v) * vStar;
}

double Vpt_3o(double p, double T)
{
    IJnData IJn[] = {{0, -12, 1.28746023979718e-35},
                     {0, -4, -7.35234770382342e-12},
                     {0, -1, 2.89078692149150e-03},
                     {2, -1, 2.44482731907223e-01},
                     {3, -10, 1.41733492030985e-24},
                     {4, -12, -3.54533853059476e-29},
                     {4, -8, -5.94539202901431e-18},
                     {4, -5, -5.85188401782779e-09},
                     {4, -4, 2.01377325411803e-06},
                     {4, -1, 1.38647388209306},
                     {5, -4, -1.73959365084772e-05},
                     {5, -3, 1.37680878349369e-03},
                     {6, -8, 8.14897605805513e-15},
                     {7, -12, 4.25596631351839e-26},
                     {8, -10, -3.87449113787755e-18},
                     {8, -8, 1.39814747930240e-13},
                     {8, -4, -1.71849638951521e-03},
                     {10, -12, 6.41890529513296e-22},
                     {10, -8, 1.18960578072018e-11},
                     {14, -12, -1.55282762571611e-18},
                     {14, -8, 2.33907907347507e-08},
                     {20, -12, -1.74093247766213e-13},
                     {20, -10, 3.77682649089149e-09},
                     {24, -12, -5.16720236575302e-11}};
    para d =
        {
            .PS = 23.0,
            .TS = 650.0,
            .VS = 0.0034,
            .a = 0.974,
            .b = 0.996,
            .c = 0.5,
            .d = 1.0,
            .e = 1.0,
        };
    return pT2v_sum(p, T, 24, IJn, d);
}

double Vpt_3p(double p, double T)
{
    IJnData IJn[] = {{0, -1, -9.82825342010366e-05},
                     {0, 0, 1.05145700850612},
                     {0, 1, 1.16033094095084e+02},
                     {0, 2, 3.24664750281543e+03},
                     {1, 1, -1.23592348610137e+03},
                     {2, -1, -5.61403450013495e-02},
                     {3, -3, 8.56677401640869e-08},
                     {3, 0, 2.36313425393924e+02},
                     {4, -2, 9.72503292350109e-03},
                     {6, -2, -1.03001994531927},
                     {7, -5, -1.49653706199162e-09},
                     {7, -4, -2.15743778861592e-05},
                     {8, -2, -8.34452198291445},
                     {10, -3, 5.86602660564988e-01},
                     {12, -12, 3.43480022104968e-26},
                     {12, -6, 8.16256095947021e-06},
                     {12, -5, 2.94985697916798e-03},
                     {14, -10, 7.11730466276584e-17},
                     {14, -8, 4.00954763806941e-10},
                     {14, -3, 1.07766027032853e+01},
                     {16, -8, -4.09449599138182e-07},
                     {18, -8, -7.29121307758902e-06},
                     {20, -10, 6.77107970938909e-09},
                     {22, -10, 6.02745973022975e-08},
                     {24, -12, -3.82323011855257e-11},
                     {24, -8, 1.79946628317437e-03},
                     {36, -12, -3.45042834640005e-04}};
    para d =
        {
            .PS = 23.0,
            .TS = 650.0,
            .VS = 0.0041,
            .a = 0.972,
            .b = 0.997,
            .c = 0.5,
            .d = 1.0,
            .e = 1.0,
        };
    return pT2v_sum(p, T, 27, IJn, d);
}

double Vpt_3q(double p, double T)
{
    IJnData IJn[] = {{-12, 10, -8.20433843259950e+04},
                     {-12, 12, 4.73271518461586e+10},
                     {-10, 6, -8.05950021005413e-02},
                     {-10, 7, 3.28600025435980e+01},
                     {-10, 8, -3.56617029982490e+03},
                     {-10, 10, -1.72985781433335e+09},
                     {-8, 8, 3.51769232729192e+07},
                     {-6, 6, -7.75489259985144e+05},
                     {-5, 2, 7.10346691966018e-05},
                     {-5, 5, 9.93499883820274e+04},
                     {-4, 3, -6.42094171904570e-01},
                     {-4, 4, -6.12842816820083e+03},
                     {-3, 3, 2.32808472983776e+02},
                     {-2, 0, -1.42808220416837e-05},
                     {-2, 1, -6.43596060678456e-03},
                     {-2, 2, -4.28577227475614},
                     {-2, 4, 2.25689939161918e+03},
                     {-1, 0, 1.00355651721510e-03},
                     {-1, 1, 3.33491455143516e-01},
                     {-1, 2, 1.09697576888873},
                     {0, 0, 9.61917379376452e-01},
                     {1, 0, -8.38165632204598e-02},
                     {1, 1, 2.47795908411492},
                     {1, 3, -3.19114969006533e+03}};
    para d =
        {
            .PS = 23.0,
            .TS = 650.0,
            .VS = 0.0022,
            .a = 0.848,
            .b = 0.983,
            .c = 1.0,
            .d = 1.0,
            .e = 4.0,
        };
    return pT2v_sum(p, T, 24, IJn, d);
}

double Vpt_3r(double p, double T)
{
    IJnData IJn[] = {{-8, 6, 1.44165955660863e-03},
                     {-8, 14, -7.01438599628258e+12},
                     {-3, -3, -8.30946716459219e-17},
                     {-3, 3, 2.61975135368109e-01},
                     {-3, 4, 3.93097214706245e+02},
                     {-3, 5, -1.04334030654021e+04},
                     {-3, 8, 4.90112654154211e+08},
                     {0, -1, -1.47104222772069e-04},
                     {0, 0, 1.03602748043408},
                     {0, 1, 3.05308890065089},
                     {0, 5, -3.99745276971264e+06},
                     {3, -6, 5.69233719593750e-12},
                     {3, -2, -4.64923504407778e-02},
                     {8, -12, -5.35400396512906e-18},
                     {8, -10, 3.99988795693162e-13},
                     {8, -8, -5.36479560201811e-07},
                     {8, -5, 1.59536722411202e-02},
                     {10, -12, 2.70303248860217e-15},
                     {10, -10, 2.44247453858506e-08},
                     {10, -8, -9.83430636716454e-06},
                     {10, -6, 6.63513144224454e-02},
                     {10, -5, -9.93456957845006},
                     {10, -4, 5.46491323528491e+02},
                     {10, -3, -1.43365406393758e+04},
                     {10, -2, 1.50764974125511e+05},
                     {12, -12, -3.37209709340105e-10},
                     {14, -12, 3.77501980025469e-09}};
    para d =
        {
            .PS = 23.0,
            .TS = 650.0,
            .VS = 0.0054,
            .a = 0.874,
            .b = 0.982,
            .c = 1.0,
            .d = 1.0,
            .e = 1.0,
        };
    return pT2v_sum(p, T, 27, IJn, d);
}

double Vpt_3s(double p, double T)
{
    IJnData IJn[] = {{-12, 20, -5.32466612140254e+22},
                     {-12, 24, 1.00415480000824e+31},
                     {-10, 22, -1.91540001821367e+29},
                     {-8, 14, 1.05618377808847e+16},
                     {-6, 36, 2.02281884477061e+58},
                     {-5, 8, 8.84585472596134e+07},
                     {-5, 16, 1.66540181638363e+22},
                     {-4, 6, -3.13563197669111e+05},
                     {-4, 32, -1.85662327545324e+53},
                     {-3, 3, -6.24942093918942e-02},
                     {-3, 8, -5.04160724132590e+09},
                     {-2, 4, 1.87514491833092e+04},
                     {-1, 1, 1.21399979993217e-03},
                     {-1, 2, 1.88317043049455},
                     {-1, 3, -1.67073503962060e+03},
                     {0, 0, 9.65961650599775e-01},
                     {0, 1, 2.94885696802488},
                     {0, 4, -6.53915627346115e+04},
                     {0, 28, 6.04012200163444e+49},
                     {1, 0, -1.98339358557937e-01},
                     {1, 32, -1.75984090163501e+57},
                     {3, 0, 3.56314881403987},
                     {3, 1, -5.75991255144384e+02},
                     {3, 2, 4.56213415338071e+04},
                     {4, 3, -1.09174044987829e+07},
                     {4, 18, 4.37796099975134e+33},
                     {4, 24, -6.16552611135792e+45},
                     {5, 4, 1.93568768917797e+09},
                     {14, 24, 9.50898170425042e+53}};
    para d =
        {
            .PS = 21.0,
            .TS = 640.0,
            .VS = 0.0022,
            .a = 0.886,
            .b = 0.990,
            .c = 1.0,
            .d = 1.0,
            .e = 4.0,
        };
    return pT2v_sum(p, T, 29, IJn, d);
}

double Vpt_3t(double p, double T)
{
    IJnData IJn[] = {{0, 0, 1.55287249586268},
                     {0, 1, 6.64235115009031},
                     {0, 4, -2.89366236727210e+03},
                     {0, 12, -3.85923202309848e+12},
                     {1, 0, -2.91002915783761},
                     {1, 10, -8.29088246858083e+11},
                     {2, 0, 1.76814899675218},
                     {2, 6, -5.34686695713469e+08},
                     {2, 14, 1.60464608687834e+17},
                     {3, 3, 1.96435366560186e+05},
                     {3, 8, 1.56637427541729e+12},
                     {4, 0, -1.78154560260006},
                     {4, 10, -2.29746237623692e+15},
                     {7, 3, 3.85659001648006e+07},
                     {7, 4, 1.10554446790543e+09},
                     {7, 7, -6.77073830687349e+13},
                     {7, 20, -3.27910592086523e+30},
                     {7, 36, -3.41552040860644e+50},
                     {10, 10, -5.27251339709047e+20},
                     {10, 12, 2.45375640937055e+23},
                     {10, 14, -1.68776617209269e+26},
                     {10, 16, 3.58958955867578e+28},
                     {10, 22, -6.56475280339411e+35},
                     {18, 18, 3.55286045512301e+38},
                     {20, 32, 5.69021454413270e+57},
                     {22, 22, -7.00584546433113e+47},
                     {22, 36, -7.05772623326374e+64},
                     {24, 24, 1.66861176200148e+52},
                     {28, 28, -3.00475129680486e+60},
                     {32, 22, -6.68481295196808e+50},
                     {32, 32, 4.28432338620678e+68},
                     {32, 36, -4.44227367758304e+71},
                     {36, 36, -2.81396013562745e+76}};
    para d =
        {
            .PS = 20.0,
            .TS = 650.0,
            .VS = 0.0088,
            .a = 0.803,
            .b = 1.02,
            .c = 1.0,
            .d = 1.0,
            .e = 1.0,
        };
    return pT2v_sum(p, T, 33, IJn, d);
}

double Vpt_3u(double p, double T)
{
    IJnData IJn[] = {{-12, 14, 1.22088349258355e+17},
                     {-10, 10, 1.04216468608488e+09},
                     {-10, 12, -8.82666931564652e+15},
                     {-10, 14, 2.59929510849499e+19},
                     {-8, 10, 2.22612779142211e+14},
                     {-8, 12, -8.78473585050085e+17},
                     {-8, 14, -3.14432577551552e+21},
                     {-6, 8, -2.16934916996285e+12},
                     {-6, 12, 1.59079648196849e+20},
                     {-5, 4, -3.39567617303423e+02},
                     {-5, 8, 8.84387651337836e+12},
                     {-5, 12, -8.43405926846418e+20},
                     {-3, 2, 1.14178193518022e+01},
                     {-1, -1, -1.22708229235641e-04},
                     {-1, 1, -1.06201671767107e+02},
                     {-1, 12, 9.03443213959313e+24},
                     {-1, 14, -6.93996270370852e+27},
                     {0, -3, 6.48916718965575e-09},
                     {0, 1, 7.18957567127851e+03},
                     {1, -2, 1.05581745346187e-03},
                     {2, 5, -6.51903203602581e+14},
                     {2, 10, -1.60116813274676e+24},
                     {3, -5, -5.10254294237837e-09},
                     {5, -4, -1.52355388953402e-01},
                     {5, 2, 6.77143292290144e+11},
                     {5, 3, 2.76378438378930e+14},
                     {6, -5, 1.16862983141686e-02},
                     {6, 2, -3.01426947980171e+13},
                     {8, -8, 1.69719813884840e-08},
                     {8, 8, 1.04674840020929e+26},
                     {10, -4, -1.08016904560140e+04},
                     {12, -12, -9.90623601934295e-13},
                     {12, -4, 5.36116483602738e+06},
                     {12, 4, 2.26145963747881e+21},
                     {14, -12, -4.88731565776210e-10},
                     {14, -10, 1.51001548880670e-05},
                     {14, -6, -2.27700464643920e+04},
                     {14, 6, -7.81754507698846e+27}};

    para d =
        {
            .PS = 23.0,
            .TS = 650.0,
            .VS = 0.0026,
            .a = 0.902,
            .b = 0.988,
            .c = 1.0,
            .d = 1.0,
            .e = 1.0,
        };
    return pT2v_sum(p, T, 38, IJn, d);
}

double Vpt_3v(double p, double T)
{
    IJnData IJn[] = {{-10, -8, -4.15652812061591e-55},
                     {-8, -12, 1.77441742924043e-61},
                     {-6, -12, -3.57078668203377e-55},
                     {-6, -3, 3.59252213604114e-26},
                     {-6, 5, -2.59123736380269e+01},
                     {-6, 6, 5.94619766193460e+04},
                     {-6, 8, -6.24184007103158e+10},
                     {-6, 10, 3.13080299915944e+16},
                     {-5, 1, 1.05006446192036e-09},
                     {-5, 2, -1.92824336984852e-06},
                     {-5, 6, 6.54144373749937e+05},
                     {-5, 8, 5.13117462865044e+12},
                     {-5, 10, -6.97595750347391e+18},
                     {-5, 14, -1.03977184454767e+28},
                     {-4, -12, 1.19563135540666e-48},
                     {-4, -10, -4.36677034051655e-42},
                     {-4, -6, 9.26990036530639e-30},
                     {-4, 10, 5.87793105620748e+20},
                     {-3, -3, 2.80375725094731e-18},
                     {-3, 10, -1.92359972440634e+22},
                     {-3, 12, 7.42705723302738e+26},
                     {-2, 2, -5.17429682450605e+01},
                     {-2, 4, 8.20612048645469e+06},
                     {-1, -2, -1.88214882341448e-09},
                     {-1, 0, 1.84587261114837e-02},
                     {0, -2, -1.35830407782663e-06},
                     {0, 6, -7.23681885626348e+16},
                     {0, 10, -2.23449194054124e+26},
                     {1, -12, -1.11526741826431e-35},
                     {1, -10, 2.76032601145151e-29},
                     {3, 3, 1.34856491567853e+14},
                     {4, -6, 6.52440293345860e-10},
                     {4, 3, 5.10655119774360e+16},
                     {4, 10, -4.68138358908732e+31},
                     {5, 2, -7.60667491183279e+15},
                     {8, -12, -4.17247986986821e-19},
                     {10, -2, 3.12545677756104e+13},
                     {12, -3, -1.00375333864186e+14},
                     {14, 1, 2.47761392329058e+26}};

    para d =
        {
            .PS = 23.0,
            .TS = 650.0,
            .VS = 0.0031,
            .a = 0.960,
            .b = 0.995,
            .c = 1.0,
            .d = 1.0,
            .e = 1.0,
        };
    return pT2v_sum(p, T, 39, IJn, d);
}

double Vpt_3w(double p, double T)
{
    IJnData IJn[] = {{-12, 8, -5.86219133817016e-08},
                     {-12, 14, -8.94460355005526e+10},
                     {-10, -1, 5.31168037519774e-31},
                     {-10, 8, 1.09892402329239e-01},
                     {-8, 6, -5.75368389425212e-02},
                     {-8, 8, 2.28276853990249e+04},
                     {-8, 14, -1.58548609655002e+18},
                     {-6, -4, 3.29865748576503e-28},
                     {-6, -3, -6.34987981190669e-25},
                     {-6, 2, 6.15762068640611e-09},
                     {-6, 8, -9.61109240985747e+07},
                     {-5, -10, -4.06274286652625e-45},
                     {-4, -1, -4.71103725498077e-13},
                     {-4, 3, 7.25937724828145e-01},
                     {-3, -10, 1.87768525763682e-39},
                     {-3, 3, -1.03308436323771e+03},
                     {-2, 1, -6.62552816342168e-02},
                     {-2, 2, 5.79514041765710e+02},
                     {-1, -8, 2.37416732616644e-27},
                     {-1, -4, 2.71700235739893e-15},
                     {-1, 1, -9.07886213483600e+01},
                     {0, -12, -1.71242509570207e-37},
                     {0, 1, 1.56792067854621e+02},
                     {1, -1, 9.23261357901470e-01},
                     {2, -1, -5.97865988422577},
                     {2, 2, 3.21988767636389e+06},
                     {3, -12, -3.99441390042203e-30},
                     {3, -5, 4.93429086046981e-08},
                     {5, -10, 8.12036983370565e-20},
                     {5, -8, -2.07610284654137e-12},
                     {5, -6, -3.40821291419719e-07},
                     {8, -12, 5.42000573372233e-18},
                     {8, -10, -8.56711586510214e-13},
                     {10, -12, 2.66170454405981e-14},
                     {10, -8, 8.58133791857099e-06}};
    para d =
        {
            .PS = 23.0,
            .TS = 650.0,
            .VS = 0.0039,
            .a = 0.959,
            .b = 0.995,
            .c = 1.0,
            .d = 1.0,
            .e = 4.0,
        };
    return pT2v_sum(p, T, 35, IJn, d);
}

double Vpt_3x(double p, double T)
{
    IJnData IJn[] = {{-8, 14, 3.77373741298151e+18},
                     {-6, 10, -5.07100883722913e+12},
                     {-5, 10, -1.03363225598860e+15},
                     {-4, 1, 1.84790814320773e-06},
                     {-4, 2, -9.24729378390945e-04},
                     {-4, 14, -4.25999562292738e+23},
                     {-3, -2, -4.62307771873973e-13},
                     {-3, 12, 1.07319065855767e+21},
                     {-1, 5, 6.48662492280682e+10},
                     {0, 0, 2.44200600688281},
                     {0, 4, -8.51535733484258e+09},
                     {0, 10, 1.69894481433592e+21},
                     {1, -10, 2.15780222509020e-27},
                     {1, -1, -3.20850551367334e-01},
                     {2, 6, -3.82642448458610e+16},
                     {3, -12, -2.75386077674421e-29},
                     {3, 0, -5.63199253391666e+05},
                     {3, 8, -3.26068646279314e+20},
                     {4, 3, 3.97949001553184e+13},
                     {5, -6, 1.00824008584757e-07},
                     {5, -2, 1.62234569738433e+04},
                     {5, 1, -4.32355225319745e+10},
                     {6, 1, -5.92874245598610e+11},
                     {8, -6, 1.33061647281106},
                     {8, -3, 1.57338197797544e+06},
                     {8, 1, 2.58189614270853e+13},
                     {8, 8, 2.62413209706358e+24},
                     {10, -8, -9.20011937431142e-02},
                     {12, -10, 2.20213765905426e-03},
                     {12, -8, -1.10433759109547e+01},
                     {12, -5, 8.47004870612087e+06},
                     {12, -4, -5.92910695762536e+08},
                     {14, -12, -1.83027173269660e-05},
                     {14, -10, 1.81339603516302e-01},
                     {14, -8, -1.19228759669889e+03},
                     {14, -6, 4.30867658061468e+06}};
    para d =
        {
            .PS = 23.0,
            .TS = 650.0,
            .VS = 0.0049,
            .a = 0.910,
            .b = 0.988,
            .c = 1.0,
            .d = 1.0,
            .e = 1.0,
        };
    return pT2v_sum(p, T, 36, IJn, d);
}

double Vpt_3y(double p, double T)
{
    IJnData IJn[] = {{0, -3, -5.25597995024633e-10},
                     {0, 1, 5.83441305228407e+03},
                     {0, 5, -1.34778968457925e+16},
                     {0, 8, 1.18973500934212e+25},
                     {1, 8, -1.59096490904708e+26},
                     {2, -4, -3.15839902302021e-07},
                     {2, -1, 4.96212197158239e+02},
                     {2, 4, 3.27777227273171e+18},
                     {2, 5, -5.27114657850696e+21},
                     {3, -8, 2.10017506281863e-17},
                     {3, 4, 7.05106224399834e+20},
                     {3, 8, -2.66713136106469e+30},
                     {4, -6, -1.45370512554562e-08},
                     {4, 6, 1.49333917053130e+27},
                     {5, -2, -1.49795620287641e+07},
                     {5, 1, -3.81881906271100e+15},
                     {8, -8, 7.24660165585797e-05},
                     {8, -2, -9.37808169550193e+13},
                     {10, -5, 5.14411468376383e+09},
                     {12, -8, -8.28198594040141e+04}};

    para d =
        {
            .PS = 22.0,
            .TS = 650.0,
            .VS = 0.0031,
            .a = 0.996,
            .b = 0.994,
            .c = 1.0,
            .d = 1.0,
            .e = 4.0,
        };
    return pT2v_sum(p, T, 20, IJn, d);
}

double Vpt_3z(double p, double T)
{
    IJnData IJn[] = {{-8, 3, 2.4400789229065e-11},
                     {-6, 6, -4.6305743033124e+06},
                     {-5, 6, 7.2880327477771e+09},
                     {-5, 8, 3.2777630285886e+15},
                     {-4, 5, -1.1059817011841e+09},
                     {-4, 6, -3.2389991572996e+12},
                     {-4, 8, 9.2381400702325e+15},
                     {-3, -2, 8.4225008041371e-13},
                     {-3, 5, 6.6322143624551e+11},
                     {-3, 6, -1.6717018667214e+14},
                     {-2, 2, 2.5374935870139e+03},
                     {-1, -6, -8.1973155961052e-21},
                     {0, 3, 3.2838058789066e+11},
                     {1, 1, -6.2500479117154e+07},
                     {2, 6, 8.0319795746202e+20},
                     {3, -6, -2.0439701133835e-11},
                     {3, -2, -3.7839104705594e+03},
                     {6, -6, 9.7287654593862e-03},
                     {6, -5, 1.5435572168146e+01},
                     {6, -4, -3.7396286292864e+03},
                     {6, -1, -6.8285901137457e+10},
                     {8, -8, -2.4848801561454e-04},
                     {8, -4, 3.9453604949707e+06}};
    para d =
        {
            .PS = 22.0,
            .TS = 650.0,
            .VS = 0.0038,
            .a = 0.993,
            .b = 0.994,
            .c = 1.0,
            .d = 1.0,
            .e = 4.0,
        };
    return pT2v_sum(p, T, 23, IJn, d);
}

double Vpt_3subreg(double p, double t, char subRegion)
{
    double v;
    switch (subRegion)
    {
    case 'a':
        v = Vpt_3a(p, t);
        break;
    case 'b':
        v = Vpt_3b(p, t);
        break;
    case 'c':
        v = Vpt_3c(p, t);
        break;
    case 'd':
        v = Vpt_3d(p, t);
        break;
    case 'e':
        v = Vpt_3e(p, t);
        break;
    case 'f':
        v = Vpt_3f(p, t);
        break;
    case 'g':
        v = Vpt_3g(p, t);
        break;
    case 'h':
        v = Vpt_3h(p, t);
        break;
    case 'i':
        v = Vpt_3i(p, t);
        break;
    case 'j':
        v = Vpt_3j(p, t);
        break;
    case 'k':
        v = Vpt_3k(p, t);
        break;
    case 'l':
        v = Vpt_3l(p, t);
        break;
    case 'm':
        v = Vpt_3m(p, t);
        break;
    case 'n':
        v = Vpt_3n(p, t);
        break;
    case 'o':
        v = Vpt_3o(p, t);
        break;
    case 'p':
        v = Vpt_3p(p, t);
        break;
    case 'q':
        v = Vpt_3q(p, t);
        break;
    case 'r':
        v = Vpt_3r(p, t);
        break;
    case 's':
        v = Vpt_3s(p, t);
        break;
    case 't':
        v = Vpt_3t(p, t);
        break;
    case 'u':
        v = Vpt_3u(p, t);
        break;
    case 'v':
        v = Vpt_3v(p, t);
        break;
    case 'w':
        v = Vpt_3w(p, t);
        break;
    case 'x':
        v = Vpt_3x(p, t);
        break;
    case 'y':
        v = Vpt_3y(p, t);
        break;
    case 'z':
        v = Vpt_3z(p, t);
        break;
    default:
        v = -1000.0;
        break;
    }
    return v;
}