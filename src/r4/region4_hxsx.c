/*
 Region 4 - The extended input pairs: (h,x),(s,x)
 x: Steam quality
*/
#include <math.h>
#include <stdio.h>
#include "region4.h"

// function for getting the steam quality,residuals: x(T,y)-x
double Ty2x_residuals(double T, double x, double y, int y_id)
{
    double sw = T2SatWater(T, y_id);
    double ss = T2SatSteam(T, y_id);
    return (y - sw) / (ss - sw) - x;
}

/// Bisection for the root : Ty2x_residuals(T,x, y, y_id)=0
double bisection_reg4(double x, double y, int y_id, double Tl, double Tr, double tol, int maxiter)
{
    double T = 0.0;
    double fl = Ty2x_residuals(Tl, x, y, y_id); // residual for left  bound
    double fr = Ty2x_residuals(Tr, x, y, y_id); // resdiual for right bound
    double f = 0.0;
    int numIters = 0;

    for (int i = 0; i < maxiter; i++)
    {
        numIters += 1;
        // get midpoint
        T = 0.5 * (Tl + Tr);
        // evaluate resdiual at midpoint
        f = Ty2x_residuals(T, x, y, y_id);
        //  check for convergence
        if (fabs(f) < tol)
        {
            break;
        };

        // reset the bounds
        if (f * fl < 0.0)
        {
            // move right bound info to mid
            Tr = T;
            fr = f;
        }
        else
        {
            // move left bound info to mid
            Tl = T;
            fl = f;
        }
    };
    return T;
}
