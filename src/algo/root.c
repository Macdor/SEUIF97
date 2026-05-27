//-----------------------------------------------------------------------------
// Root Finding Methods For IAPWS - IF97
// Ref: Numerical Reciples - Chapter 9 :Pages:347 ~ 368
// Last updated: 2004.01.03 By Maohua Cheng
//-----------------------------------------------------------------------------

#include <math.h>
#include <stdio.h>
#include "algorithm.h"

#define EPS 3.0e-8
//----------------------------------------------------------------------------
// SECANT METHOD : Ch.9.2: Pages 357, 
// Using the secant method, find the root of a func throught lie between x1 and
// x2. The root returned as rtsec, is refined until its accuracy is ABS(xacc)
//----------------------------------------------------------------------------

double rtsec2(callfunc func, double cVar1, double fr, double x1, double x2, double fl,double f, double xacc, int iMAX)
{
  double xl,rts,swap,dx;
  // pick the bound with the smaller function value as the most recent guess
  if (fabs(fl) < fabs(f)) {
        rts=x1;
        xl=x2;
        swap=fl;
        fl=f;
        f=swap;
  } else {
        xl=x1;
        rts=x2;
  }
  //scant loop
  int i=0;
  if ((f-fl)!=0.0)
  {
      do
      {
          dx=(xl-rts)*f/(f-fl); // increment with respect to latest value
          xl=rts;
          fl=f;
          rts += dx;
          // TODO： 可将解的上下限作为参数带进来，保证迭代过程的解不超上下限
          //rts must bounded in region X
          f=fr-(*func)(cVar1,rts);
          i++;
      }
    //  while ((fabs(dx) > xacc)&&(i<iMAX)); // Covergence
    //  while ((fabs(dx) > xacc)&&(i<iMAX)&&(f!=0.0)); // Covergence
      while ( fabs(dx) > xacc && i<iMAX &&
              f!=0.0 && (f-fl)!=0.0 ); // Covergence

     // printf("rtsec2 i= %d \n",i);
  }
  
  return rts;
}


double rtsec1(callfunc func, double cVar2, double fr, double x1, double x2, double fl, double f, double xacc,int iMAX)
{
  double xl,rts,swap,dx;
  // pick the bound with the smaller function value as the most recent guess
  if (fabs(fl) < fabs(f)) { 
      rts=x1;
      xl=x2;
      swap=fl;
      fl=f;
      f=swap;
  }  else {
      xl=x1;
      rts=x2;
  }
  //scant loop
  int i=0;
  if ((f-fl)!=0.0)
  {
     do
     {
        dx=(xl-rts)*f/(f-fl); // increment with respect to latest value
        xl=rts;
        fl=f;
        rts += dx;
        // TODO： 可将解的上下限作为参数带进来，保证迭代过程的解不超上下限
        //rts may out-bounded in region X
        if (rts<=0) rts=0.000001;
        //if (rts>100) rts=100;
      
        f=fr-(*func)(rts,cVar2);
        i++;
     
     }
    //  while ((fabs(dx) > xacc)&&(i<iMAX)); // Covergence
     while ((fabs(dx) > xacc)&&(i<iMAX)&&(f!=0.0)&&((f-fl)!=0.0)); // Covergence
    
     //printf("rtsec1 i= %d \n",i);
  
  };
  
  return (rts);
}

#undef EPS


