/*
1. Windows with MinGW-W64：libseuif97.dll
   Add SEUIF97 the paths of lib and header file to  the  environment variables of GCC/MinGW-W64  :
       C:/Windows/system/libseuif97.dll
       ./include/seuif97.h

   Build:
      gcc -o demo.exe demo.c -I./include  -LC:/Windows/system/  -lseuif97

2. Linux: ibseuif97.so
   Add SEUIF97 the paths of lib and header file to  the  environment variables of GCC/MinGW-W64  :
      ./usr/lib/libseuif97.so
      ./include/seuif97.h
   Build:
     gcc -o demo demo.c  -I./include  -L/usr/lib/ -lseuif97 -lm

3. Windows with MSVC: seuif97.dll
   Using cmake to to build the  seuif97.dll with MSVC under Developer Command Prompt for VS
      ./seuif97.dll
      ./lib/seuif97.lib
      ./include/seuif97.h

   Build:
      cl /Fedemo.exe /Fo./obj/demo.obj demo.c -I./include/ ./lib/seuif97.lib

Author: Cheng Maohua
*/

#include <stdio.h>
#include <stdlib.h>
#include "seuif97.h"

int main(void)
{

    double p = 16.13;
    double t = 535;
    double h, s, v;

    h = pt(p, t, 4);
    s = pt(p, t, 5);
    v = pt(p, t, 3);
    printf("(p,t)(%.2f,%.2f) h= %.2f, s= %.4f, v= %.4f\n", p, t, h, s, v);
    return EXIT_SUCCESS;
}
