/*
1. Windows with MinGW-W64：libseuif97.dll
   Add SEUIF97 the paths of lib and header file to  the  environment variables of GCC/MinGW-W64  :
       C:/Windows/system/libseuif97.dll
       ./include/seuif97.h
   Build:
      gcc -o demo.exe demo.c -I./include  -LC:/Windows/system/  -lseuif97

2. Linux: libseuif97.so
   Add SEUIF97 the paths of lib and header file to  the  environment variables of GCC/MinGW-W64  :
      ./usr/lib/libseuif97.so
      ./include/seuif97.h
   Build:
     gcc -o demo demo.c  -I./include  -L/usr/lib/ -lseuif97 -lm

 3. Windows with MSVC using libseuif97.dll
      C:/Windows/system/libseuif97.dll
      64:  ../../shared_lib/Windows/x64/libseuif97.lib
      32:  ../../shared_lib/Windows/x86/libseuif97.lib
      ./include/seuif97.h

   Build:
   64bits：x64 Native Tools Command Prompt for VS
      cl /Fedemo.exe demo.c -I./include/ ../../shared_lib/Windows/x64/libseuif97.lib
   32 bits:Developer Command Prompt for VS
      cl /Fedemo.exe demo.c -I./include/ ../../shared_lib/Windows/x86/libseuif97.lib

4. Windows with MSVC: seuif97.dll
   Using cmake to to build the seuif97.dll with msvc
      C:\Windows\System32\seuif97.dll //64bits under x64 Native Tools Command Prompt for VS
      C:\Windows\SysWOW64\seuif97.dll //32bits under Developer Command Prompt for VS
      lib:  ../../build/seuif97.lib
      ./include/seuif97.h

   Build:
     64 bits under x64 Native Tools Command Prompt for VS 
        cl /Fedemo.exe demo.c -I./include/ ../../build/seuif97.lib
     32 bits under Developer Command Prompt for VS
        cl /Fedemo.exe demo.c -I./include/ ../../buildc/seuif97.lib
     
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
