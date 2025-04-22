/*
Add  SEUIF97 the paths of lib and header file to  the  environment variables of GCC/MinGW-W64  :
    1 Windows: C:/Windows/system/libseuif97.dll
    2 Linux: ./usr/lib/libseuif97.so
    3 ./include/seuif97.h

Build:
Windows with MinGW-W64：
   gcc -o demo.exe demo.c -I./include  -LC:/Windows/system/  -lseuif97

Linux：
   gcc -o demo demo.c  -I./include  -L/usr/lib/ -lseuif97 -lm

Run:
./demo

Windows  MSVC: seuif97.dll
   cl /Fedemo.exe /Fo./obj/demo.obj demo.c -I./include/  ./lib/seuif97.lib

Author: Cheng Maohua
*/

#include <stdio.h>
#include <stdlib.h>
#include "../../src/common/seuif97.h"

int main(void)
{

    double in1 = 0.1;
    double in2 = 200;
    double out;

    out = ph(in1, in2, 31);
    printf("(in1,in2)(%.5f,%.5f) r= %.10f", in1, in2, out);
    return EXIT_SUCCESS;
}
