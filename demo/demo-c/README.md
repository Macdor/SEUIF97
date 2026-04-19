# Using SEUIF97 with C/C++ 

## SEUIF97 

1. The dynamic library with MinGW-W64/GCC: libseuif97.dll
   * Windows `C:/Windows/system/libseuif97.dll`
   * Linux `./usr/lib/libseuif97.so`

2. The dynamic library with MSVC: seuif97.dll
   *  Using cmake to to build the  seuif97.dll with MSVC
      * ./seuif97.dll
      * ./lib/seuif97.lib

3. header file: `./include/seuif97.h`

## Examples

* [demo.c](./demo.c)

* [H-S(Mollier) Diagram of Steam Turbine Expansion](./Turbine_H-S.cpp)
