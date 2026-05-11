# Using SEUIF97 with C/C++ 

## SEUIF97 

1. The dynamic library with MinGW-W64/GCC: `libseuif97`
   * Windows `C:/Windows/system/libseuif97.dll`
   * Linux `./usr/lib/libseuif97.so`

2. Windows with MSVC using `libseuif97.dll`
   * DLL: `C:/Windows/system/libseuif97.dll`
   * Lib:  
      *  64:  `../../shared_lib/Windows/x64/libseuif97.lib`
      *  32:  `../../shared_lib/Windows/x86/libseuif97.lib`

3. The dynamic library with MSVC: `seuif97.dll`
   Using `cmake` to build the  `seuif97.dll` with MSVC
   * DLL: copy `./build/seuif97.dll` to
       * `./` or `C:\Windows\SysWOW64\seuif97.dll`, (32bits  build under Developer Command Prompt for VS)
       * `./` or  `C:\Windows\System32\seuif97.dll`, (64bits  build under **x64** Native Tools Command Prompt for VS)
   * Lib: `../../build/seuif97.lib`
 
4. header file: `./include/seuif97.h`

## Examples

* [demo.c](./demo.c)

* [H-S(Mollier) Diagram of Steam Turbine Expansion](./Turbine_H-S.cpp)
