
#  Building the Shared Library

## Make

### GCC at Linux

```base
make
```

### MinGW-GCC at Windows

64bit __stdcall  dll(VBA)

```makefile
#CFLAGS=-O3 -DBUILD_DLL 
```

32bit  __stdcall DLL(VBA)

```makefile
#CFLAGS=-O3 -DBUILD_DLL -m32
```

## CMake

### GCC at Linux

```bash
cmake -B ./build/   
cmake --build ./build/ --config Release
```

###  MinGW-GCC at Windows 

For 64-bit __stdcall DLL 
```bash
cmake -B ./build -G "MinGW Makefiles"  -DCMAKE_BUILD_TYPE=Release
cmake --build ./build/ --config Release
```
###  Windows: Microsoft Visual Studio c++ 

For 64-bit __stdcall DLL 

```bash
cmake -B ./build/  -DBUILD_DLL=ON -G "Visual Studio 18 2026" -A x64
cmake --build ./build/ --config Release
```

For 32-bit __stdcall DLL:
```bash
cmake -B ./build/  -DBUILD_DLL=ON -G "Visual Studio 18 2026" -A Win32
cmake --build ./build/ --config Release
```
