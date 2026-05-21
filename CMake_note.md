
## CMake

## Linux

```bash
cmake -B ./build/   
cmake --build ./build/ --config Release
```

##  Windows: MinGW 

```bash
cmake -B ./build -G "MinGW Makefiles"  -DCMAKE_BUILD_TYPE=Release
cmake --build ./build/ --config Release
```
##  Windows: Microsoft Visual Studio c++ 

For 64-bit __stdcall DLL 

```bash
cmake -B ./build/  -DBUILD_DLL=ON -DSTDCALL_LIBRARY=ON-G "Visual Studio 18 2026" -A x64
cmake --build ./build/ --config Release
```

For 64-bit __cdecl DLL:
```bash
cmake -B ./build/  -DBUILD_DLL=ON -DCDECL_LIBRARY=ON -G "Visual Studio 18 2026" -A x64
cmake --build ./build/ --config Release
```

For 32-bit __stdcall DLL:
```bash
cmake -B ./build/  -DBUILD_DLL=ON -DSTDCALL_LIBRARY=ON -G "Visual Studio 18 2026" -A Win32
cmake --build ./build/ --config Release
```

For 32-bit __cdecl DLL:
```bash
cmake -B ./build/  -DBUILD_DLL=ON -DCDECL_LIBRARY=ON -G "Visual Studio 18 2026" -A Win32
cmake --build ./build/ --config Release
```