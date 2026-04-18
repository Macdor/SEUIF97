# SEUIF97

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8242452.svg)](https://doi.org/10.5281/zenodo.8242452)

This is the C implementation of the high-speed IAPWS-IF97 package **SEUIF97**. It is suitable for computation-intensive calculations，such as the simulation of non-stationary processes, on-line process monitoring and optimizations.
 
Through the high-speed library, IAPWS-IF97 results are generated with high accuracy and at speeds significantly faster than the **repeated squaring method** or the standard C library's `math.pow()`

**SEUIF97** is also significantly faster than approximate equations and fast-calculation algorithms for water and steam properties.

## Key Acceleration Methods

* Shortest Addition Chain Algorithm: Used for the rapid computation of integer powers via optimal multiplication sequences.
* Power Caching for Multi-Polynomial Evaluation: Precomputes powers to avoid redundant calculations, directly retrieving them to boost performance across multiple polynomials.
* Recurrence Method for Multi-Polynomial Evaluation：Utilizes the derivative relationship to compute a base polynomial, then derives the rest via base scaling  to maximize efficiency

## The shared library

**Building the shared library**

You can build the library using either `make` or `cmake`:

* make

```bash
make
```

* cmake

```bash
cmake -B ./build/   
cmake --build ./build/ --config Release
```

**The compiled shared libraries**

The compiled shared libraries using GCC are available in the [/shared_lib/](./shared_lib/) directory:

* [Windows(x86/64)](./shared_lib/Windows): `libseuif97.dll` 

* [Linux(x64)](./shared_lib/Linux): `libseuif97.so`

**Interfaces and examples**

Interfaces and examples are provided in the  [/demo/](./demo) directory, supporting a wide range of languages and environments

*  C/C++, Python, C#, Java, Excel VBA, MATLAB, Rust, Fortran, Pascal, Golang, Modelica

## Functions of the SEUIF97 Shared Library

The SEUIF97 library provides comprehensive functions for calculating water and steam properties, as well as the thermodynamic processes of steam turbines.

**Water and Steam Properties**

Using SEUIF97, you can define the state of steam using various pairs of known properties to determine any desired output property, including in the [30 properties](#properties-in-libseuif97).

The following 12 input pairs are implemented: 

```
(p,t) (p,h) (p,s) (p,v) 

(t,h) (t,s) (t,v) 

(p,x) (t,x) (h,x) (s,x) 

(h,s)  
```

The type of property functions are provided in the package

```c 
  ??(in1,in2,o_id)
```

* the first,second input parameters : the input property pairs
* the third input parameters: the ID of the calculated property - [o_id](#properties)
* the return: the calculated property value of `o_id`

**Thermodynamic Process of Steam Turbine**
   
* Isentropic Enthalpy Drop：ishd(pi,ti,pe)

  ```txt
    pi - inlet pressure(MPa); ti -inlet temperature(°C)
    pe - outlet pressure(MPa)
  ```
*  Isentropic Efficiency(0~100)： ief(pi,ti,pe,te) (superheated steam zone)
   ```txt
       pi - inlet pressure(MPa);  ti - inlet temperature(°C)
       pe - outlet pressure(MPa); te - outlet temperature(°C)
   ```

**The Function Prototype in C**

```c
// Functions of Properties
double pt(double p, double t, int o_id);
double ph(double p, double h, int o_id);
double ps(double p, double s, int o_id);
double pv(double p, double v, int o_id);

double th(double t, double h, int o_id);
double ts(double t, double s, int o_id);
double tv(double t, double v,  int o_id);

double hs(double h, double s,  int o_id);

double px(double p, double x,  int o_id);
double tx(double t, double x,  int o_id);
double hx(double h, double x,  int o_id);
double sx(double s, double x,  int o_id);

//The Functions for Thermodynamic Process of Steam Turbine
double ishd(double pi, double ti, double pe);
double ief(double pi, double ti, double pe, double te);
```

## Using SEUIF97 library 

To use the shared library in your preferred programming language, follow these steps:

1. Put the shared library in the default `Lib` path of OS or the programming language
   * **Windows(x86/64)** 
     * copy `libseuif97.dll` in the [Windows/x86](./shared_lib/Windows/x86) or [Windows/x64](./shared_lib/Windows/x64) folder to a default path of Windows32/64's DLL:  `C:\Windows\System`
   
   * **Linux(x64)** 
     * copy `libseuif97.so` in the [Linux/x64](./shared_lib/Linux/x64) folder to a default path of Linux shared lib : `/usr/lib`/
   
2. Add the specific API file for your programming language to its corresponding`API` path.

  * **Note**: API paths vary by language. Please refer to the specific  [Examples](./demo/) provided for your programming language to locate the correct path and setup instructions.

## Interfaces and Examples

* [./demo/](./demo)

|  Language                           |    Interface                                     | 
|:-----------------------------------:|:-------------------------------------------------|
| [C/C++](./demo/demo-c)              |  [seuif97.h](./demo/demo-c/include/seuif97.h)    |
| [Python](./demo/demo-python)        |  [seuif97.py](./demo/demo-python/seuif97.py)     |
| [C#](./demo/demo-csharp)            |  [seuif97.cs](./demo/demo-csharp/seuif97.cs)     |
| [Excel VBA](./demo/ExcelVBA)        |  [seuif97.bas](./demo/ExcelVBA/seuif97.bas)      |
| [Java](./demo/demo-java)            |  [seuif97.java](./demo/demo-java/seuif97.java)   |
| [MATLAB64](./demo/MATLAB64)         |  [seuif97.m](./demo/MATLAB64/seuif97/seuif97.m)  |
| [Rust](./demo/demo-rust)            |  [seuif97.rs](./demo/demo-rust/src/seuif97.rs)   |
| [Fortran](./demo/demo-Fortran)      |  [seuif97.f08](./demo/demo-fortran/seuif97.f08)  |
| [Pascal](./demo/demo-pascal)        |  [seuif97.pas](./demo/demo-pascal/seuif97.pas)   |
| [Modelica](./demo/demo-modelica)    |  [seuif97.mo](./demo/demo-modelica/demomodelica/seuif97.mo) |
| [Golang](./demo-go)                |  **Example** [demo.go](./demo/demo-go/demo.go) |

You can modify the provided interfaces to match your own API needs.

**The Selected Examples**

**Python**

* [T-S Diagram](./demo/demo-python/Diagram_T-S.py)

* [H-S Diagram](./demo/demo-python/Diagram_H-S.py)

* [H-S Diagram of Steam Turbine Expansion](./demo/demo-python/Turbine_H-S.py)

**C++**

* [H-S Diagram of Steam Turbine Expansion](./demo/demo-c/Turbine_H-S.cpp)


## Properties 

| Property                              |    Unit     | Symbol | o_id       |
| ------------------------------------- | :---------: |:-----: |:---------: |
| Pressure                              |     MPa     |      p |          0 |
| Temperature                           |     °C      |      t |          1 |
| Density                               |   kg/m^3    |      ρ |          2 |
| Specific Volume                       |   m^3/kg    |      v |          3 |
| Specific enthalpy                     |    kJ/kg    |      h |          4 |
| Specific entropy                      |  kJ/(kg·K)  |      s |          5 |
| Specific exergy                       |    kJ/kg    |      e |          6 |
| Specific internal energy              |    kJ/kg    |      u |          7 |
| Specific isobaric heat capacity       |  kJ/(kg·K)  |     cp |          8 |
| Specific isochoric heat capacity      |  kJ/(kg·K)  |     cv |          9 |
| Speed of sound                        |     m/s     |      w |         10 |
| Isentropic exponent                   |             |      k |         11 |
| Specific Helmholtz free energy        |    kJ/kg    |      f |         12 |
| Specific Gibbs free energy            |    kJ/kg    |      g |         13 |
| Compressibility factor                |             |      z |         14 |
| Steam quality                         |             |      x |         15 |
| Region                                |             |      r |         16 |
| Isobaric cubic expansion coefficient  |     1/K     |    ɑv  |         17 |
| Isothermal compressibility            |    1/MPa    |     kT |         18 |
| Partial derivative (∂V/∂T)p           |  m3/(kg·K)  |(∂V/∂T)p|         19 |
| Partial derivative (∂V/∂P)T           | m3/(kg·MPa) |(∂V/∂P)T|         20 |
| Partial derivative (∂p/∂t)v           |    MPa/K    |(∂p/∂t)v|         21 |
| Isothermal throttling coefficient     | kJ/(kg·MPa) |    δt  |         22 |
| Joule-Thomson coefficient             |    K/MPa    |     μ  |         23 |
| Dynamic viscosity                     |  kg/(m·s)   |     η  |         24 |
| Kinematic viscosity                   |    m^2/s    |     ν  |         25 |
| Thermal conductivity                  |   W/(m.K)   |     λ  |         26 |
| Thermal diffusivity                   |   um^2/s    |     a  |         27 |
| Prandtl number                        |             |     Pr |         28 |
| Surface tension                       |     N/m     |     σ  |         29 |

## Publications

* 王培红,贾俊颖,程懋华. 水和水蒸汽热力性质IAPWS-IF97公式的通用计算模型[J]. 动力工程,2001,21(6)：1564-1567 [[ pdf ]](./doc/水和水蒸汽热力性质IAPWS-IF97公式的通用计算模型.pdf)

* 芮嘉敏,孙振业,程懋华. 基于最短加法链状态空间树的IAPWS-IF97快速计算方法[J]. 汽轮机技术,2017,59(4):245-247 [[ pdf ]](./doc/基于最短加法链状态空间树的IAPWS-IF97快速计算方法.pdf)
 

## Cite as

* Cheng Maohua. (2023). The C implementation of the high-speed IAPWS-IF97 package: SEUIF97 (1.2.0). Zenodo. https://doi.org/10.5281/zenodo.8242452