/*
The header file of SEUIF97

   Author:   Cheng Maohua
   Email:    cmh@seu.edu.cn
*/

#ifndef SEUIF97_H

#define SEUIF97_H

#ifdef __cplusplus
extern "C"
{
#endif

// Кроссплатформенные макросы для экспорта/импорта
#if defined(_WIN32) || defined(_WIN64)
    #ifdef LIB_EXPORT
        #define IF97_API __declspec(dllexport)
    #elif defined(LIB_IMPORT)
        #define IF97_API __declspec(dllimport)
    #else
        #define IF97_API
    #endif
#elif defined(__linux__) || defined(__ANDROID__)
    // Для Linux/Android используем GCC атрибуты
    #ifdef LIB_EXPORT
        #define IF97_API __attribute__((visibility("default")))
    #else
        #define IF97_API
    #endif
#else
    #define IF97_API
#endif

    IF97_API double pt(double p, double t, int o_id);
    IF97_API double ph(double p, double h, int o_id);
    IF97_API double ps(double p, double s, int o_id);
    IF97_API double pv(double p, double v, int o_id);
    IF97_API double th(double t, double h, int o_id);

    IF97_API double ts(double t, double s, int o_id);
    IF97_API double tv(double t, double v, int o_id);
    IF97_API double hs(double h, double s, int o_id);
    IF97_API double px(double p, double x, int o_id);
    IF97_API double tx(double t, double x, int o_id);
    //
    IF97_API double hx(double h, double x, int o_id);
    IF97_API double sx(double s, double x, int o_id);

    //  Functions of Thermodynamic Process
    IF97_API double ishd(double pi, double ti, double pe);
    IF97_API double ief(double pi, double ti, double pe, double te);
 
#ifdef __cplusplus
}
#endif

#endif
