/*
 The APT of region 2

*/
#include "region2.h"

double pT_reg2(double p, double T, int o_id)
// o_id: output propertry
{
    double value = 0.0;
    switch (o_id)
    {
    case OT:
        value = T;
        break;
    case OP:
        value = p;
        break;
    case OX:
        value = 1.0;
        break;
    case OR:
        value = 2.0;
        break;
    default:
        value = pT_reg(p, T, o_id, pT_thermal_reg2, pT_ext_reg2);
        break;
    }
    return value;
}

double ph_reg2(double p, double h, int o_id)
// o_id: output propertry
{
    if (o_id == OP)
        return p;
    if (o_id == OH)
        return h;
    if (o_id == OX)
        return 1.0;
    if (o_id == OR)
        return 2.0;

    double T;
    T = ph2T_reg2(p, h);
    return pT_reg2(p, T, o_id);
}

double ps_reg2(double p, double s, int o_id)
// o_id: output propertry
{
    if (o_id == OP)
        return p;
    if (o_id == OS)
        return s;
    if (o_id == OX)
        return 1.0;
    if (o_id == OR)
        return 2.0;

    double T;
    T = ps2T_reg2(p, s);
    return pT_reg2(p, T, o_id);
}

double hs_reg2(double h, double s, int o_id)
// o_id: output propertry
{
    if (o_id == OH)
        return h;
    if (o_id == OS)
        return s;
    if (o_id == OX)
        return 1.0;
    if (o_id == OR)
        return 2.0;

    double p;
    p = hs2p_reg2(h, s);
    return ph_reg2(p, h, o_id);
}

//----------------------------------------
// the extend input pairs (p,v) T
//------------------------------------------------------------------
double pv_reg2(double p, double v, int o_id)
// o_id: output propertry
{
    if (o_id == OP)
        return p;
    if (o_id == OV)
        return v;
    if (o_id == OX)
        return 1.0;
    if (o_id == OR)
        return 2.0;

    double T;
    T = pv2T_reg2(p, v);
    return pT_reg2(p, T, o_id);
}

double Tv_reg2(double T, double v, int o_id)
// o_id: output propertry
{
    if (o_id == OT)
        return T;
    if (o_id == OV)
        return v;
    if (o_id == OD)
        return 1/v;
    if (o_id == OX)
        return 1.0;
    if (o_id == OR)
        return 2.0;

    double p;
    p = Tv2p_reg2(T, v);
    return pT_reg2(p, T, o_id);
}

// (T,s)
double Ts_reg2(double T, double s, int o_id)
// o_id: output propertry
{
    if (o_id == OT)
        return T;
    if (o_id == OS)
        return s;
    if (o_id == OX)
        return 1.0;
    if (o_id == OR)
        return 2.0;

    double p;
    p = Ts2p_reg2(T, s);
    return pT_reg2(p, T, o_id);
}

// (T,h)
double Th_reg2(double T, double h, int o_id)
// o_id: output propertry
{
    if (o_id == OT)
        return T;
    if (o_id == OH)
        return h;
    if (o_id == OX)
        return 1.0;
    if (o_id == OR)
        return 2.0;

    double p = Th2p_reg2(T, h);
    if (o_id == OP)
        return p;

    return pT_reg2(p, T, o_id);
}

double pT_thermal_reg2(double p, double T, int o_id)
// o_id: output propertry
{
    double value;
    switch (o_id)
    {
    case OT:
        value = T;
        break;
    case OP:
        value = p;
        break;
    case OX:
        value = 1.0;
        break;
    case OR:
        value = 2.0;
        break;
    case OV:
        value = pT2v_reg2(p, T);
        break;
    case OD:
        value = 1.0 / pT2v_reg2(p, T);
        break;
    case OH:
        value = pT2h_reg2(p, T);
        break;
    case OS:
        value = pT2s_reg2(p, T);
        break;
    case OU:
        value = pT2u_reg2(p, T);
        break;
    case OCV:
        value = pT2cv_reg2(p, T);
        break;
    case OCP:
        value = pT2cp_reg2(p, T);
        break;
    case OW:
        value = pT2w_reg2(p, T);
        break;
    default:
        value = INVALID_OUTID;
        break;
    }
    return value;
}

// Region1: the extended properties
double pT_ext_reg2(double p, double T, int o_id)
{
    double r = 0.0;
    switch (o_id)
    {
    case OT:
        r = T;
        break;
    case OP:
        r = p;
        break;
    case OX:
        r = 1.0;
        break;
    case OR:
        r = 2.0;
        break;
    case OKISE:
        r = pT2k_reg2(p, T);
        break;
    case OIPCEC:
        r = pT2ipcec_reg2(p, T);
        break;
    case OKT:
        r = pT2kt_reg2(p, T);
        break;
    case OG:
        r = pT2g_reg2(p, T);
        break;
    case OF:
        r = pT2f_reg2(p, T);
        break;
    case OJTC:
        r = pT2joule_reg2(p, T);
        break;
        break;
    case OIJTC:
        r = pT2ijoule_reg2(p, T);
        break;
    case OZ:
        r = pT2z_reg2(p, T);
        break;
    case OE:
        r = pT2e_reg2(p, T);
        break;
    case ODPDT:
        r = pT2dpdtcv_reg2(p, T);
        break;
    case ODVDP:
        r = pT2dvdpct_reg2(p, T);
        break;
    case ODVDT:
        r = pT2dvdtcp_reg2(p, T);
        break;
    case ODVDHCP:
    {
        double dvdtcp = pT2dvdtcp_reg2(p, T);
        double Cp = pT2cp_reg2(p, T);
        r = dvdtcp/Cp;
        break;
    }
    case ODVDPCH:
    {
        double dvdpct = pT2dvdpct_reg2(p, T);
        double dvdtcp = pT2dvdtcp_reg2(p, T);
        double dhdpct = pT2ijoule_reg2(p, T);
        double Cp = pT2cp_reg2(p, T);
        r = dvdpct - dvdtcp*dhdpct/Cp;
        break;
    }
    /*
       OPC => pT2pc_reg2(p, T),
       OBETAP => pT2batap_reg2(p, T),
       OFI => pT2fi_reg2(p, T),
       OFU => pT2fu_reg2(p, T),
       OALFAP => pT2alfap_reg2(p, T),*/
    default:
        r = INVALID_OUTID;
    }
    return r;
}