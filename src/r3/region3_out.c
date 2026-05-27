/*
 Tha API of region3

*/
#include "region3.h"

double Td_reg3(double T, double d, int o_id)
// o_id: output propertry
{
    double value;
    switch (o_id)
    {
    case OX:
        value = 1.0;
        break;
    case OR:
        value = 3.0;
        break;
    case OT:
        value = T;
        break;
    case OD:
        value = d;
        break;
    case OV:
        value = 1.0 / d;
        break;
    case OP:
    case OH:
    case OS:
    case OU:
    case OCV:
    case OCP:
    case OW:
        value = Td_thermal_reg3(T, d, o_id);
        break;
    case ODV:
    case OKV:
    case OTC:
    case OTD:
    case OPR:
    case OST:
        value = Td_transport_reg3(T, d, o_id);
        break;
    default:
        value = Td_ext_reg3(T, d, o_id);
        break;
    }
    return value;
}

double Td_transport_reg3(double T, double d, int o_id)
{
    double value;
    double cp;
    double tc;
    double dv;
    switch (o_id)
    {
    case OX:
        value = 1.0;
        break;
    case OR:
        value = 3.0;
        break;
    case OT:
        value = T;
        break;
    case OD:
        value = d;
        break;
    case OV:
        value = 1.0 / d;
        break;
    case ODV:
    {
        value = viscosity(d, T);
        break;
    }
    case OKV:
    {
        value = kinematic_viscosity(d, T);
        break;
    }
    case OTC:
    {
        value = thCond(d, T);
        break;
    }
    case OTD:
    {
        cp = Td2cp_reg3(T, d);
        tc = thCond(d, T);
        value = thermal_diffusivity(tc, cp, d);
        break;
    }
    case OPR:
    {
        dv = viscosity(d, T);
        cp = Td2cp_reg3(T, d);
        tc = thCond(d, T);
        value = prandtl_number(dv, cp, tc);
        break;
    }
    case OST:
    {
        value = tension(T);
        break;
    }
    default:
    {
        value = INVALID_OUTID;
        break;
    }
    }
    return value;
}

double pT_reg3(double p, double T, int o_id)
// o_id: output propertry
{
    if (o_id == OP)
        return p;
    if (o_id == OT)
        return T;
    if (o_id == OX)
        return 1.0;
    if (o_id == OR)
        return 3.0;
    double d;
    d = 1.0 / pT2v_reg3(p, T);
    return Td_reg3(T, d, o_id);
}

double ph_reg3(double p, double h, int o_id)
// o_id: output propertry
{
    if (o_id == OP)
        return p;
    if (o_id == OH)
        return h;
    if (o_id == OX)
        return 1.0;
    if (o_id == OR)
        return 3.0;

    double d, T;
    d = 1.0 / ph2v_reg3(p, h);
    T = ph2T_reg3(p, h);
    return Td_reg3(T, d, o_id);
}

double ps_reg3(double p, double s, int o_id)
// o_id: output propertry
{
    if (o_id == OP)
        return p;
    if (o_id == OS)
        return s;
    if (o_id == OX)
        return 1.0;
    if (o_id == OR)
        return 3.0;
    
    double d, T;
    d = 1.0 / ps2v_reg3(p, s);
    T = ps2T_reg3(p, s);
    return Td_reg3(T, d, o_id);
}

double hs_reg3(double h, double s, int o_id)
// o_id: output propertry
{
    if (o_id == OH)
        return h;
    if (o_id == OS)
        return s;
    if (o_id == OX)
        return 1.0;
    if (o_id == OR)
        return 3.0;

    double p, T, d;
    p = hs2p_reg3(h, s);

    T = ph2T_reg3(p, h);
    d = 1.0 / ph2v_reg3(p, h);

    return Td_reg3(T, d, o_id);
}

// IAPWS-IF97 Region3: The extended input pair
// * (p,v) (t,v), (t,h),(t,s)

///  Region3:  (p,v) T
double pv_reg3(double p, double v, int o_id)
{
    if (o_id == OP)
        return p;
    if (o_id == OV)
        return v;
    if (o_id == OX)
        return 1.0;
    if (o_id == OR)
        return 3.0;
    double d = 1.0 / v;
    if (o_id == OD)
    {
        return d;
    };
    double T = pv2T_reg3(p, v);
    if (o_id == OT)
    {
        return T;
    };
    return Td_reg3(T, d, o_id);
}

///  Region3:  (T,v)
double Tv_reg3(double T, double v, int o_id)
{
    if (o_id == OT)
        return T;
    if (o_id == OV)
        return v;
    if (o_id == OX)
        return 1.0;
    if (o_id == OR)
        return 3.0;
    double d = 1.0 / v;
    if (o_id == OD)
    {
        return d;
    };
    return Td_reg3(T, d, o_id);
}

//  Region3:  (T,h)
double Th_reg3(double T, double h, int o_id)
{
    if (o_id == OT)
        return T;
    if (o_id == OH)
        return h;
    if (o_id == OX)
        return 1.0;
    if (o_id == OR)
        return 3.0;
    double d = Th2d_reg3(T, h);
    if (o_id == OD)
    {
        return d;
    };
    if (o_id == OV)
    {
        return 1.0 / d;
    };
    return Td_reg3(T, d, o_id);
}

//  Region3:  (T,s)
double Ts_reg3(double T, double s, int o_id)
{
    if (o_id == OT)
        return T;
    if (o_id == OS)
        return s;
    if (o_id == OX)
        return 1.0;
    if (o_id == OR)
        return 3.0;
    double d = Ts2d_reg3(T, s);
    if (o_id == OD)
    {
        return d;
    };
    if (o_id == OV)
    {
        return 1.0 / d;
    };
    return Td_reg3(T, d, o_id);
}

double Td_thermal_reg3(double T, double d, int o_id)
// o_id: output propertry
{
    double value;
    switch (o_id)
    {
    case OT:
        value = T;
        break;
    case OX:
        value = 1.0;
        break;
    case OR:
        value = 3.0;
        break;
    case OD:
        value = d;
        break;
    case OV:
        value = 1 / d;
        break;
    case OP:
        value = Td2p_reg3(T, d);
        break;
    case OH:
        value = Td2h_reg3(T, d);
        break;
    case OS:
        value = Td2s_reg3(T, d);
        break;
    case OU:
        value = Td2u_reg3(T, d);
        break;
    case OCV:
        value = Td2cv_reg3(T, d);
        break;
    case OCP:
        value = Td2cp_reg3(T, d);
        break;
    case OW:
        value = Td2w_reg3(T, d);
        break;
    default:
        value = INVALID_OUTID;
        break;
    }
    return value;
}

double Td_ext_reg3(double T, double d, int o_id)
{
    double r = 0.0;
    switch (o_id)
    {
    case OT:
        r = T;
        break;
    case OX:
        r = 1.0;
        break;
    case OR:
        r = 3.0;
        break;
    case OD:
        r = d;
        break;
    case OV:
        r = 1 / d;
        break;
    case OE:
        r = Td2e_reg3(T, d);
        break;
    case OKISE:
        r = Td2k_reg3(T, d);
        break;
    case OF:
        r = Td2f_reg3(T, d);
        break;
    case OG:
        r = Td2g_reg3(T, d);
        break;
    case OZ:
        r = Td2z_reg3(T, d);
        break;
    case OKT:
        r = Td2kt_reg3(T, d);
        break;
    case OIPCEC:
        r = Td2ipcec_reg3(T, d);
        break;
    case ODVDP:
        r = Td2dvdpct_reg3(T, d);
        break;

    case ODVDT:
        r = Td2dvdtcp_reg3(T, d);
        break;
    case ODPDT:
        r = Td2dpdtcv_reg3(T, d);
        break;
    case OJTC:
        r = Td2joule_reg3(T, d);
        break;
    case OIJTC:
        r = Td2ijoule_reg3(T, d);
        break;
    case ODVDHCP:
        r = Td2dvdhcp_reg3(T, d);
        break;
    case ODVDPCH:
        r = Td2dvdpch_reg3(T, d);
        break;
    default:
        r = INVALID_OUTID;
    }
    return r;
}