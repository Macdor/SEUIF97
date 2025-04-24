#include <wchar.h>
#include "mex.h"
#include "matrix.h"
#include "seuif97.h"

/* The computational routine */
static void calc_seuif97(double* out, mxChar* funChar, double *in1, double *in2, int *o_id)
{
    bool flag = false;
    if (wcscmp(funChar, L"pt")==0) *out = pt(*in1, *in2, *o_id); flag = true;
    if (wcscmp(funChar, L"ph")==0) *out = ph(*in1, *in2, *o_id); flag = true;
    if (wcscmp(funChar, L"ps")==0) *out = ps(*in1, *in2, *o_id); flag = true;
    if (wcscmp(funChar, L"pv")==0) *out = pv(*in1, *in2, *o_id); flag = true;
    if (wcscmp(funChar, L"th")==0) *out = th(*in1, *in2, *o_id); flag = true;
    if (wcscmp(funChar, L"ts")==0) *out = ts(*in1, *in2, *o_id); flag = true;
    if (wcscmp(funChar, L"tv")==0) *out = tv(*in1, *in2, *o_id); flag = true;
    if (wcscmp(funChar, L"hs")==0) *out = hs(*in1, *in2, *o_id); flag = true;
    if (wcscmp(funChar, L"px")==0) *out = px(*in1, *in2, *o_id); flag = true;
    if (wcscmp(funChar, L"tx")==0) *out = tx(*in1, *in2, *o_id); flag = true;
    if (wcscmp(funChar, L"hx")==0) *out = hx(*in1, *in2, *o_id); flag = true;
    if (wcscmp(funChar, L"sx")==0) *out = sx(*in1, *in2, *o_id); flag = true;
    if (!flag) *out = -1.0;
}

/* The gateway function */
void mexFunction(int nlhs, mxArray* plhs[],
    int nrhs, const mxArray* prhs[])
{

    if (nrhs != 4) {
        mexErrMsgIdAndTxt("MyToolbox:mex_seuif97:nrhs", "Four inputs required.");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("MyToolbox:mex_seuif97:nlhs", "One output required.");
    }
    if (!mxIsChar(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:mex_seuif97:notString", "First input must be a function string.");
    }
    if (!mxIsNumeric(prhs[1]) ||
        mxIsComplex(prhs[1]) ||
        mxGetNumberOfElements(prhs[1]) != 1) {
        mexErrMsgIdAndTxt("MyToolbox:mex_seuif97:notRealScalar", "Second input must be real scalar.");
    }
    if (!mxIsNumeric(prhs[2]) ||
        mxIsComplex(prhs[2]) ||
        mxGetNumberOfElements(prhs[2]) != 1) {
        mexErrMsgIdAndTxt("MyToolbox:mex_seuif97:notRealScalar", "Third input must be real scalar.");
    }
    if (!mxIsInt32(prhs[3]) ||
        mxIsComplex(prhs[3]) ||
        mxGetNumberOfElements(prhs[3]) != 1) {
        mexErrMsgIdAndTxt("MyToolbox:mex_seuif97:notRealIntScalar", "Fourth input must be real int scalar.");
    }
    mxChar* strfun = mxGetChars(prhs[0]);
    mxDouble* in1 = mxGetDoubles(prhs[1]);
    mxDouble* in2 = mxGetDoubles(prhs[2]);
    mxInt32* o_id = mxGetInt32s(prhs[3]);

    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    mxDouble* out = mxGetDoubles(plhs[0]);

    calc_seuif97(out, strfun, in1, in2, o_id);
}