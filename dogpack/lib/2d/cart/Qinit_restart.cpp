
class dTensorBC4;
double ReadStateArray(const char* framedir, const char* varname,
    dTensorBC4& q, int nstart);
//
//
// deprecated
#if 0
void QinitRestart(
    int nstart,
    const char* varname,
    dTensorBC4& q,
    const char* outputdir)
{
    ReadStateArray(outputdir, varname, q, nstart);
}
#endif

// deprecated
void Qinit_restart(
    int nstart,
    dTensorBC4& q,
    const char* outputdir)
{
    ReadStateArray(outputdir, "q", q, nstart);
}

