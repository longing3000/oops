#include "include/quantum/QuantumState.h"

////////////////////////////////////////////////////////////////////////////////
//{{{ QuantumState 

QuantumState::QuantumState()
{ //LOG(INFO) << "Default constructor: QuantumState";
}

QuantumState::~QuantumState()
{ //LOG(INFO) << "Default destructor: QuantumState";
}

#ifdef HAS_MATLAB
void QuantumState::saveVector(string filename)
{
    cx_vec m=this->getVector();
    vec m_r=real(m);
    vec m_i=imag(m);

    mxArray *pArray = mxCreateDoubleMatrix(_dimension*_dimension,1,mxCOMPLEX);

    int dim=_dimension*_dimension;
    memcpy((void *)(mxGetPr(pArray)), (void *) m_r.memptr(), dim*sizeof(double));
    memcpy((void *)(mxGetPi(pArray)), (void *) m_i.memptr(), dim*sizeof(double));

    string dpg_filename =DEBUG_PATH + filename + ".mat";
    cout << dpg_filename << " is exported for debug!" << endl;
    MATFile *mFile = matOpen(dpg_filename.c_str(), "w");
    matPutVariableAsGlobal(mFile, filename.c_str(), pArray);
    matClose(mFile);
    
    mxDestroyArray(pArray);
}
#else
void QuantumState::saveVector(string filename)
{
    cout << "MATLAB not installed." << endl;
}
#endif
//}}}
////////////////////////////////////////////////////////////////////////////////
