#include "include/quantum/MixedState.h"

////////////////////////////////////////////////////////////////////////////////
//{{{ MixedState 

MixedState::MixedState()
{ //LOG(INFO) << "Default constructor: MixedState";
}

MixedState::~MixedState()
{ //LOG(INFO) <<  "Default destructor: MixedState";
}
//}}}
/////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ DensityOperator
DensityOperator::DensityOperator()
{ //LOG(INFO) << "Default constructor: DensityOperator";
}

DensityOperator::~DensityOperator()
{ //LOG(INFO) << "Default destructor: DensityOperator";
}

DensityOperator::DensityOperator(const vector<cSPIN>& spin_list)
{
    _op = HilbertSpaceOperator(spin_list);
    _dimension = _op.getDimension();
}

#ifdef HAS_MATLAB
void DensityOperator::saveMatrix(string filename)
{
    cx_mat m=this->getMatrix();
    mat m_r= real(m).t();
    mat m_i=-imag(m).t();
    
    mxArray *pArray = mxCreateDoubleMatrix(_dimension,_dimension,mxCOMPLEX);

    int dim=_dimension*_dimension;
    memcpy((void *)(mxGetPr(pArray)), (void *) m_r.memptr(), dim*sizeof(double));
    memcpy((void *)(mxGetPi(pArray)), (void *) m_i.memptr(), dim*sizeof(double));

    string dpg_filename = DEBUG_PATH + filename + ".mat";
    cout << dpg_filename << " is exported for debug!" << endl;
    MATFile *mFile = matOpen(dpg_filename.c_str(), "w");
    matPutVariableAsGlobal(mFile, filename.c_str(), pArray);
    matClose(mFile);

    mxDestroyArray(pArray);      
}
#else
void DensityOperator::saveMatrix(string filename)
{
    cout << "MATLAB not installed." << endl;
}
#endif
//}}}
////////////////////////////////////////////////////////////////////////////////
