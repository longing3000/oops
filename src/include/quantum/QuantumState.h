#ifndef QUANTUMSTATE_H
#define QUANTUMSTATE_H

#ifdef HAS_MATLAB
#include <mat.h>
#endif

#include <armadillo>
#include <string>

using namespace std;
using namespace arma;
extern string DEBUG_PATH;
/// \addtogroup Quantum 
/// @{

/// \defgroup QuantumState QuantumState
/// @{


////////////////////////////////////////////////////////////////////////////////
//{{{ QuantumState
class QuantumState
{
public:
    QuantumState();
    ~QuantumState();

    cx_vec getVector() const {return _vector;};
    void saveVector(string filename);
    size_t    getDimension() const {return _dimension;};
protected:
    bool _is_pure;
    size_t _dimension;
    cx_vec _vector;
private:
};
//}}}
////////////////////////////////////////////////////////////////////////////////

/// @}
#endif
