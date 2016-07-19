#include "include/app/app.h"
#include "include/app/cce.h"
#include "include/math/MatExp.h"
#include <complex>
#include "expv/include/expv.h"

cx_mat MAT;
cx_vec VEC;
vec    TIME_LIST;
SumKronProd SKP;
cx_double PREFACTOR;


void prepare_data(string filename);

int  main(int argc, char* argv[])
{
    string filename = "./dat_example/input/C13Bath/RoyCoord.xyz";
    prepare_data(filename);

    return 0;
}

void prepare_data(string filename)
{/*{{{*/
    // create defect center
    double x =  0.0, y = 0.0, z = 0.0;
    vec coord; coord << x << y << z;
    NVCenter nv(NVCenter::N14, coord);//creat a defect center call NVCenter 
    
    double magBx = 0.0,  magBy =  0.0, magBz = 0.2;
    nv.set_magB(magBx, magBy, magBz);
    nv.make_espin_hamiltonian();//nv state

    cout << nv.get_eigen_state(0) << endl;
    cout << nv.get_eigen_state(1) << endl;

    cSPIN espin=nv.get_espin();
    PureState st0(nv.get_eigen_state(0));
    PureState st1(nv.get_eigen_state(1));
    espin.set_coordinate(coord);
    cout << "espin coordinate = " << espin.get_coordinate() << endl;

    // create bath spins
    cout << "start to generate." << endl;
    cSpinSourceFromFile spin_file(filename);
    cSpinCollection spins(&spin_file);
    spins.make();
    cout << "generate success." << endl;

    vector<cSPIN> sl = spins.getSpinList();

    cout << sl[0].get_coordinate() << sl[0].get_gamma() << endl;
    cout << sl[1].get_coordinate() << sl[1].get_gamma() << endl;


    vec magB; magB << magBx << magBy << magBz;
    double amplitude=10000.0;
    double phase=1.5;double omega=sl[0].get_gamma()*magBz;
    vec field_axis;field_axis << 1.0 << 0.0 << 0.0;
    RWASpinZeemanInteraction zee(sl, magB,omega);//Zeeman interaction of bath spins
    RWASpinDipolarInteraction dip(sl);//Dipolar interaction of bath spins
    RWADipolarField hf_field0(sl, espin, st0);//hyperfine interaction from spins state1
    RWADipolarField hf_field1(sl, espin, st1);//hyperfine interaction from center spins state2
    ExternalField ex_field(sl,amplitude,phase,field_axis); 

    Hamiltonian hami0(sl);
    hami0.addInteraction(zee);
    hami0.addInteraction(dip);
    hami0.addInteraction(hf_field0);
    hami0.addInteraction(ex_field);
    hami0.make();//zee+dip+hf_field

    Hamiltonian hami1(sl);
    hami1.addInteraction(zee);
    hami1.addInteraction(dip);
    hami1.addInteraction(hf_field1);
    hami1.addInteraction(ex_field);
    hami1.make();//zee+dip-hf_field


    //cout << hami0.getMatrix() << endl;
    //cout << hami1.getMatrix() << endl;
    hami0.saveMatrix("hami0");
    hami1.saveMatrix("hami1");

    Liouvillian lv0(hami0, SHARP);
    Liouvillian lv1(hami1, FLAT);
    Liouvillian lvH = lv0 - lv1;//Liouvillian


//    double rate =0.0*2.0*datum::pi*1e4;//dephasing rate
//    vec axis; axis << 1.0 << 1.0 << 1.0;
//    SpinDephasing dephasing(sl, rate, normalise(axis));//set dephaing for each spins
//    LiouvilleSpaceOperator dephaseOperator(sl);
//    dephaseOperator.addInteraction(dephasing);
//    dephaseOperator.make();//generate dephasing Liouvillian
//
//    QuantumOperator lv = lvH + dephaseOperator;
    QuantumOperator lv = lvH;
    lv.saveMatrix("lv");//generate all the Liouvillian


    vec _bath_polarization = zeros<vec>(3);//prepare for initial state of bath spins
    SpinPolarization p(sl, _bath_polarization);//this means no polarization of the bath spins, which means that the initial state of bath spins is identity matrix
    DensityOperator ds(sl);
    ds.addStateComponent(p);
    ds.make();//generate initial state in Hilbert space 
    ds.makeVector();//generate initial state in Liouville space
    ds.saveVector("state_vec");//saveVector() function is wrong!!!
    ds.saveMatrix("state_mat");//saveMatrix() function is wrong!!!

}/*}}}*/
