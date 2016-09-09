#include "include/app/app.h"
#include "include/app/DefectCenter.h"
#include <fstream>
#include <armadillo>
#include <assert.h>

int  main()
{
    //set nv
    vec nv_coord;nv_coord << 0.0 << 0.0 << -20.0;
    vec magB;magB << 0.1 << 0.1 << 0.1;
    NVCenter nv(NVCenter::N14,nv_coord);
    nv.set_magB(magB);
    nv.make_espin_hamiltonian();
    cSPIN nv_espin=nv.get_espin();
    int state_idx0=0;
    int state_idx1=1;
    cx_vec state0=nv.get_eigen_state(state_idx0);
    cx_vec state1=nv.get_eigen_state(state_idx1);
    
    cout << "state 0:" << nv.get_eigen_value(0) << endl;
    cout << state0 << endl;
    cout << "state 1:" << nv.get_eigen_value(1) << endl;
    cout << state1 << endl;
    cout << "state 2:" << nv.get_eigen_value(2) << endl;
    cout << nv.get_eigen_state(2) << endl;

    //set bath spins
    string filename="/home/david/code/oops/src/coord.xyz";
    cSpinSourceFromFile spin_file(filename);
    cSpinCollection bath_spins(&spin_file);
    bath_spins.make();
    vector<cSPIN> spin_list=bath_spins.getSpinList();


    //{{{check Hamiltonian and evolution
    //check single spin Hamiltonian and evolution
    vector<cSPIN> spin2;
    spin2.push_back(spin_list[0]);
    //SpinDipolarInteraction dip(spin2);
    SpinZeemanInteraction zee(spin2,magB);
    DipolarField hf_field0(spin2,nv_espin,state0);
    DipolarField hf_field1(spin2,nv_espin,state1);

    Hamiltonian hami0(spin2);
    //hami0.addInteraction(dip);
    hami0.addInteraction(zee);
    hami0.addInteraction(hf_field0);
    hami0.make();
    Hamiltonian hami1(spin2);
    //hami1.addInteraction(dip);
    hami1.addInteraction(zee);
    hami1.addInteraction(hf_field1);
    hami1.make();
    Liouvillian lv0(hami0,SHARP);
    Liouvillian lv1(hami1,FLAT);
    lv0.saveMatrix("lv0_1");
    lv1.saveMatrix("lv1_1");


    double rate=2.0*datum::pi;
    vec dephase_axis;dephase_axis << 1.0 << 1.0 << 1.0;
    SpinDephasing dephasing(spin2,rate,dephase_axis);
    LiouvilleSpaceOperator dephaseOperator(spin2);
    dephaseOperator.addInteraction(dephasing);
    dephaseOperator.make();
    dephaseOperator.saveMatrix("dephaseOperator1");

    vec bath_polarization=zeros<vec>(3);
    SpinPolarization p(spin2,bath_polarization);
    DensityOperator ds(spin2);
    ds.addStateComponent(p);
    ds.make();
    ds.makeVector();
    ds.saveVector("state_vec1");

    //}}}
}
