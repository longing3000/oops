#include "include/app/app.h"
#include "include/app/DefectCenter.h"
#include <fstream>
#include <armadillo>
#include <assert.h>

int  main()
{
    //set nv
    vec nv_coord;nv_coord << 0.0 << 0.0 << 0.0;
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
    string filename="/home/david/code/test/dat/input/C13Bath/RoyCoord.xyz";
    cSpinSourceFromFile spin_file(filename);
    cSpinCollection bath_spins(&spin_file);
    bath_spins.make();
    vector<cSPIN> spin_list=bath_spins.getSpinList();


    //{{{check dephasing Hamiltonian and evolution
    vector<cSPIN> spin2;
    spin2.push_back(spin_list[0]);spin2.push_back(spin_list[1]);
    SpinDipolarInteraction dip(spin2);
    SpinZeemanInteraction zee(spin2,magB);
    DipolarField hf_field0(spin2,nv_espin,state0);
    DipolarField hf_field1(spin2,nv_espin,state1);
    
    //check Liouvillian
    Hamiltonian dip_H(spin2);
    dip_H.addInteraction(dip);
    dip_H.make();
    Liouvillian dip_H_sharp(dip_H,SHARP);
    Liouvillian dip_H_flat(dip_H,FLAT);
    dip_H_sharp.saveMatrix("sharp");
    dip_H_flat.saveMatrix("flat");

    Hamiltonian hami0(spin2);
    hami0.addInteraction(dip);
    hami0.addInteraction(zee);
    hami0.addInteraction(hf_field0);
    hami0.make();
    Hamiltonian hami1(spin2);
    hami1.addInteraction(dip);
    hami1.addInteraction(zee);
    hami1.addInteraction(hf_field1);
    hami1.make();
    double rate=1.0*2.0*datum::pi*1e4;
    vec dephase_axis;dephase_axis << 0.0 << 0.0 << 1.0;
    SpinDephasing dephasing(spin2,rate,dephase_axis);
    LiouvilleSpaceOperator dephaseOperator(spin2);
    dephaseOperator.addInteraction(dephasing);
    dephaseOperator.make();
    Liouvillian lv0(hami0,SHARP);
    Liouvillian lv1(hami1,FLAT);
    Liouvillian lv_d=lv0-lv1+dephaseOperator;
    Liouvillian lv=lv0-lv1;
    lv.saveMatrix("lv");
    hami0.saveMatrix("hami0");
    hami1.saveMatrix("hami1");
    lv_d.saveMatrix("lv_d");
    dephaseOperator.saveMatrix("dephase");

    vec bath_polarization=zeros<vec>(3);
    SpinPolarization p(spin2,bath_polarization);
    DensityOperator ds(spin2);
    ds.addStateComponent(p);
    ds.make();
    ds.makeVector();
    ds.saveVector("state_vec");

    //}}}
    
    //{{{hy perfine,dip
    string hype_file="/home/david/code/test/hyperfine_coeff.dat";
    string dip_file="/home/david/code/test/dip_coeff.dat";
    ofstream foutput1(dip_file.c_str());
    ofstream foutput(hype_file.c_str());
    if(!foutput) assert(0);
    if(!foutput1) assert(0);
    

    foutput << state_idx0 << "  " << state_idx1 << endl;

    //get coeff of hyperfine
    for(int i=0; i<spin_list.size(); ++i)
    {
        vec dip_field1=dipole_field(spin_list[i],nv_espin,state0);
        vec dip_field2=dipole_field(spin_list[i],nv_espin,state1);
        cout << dip_field1[2] << "...." <<  dip_field2[2] << endl; 
        foutput << dip_field1[2] << "  " << dip_field2[2] << endl;
   }
    
    foutput.close();

    //get coeff of dipolar
    for(int i=0; i<spin_list.size(); ++i)
        for(int j=i+1; j<spin_list.size(); ++j)
        {
            vec dip_coeff=dipole(spin_list[i],spin_list[j]);
            dip_coeff *= 100000.0;
            foutput1 << dip_coeff[0] << "   " << dip_coeff[1] << "  " << dip_coeff[2] << "  "
                     << dip_coeff[3] << "   " << dip_coeff[4] << "  " << dip_coeff[5] << "  "
                     << dip_coeff[6] << "   " << dip_coeff[7] << "  " << dip_coeff[8] << endl;
        }
    foutput1.close();
    //}}}

}
