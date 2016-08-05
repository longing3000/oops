#include "include/app/app.h"
#include "include/app/DefectCenter.h"
#include <fstream>
#include <armadillo>
#include <assert.h>

int  main()
{
    string save_file="/home/david/code/test/hyperfine_coeff.dat";
    ofstream foutput(save_file.c_str());
    if(!foutput) assert(0);

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

    foutput << state_idx0 << "  " << state_idx1 << endl;

    //set bath spins
    string filename="/home/david/code/test/dat/input/C13Bath/RoyCoord.xyz";
    cSpinSourceFromFile spin_file(filename);
    cSpinCollection bath_spins(&spin_file);
    bath_spins.make();
    vector<cSPIN> spin_list=bath_spins.getSpinList();

    //get coeff of hyperfine
    for(int i=0; i<spin_list.size(); ++i)
    {
        vec dip_field1=dipole_field(spin_list[i],nv_espin,state0);
        vec dip_field2=dipole_field(spin_list[i],nv_espin,state1);
        cout << dip_field1[2] << "...." << dip_field2[2] << endl; 
        foutput << dip_field1[2] << "  " << dip_field2[2] << endl;
   }
    
    foutput.close();
}

