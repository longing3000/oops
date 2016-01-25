#include <vector>
#include <string>
#include <iostream>
#include "include/spin/Spin.h"
#include "include/spin/SpinData.h"
#include "include/spin/SpinCollection.h"
#include "include/spin/SpinSource.h"
#include "include/spin/SpinCluster.h"
#include "include/spin/SpinClusterAlgorithm.h"
#include "include/spin/SpinInteraction.h"
#include "include/kron/KronProd.h"

#include "include/easylogging++.h"
#include "include/misc/misc.h"
#include "include/quantum/HilbertSpaceOperator.h"
#include "include/quantum/LiouvilleSpaceOperator.h"
#include "include/quantum/MixedState.h"
#include "include/quantum/PureState.h"
#include "include/spin/SpinState.h"
#include "include/quantum/QuantumEvolutionAlgorithm.h"
#include "include/quantum/QuantumEvolution.h"

_INITIALIZE_EASYLOGGINGPP

using namespace std;
using namespace arma;

cSPINDATA SPIN_DATABASE=cSPINDATA();

cSPIN            create_e_spin();
cSpinCollection  create_bath_spins_from_file();
cSpinCluster     create_spin_clusters(const cSpinCollection& sc);
Hamiltonian      create_spin_hamiltonian(const cSPIN& espin, const int spin_state, const vector<cSPIN>& spin_list);
Liouvillian      create_spin_liouvillian(const Hamiltonian& hami0, const Hamiltonian hami1);
DensityOperator  create_spin_density_state(const vector<cSPIN>& spin_list);

int  main(int argc, char* argv[])
{
    _START_EASYLOGGINGPP(argc, argv);
    easyloggingpp::Configurations confFromFile("../src/logs/log.conf");  // Load configuration from file
    easyloggingpp::Loggers::reconfigureAllLoggers(confFromFile); // Re-configures all the loggers to current configuration file
    LOG(INFO) << "################################################### Program begins ###################################################"; 

    cSPIN espin = create_e_spin();

    cSpinCollection spin_collection = create_bath_spins_from_file();

    cSpinCluster spin_clusters = create_spin_clusters(spin_collection);

    vector<cSPIN> spin_list = spin_clusters.getCluster(2, 4);

    int spin_up = 0, spin_down = 1;
    Hamiltonian hami0 = create_spin_hamiltonian(espin, spin_up, spin_list);
    Hamiltonian hami1 = create_spin_hamiltonian(espin, spin_down, spin_list);

    Liouvillian lv = create_spin_liouvillian(hami0, hami1);

    DensityOperator ds = create_spin_density_state(spin_list);

    SimpleFullMatrixVectorEvolution kernel(lv, ds);
    kernel.setTimeSequence( linspace<vec>(0.0, 1.0, 101) );

    QuantumEvolution dynamics(&kernel);
    dynamics.run();
}

////////////////////////////////////////////////////////////////////////////////
//{{{ Create an electrion spin
cSPIN create_e_spin()
{
    string isotope="E";

    double coord[] = {0.0, 0.0, 0.0};
    vector<double> coordinate(coord, coord+3);

    cSPIN espin=cSPIN(coordinate, isotope);
    return espin;
}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ Create bath spin list from xyz file
cSpinCollection create_bath_spins_from_file()
{
    cSpinSourceFromFile spin_file("../bin/RoyCoord.xyz");
    cSpinCollection sc(&spin_file);
    sc.make();
    return sc;
}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ Create spin clusters from a given spin list
cSpinCluster create_spin_clusters(const cSpinCollection& sc)
{
    sp_mat c=sc.getConnectionMatrix(8.0);

    cDepthFirstPathTracing dfpt(c, 3);
    cSpinCluster cluster(sc, &dfpt);

    cluster.make();
    return cluster;
}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ Create spin Hamiltonian for a given cluster
Hamiltonian create_spin_hamiltonian(const cSPIN& espin, const int spin_state, const vector<cSPIN>& spin_list)
{
    SpinDipolarInteraction dip(spin_list);

    vec magB; 
    magB << 3.0e-4 << 2.0e-4 << 1.0e-4;
    SpinZeemanInteraction zee(spin_list, magB);

    PureState center_spin_state(espin); 
    center_spin_state.setComponent(spin_state, 1.0);
    DipolarField hf_field(spin_list, espin, center_spin_state);

    Hamiltonian hami(spin_list);
    hami.addInteraction(dip);
    hami.addInteraction(zee);
    hami.addInteraction(hf_field);
    hami.make();
    return hami;
}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ Create Liouvillian operator from given Hamiltonians
Liouvillian create_spin_liouvillian(const Hamiltonian& hami0, const Hamiltonian hami1)
{
    Liouvillian lv0(hami0, SHARP);
    Liouvillian lv1(hami1, FLAT);

    Liouvillian lv = lv0 - lv1;
    cx_mat lvMat = lv.getMatrix();
    cout << lvMat << endl;
    return lv;
}
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{ Create spin density matrix
DensityOperator create_spin_density_state(const vector<cSPIN>& spin_list)
{
    vec pol = zeros<vec>(3);
    SpinPolarization p(spin_list, pol);

    DensityOperator ds(spin_list);
    ds.addStateComponent(p);
    ds.make();
    ds.makeVector();
    cout << ds.getVector() << endl;
    cout << ds.getMatrix() << endl;
    return ds;
}
//}}}
////////////////////////////////////////////////////////////////////////////////


