#include "include/app/ensemble_cce_exf.h"


////////////////////////////////////////////////////////////////////////////////
//{{{  EXFEnsembleCCE
void EXFEnsembleCCE::set_parameters()
{/*{{{*/
    string input_filename  = _para["input"].as<string>();
    string output_filename = _para["output"].as<string>();
    string input_field_filename  = _para["input_data"].as<string>();

    _state_idx0            = _para["state0"].as<int>();
    _state_idx1            = _para["state1"].as<int>();

    _max_order             = _para["cce"].as<int>();
    _cut_off_dist          = _para["cutoff"].as<double>();
    _bath_polarization     = vec( _para["polarization"].as<string>() );
    _bath_dephasing_rate   = _para["dephasing_rate"].as<double>();
    _bath_dephasing_axis   = vec( _para["dephasing_axis"].as<string>() );

    _t0                    = _para["start"].as<double>(); 
    _t1                    = _para["finish"].as<double>();

    _field_axis            = vec( _para["field axis"].as<string>() );//two confusion, the vector?, as<double>?
    _omega                 =_para["omega"].as<double>();
    _magB                  = vec( _para["magnetic_field"].as<string>() );
    
    _bath_spin_filename = INPUT_PATH + input_filename;
    _result_filename    = OUTPUT_PATH + output_filename + ".mat";
    _external_field_filename = INPUT_PATH + input_field_filename;
}/*}}}*/

vec EXFEnsembleCCE::cluster_evolution(int cce_order, int index)
{/*{{{*/
    vector<cSPIN> spin_list = _my_clusters.getCluster(cce_order, index);
    
    vector<QuantumOperator> lv_list;
    vector<double> time_segment;
    
    for(int i =0; i<_amplitude_list.size(); ++i)
    {
        Hamiltonian hami0 = create_spin_hamiltonian(_center_spin, _state_pair.first, spin_list,_amplitude_list[i],_phase_list[i],_field_axis,_omega);
        Hamiltonian hami1 = create_spin_hamiltonian(_center_spin, _state_pair.second, spin_list,_amplitude_list[i],_phase_list[i],_field_axis,_omega);
        LiouvilleSpaceOperator dephase = create_incoherent_operator(spin_list);
        QuantumOperator lvA = create_spin_liouvillian(hami0, hami1) + dephase;
        lv_list.push_back(lvA);
        time_segment.push_back(_time_list[i+1]-_time_list[i]);
    }
    DensityOperator ds = create_spin_density_state(spin_list);//no polarization
    NEQPiecewiseFullMatrixVectorEvolution kernel(lv_list, time_segment, ds);
    kernel.setTimeSequence( _t0, _t1, _nTime);

    ClusterCoherenceEvolution dynamics(&kernel);
    dynamics.run();

    return calc_observables(&kernel);
}/*}}}*/

Hamiltonian EXFEnsembleCCE::create_spin_hamiltonian(const cSPIN& espin, const PureState& center_spin_state, const vector<cSPIN>& spin_list, const double& amplitude, const double& phase, const vec& axis, const double& omega)
{/*{{{*/
    RWASpinDipolarInteraction dip(spin_list);

    RWASpinZeemanInteraction zee(spin_list, _magB,_omega);

    RWADipolarField hf_field(spin_list, espin, center_spin_state);
    
    ExternalField ex_field(spin_list, amplitude, phase,axis);

    Hamiltonian hami(spin_list);
    hami.addInteraction(dip);
    hami.addInteraction(zee);
    hami.addInteraction(hf_field);
    hami.addInteraction(ex_field);
    hami.make();
    return hami;
}/*}}}*/

LiouvilleSpaceOperator EXFEnsembleCCE::create_incoherent_operator(const vector<cSPIN>& spin_list)
{/*{{{*/
    double rate = 2.0*datum::pi*_bath_dephasing_rate;
    vec axis = normalise(_bath_dephasing_axis);

    SpinDephasing dephasing(spin_list, rate, axis);
    LiouvilleSpaceOperator dephaseOperator(spin_list);
    dephaseOperator.addInteraction(dephasing);
    dephaseOperator.make();
    return dephaseOperator;

}/*}}}*/
Liouvillian EXFEnsembleCCE::create_spin_liouvillian(const Hamiltonian& hami0, const Hamiltonian hami1)
{/*{{{*/
    Liouvillian lv0(hami0, SHARP);
    Liouvillian lv1(hami1, FLAT);
    Liouvillian lv = lv0 - lv1;
    return lv;
}/*}}}*/

DensityOperator EXFEnsembleCCE::create_spin_density_state(const vector<cSPIN>& spin_list)
{/*{{{*/
    SpinPolarization p(spin_list, _bath_polarization);

    DensityOperator ds(spin_list);
    ds.addStateComponent(p);
    ds.make();
    ds.makeVector();
    return ds;
}/*}}}*/

vec EXFEnsembleCCE::calc_observables(QuantumEvolutionAlgorithm* kernel)
{/*{{{*/
    int dim = kernel->getMatrixDim();
    cx_vec init_st = kernel->getInitalState();
    vector<cx_vec>  state = kernel->getResult();
    vec res = ones<vec>(_nTime);
    for(int i=0; i<_nTime; ++i)
        res(i) = real( cdot( dim*init_st, state[i]) );
    return res;
}/*}}}*/
//}}}
////////////////////////////////////////////////////////////////////////////////





