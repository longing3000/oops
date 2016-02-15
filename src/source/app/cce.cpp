#include "include/app/cce.h"

////////////////////////////////////////////////////////////////////////////////
//{{{  CCE
CCE::CCE(int my_rank, int worker_num, const string& config_file)
{ 
    _my_rank = my_rank;
    _worker_num = worker_num;
    
    char cfg_path[500];
    strcpy(cfg_path, PROJECT_PATH);
    strcat(cfg_path, "/dat/config/");
    strcat(cfg_path, config_file.c_str() );
    _cfg = ConfigXML(cfg_path);

    if(my_rank == 0)
        _cfg.printParameters();
}

void CCE::run()
{
    set_parameters();
    create_center_spin();
    create_bath_spins();
    create_spin_clusters();

    run_each_clusters();
    post_treatment();

}

cSPIN CCE::create_center_spin()
{
    _center_spin=cSPIN(_center_spin_coord, _center_spin_isotope);
    return _center_spin;
}

cSpinCollection CCE::create_bath_spins()
{
    cSpinSourceFromFile spin_file(_bath_spin_filename);
    _bath_spins = cSpinCollection(&spin_file);
    _bath_spins.make();

    if(_my_rank == 0)
        cout << _bath_spins.getSpinNum() << " spins are read from file: " << _bath_spin_filename << endl << endl;
    return _bath_spins;
}

void CCE::create_spin_clusters()
{
    if(_my_rank == 0)
    {
        sp_mat c=_bath_spins.getConnectionMatrix(_cut_off_dist);
        cDepthFirstPathTracing dfpt(c, _max_order);
        _spin_clusters=cSpinCluster(_bath_spins, &dfpt);
        _spin_clusters.make();
    }

    job_distribution();
}

void CCE::job_distribution()
{/*{{{*/
    uvec clstLength;     vector<umat> clstMat;
    if(_my_rank == 0)
    {
        _spin_clusters.MPI_partition(_worker_num);
        
        clstLength = _spin_clusters.getMPI_ClusterLength(0);
        clstMat = _spin_clusters.getMPI_Cluster(0);
        
        for(int i=1; i<_worker_num; ++i)
        {
            uvec clstNum = _spin_clusters.getMPI_ClusterLength(i);
            MPI_Send(clstNum.memptr(), _max_order, MPI_UNSIGNED, i, 0, MPI_COMM_WORLD);
            
            vector<umat> clstMatList = _spin_clusters.getMPI_Cluster(i);
            for(int j=0; j<_max_order; ++j)
            {
                umat clstMat_j = clstMatList[j];
                MPI_Send(clstMat_j.memptr(), (j+1)*clstNum(j), MPI_UNSIGNED, i, j+1, MPI_COMM_WORLD);
            }
        }
    }
    else
    {
        unsigned int * clstLengthData = new unsigned int [_max_order];
        MPI_Recv(clstLengthData, _max_order, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        uvec tempV(clstLengthData, _max_order);
        clstLength = tempV;
        delete [] clstLengthData;
        
        for(int j=0; j<_max_order;++j)
        {
            unsigned int * clstMatData = new unsigned int [(j+1)*clstLength(j)];
            MPI_Recv(clstMatData, (j+1)*clstLength(j), MPI_UNSIGNED, 0, j+1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            umat tempM(clstMatData, clstLength(j), j+1);
            clstMat.push_back(tempM);
            delete [] clstMatData;
        }
    }
    _my_clusters = cSpinCluster(_bath_spins, clstLength, clstMat);
}/*}}}*/

void CCE::run_each_clusters()
{
    for(int cce_order = 0; cce_order < _max_order; ++cce_order)
    {
        cout << "my_rank = " << _my_rank << ": " << "calculating order = " << cce_order << endl;
        size_t clst_num = _my_clusters.getClusterNum(cce_order);
        
        mat resMat(_nTime, clst_num, fill::ones);
        for(int i = 0; i < clst_num; ++i)
        {
            cout << "my_rank = " << _my_rank << ": " << i << "/" << clst_num << endl;
            resMat.col(i) = cluster_evolution(cce_order, i);
        }
        
        DataGathering(resMat, cce_order, clst_num);
    }
}
void CCE::DataGathering(const mat& resMat, int cce_order, int clst_num)
{/*{{{*/

    if(_my_rank != 0)
        MPI_Send(resMat.memptr(), _nTime*clst_num, MPI_DOUBLE, 0, 100+_my_rank, MPI_COMM_WORLD);
    else
    {
        double * cce_evolve_data= new double [_nTime * _spin_clusters.getClusterNum(cce_order)];
        memcpy(cce_evolve_data, resMat.memptr(), _nTime*clst_num*sizeof(double));
        
        size_t prev_clst_num = clst_num;
        for(int source = 1; source < _worker_num; ++source)
        {
            pair<size_t, size_t> pos = _spin_clusters.getMPI_ClusterSize(cce_order, source);
            MPI_Recv(cce_evolve_data + _nTime*pos.first, _nTime*(pos.second - pos.first), MPI_DOUBLE, source, 100+source, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        
        mat res_i(cce_evolve_data, _nTime, _spin_clusters.getClusterNum(cce_order) );
        _cce_evovle_result.push_back(res_i);

        delete [] cce_evolve_data;
    }
}/*}}}*/

void CCE::post_treatment()
{
    if(_my_rank == 0)
    {
        cce_coherence_reduction();
        compuate_final_coherence();
        export_mat_file();
    }
}

void CCE::cce_coherence_reduction()
{/*{{{*/
    for(int cce_order = 0; cce_order<_max_order; ++cce_order)
    {
        mat tilder_mat =  ones( size(_cce_evovle_result[cce_order]) );
        for(int j = 0; j<_cce_evovle_result[cce_order].n_cols; ++j)
        {
            vec res_j = _cce_evovle_result[cce_order].col(j);
            
            set<CluserPostion > sub_pos = _spin_clusters.getSubClusters(cce_order, j);
            for(set<CluserPostion >::iterator it=sub_pos.begin(); it!=sub_pos.end(); ++it)
            {
                vec sub_res = _cce_evovle_result_tilder[it->first].col(it->second);
                res_j = res_j / sub_res;
            }
            tilder_mat.col(j) = res_j;
        }
        _cce_evovle_result_tilder.push_back( tilder_mat );
    }

}/*}}}*/

void CCE::compuate_final_coherence()
{/*{{{*/
    _final_result = mat(_nTime, _max_order, fill::zeros);
    _final_result_each_order = mat(_nTime, _max_order, fill::zeros);

    vec final_res_vec = ones<vec> (_nTime);
    for(int cce_order = 0; cce_order<_max_order; ++cce_order)
    {
        vec res_vec = ones<vec> (_nTime);
        for(int j=0; j<_cce_evovle_result_tilder[cce_order].n_cols; ++j)
            res_vec = res_vec % _cce_evovle_result_tilder[cce_order].col(j);
        _final_result_each_order.col(cce_order) = res_vec;
        
        final_res_vec = final_res_vec % res_vec;
        _final_result.col(cce_order)= final_res_vec;
    }
}/*}}}*/

void CCE::export_mat_file() 
{/*{{{*/
#ifdef HAS_MATLAB
    cout << "begin post_treatement ... storing cce_data to file: " << _result_filename << endl;
    MATFile *mFile = matOpen(_result_filename, "w");
    for(int i=0; i<_max_order; ++i)
    {
        char i_str [10];
        sprintf(i_str, "%d", i);
        string idx_str = i_str;
        string label = "CCE" + idx_str;
        string label1 = "CCE" + idx_str+"_tilder";
        
        size_t nClst = _spin_clusters.getClusterNum(i);
        mxArray *pArray = mxCreateDoubleMatrix(_nTime, nClst, mxREAL);
        mxArray *pArray1 = mxCreateDoubleMatrix(_nTime, nClst, mxREAL);
        
        size_t length= _nTime * nClst;
        memcpy((void *)(mxGetPr(pArray)), (void *) _cce_evovle_result[i].memptr(), length*sizeof(double));
        memcpy((void *)(mxGetPr(pArray1)), (void *) _cce_evovle_result_tilder[i].memptr(), length*sizeof(double));
        
        matPutVariableAsGlobal(mFile, label.c_str(), pArray);
        matPutVariableAsGlobal(mFile, label1.c_str(), pArray1);
        
        mxDestroyArray(pArray);
        mxDestroyArray(pArray1);
    }

    mxArray *pRes = mxCreateDoubleMatrix(_nTime, _max_order, mxREAL);
    mxArray *pRes1 = mxCreateDoubleMatrix(_nTime, _max_order, mxREAL);
    size_t length= _nTime*_max_order;
    memcpy((void *)(mxGetPr(pRes)), (void *) _final_result_each_order.memptr(), length*sizeof(double));
    memcpy((void *)(mxGetPr(pRes1)), (void *) _final_result.memptr(), length*sizeof(double));
    matPutVariableAsGlobal(mFile, "final_result_each_order", pRes);
    matPutVariableAsGlobal(mFile, "final_result", pRes1);
    mxDestroyArray(pRes);
    mxDestroyArray(pRes1);
    matClose(mFile);
#endif
}/*}}}*/
//}}}
////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////
//{{{  EnsembleCCE
EnsembleCCE::EnsembleCCE(int my_rank, int worker_num, const string& config_file)
{ /*{{{*/
    _my_rank = my_rank;
    _worker_num = worker_num;

    char cfg_path[500];
    strcpy(cfg_path, PROJECT_PATH);
    strcat(cfg_path, "/dat/config/");
    strcat(cfg_path, config_file.c_str() );
    _cfg = ConfigXML(cfg_path);

    if(my_rank == 0)
        _cfg.printParameters();
}/*}}}*/

void EnsembleCCE::set_parameters()
{/*{{{*/
    strcpy(_bath_spin_filename, PROJECT_PATH); 
    strcpy(_result_filename, PROJECT_PATH); 
    strcat(_bath_spin_filename, "/dat/input/RoyCoord.xyz");
    strcat(_result_filename, "/dat/output/cce_res.mat");

    double x = _cfg.getDoubleParameter("CenterSpin", "coordinateX");
    double y = _cfg.getDoubleParameter("CenterSpin", "coordinateY");
    double z = _cfg.getDoubleParameter("CenterSpin", "coordinateZ");
    _center_spin_coord << x << y << z;

    _center_spin_isotope = _cfg.getStringParameter("CenterSpin", "isotope");
    _cut_off_dist        = _cfg.getDoubleParameter("SpinBath", "cut_off_dist");
    _max_order           = _cfg.getIntParameter("CCE", "max_order");
    _nTime               = _cfg.getIntParameter("Dynamics", "nTime");
    _t0                  = _cfg.getDoubleParameter("Dynamics", "t0"); 
    _t1                  = _cfg.getDoubleParameter("Dynamics", "t1"); 
    _pulse_name          = _cfg.getStringParameter("Condition", "pulse_name");
    _pulse_num           = _cfg.getIntParameter("Condition", "pulse_number");

    double magBx = _cfg.getDoubleParameter("Condition", "magnetic_fieldX");
    double magBy = _cfg.getDoubleParameter("Condition", "magnetic_fieldY");
    double magBz = _cfg.getDoubleParameter("Condition", "magnetic_fieldZ");
    _magB << magBx << magBy << magBz; 

}/*}}}*/

vec EnsembleCCE::cluster_evolution(int cce_order, int index)
{
    vector<cSPIN> spin_list = _my_clusters.getCluster(cce_order, index);

    int spin_up = 0, spin_down = 1;
    Hamiltonian hami0 = create_spin_hamiltonian(_center_spin, spin_up, spin_list);
    Hamiltonian hami1 = create_spin_hamiltonian(_center_spin, spin_down, spin_list);

    Liouvillian lv1 = create_spin_liouvillian(hami0, hami1);
    Liouvillian lv2 = create_spin_liouvillian(hami1, hami0);
    
    vector<QuantumOperator> lv_list = riffle((QuantumOperator) lv1, (QuantumOperator) lv2, _pulse_num);
    vector<double> time_segment = Pulse_Interval(_pulse_name, _pulse_num);

    DensityOperator ds = create_spin_density_state(spin_list);

    PiecewiseFullMatrixVectorEvolution kernel(lv_list, time_segment, ds);
    kernel.setTimeSequence( _t0, _t1, _nTime);

    ClusterCoherenceEvolution dynamics(&kernel);
    dynamics.run();

    return dynamics.calc_obs();
}

Hamiltonian EnsembleCCE::create_spin_hamiltonian(const cSPIN& espin, const int spin_state, const vector<cSPIN>& spin_list)
{
    SpinDipolarInteraction dip(spin_list);

    SpinZeemanInteraction zee(spin_list, _magB);

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

Liouvillian EnsembleCCE::create_spin_liouvillian(const Hamiltonian& hami0, const Hamiltonian hami1)
{
    Liouvillian lv0(hami0, SHARP);
    Liouvillian lv1(hami1, FLAT);
    Liouvillian lv = lv0 - lv1;
    return lv;
}

DensityOperator EnsembleCCE::create_spin_density_state(const vector<cSPIN>& spin_list)
{
    vec pol = zeros<vec>(3);
    SpinPolarization p(spin_list, pol);

    DensityOperator ds(spin_list);
    ds.addStateComponent(p);
    ds.make();
    ds.makeVector();
    return ds;
}
//}}}
////////////////////////////////////////////////////////////////////////////////