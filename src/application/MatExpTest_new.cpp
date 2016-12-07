#include "include/app/app.h"
#include "include/math/MatExp.h"
#include <complex>
#include "expv/include/expv.h"

cx_mat MAT;
cx_vec VEC;
vec    TIME_LIST;
SumKronProd SKP;
cx_double II = cx_double(0.0,1.0);

void test_corr(const int& seed, ofstream& file);
void prepare_data(string filename, int i);
cx_mat test_arma_mat();
cx_mat test_pade_mat();
cx_mat test_large_mat();
cx_mat test_large_mat_sparse();
cx_mat test_very_large_mat_CPU();
cx_mat test_very_large_mat_GPU();
cx_mat time_cal(cx_mat (*func)(), ofstream& file);

int  main(int argc, char* argv[])
{
    string filename_expm = "./dat/output/expm.txt";
    ofstream expm_file(filename_expm.c_str());
    if(!expm_file) assert(0);
    prepare_data(filename,4);
    test_corr(5,expm);


    expm_file << "|" << "        | " << "   ArmaExpM   | " << "   PadeExpM   | " << "   VecExpM    | " << "  SpVecExpM   | " << "InexplicitCPU | " << "InexplicitGPU | " << endl;
    int maxSpin=8;
    string filename = "./dat/input/C13Bath/RoyCoord.xyz";
    for (int i=3; i<maxSpin; ++i)
    {   
        cout << "prepare data for spin-" << i << " start..."<< endl;
        prepare_data(filename,i);
        expm_file << "|" << " " << i << " spin | ";
    	
        cx_mat res_arma = time_cal(&test_arma_mat, expm_file);
        cx_mat res_pade = time_cal(&test_pade_mat, expm_file);
        cx_mat res_large = time_cal(&test_large_mat, expm_file);
        cx_mat res_large_sp = time_cal(&test_large_mat_sparse, expm_file);
        cx_mat res_very_large_CPU = time_cal(&test_very_large_mat_CPU, expm_file);
        cx_mat res_very_large_GPU = time_cal(&test_very_large_mat_GPU, expm_file);
    
        cout << "calculate for spin-" << i  << " done." << endl;
        expm_file << endl;
    }
    expm_file.close();
    return 0;
}

void test_corr(const int& seed, ofstream& expm) 
/*{{{*/{
    arma_rng::set_seed(seed);
    cx_mat sigma_n = randu<cx_mat>(2,2);
    sigma_n=( sigma_n+sigma_n.t() );//Hremitian
    cx_mat sigma_x(2,2);
    sigma_x << 0.0 << 1.0 << endr
            << 1.0 << 0.0;
    cx_mat sigma_y(2,2);
    sigma_y << 0.0 << -II << endr
            << II  << 0.0;
    cx_mat sigma_z(2,2);
    sigma_z << 1.0 << 0.0 << endr
            << 0.0 << -1.0;
    cx_mat iden_i(2,2);
    iden_i << 1.0 << 0.0 << endr
           << 0.0 << 1.0;

    double cons1 = sqrt( 4*pow(real(sigma_n(0,1)),2) + 4*pow(imag(sigma_n(0,1)),2)+pow( real(sigma_n(0,0)-sigma_n(1,1)),2) ); 
    cx_double cons2 = exp( II*(sigma_n(0,0) + sigma_n(1,1))/2 );
    cx_mat exp_sigma_n = ( cons2*cos(cons1/2.0) )*iden_i + II*( ((sigma_n(0,0) - sigma_n(1,1))*cons2*sin(cons1/2.0))/cons1 )*sigma_z
                         +II*( (2.0*real(sigma_n(0,1))*cons2*sin(cons1/2.0))/cons1 )*sigma_x -II*( (2.0*imag(sigma_n(0,1))*cons2*sin(cons1/2.0))/cons1 )*sigma_y;
    //cout << "sigma_n:" << endl << sigma_n << endl;
    //cout << "exp_sigma_n:" << endl << std::setprecision(9) << exp_sigma_n << endl;
    
    MatExp Arma_exp(sigma_n, II, MatExp::ArmadilloExpMat);
    Arma_exp.run();
    cx_mat exp_Arma = Arma_exp.getResultMatrix();
    
    MatExp Pade_exp(sigma_n, II, MatExp::PadeApproximation);
    Pade_exp.run();
    cx_mat exp_Pade = Pade_exp.getResultMatrix();
    //cout << "exp_arma:" << endl << std::setprecision(9) << exp_Arma << endl;
    //cout << "exp_pade:" << endl << std::setprecision(9) << exp_Pade << endl;

    expm << "diff between exact, exp_arma, exp_pade for a two-dimensional random complex matrix." << endl;
    expm << "                   diff" << endl;
    expm << " norm(Exact-Arma): " << std::setprecision(9) << norm(exp_sigma_n-exp_Arma) << endl;
    expm << " norm(Exact-Pade): " << std::setprecision(9) << norm(exp_sigma_n-exp_Pade) << endl;
    
    cx_mat res_Arma = test_arma_mat()*VEC;
    cx_mat resVecExp = test_large_mat();
    cx_mat resSpVecExp = test_large_mat_sparse();
    cx_mat resInplicitCPU = test_very_large_mat_CPU();
    cx_mat resInplicitGPU = test_very_large_mat_GPU();
    expm << "diff between arma,VecExp,SpVecExp,InexplicitCPU,InexplicitGPU for a random 256-dimensional initial vector, the matrix includes dipolar and zeeman interaction." << endl;
    expm << "                          diff" << endl;
    expm << " norm(Arma-VecExp):       " << std::setprecision(9) << norm(res_Arma-resVecExp) << endl;
    expm << " norm(Arma-SpVecExp):     " << std::setprecision(9) << norm(res_Arma-resSpVecExp) << endl;
    expm << " norm(Arma-InplicitCPU):  " << std::setprecision(9) << norm(res_Arma-resInplicitCPU) << endl;
    expm << " norm(Arma-InplicitGPU):  " << std::setprecision(9) << norm(res_Arma-resInplicitGPU) << endl;
    expm << endl;
    expm << endl;
    expm << endl;
}/*}}}*/;

void prepare_data(string filename, int i)
{/*{{{*/
    cSpinSourceFromFile spin_file(filename);
    cSpinCollection spins(&spin_file);
    spins.make();
    
    uvec idx(i);
    for(int j=0; j<i; ++j)
        idx(j)=j;
    cClusterIndex clst(idx);
    vector<cSPIN> sl = spins.getSpinList(clst);
    vec magB;
    magB << 0.1 << 0.1 << 0.1;

    SpinDipolarInteraction dip(sl);
    SpinZeemanInteraction zee(sl, magB);
    Hamiltonian hami(sl);
    hami.addInteraction(dip);
    hami.make();

    Liouvillian lv0(hami);

    double rate = 1.0;
    vec axis; axis << 1.0 << 1.0 << 1.0;
    SpinDephasing dephasing(sl, rate, normalise(axis));
    LiouvilleSpaceOperator dephaseOperator(sl);
    dephaseOperator.addInteraction(dephasing);
    dephaseOperator.make();

    QuantumOperator lv = lv0 + dephaseOperator;

    MAT = lv.getMatrix(); 
    cout << "hamiltonian mat generated." <<endl;

    SKP =lv.getKronProdForm();

    int dim = MAT.n_cols;
    PureState psi(dim);
    cx_vec res = randu<cx_vec>(dim,1);
    psi.setVector(res);
    VEC = psi.getVector();
    cout << "vector generated." << endl;

    TIME_LIST = linspace<vec>(0.1, 0.1, 1);
}/*}}}*/

cx_mat test_arma_mat()
{ /*{{{*/
    MatExp expM(MAT, -II*TIME_LIST(9), MatExp::ArmadilloExpMat);
    expM.run();
    return expM.getResultMatrix();
}/*}}}*/

cx_mat test_pade_mat()
{ /*{{{*/
    MatExp expM(MAT, -II*TIME_LIST(9), MatExp::PadeApproximation);
    expM.run();
    return expM.getResultMatrix();
}/*}}}*/

cx_mat test_large_mat()
{/*{{{*/
    MatExpVector expM(MAT, VEC, -II, TIME_LIST);
    return expM.run();
}/*}}}*/

cx_mat test_large_mat_sparse()
{/*{{{*/
    sp_cx_mat mat_sparse = sp_cx_mat(MAT);
    MatExpVector expM(mat_sparse, VEC, -II, TIME_LIST);
    return expM.run();
}/*}}}*/

cx_mat test_very_large_mat_CPU()
{/*{{{*/
    MatExpVector expM(SKP, VEC, TIME_LIST, MatExpVector::Inexplicit);  
    return expM.run();
}/*}}}*/

cx_mat test_very_large_mat_GPU()
{/*{{{*/
    MatExpVector expM(SKP, VEC, TIME_LIST, MatExpVector::InexplicitGPU);  
    return expM.run();
}/*}}}*/

cx_mat time_cal(cx_mat (*func)(), ofstream& file)
{/*{{{*/
    cx_mat res;
    clock_t start,end;
    int max=100;
    double time=0;
    for(int j=0; j<max; ++j)
	{   
        start = clock();
        res = func();
        end = clock();
        if(((double)(end-start)/CLOCKS_PER_SEC) > 1.0)
            max=1;
        time += (double)(end-start)/CLOCKS_PER_SEC;
    }
    file << std::scientific <<std::setprecision(7) << time/max  << " | ";
    return res;
}/*}}}*/
