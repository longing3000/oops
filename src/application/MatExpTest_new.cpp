#include "include/app/app.h"
#include "include/math/MatExp.h"
#include <complex>
#include "expv/include/expv.h"

cx_mat MAT;
cx_vec VEC;
vec    TIME_LIST;
SumKronProd SKP;
cx_double PREFACTOR;

void prepare_data(string filename);
cx_mat test_arma_mat();
cx_mat test_pade_mat();
cx_mat test_large_mat();
cx_mat test_large_mat_sparse();
cx_mat test_very_large_mat_CPU();
cx_mat test_very_large_mat_GPU();

int  main(int argc, char* argv[])
{
    string filename_expm = "./dat/output/expm.txt";
    ofstream expm_file(filename_expm.c_str());
    if(!expm_file) assert(0);
    expm_file << "|" << "        | " << "   ArmaExpM   | " << "   PadeExpM   ||||| " << "   VecExpM    | " << "  SpVecExpM   | " << "InexplicitCPU | " << "InexplicitGPU | " << endl;
    
    int maxSpin=8;
    string filepath = "./dat/input/RoyCoord.xyz";
    clock_t start,end;
    for (int i=3; i<maxSpin; ++i)
    {   
        stringstream char_i;
        char_i << i;
        string filename=filepath + char_i.str();
        cout << "prepare data for spin-" << char_i.str() << " start"<< endl;
        prepare_data(filename);
        expm_file << "|" << " " << i << " spin | ";
	int max=10;
	if(i<5)
       	    max=100;

	cx_mat res_arma;
        start = clock();
	for(int j=0; j<max; ++j)
	{
	    res_arma = test_arma_mat();
        }
	end = clock();
        expm_file << std::scientific <<std::setprecision(7) << ((float)(end-start)/CLOCKS_PER_SEC)/max  << " | ";

	cx_mat res_pade;
        start = clock();
	for(int j=0; j<max; ++j)
	{
	    res_pade = test_pade_mat();
        }
	end = clock();
        expm_file << std::scientific << std::setprecision(7) << ((float)(end-start)/CLOCKS_PER_SEC)/max  << " ||||| ";
	
	cx_mat res_large;
        start = clock();
	for(int j=0; j<max; ++j)
	{
	    res_large = test_large_mat();
        }
	end = clock();
        expm_file << std::scientific << std::setprecision(7) << ((float)(end-start)/CLOCKS_PER_SEC)/max  << " | ";

	cx_mat res_large_sp;
        start = clock();
	for(int j=0; j<max; ++j)
	{
            res_large_sp = test_large_mat_sparse();
 	}
        end = clock();
        expm_file << std::scientific << std::setprecision(7) << ((float)(end-start)/CLOCKS_PER_SEC)/max  << " | ";
	
	cx_mat res_very_large_CPU;
        start = clock();
	for(int j=0; j<max; ++j)
	{
            res_very_large_CPU = test_very_large_mat_CPU();
 	}
        end = clock();
        expm_file << std::scientific << std::setprecision(7) << ((float)(end-start)/CLOCKS_PER_SEC)/max  << " | ";

	cx_mat res_very_large_GPU;
        start = clock();
	for(int j=0; j<max; ++j)
	{
            res_very_large_GPU = test_very_large_mat_GPU();
 	}
        end = clock();
        expm_file << std::scientific << std::setprecision(7) << ((float)(end-start)/CLOCKS_PER_SEC)/max  << " | ";

        if(i==4)
        {
            cout << "the differences of two kind of ExpM for spin-4 are:" << endl;
            cout << "diff 0 = " << norm(res_pade - res_arma) << endl;
            cout << "the differences of four kind of VecExpM for spin-4 are:" << endl;
            cout << "diff 1 = " << norm(res_large_sp - res_large) << endl;
            cout << "diff 2 = " << norm(res_very_large_CPU - res_large) << endl;
            cout << "diff 3 = " << norm(res_very_large_GPU - res_large) << endl;
            cout << "diff 4 = " << norm(res_very_large_GPU - res_very_large_CPU) << endl;
        }

        cout << "calculate for spin-" << char_i.str()  << " done" << endl;
        expm_file << endl;
    }
    cout << "!!!!!pay attention to a little difference in function test_large_mat_sparse!!!!!" << endl;
    expm_file.close();
    return 0;
}

void prepare_data(string filename)
{/*{{{*/
    cSpinSourceFromFile spin_file(filename);
    cSpinCollection spins(&spin_file);
    spins.make();

    vector<cSPIN> sl = spins.getSpinList();
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

    PREFACTOR = cx_double(0.0, -1.0);
    MAT = lv.getMatrix(); 
    cout << "hamiltonian mat generated." <<endl;

    SKP =lv.getKronProdForm();

    int dim = MAT.n_cols;
    PureState psi(dim);
    psi.setComponent(0, 1.0);
    VEC = psi.getVector();
    cout << "vector generated." <<endl;

    TIME_LIST = linspace<vec>(0.01, 0.1, 10);
}/*}}}*/

cx_mat test_arma_mat()
{ /*{{{*/
    MatExp expM(MAT, PREFACTOR*TIME_LIST(9), MatExp::ArmadilloExpMat);
    expM.run();
    return expM.getResultMatrix();
}/*}}}*/

cx_mat test_pade_mat()
{ /*{{{*/
    MatExp expM(MAT, PREFACTOR*TIME_LIST(9), MatExp::PadeApproximation);
    expM.run();
    return expM.getResultMatrix();
}/*}}}*/

cx_mat test_large_mat()
{/*{{{*/
    MatExpVector expM(MAT, VEC, PREFACTOR, TIME_LIST);
    return expM.run();
}/*}}}*/

cx_mat test_large_mat_sparse()
{/*{{{*/
    sp_cx_mat mat_sparse = sp_cx_mat(MAT);
    MatExpVector expM(mat_sparse, VEC, PREFACTOR, TIME_LIST);
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
