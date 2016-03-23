#include "include/app/cce.h"
#include "include/misc/xmlreader.h"
#include <cstdlib>
#include "include/misc/lattice.h"
#include "include/spin/SpinClusterFromLattice.h"

_INITIALIZE_EASYLOGGINGPP

string PROJECT_PATH;
string LOG_PATH;
string INPUT_PATH;
string OUTPUT_PATH;
string CONFIG_PATH;
string DEBUG_PATH;

cSPINDATA SPIN_DATABASE=cSPINDATA();
ConfigXML set_parameters(const string& xml_file_name);

int  main(int argc, char* argv[])
{
    int dim = 2;
    vec base1, base2;
    base1 << 1.0 << 0.0 << 0.0; base2 << 0.0 << 1.0 << 0.0;
    vector<vec> bases; bases.push_back(base1); bases.push_back(base2);
    int atom_num = 2;
    vec coord1, coord2, coord3;
    coord1 << 0.0 << 0.0 << 0.0;
    coord2 << 0.5*3.57 << 0.5*3.57 << 0.0;
    coord3 << 0.15 << 0.1 << 0.0;
    vector<vec> pos; pos.push_back(coord1); pos.push_back(coord2);// pos.push_back(coord3);
    vector<double> latt_const;
    latt_const.push_back(3.57); latt_const.push_back(3.57);
    vector<string> iso; 
    iso.push_back("13C"); iso.push_back("13C");// iso.push_back("13C"); 
    Lattice latt(dim, bases, latt_const, atom_num, pos, iso);
    imat range; range << -20 << 21 << endr << -20 << 21;
    latt.setRange(range);

    cSpinSourceFromLattice spin_on_lattice(latt, range);
    cSpinCollection _bath_spins(&spin_on_lattice);
    _bath_spins.make();

    int maxOrder = 6;
    sp_mat c=_bath_spins.getConnectionMatrix(4.0);
    imat root_range; root_range << -6 << 7 << endr << -6 << 7;
    cUniformBathOnLattice bath_on_lattice(c, maxOrder, _bath_spins, latt, root_range);
    cSpinCluster _spin_clusters(_bath_spins, &bath_on_lattice);
    _spin_clusters.make();
    
    _spin_clusters.enable_sub_cluster_position();
    cout << _spin_clusters << endl;

    _spin_clusters.MPI_partition(200);
    
    return 0;



    ConfigXML cfg = set_parameters("SingleSampleCCE.xml");

    string log_file = LOG_PATH + cfg.getStringParameter("Data", "log_file");
    _START_EASYLOGGINGPP(argc, argv);
    easyloggingpp::Configurations confFromFile(log_file.c_str());
    easyloggingpp::Loggers::reconfigureAllLoggers(confFromFile);

    
    int worker_num(0), my_rank(0);
    int mpi_status = MPI_Init(&argc, &argv);
    assert (mpi_status == MPI_SUCCESS);

    MPI_Comm_size(MPI_COMM_WORLD, &worker_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    LOG(INFO) << "my_rank = " << my_rank << "  vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv Program begins vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv"; 

    // create defect center
    double x = cfg.getDoubleParameter("CenterSpin",  "coordinate_x");
    double y = cfg.getDoubleParameter("CenterSpin",  "coordinate_y");
    double z = cfg.getDoubleParameter("CenterSpin",  "coordinate_z");
    vec coord; coord << x << y << z;
    NVCenter nv(NVCenter::N14, coord);
    
    double magBx = cfg.getDoubleParameter("Condition",  "magnetic_fieldX");
    double magBy = cfg.getDoubleParameter("Condition",  "magnetic_fieldY");
    double magBz = cfg.getDoubleParameter("Condition",  "magnetic_fieldZ");
    nv.set_magB(magBx, magBy, magBz);
    nv.make_espin_hamiltonian();

    // CCE
    SingleSampleCCE sol(my_rank, worker_num, &nv, cfg);
    sol.run();

    LOG(INFO) << "my_rank = " << my_rank << "  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ Program ends ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"; 

    mpi_status = MPI_Finalize();
    assert (mpi_status == MPI_SUCCESS);
}

ConfigXML set_parameters(const string& xml_file_name)
{/*{{{*/
    char *env_path = std::getenv("CCE_PROJ_PATH");
    if(env_path!=NULL)
        PROJECT_PATH = env_path;
    else
    {
        char pwd[500];
        getcwd(pwd, sizeof(pwd));
        PROJECT_PATH = pwd;
    }

    LOG_PATH    = PROJECT_PATH + "/dat/log/";
    INPUT_PATH  = PROJECT_PATH + "/dat/input/";
    OUTPUT_PATH = PROJECT_PATH + "/dat/output/";
    CONFIG_PATH = PROJECT_PATH + "/dat/config/";
    DEBUG_PATH  = PROJECT_PATH = "/dat/debug/";

    ConfigXML cfg( CONFIG_PATH+xml_file_name );
    return cfg;
}/*}}}*/

