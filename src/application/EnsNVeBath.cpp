#include "include/app/app.h"
#include "include/app/ensemble_cce.h"

_INITIALIZE_EASYLOGGINGPP

namespace po = boost::program_options;

po::variables_map ParseCommandLineOptions(int argc, char* argv[]);
void set_parameters(const string& xml_file_name);
NVCenter create_defect_center(const po::variables_map& para);
cSpinSourceUniformRandom create_spin_source(const po::variables_map& para);
cDepthFirstPathTracing create_spin_cluster_algrithm(const po::variables_map& para, const cSpinCollection& bath_spins);

int  main(int argc, char* argv[])
{
    cout << argc << endl;
    string a("");
    for(int i=1; i<argc; ++i)
    {
        cout << i << "\t";
        cout << string(argv[i]) << endl;
        a += argv[i];
    }
    cout << a << endl;
    exit(0);
    po::variables_map para = ParseCommandLineOptions(argc, argv);

    ////////////////////////////////////////////////////////////////////////////////
    //{{{  MPI Preparation
    int worker_num(0), my_rank(0);
    int mpi_status = MPI_Init(&argc, &argv);
    assert (mpi_status == MPI_SUCCESS);

    MPI_Comm_size(MPI_COMM_WORLD, &worker_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    //}}}
    ////////////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////////////////
    //{{{ LOGGING 
    string log_file = LOG_PATH +  para["logfile"].as<string>();
    _START_EASYLOGGINGPP(argc, argv);
    easyloggingpp::Configurations confFromFile(log_file.c_str());
    easyloggingpp::Loggers::reconfigureAllLoggers(confFromFile);
    //}}}
    ////////////////////////////////////////////////////////////////////////////////

    LOG(INFO) << "my_rank = " << my_rank << "  vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv Program begins vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv"; 

    EnsembleCCE sol(my_rank, worker_num, para);

    // Step 1: make a defect center
    NVCenter nv = create_defect_center(para);  
    sol.set_defect_center(&nv);

    // Step 2: make bath spins 
    cSpinSourceUniformRandom spinUR = create_spin_source(para);
    sol.set_bath_spin(&spinUR);
    
    
    // Step 3: make clusters
    cSpinCollection bath_spins = sol.getSpinCollecion();
    cDepthFirstPathTracing   dfpt = create_spin_cluster_algrithm(para, bath_spins);
    sol.set_bath_cluster(&dfpt);

    // Step 4: run_each_cluster 
    sol.run_each_clusters();

    // Step 5: post_treatment
    sol.post_treatment();

    LOG(INFO) << "my_rank = " << my_rank << "  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ Program ends ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"; 

    ////////////////////////////////////////////////////////////////////////////////
    //{{{ MPI Finializing
    mpi_status = MPI_Finalize();
    assert (mpi_status == MPI_SUCCESS);
    //}}}
    ////////////////////////////////////////////////////////////////////////////////
}


po::variables_map ParseCommandLineOptions(int argc, char* argv[])
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
    DEBUG_PATH  = PROJECT_PATH + "/dat/debug/";


    po::variables_map para;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "Print help message")
        ("input,i",          po::value<string>()->default_value("None"),                 "Input file name")
        ("output,o",         po::value<string>()->default_value("EnsNVeBath"),           "Output .mat file of results")
        ("logfile,l",        po::value<string>()->default_value("EnsembleCCE.conf"),     "Config. file of logging")
        
        ("position,P",       po::value<string>()->default_value("0.0 0.0 0.0"),          "Central spin position")
        ("state0,a",         po::value<int>()->default_value(0),                         "Central spin state index - a")
        ("state1,b",         po::value<int>()->default_value(1),                         "Central spin state index - b")

        ("cce,c",            po::value<int>()->default_value(3),                         "CCE order")
        ("cutoff,d",         po::value<double>()->default_value(300.0),                    "Cut-off distance of bath spins")
        ("polarization,z",   po::value<string>()->default_value("0.0 0.0 0.0"),          "magnetic field vector in Tesla")
        ("dephasing_rate,r", po::value<double>()->default_value(0.0),                    "dephasing rate of bath spins")
        ("dephasing_axis,x", po::value<string>()->default_value("1.0 1.0 1.0"),          "dephasing axis of bath spins")
        ("range,g",          po::value<double>()->default_value(1000.0),               "bath range [-g, g]")
        ("number,N",         po::value<int>()->default_value(100),                     "bath spin number")
        ("seed,D",           po::value<int>()->default_value(1),                       "bath seed")
        ("isotope",          po::value<string>()->default_value("E"),                    "bath spin isotope")

        ("nTime,n",          po::value<int>()->default_value(101),                       "Number of time points")
        ("start,s",          po::value<double>()->default_value(0.0),                    "Start time (in unit of sec.)")
        ("finish,f",         po::value<double>()->default_value(0.002),                  "Finish time (in unit of sec.)")
        ("magnetic_field,B", po::value<string>()->default_value("0.1 0.1 0.1"),          "magnetic field vector in Tesla")
        ("pulse,p",          po::value<string>()->default_value("CPMG"),                 "Pulse name")
        ("pulse_num,m",      po::value<int>()->default_value(1),                         "Pulse number")
        ;

    po::store(parse_command_line(argc, argv, desc), para);
    ifstream ifs("config.cfg");
    if (ifs!=NULL)
        po::store(parse_config_file(ifs,desc),para);
    else
        cout << "No configure file found!" << endl;
    po::notify(para);    

    if (para.count("help")) {
        cout << desc;
        exit(0);
    }
    return para;
}/*}}}*/

NVCenter create_defect_center(const po::variables_map& para)
{/*{{{*/
    vec coord( para["position"].as<string>() );
    vec magB( para["magnetic_field"].as<string>() );

    NVCenter nv(NVCenter::N14, coord);
    nv.set_magB(magB);
    nv.make_espin_hamiltonian();

    return nv;
}/*}}}*/

cSpinSourceUniformRandom create_spin_source(const po::variables_map& para)
{/*{{{*/
    double   range = para["range"].as<double>();
    int        num = para["number"].as<int>();
    int       seed = para["seed"].as<int>();
    string isotope = para["isotope"].as<string>();

    cSpinSourceUniformRandom spinUR(range, num, isotope, seed);
    return spinUR;
}/*}}}*/

cDepthFirstPathTracing create_spin_cluster_algrithm(const po::variables_map& para, const cSpinCollection& bath_spins)
{/*{{{*/
    double cut_off_dist = para["cutoff"].as<double>();
    int    max_order    = para["cce"].as<int>();
    sp_mat c = bath_spins.getConnectionMatrix(cut_off_dist);
    cDepthFirstPathTracing dfpt(c, max_order);
    return dfpt;
}/*}}}*/
