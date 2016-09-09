#include "include/app/app.h"
#include "include/app/ensemble_cce.h"
#include "include/my_debug.h"

_INITIALIZE_EASYLOGGINGPP

namespace po = boost::program_options;

po::variables_map ParseCommandLineOptions(int argc, char* argv[]);
void set_parameters(const string& xml_file_name);
NVCenter create_defect_center(const po::variables_map& para);
cSpinSourceFromLattice create_spin_source(const po::variables_map& para);
cUniformBathOnLattice create_spin_cluster_algrithm(const po::variables_map& para, const cSpinCollection& bath_spins);
void output_cluster(const po::variables_map& para, const cSpinCluster& spin_clusters);

int  main(int argc, char* argv[])
{
    global_log_file.open((OUTPUT_PATH+"dat/output/log_file.txt").c_str(),std::ios_base::app);
    global_log_file << "program begins" << endl;

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
    
    global_start=clock();
    EnsembleCCE sol(my_rank, worker_num, para);

    // Step 1: make a defect center
    NVCenter nv = create_defect_center(para);  
    sol.set_defect_center(&nv);
    
    // Step 2: make bath spins 
    cSpinSourceFromLattice spin_from_latt = create_spin_source(para);
    sol.set_bath_spin(&spin_from_latt);

    // Step 3: make clusters
    cSpinCollection bath_spins = sol.getSpinCollecion();
    cUniformBathOnLattice bath_on_lattice = create_spin_cluster_algrithm(para, bath_spins);
    global_start=clock();
    sol.set_bath_cluster(&bath_on_lattice);
    global_end=clock();
    global_log_file << "generate all cluster up to order " << para["cce"].as<int>() << ", it consumes:" << (double)(global_end-global_start)/CLOCKS_PER_SEC << "s"  << endl;
    output_cluster(para,sol.getSpinClusters());

    // Step 4: run_each_cluster
    global_start=clock();
    sol.run_each_clusters();
    global_end=clock();
    global_log_file << "cluster evolution total time:" <<  (double)(global_end-global_start)/CLOCKS_PER_SEC << "s" << endl;

    // Step 5: post_treatment
    global_start=clock();
    sol.post_treatment();
    global_end=clock();
    global_log_file << "posttreatment total time:" << (double)(global_end-global_start)/CLOCKS_PER_SEC << "s" << endl;
    global_log_file << "program ends." << endl;
    global_log_file.close();

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

    ////////////////////////////////////////////////////////////////////////////////
    //{{{ record command line options
    string output_filename("EnsNVLattice2D");
    string command_opt("");
    for(int i=1; i<argc; ++i)
        command_opt += argv[i];
    cout << "Command line options recieved: " << command_opt << endl;
    output_filename += command_opt;
    //}}}
    ////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////
    //{{{ Set path
    char *env_path = std::getenv("CCE_PROJ_PATH");
    if(env_path!=NULL)
    {
        PROJECT_PATH = env_path;
        cout << "PROJECT_PATH got from environment varialbe: ";
    }
    else
    {
        char pwd[500];
        getcwd(pwd, sizeof(pwd));
        PROJECT_PATH = pwd;
        cout << "PROJECT_PATH set as present working directory (pwd): ";
    }
    cout << PROJECT_PATH << endl;

    LOG_PATH    = PROJECT_PATH + "/dat/log/";
    INPUT_PATH  = PROJECT_PATH + "/dat/input/";
    OUTPUT_PATH = PROJECT_PATH + "/dat/output/";
    CONFIG_PATH = PROJECT_PATH + "/dat/config/";
    DEBUG_PATH  = PROJECT_PATH + "/dat/debug/";
    //}}}
    ////////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////////
    //{{{ Pars e command line options
    po::variables_map para;
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "Print help message")
        ("initialize,I", "Initialize a dat folder")
 
        ("input,i",          po::value<string>()->default_value("None"),                 "Input file name (not used)")
        ("output,o",         po::value<string>()->default_value(output_filename),        "Output .mat file of results")
        ("logfile,l",        po::value<string>()->default_value("EnsembleCCE.conf"),     "Config. file of logging")
        
        ("position,P",       po::value<string>()->default_value("0.0 0.0 -20.0"),          "Central spin position")
        ("state0,a",         po::value<int>()->default_value(0),                         "Central spin state index - a")
        ("state1,b",         po::value<int>()->default_value(1),                         "Central spin state index - b")

        ("cce,c",            po::value<int>()->default_value(6),                         "CCE order")
        ("cutoff,d",         po::value<double>()->default_value(4.0),                    "Cut-off distance of bath spins")
        ("polarization,z",   po::value<string>()->default_value("0.0 0.0 0.0"),          "bath spin polarization")
        ("dephasing_rate,r", po::value<double>()->default_value(10000.0),                "dephasing rate of bath spins")
        ("dephasing_axis,x", po::value<string>()->default_value("1.0 1.0 1.0"),          "dephasing axis of bath spins")

        ("latt_const,L",     po::value<double>()->default_value(3.57),                   "lattice constant (default diamond)")
        ("isotope",          po::value<string>()->default_value("13C"),                  "bath spin isotope")
        ("range",            po::value<int>()->default_value(20),                        "full range of bath spins")
        ("root_range",       po::value<int>()->default_value(8),                         "root range of bath spins")

        ("nTime,n",          po::value<int>()->default_value(101),                       "Number of time points")
        ("start,s",          po::value<double>()->default_value(0.0),                    "Start time (in unit of sec.)")
        ("finish,f",         po::value<double>()->default_value(0.005),                  "Finish time (in unit of sec.)")
        ("magnetic_field,B", po::value<string>()->default_value("0.1 0.1 0.1"),          "magnetic field vector in Tesla")
        ("pulse,p",          po::value<string>()->default_value("CPMG"),                 "Pulse name")
        ("pulse_num,m",      po::value<int>()->default_value(1),                         "Pulse number")
        ;

    po::store(parse_command_line(argc, argv, desc), para);
    ifstream ifs("config.cfg");
    if (ifs!=NULL)
        po::store(parse_config_file(ifs,desc),para);
    else
        cout << "No configure file found." << endl;
    po::notify(para);    

    PrintVariableMap(para);

    if (para.count("help")) {
        cout << desc;
        exit(0);
    }
    if (para.count("initialize")) {
        string path = getenv("OOPS_PATH");
        string cmd = "cp -r " + path + "/src/dat_example dat" ;
        cout << cmd << endl;
        system( cmd.c_str() );
        cout << "default dat folder initialized." << endl;
        exit(0);
    }
    //}}}
    ////////////////////////////////////////////////////////////////////////////////
    return para;
}/*}}}*/

NVCenter create_defect_center(const po::variables_map& para)
{/*{{{*/
    vec coord( para["position"].as<string>() );
    vec magB( para["magnetic_field"].as<string>() );

    NVCenter nv(NVCenter::N14, coord);
    nv.set_magB(magB);
    nv.make_espin_hamiltonian();
    
    global_log_file << "generate NV center at (" << coord(0) << "   " << coord(1) << "  " << coord(2) << ") with magnetic field (" << magB(0) << "   " << magB(1) << "   " << magB(2) << ")." << endl;
    return nv;
}/*}}}*/

cSpinSourceFromLattice create_spin_source(const po::variables_map& para)
{/*{{{*/
    double lattice_const = para["latt_const"].as<double>();
    double cut_off= para["cutoff"].as<double>();
    string isotope = para["isotope"].as<string>();
    int range_i = para["range"].as<int>();

    TwoDimFaceCenterLattice latt(lattice_const, isotope);
    latt.setRange(range_i);

    cSpinSourceFromLattice spin_on_lattice(latt);
    return spin_on_lattice;
}/*}}}*/

cUniformBathOnLattice create_spin_cluster_algrithm(const po::variables_map& para, const cSpinCollection& bath_spins)
{/*{{{*/
    double lattice_const = para["latt_const"].as<double>();
    double cut_off       = para["cutoff"].as<double>();
    string isotope       = para["isotope"].as<string>();
    int    maxOrder      = para["cce"].as<int>();
    int    range_i       = para["range"].as<int>();
    int    root_range_i  = para["root_range"].as<int>();

    TwoDimFaceCenterLattice latt(lattice_const, isotope);
    latt.setRange(range_i);

    //cout coordnate
    //TwoDimFaceCenterLattice latt_coord(lattice_const,isotope);
    //latt_coord.setRange(root_range_i);
    //latt_coord.save_to_file(OUTPUT_PATH+"coord.xyz");
 
    global_log_file << "generate " << (2*root_range_i+1)*(2*root_range_i+1)*2 << " bath spins from lattice." << endl;

    sp_mat c=bath_spins.getConnectionMatrix(cut_off);
    cUniformBathOnLattice bath_on_lattice(c, maxOrder, bath_spins, latt, root_range_i);
    return  bath_on_lattice;
}/*}}}*/

void output_cluster(const po::variables_map& para, const cSpinCluster& spin_clusters)
{/*{{{*/
   int maxOrder = para["cce"].as<int>();
   global_log_file << "up to now";
   for(int i=0; i<maxOrder; ++i)
       global_log_file << ",generate " << spin_clusters.getClusterNum(i) << " CCE-" << i << " clusters";
   global_log_file << endl;
}/*}}}*/
