/*! \file */


#include <simulation/nmspc_settings.h>
#include <general/nmspc_random_numbers.h>
#include <algorithms/class_algorithm_launcher.h>
#include <io/class_settings_reader.h>
#include <io/nmspc_logger.h>

#include <iostream>
#include <h5pp/h5pp.h>
#include <experimental/filesystem>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef OpenBLAS_AVAILABLE
#include <cblas.h>
#include <openblas_config.h>
#endif


#ifdef MKL_AVAILABLE
#define MKL_Complex8 std::complex<float>
#define MKL_Complex16 std::complex<double>
#include <mkl_service.h>
#include <mkl.h>
#endif

#include <thread>

/*!
    \brief  Main function. Sets simulation parameters and excecutes the desired algorithms.
    \return an integer 0 upon exit success
*/

int main(int argc, char* argv[]) {
    auto log = Logger::setLogger("DMRG",0);




    //print all given parameters
    //Load input and output files from command line. If none were given use defaults.
    //Normally an output filename is given in the input file. But it can also be given from command line.
    std::string inputfile  = "input.cfg";
    std::string outputfile = "output.h5";
    int seed_init = -1; //Only accept non-negative seeds
    int i = 0;
    std::vector<std::string> allArgs(argv+1, argv + argc);
    for (auto &arg_word : allArgs){
        std::istringstream iss(arg_word);
        std::string arg;
        while(iss >> arg){
            log->info("Input argument {} : {}",i++,arg);
            if (arg.find(".cfg") != std::string::npos) {inputfile  = arg;continue;}
            if (arg.find(".h5")  != std::string::npos) {outputfile = arg;continue;}
            if (arg.find_first_not_of( "0123456789" ) == std::string::npos and seed_init < 0){seed_init = std::stoi(arg); continue;}
        }
    }


    class_settings_reader indata(inputfile);
    if(indata.found_file){
        settings::load_from_file(indata);
    }else{
        try{
            auto h5ppFile = std::make_shared<h5pp::File> (outputfile,h5pp::AccessMode::READONLY,h5pp::CreateMode::OPEN);
            log->info("Loading settings from existing file [{}]", h5ppFile->getFilePath());
            settings::load_from_hdf5(*h5ppFile);
        }catch(std::exception &ex){
            log->info("Couldn't find an inputfile or previous outputfile to load settings: {}", outputfile,ex.what() );
            log->info("Running defaults");
        }
    }
    if (outputfile != "output.h5"){
        log->info("Replacing output filename {} --> {}",settings::hdf5::output_filename, outputfile);
        settings::hdf5::output_filename = outputfile;
    }
    if (seed_init >= 0){
        log->info("Replacing seed_init {} --> {}", settings::model::seed_init, seed_init);
        settings::model::seed_init = seed_init;
        //Append the seed_init to the output filename
        namespace fs = std::experimental::filesystem;
        fs::path oldFileName = settings::hdf5::output_filename;
        fs::path newFileName = settings::hdf5::output_filename;
        newFileName.replace_filename(oldFileName.stem().string() + "_" + std::to_string(seed_init) + oldFileName.extension().string() );
        settings::hdf5::output_filename = newFileName.string();
        log->info("Appending seed_init to output filename: [{}] --> [{}]",oldFileName.string(), newFileName.string());
    }


    //Set the number of threads to be used


    #ifdef _OPENMP
        if(settings::threading::num_threads_omp   <= 0) { settings::threading::num_threads_omp   = std::thread::hardware_concurrency(); }
        if(settings::threading::num_threads_eigen <= 0) { settings::threading::num_threads_eigen = std::thread::hardware_concurrency(); }
        if(settings::threading::num_threads_blas  <= 0) { settings::threading::num_threads_blas  = std::thread::hardware_concurrency(); }

        omp_set_num_threads(settings::threading::num_threads_omp);
        Eigen::setNbThreads(settings::threading::num_threads_eigen);
        log->info("Using Eigen  with {} threads",Eigen::nbThreads());
        log->info("Using OpenMP with {} threads",omp_get_max_threads());

        #ifdef OpenBLAS_AVAILABLE
                openblas_set_num_threads(settings::threading::num_threads_blas);
                std::cout << OPENBLAS_VERSION
                          << " compiled with parallel mode " << openblas_get_parallel()
                          << " for target " << openblas_get_corename()
                          << " with config " << openblas_get_config()
                          << " with multithread threshold " << OPENBLAS_GEMM_MULTITHREAD_THRESHOLD
                          << ". Running with " << openblas_get_num_threads() << " thread(s)" << std::endl;
        #endif

        #ifdef MKL_AVAILABLE
            if(settings::threading::num_threads_blas   <= 0){ settings::threading::num_threads_blas   = std::thread::hardware_concurrency(); }
            mkl_set_num_threads(settings::threading::num_threads_blas);
            log->info("Using Intel MKL with {} threads", mkl_get_max_threads());
        #endif

    #endif










    // Seed only this once (This also takes care of srand used by Eigen
    rn::seed(settings::model::seed_init);


    //Initialize the algorithm class
    //This class stores simulation data automatically to a file specified in the input file
    class_algorithm_launcher launcher;

    //Run the algorithms
    launcher.run_algorithms();

    return 0;
}



