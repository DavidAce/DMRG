/*! \file */


#include <sim_parameters/nmspc_sim_settings.h>
#include <sim_parameters/nmspc_model.h>
#include <algorithms/class_algorithm_launcher.h>
#include <io/class_settings_reader.h>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <h5pp/h5pp.h>
#include <experimental/filesystem>


#ifdef OpenBLAS_AVAILABLE
#include <cblas.h>
#endif


#ifdef OpenMP_AVAILABLE
#include <omp.h>
#endif

#ifdef MKL_AVAILABLE
#define MKL_Complex8 std::complex<float>
#define MKL_Complex16 std::complex<double>
#include <mkl_service.h>
#include <mkl.h>
#endif



/*!
    \brief  Main function. Sets simulation parameters and excecutes the desired algorithms.
    \return an integer 0 upon exit success
*/

int main(int argc, char* argv[]) {

    spdlog::set_default_logger(spdlog::stdout_color_mt("DMRG"));
    spdlog::set_pattern("[%Y-%m-%d %H:%M:%S][%n]%^[%=8l]%$ %v");
    spdlog::set_level(spdlog::level::trace);

    int openblas_num_threads = 1;
    #ifdef OpenBLAS_AVAILABLE
        openblas_set_num_threads(openblas_num_threads);
    std::cout << "OpenBLAS compiled with mode " << openblas_get_parallel()
              << " for target " << openblas_get_corename()
              << " with config " << openblas_get_config()
              << ". Running with " << openblas_get_num_threads() << " thread(s)" << std::endl;
    #endif

    #ifdef OpenMP_AVAILABLE
        Eigen::initParallel();
//        omp_set_num_threads(OpenMP_NUM_THREADS);
//        omp_set_dynamic(0);
//        Eigen::setNbThreads(OpenMP_NUM_THREADS);
//        Eigen::setNbThreads(0);
        spdlog::info("Using Eigen  with {} threads",Eigen::nbThreads( )  );
        spdlog::info("Using OpenMP with {} threads",omp_get_max_threads());
        #ifdef MKL_AVAILABLE
//            mkl_set_num_threads(OpenMP_NUM_THREADS);
            spdlog::info("Using Intel MKL with {} threads", mkl_get_max_threads());
        #endif
    #endif





    //Print all given parameters
    //Load input and output files from command line. If none were given use defaults.
    //Normally an output filename is given in the input file. But it can also be given from command line.
    std::string inputfile  = "input.cfg";
    std::string outputfile = "output.h5";
    int seed = -1; //Only accept non-negative seeds
    for (int i=0; i < argc; i++){
        std::string arg_string = std::string(argv[i]);
        spdlog::info("Input argument {} : {}",i,arg_string);
        if (arg_string.find(".cfg") != std::string::npos) {inputfile  = arg_string;}
        if (arg_string.find(".h5")  != std::string::npos) {outputfile = arg_string;}
        if (arg_string.find_first_not_of( "0123456789" ) == std::string::npos){seed = std::stoi(arg_string);}
    }
    class_settings_reader indata(inputfile);
    if(indata.found_file){
        settings::load_from_file(indata);
    }else{
        try{
            auto h5ppFile = std::make_shared<h5pp::File> (outputfile,h5pp::AccessMode::READONLY,h5pp::CreateMode::OPEN);
            spdlog::info("Loading settings from existing file [{}]", h5ppFile->getFilePath());
            settings::load_from_hdf5(*h5ppFile);
        }catch(std::exception &ex){
            spdlog::info("Couldn't find an inputfile or previous outputfile to load settings: {}", outputfile,ex.what() );
            spdlog::info("Running defaults");
        }
    }
    if (outputfile != "output.h5"){
        spdlog::info("Replacing output filename {} --> {}",settings::hdf5::output_filename, outputfile);
        settings::hdf5::output_filename = outputfile;
    }
    if (seed >= 0){
        spdlog::info("Replacing seed {} --> {}",settings::model::seed, seed);
        settings::model::seed = seed;
        //Append the seed to the output filename
        namespace fs = std::experimental::filesystem;
        fs::path oldFileName = settings::hdf5::output_filename;
        fs::path newFileName = settings::hdf5::output_filename;
        newFileName.replace_filename(oldFileName.stem().string() + "_" + std::to_string(seed) + oldFileName.extension().string() );
        settings::hdf5::output_filename = newFileName.string();
        spdlog::info("Appending seed to output filename: [{}] --> [{}]",oldFileName.string(), newFileName.string());
    }



    //Initialize the algorithm class
    //This class stores simulation data automatically to a file specified in the input file
    class_algorithm_launcher launcher;

    //Run the algorithms
    launcher.run_algorithms();


//    std::string tar_filename = fs::path(launcher.hdf5_path).filename().string() + ".tar.gz";
//    fs::path tar_path        = fs::path(launcher.hdf5_path).parent_path().string() + "/"+ tar_filename;
//    gzFile tar_file = gzopen(tar_path.c_str(), "wb");
//    if(tar_file == nullptr)
//    {
//        throw std::runtime_error("Could not open gzip file for writing: " );
//    }

    /* Open file to compress */

//    std::ifstream h5file (launcher.hdf5_path);
//    std::vector<char> buffer (1024,0); //reads only the first 1024 bytes
//
//    while(h5file.read(buffer.data(), buffer.size())) {
//        std::streamsize bytes_loaded=h5file.gcount();
//        int bytes_written = gzwrite(tar_file, buffer.data(), bytes_loaded);
//        if (bytes_written == 0)
//        {
//            spdlog::warn("No bytes written!");
//            break;
//        }
//        ///do with buffer
//    }
//    gzclose(tar_file);
//    h5file.close();





    return 0;
}



