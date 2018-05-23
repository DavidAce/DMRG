/*! \file */


#include <sim_parameters/nmspc_sim_settings.h>
#include <sim_parameters/nmspc_model.h>
#include <algorithms/class_algorithm_launcher.h>
#include <gitversion.h>
#include <IO/class_file_reader.h>

/*!
    \brief  Main function. Sets simulation parameters and excecutes the desired algorithms.
    \return an integer 0 upon exit success
*/

int main(int argc, char* argv[]) {
    // Print current Git status
    std::cout << "Git Branch: " + GIT::BRANCH +
            " | Commit hash: "  + GIT::COMMIT_HASH +
            " | Revision: "     + GIT::REVISION << std::endl << std::flush;

    //Print all given parameters
    //Load input and output files from command line. If none were given use defaults.
    //Normally an output filename is given in the input file. But it can also be given from command line.
    std::string inputfile  = "input.cfg";
    std::string outputfile;
    bool outputfile_given  = false;
    for (int i=0; i < argc; i++){
        std::string arg_string = std::string(argv[i]);
        std::cout <<"Input argument [" << i << "] : " << arg_string << std::endl;
        if (arg_string.find(".cfg") != std::string::npos) {inputfile  = arg_string;}
        if (arg_string.find(".h5")  != std::string::npos) {outputfile = arg_string;outputfile_given=true;}
    }
    class_file_reader indata(inputfile);
    settings::load_from_file(indata);

    //If an output filename was given explicitly, overwrite the default , if it was given explicitly in command line.
    settings::hdf5::output_filename = outputfile_given ? outputfile : settings::hdf5::output_filename;

    //Initialize the algorithm class
    //This class stores simulationdata automatically to a file specified in the input file
    class_algorithm_launcher launcher;

    //Run the algorithms
    launcher.run_algorithms();


    return 0;
}



