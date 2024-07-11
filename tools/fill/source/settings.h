#pragma once
#include <string>
namespace settings {
    inline std::string prefix;
    inline std::string filter;
    inline std::string outdir;
    inline std::string output_dirname = "output";
    inline std::string config_dirname = "config";
    inline std::string config_default = CMAKE_SOURCE_DIR "/input/input-xdmrg-ising-majorana.cfg";
}