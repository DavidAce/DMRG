//
// Created by david on 2019-12-09.
//

#ifndef DMRG_NMSPC_FILESYSTEM_H
#define DMRG_NMSPC_FILESYSTEM_H

// Include filesystem or experimental/filesystem
#if __has_include(<filesystem>)
#include <filesystem>
    namespace tools{
        namespace fs =  std::filesystem;
    }
#elif __has_include(<experimental/filesystem>)
#include <experimental/filesystem>
namespace tools {
    namespace fs = std::experimental::filesystem;
}
#else
#error Could not find <filesystem> or <experimental/filesystem>
#endif


#endif //DMRG_NMSPC_FILESYSTEM_H
