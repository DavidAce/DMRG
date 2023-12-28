#pragma once
#include "fwd_decl.h"
#include <memory>
#include <unsupported/Eigen/CXX11/ThreadPool>
namespace tenx {
    namespace threads {
        inline unsigned int num_threads = 1;
#if defined(EIGEN_USE_THREADS)
        extern std::unique_ptr<Eigen::ThreadPool>       tp;
        extern std::unique_ptr<Eigen::ThreadPoolDevice> dev;
        template<typename T>
        extern void setNumThreads(T num);
        Eigen::ThreadPoolDevice &getDevice();
#else
        extern std::unique_ptr<Eigen::DefaultDevice> dev;
        void                                         setNumThreads([[maybe_unused]] int num);
        Eigen::DefaultDevice                        &getDevice();
#endif

    }
}