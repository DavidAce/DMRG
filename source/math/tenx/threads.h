#pragma once
#include "fwd_decl.h"
#include <memory>
#include <type_traits>
#include <unsupported/Eigen/CXX11/ThreadPool>
#include <unsupported/Eigen/CXX11/src/Tensor/TensorDeviceDefault.h>
namespace tenx {
    namespace threads {
#if defined(EIGEN_USE_THREADS)
        namespace internal {
            struct ThreadPoolWrapper {
                private:
                std::unique_ptr<Eigen::ThreadPool> tp;
                public:
                std::unique_ptr<Eigen::ThreadPoolDevice> dev;

                ThreadPoolWrapper(int nt);
            };
            inline unsigned int num_threads = 1;
            extern std::unique_ptr<ThreadPoolWrapper> singleThreadWrapper;
            extern std::unique_ptr<ThreadPoolWrapper> multiThreadWrapper;
        }

        template<typename T>
        requires std::is_integral_v<T>
        extern void                 setNumThreads(T num) noexcept;
        extern int                  getNumThreads() noexcept;
        internal::ThreadPoolWrapper &get() noexcept;
#else
        extern std::unique_ptr<Eigen::DefaultDevice> dev;
        void                                         setNumThreads([[maybe_unused]] int num);
        Eigen::DefaultDevice                        &getDevice();
#endif

    }
}