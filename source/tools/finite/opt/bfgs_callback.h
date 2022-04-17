

#include <ceres/iteration_callback.h>
#include <complex>
#include <memory>

namespace spdlog {
    class logger;
}
namespace tools::finite::opt::internal {
    using cplx = std::complex<double>;
    using real = double;
    template<typename FunctorType>
    class CustomLogCallback : public ceres::IterationCallback {
        public:
        std::shared_ptr<spdlog::logger> log;
        const FunctorType              &functor;
        // The default values below are modified during construction when settings::debug == true
        int    freq_log_iter = 10; // Permit logs after this many iterations
        int    last_log_iter = 0;
        int    init_log_iter = 0; // Wait at least this long before first log
        double freq_log_time = 5; // Permit logs after this many seconds
        double last_log_time = 0;
        double init_log_time = 0; // Wait at least this long before first log
        size_t last_count    = 0;

        explicit CustomLogCallback(const FunctorType &functor_);
        ceres::CallbackReturnType operator()(const ceres::IterationSummary &summary);
    };
}