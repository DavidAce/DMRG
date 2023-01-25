#include "../rnd.h"
#include <algorithm>
#include <cstdio>
#include <omp.h>
#include <random>
#include <stdexcept>

#if defined(NDEBUG)
namespace rnd {
    static constexpr bool debug = false;
}
#else
    #include <omp.h>
namespace rnd {
    static constexpr bool debug = true;
}
#endif

namespace rnd {

    constexpr std::string_view enum2sv(const dist &d) {
        switch(d) {
            case dist::uniform: return "uniform";
            case dist::normal: return "normal";
            case dist::lognormal: return "lognormal";
        }
    }
    constexpr dist sv2enum(std::string_view d) {
        if(d == "uniform") return dist::uniform;
        if(d == "normal") return dist::normal;
        if(d == "lognormal") return dist::lognormal;
        throw std::runtime_error("rnd: unrecognized distribution: " + std::string(d));
    }

    namespace internal {
        // Make a random number engine
        //        inline pcg64 rng;
        // Commonly used distributions
        inline std::uniform_int_distribution<int>     rand_int_01(0, 1);
        inline std::uniform_real_distribution<double> rand_double_01(0.0, 1.0);
        inline std::uniform_real_distribution<double> rand_double_0_2pi(0, 2.0 * M_PI);
        inline std::normal_distribution<double>       normal_double_01(0.0, 1.0);
    }

    void seed(std::optional<long> n) {
        if(n.has_value() and n.value() >= 0) {
            auto given_seed = (unsigned long) n.value();
            std::printf("Seeding: %ld\n", given_seed);
            pcg_extras::seed_seq_from<pcg64> seq(given_seed);
            //            std::seed_seq seq{given_seed};
            internal::rng.seed(seq);
        } else {
            std::printf("Seeding: std::random_device\n");
            pcg_extras::seed_seq_from<std::random_device> seed_source;
            internal::rng.seed(seed_source);
        }
        std::srand(static_cast<unsigned>(internal::rng()));
    }

    int uniform_integer_01() {
        if constexpr(debug)
            if(omp_get_num_threads() > 1) throw std::runtime_error("rnd::uniform_integer_01 is not thread safe!");
        return internal::rand_int_01(internal::rng);
    }

    double uniform_double_01() {
        if constexpr(debug)
            if(omp_get_num_threads() > 1) throw std::runtime_error("rnd::uniform_double_01 is not thread safe!");
        return internal::rand_double_01(internal::rng);
    }

    template<typename T>
    T uniform_integer_box(T min, T max) {
        if constexpr(debug)
            if(omp_get_num_threads() > 1) throw std::runtime_error("rnd::uniform_integer_box is not thread safe!");
        std::uniform_int_distribution<T> rand_int(std::min(min, max), std::max(min, max));
        return rand_int(internal::rng);
    }
    template int      uniform_integer_box(int min, int max);
    template unsigned uniform_integer_box(unsigned min, unsigned max);
    template long     uniform_integer_box(long min, long max);
    template size_t   uniform_integer_box(size_t min, size_t max);

    double uniform_double_box(double min, double max) {
        if constexpr(debug)
            if(omp_get_num_threads() > 1) throw std::runtime_error("rnd::uniform_double_box is not thread safe!");
        std::uniform_real_distribution<> rand_real(std::min(min, max), std::max(min, max));
        return rand_real(internal::rng);
    }
    double uniform_double_box(double halfwidth) {
        if constexpr(debug)
            if(omp_get_num_threads() > 1) throw std::runtime_error("rnd::uniform_double_box is not thread safe!");
        std::uniform_real_distribution<> rand_real(-halfwidth, halfwidth);
        return rand_real(internal::rng);
    }

    std::complex<double> uniform_complex_in_unit_circle() { return std::polar(uniform_double_01(), internal::rand_double_0_2pi(internal::rng)); }

    std::complex<double> uniform_complex_on_unit_circle() { return std::polar(1.0, internal::rand_double_0_2pi(internal::rng)); }

    std::complex<double> uniform_complex_box(double real_min, double real_max, double imag_min, double imag_max) {
        return {uniform_double_box(real_min, real_max), uniform_double_box(imag_min, imag_max)};
    }

    template<typename T>
    std::vector<T> uniform_unit_n_sphere(size_t n) {
        std::vector<T> arr;
        double         norm = 0.0;
        for(size_t i = 0; i < n; i++) {
            if constexpr(std::is_same<T, std::complex<double>>::value) {
                double re   = internal::normal_double_01(internal::rng);
                double im   = internal::normal_double_01(internal::rng);
                T      cplx = T(1.0, 0.0) * re + T(0.0, 1.0) * im;
                arr.push_back(cplx);
                norm += re * re + im * im;
            } else {
                arr.push_back(internal::normal_double_01(internal::rng));
                norm += std::abs(arr[i] * arr[i]);
            }
        }

        norm = std::sqrt(norm);
        for(size_t i = 0; i < n; i++) { arr[i] /= norm; }
        return arr;
    }
    template std::vector<double>               uniform_unit_n_sphere(size_t n);
    template std::vector<std::complex<double>> uniform_unit_n_sphere(size_t n);

    std::complex<double> uniform_complex_slice(double radius_max, double angle_min, double angle_max) {
        return std::polar(uniform_double_box(0, radius_max), uniform_double_box(angle_min, angle_max));
    }

    double normal(const double mean, const double std) {
        if constexpr(debug)
            if(omp_get_num_threads() > 1) throw std::runtime_error("rnd::normal is not thread safe!");
        std::normal_distribution<double> distribution(mean, std);
        return distribution(internal::rng);
    }

    double log_normal(const double mean, const double std) {
        if constexpr(debug)
            if(omp_get_num_threads() > 1) throw std::runtime_error("rnd::log_normal is not thread safe!");
        std::lognormal_distribution<double> distribution(mean, std);
        return distribution(internal::rng);
    }

    std::vector<int> random_with_replacement(const std::vector<int> &in) {
        std::vector<int> boot;
        boot.reserve(in.size());
        for(size_t i = 0; i < in.size(); i++) boot.emplace_back(in[uniform_integer_box(0ul, in.size() - 1)]);
        return boot;
    }
    std::vector<int> random_with_replacement(const std::vector<int> &in, const size_t n) {
        if(n > in.size()) throw std::logic_error("random_with_replacement: n too large");
        std::vector<int> boot;
        boot.reserve(n);
        for(size_t i = 0; i < n; i++) { boot.emplace_back(in[uniform_integer_box(0ul, in.size() - 1)]); }
        return boot;
    }

    double gaussian_truncated(const double lowerLimit, const double upperLimit, const double mean, const double std) {
        std::normal_distribution<double> distribution(mean, std);
        double                           ul = fmax(lowerLimit, upperLimit);
        double                           ll = fmin(lowerLimit, upperLimit);
        double                           number;
        while(true) {
            number = distribution(internal::rng);
            if(number >= ll && number <= ul) { return number; }
        }
    }

    template<typename T>
    void shuffle(T &list) {
        std::shuffle(std::begin(list), std::end(list), internal::rng);
    }
    template void shuffle(std::vector<int> &list);
    template void shuffle(std::vector<unsigned> &list);
    template void shuffle(std::vector<long> &list);
    template void shuffle(std::vector<size_t> &list);
    template void shuffle(std::vector<double> &list);

    template<typename Distribution>
    std::vector<double> random(Distribution &&d, size_t num) {
        auto rndvec = std::vector<double>(num);
        for(size_t i = 0; i < num; ++i) rndvec[i] = d(internal::rng);
        return rndvec;
    }

    double random(dist d, double mean, double width) {
        switch(d) {
            case dist::uniform: return uniform_double_box(mean - width / 2, mean + width / 2);
            case dist::normal: return normal(mean, width);
            case dist::lognormal: return log_normal(mean, width);
            default: throw std::runtime_error("Invalid distribution");
        }
    }
    double random(std::string_view distribution, double mean, double width) { return random(sv2enum(distribution), mean, width); }

    std::vector<double> random(dist d, double mean, double width, size_t num) {
        switch(d) {
            case dist::uniform: return random(std::uniform_real_distribution<double>(mean - width / 2, mean + width / 2), num);
            case dist::normal: return random(std::normal_distribution<double>(mean, width), num);
            case dist::lognormal: return random(std::lognormal_distribution<double>(mean, width), num);
            default: throw std::runtime_error("Invalid distribution");
        }
    }
    std::vector<double> random(dist d, double mean, double width, const std::vector<double> &weights) {
        auto rndvec = random(d, mean, width, weights.size());
        for(size_t i = 0; i < weights.size(); ++i) rndvec[i] *= weights[i];
        return rndvec;
    }

    std::vector<double> random(std::string_view distribution, double mean, double width, size_t num) { return random(sv2enum(distribution), mean, width, num); }
    std::vector<double> random(std::string_view distribution, double mean, double width, const std::vector<double> &weights) {
        return random(sv2enum(distribution), mean, width, weights);
    }

}