//
// Created by david on 2016-07-24.
//

#include "math/rnd.h"
#include <iostream>

using namespace Eigen;
using namespace std;
namespace rnd {
    void seed(std::optional<long> n){
        if(n.has_value() and n.value() >= 0){
            auto given_seed = (unsigned long) n.value();
            std::cout << "Seeding : " << given_seed << std::endl;
            pcg_extras::seed_seq_from<pcg64> seq (given_seed);
//            std::seed_seq seq{given_seed};
            internal::rng.seed(seq);
        }else{
            std::cout << "Seeding : std::random_device" << std::endl;
            pcg_extras::seed_seq_from<std::random_device> seed_source;
            internal::rng.seed(seed_source);
        }
        std::srand(static_cast<unsigned>(internal::rng()));
    }

    int uniform_integer_01(){
        return internal::rand_int_01(internal::rng);
    }

    double uniform_double_01(){
        return internal::rand_double_01(internal::rng);
    }

    double uniform_double_box(const double min, const double max){
        std::uniform_real_distribution<>  rand_real(std::min(min,max),std::max(min,max));
        return rand_real(internal::rng);
    }


     std::complex<double> uniform_complex_1(){
        return std::polar(1.0,internal::rand_double_0_2pi(internal::rng));
    }


    double normal(const double mean, const double std){
        std::normal_distribution<double>  distribution(mean, std);
        return distribution(internal::rng);
    }


     double log_normal(const double mean, const double std){
        std::lognormal_distribution<double>  distribution(mean, std);
        return distribution(internal::rng);
    }

    ArrayXd random_with_replacement(const ArrayXd & in){
        vector<double> boot;
        for (int i = 0; i < in.derived().size(); i++){
            boot.push_back(in(uniform_integer_box(0, (int) (in.size() - 1))));
        }
        return Eigen::Map<ArrayXd>(boot.data(), static_cast<Index>(boot.size()));
    }
    ArrayXd random_with_replacement(const ArrayXd & in, const int n){
        if (n > in.size()){cout << "n too large" << endl; exit(1);}
        vector<double> boot;
        for (int i = 0; i < n; i++){
            boot.push_back(in(uniform_integer_box(0, (int) (in.size() - 1))));
        }
        return Eigen::Map<ArrayXd>(boot.data(), static_cast<Index>(boot.size()));
    }

    double gaussian_truncated(const double lowerLimit, const double upperLimit, const double mean, const double std) {
        std::normal_distribution<double> distribution(mean,std);
        double ul = fmax(lowerLimit, upperLimit);
        double ll = fmin(lowerLimit, upperLimit);
        double number;
        while (true) {
            number = distribution(internal::rng);
            if (number >= ll && number <= ul) {
                return number;
            }
        }
    }
}