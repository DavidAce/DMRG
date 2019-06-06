//
// Created by david on 2016-07-24.
//

#include "nmspc_random_numbers.h"

using namespace Eigen;
using namespace std;
namespace rn{
    std::mt19937 rng;
    void seed(unsigned long n){
        rng.seed(n);
        std::srand((unsigned int) n);
    }

    int uniform_integer_1(){
        std::uniform_int_distribution<>  rand_int(0, 1);
        return rand_int(rng);
    }

    int uniform_integer(const int min, const int max){
        std::uniform_int_distribution<>  rand_int(min,max);
        return rand_int(rng);
    }

    double uniform_double_1(){
        std::uniform_real_distribution<>  rand_real(0,1);
        return rand_real(rng);
    }

    double uniform_double(const double min, const double max){
        std::uniform_real_distribution<>  rand_real(std::min(min,max),std::max(min,max));
        return rand_real(rng);
    }


     std::complex<double> uniform_complex_1(){
        std::uniform_real_distribution<>  rand_real(0,2.0*M_PI);
        return std::polar(1.0,rand_real(rng));
    }



     double log_normal(const double mean, const double std){
        std::lognormal_distribution<double>  distribution(mean, std);
        return distribution(rng);
    }

    ArrayXd random_with_replacement(const ArrayXd & in){
        vector<double> boot;
        for (int i = 0; i < in.derived().size(); i++){
            boot.push_back(in(uniform_integer(0,(int)(in.size()-1))));
        }
        return Eigen::Map<ArrayXd>(boot.data(),boot.size());
    }
    ArrayXd random_with_replacement(const ArrayXd & in, const int n){
        if (n > in.size()){cout << "n too large" << endl; exit(1);}
        vector<double> boot;
        for (int i = 0; i < n; i++){
            boot.push_back(in(uniform_integer(0,(int)(in.size()-1))));
        }
        return Eigen::Map<ArrayXd>(boot.data(),boot.size());
    }

    double gaussian_truncated(const double lowerLimit, const double upperLimit, const double mean, const double std) {
        std::normal_distribution<double> distribution(mean,std);
        double ul = fmax(lowerLimit, upperLimit);
        double ll = fmin(lowerLimit, upperLimit);
        double number;
        while (true) {
            number = distribution(rng);
            if (number >= ll && number <= ul) {
                return number;
            }
        }
    }
}