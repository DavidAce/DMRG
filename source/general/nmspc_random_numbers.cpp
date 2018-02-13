//
// Created by david on 2016-07-24.
//

#include "nmspc_random_numbers.h"
using namespace Eigen;
using namespace std;
namespace rn{
    std::mt19937 rng;

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