//
// Created by david on 8/17/17.
//

#ifndef FINITE_DMRG_EIGEN_N_DATA_CONTAINERS_H
#define FINITE_DMRG_EIGEN_N_DATA_CONTAINERS_H
#include <vector>

struct dmrg_data_container{
    std::vector<long>   chi;
    std::vector<double> wall_time;
    std::vector<double> energy;
    std::vector<double> entropy;
    std::vector<double> trunc_error;

    void clear(){
        chi.clear();
        wall_time.clear();
        energy.clear();
        entropy.clear();
        trunc_error.clear();
    }
    void push_back(long chi_, double wall_time_, double energy_, double entropy_, double trunc_error_){
        chi.push_back(chi_);
        wall_time.push_back(wall_time_);
        energy.push_back(energy_);
        entropy.push_back(entropy_);
        trunc_error.push_back(trunc_error_);
    }
};

struct tebd_data_container{
    std::vector<long>   chi;
    std::vector<double> time_step;
    std::vector<double> phys_time;
    std::vector<double> wall_time;
    std::vector<double> energy;
    std::vector<double> entropy;
    std::vector<double> trunc_error;

    void clear(){
        chi.clear();
        time_step.clear();
        phys_time.clear();
        wall_time.clear();
        energy.clear();
        entropy.clear();
        trunc_error.clear();
    }
    void push_back(long chi_, double time_step_, double phys_time_, double wall_time_, double energy_, double entropy_, double trunc_error_){
        chi.push_back(chi_);
        time_step.push_back(time_step_);
        phys_time.push_back(phys_time_);
        wall_time.push_back(wall_time_);
        energy.push_back(energy_);
        entropy.push_back(entropy_);
        trunc_error.push_back(trunc_error_);
    }


};
#endif //FINITE_DMRG_EIGEN_N_DATA_CONTAINERS_H
