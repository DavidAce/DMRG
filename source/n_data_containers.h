//
// Created by david on 8/17/17.
//

#ifndef DMRG_N_DATA_CONTAINERS_H
#define DMRG_N_DATA_CONTAINERS_H
#include <vector>
#include <class_superblock.h>

class output_data_container{
public:
    class_hdf5_dataset_buffer<long>   chi_cur;
    class_hdf5_dataset_buffer<long>   chi_max;
    class_hdf5_dataset_buffer<double> energy;
    class_hdf5_dataset_buffer<double> entropy;
    class_hdf5_dataset_buffer<double> variance;
    class_hdf5_dataset_buffer<double> trerror;
    class_hdf5_dataset_buffer<double> e_error;
    class_hdf5_dataset_buffer<double> time_step;
    class_hdf5_dataset_buffer<double> phys_time;
    class_hdf5_dataset_buffer<double> wall_time;
    class_hdf5_dataset_buffer<int>    length;

    output_data_container( class_hdf5 &hdf5_,
                           const std::string &group_name,
                           const long iteration
                            ):
            chi_cur    (hdf5_, group_name, iteration, "chi_cur"    ),
            chi_max    (hdf5_, group_name, iteration, "chi_max"    ),
            energy     (hdf5_, group_name, iteration, "energy"     ),
            entropy    (hdf5_, group_name, iteration, "entropy"    ),
            variance   (hdf5_, group_name, iteration, "variance"    ),
            trerror    (hdf5_, group_name, iteration, "truncation_error"    ),
            e_error    (hdf5_, group_name, iteration, "energy_error"    ),
            time_step  (hdf5_, group_name, iteration, "time_step"   ),
            phys_time  (hdf5_, group_name, iteration, "phys_time"   ),
            wall_time  (hdf5_, group_name, iteration, "wall_time"   ),
            length     (hdf5_, group_name, iteration, "length"   )
    {
    }


    void push_back(class_superblock &superblock){
        chi_cur  .push_back(superblock.chi);
        chi_max  .push_back(superblock.chi_max);
        energy   .push_back(superblock.get_energy());
        entropy  .push_back(superblock.get_entropy());
        variance .push_back(superblock.get_variance());
        trerror  .push_back(superblock.truncation_error);
        length   .push_back(superblock.chain_length);
    }


    void push_back(double time_step_,
                   double phys_time_,
                   double wall_time_)
    {
        time_step.push_back(time_step_);
        phys_time.push_back(phys_time_);
        wall_time.push_back(wall_time_);
    }

    void push_back(double energy_error_){
        e_error.push_back(energy_error_);
    }

    void push_back(long chi_cur_,
                   long chi_max_,
                   double energy_,
                   double entropy_,
                   double variance_,
                   double trunc_error_,
                   double time_step_,
                   double phys_time_,
                   double wall_time_)
    {
        chi_cur  .push_back(chi_cur_);
        chi_max  .push_back(chi_max_);
        energy   .push_back(energy_);
        entropy  .push_back(entropy_);
        variance .push_back(variance_);
        trerror  .push_back(trunc_error_);
        time_step.push_back(time_step_);
        phys_time.push_back(phys_time_);
        wall_time.push_back(wall_time_);
    }
};


#endif //DMRG_N_DATA_CONTAINERS_H
