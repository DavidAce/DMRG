//
// Created by david on 8/17/17.
//

#ifndef DMRG_N_DATA_CONTAINERS_H
#define DMRG_N_DATA_CONTAINERS_H
#include <vector>
#include <general/class_hdf5.h>
#include <sim_parameters/n_model.h>
#include "mps_routines/class_superblock.h"

class output_data_container{
public:
    class_hdf5_dataset_buffer<long>   chi;
    class_hdf5_dataset_buffer<long>   chi_max;
    class_hdf5_dataset_buffer<double> energy;
    class_hdf5_dataset_buffer<double> energy_ex; //Exact energy
    class_hdf5_dataset_buffer<double> entropy;
    class_hdf5_dataset_buffer<double> variance;
    class_hdf5_dataset_buffer<double> trerror;
    class_hdf5_dataset_buffer<double> e_error;
    class_hdf5_dataset_buffer<int>    iter_step;
    class_hdf5_dataset_buffer<double> time_step;
    class_hdf5_dataset_buffer<double> phys_time;
    class_hdf5_dataset_buffer<double> wall_time;
    class_hdf5_dataset_buffer<int>    length;

    output_data_container( class_hdf5 &hdf5_,
                           const std::string &group_name,
                           const long iteration
                            ):
            chi        (hdf5_, group_name, iteration, "chi"    ),
            chi_max    (hdf5_, group_name, iteration, "chi_max"    ),
            energy     (hdf5_, group_name, iteration, "energy"     ),
            energy_ex  (hdf5_, group_name, iteration, "energy_ex"  ),
            entropy    (hdf5_, group_name, iteration, "entropy"    ),
            variance   (hdf5_, group_name, iteration, "variance"   ),
            trerror    (hdf5_, group_name, iteration, "truncation_error"    ),
            e_error    (hdf5_, group_name, iteration, "energy_error"    ),
            iter_step  (hdf5_, group_name, iteration, "iter_step"  ),
            time_step  (hdf5_, group_name, iteration, "time_step"  ),
            phys_time  (hdf5_, group_name, iteration, "phys_time"  ),
            wall_time  (hdf5_, group_name, iteration, "wall_time"  ),
            length     (hdf5_, group_name, iteration, "length"     )
    {
    }


    void push_back(class_superblock &superblock){
        chi      .push_back(superblock.chi);
        chi_max  .push_back(superblock.chi_max);
        energy   .push_back(superblock.get_energy());
        energy_ex.push_back(Model::energy_exact);
        entropy  .push_back(superblock.get_entropy());
        variance .push_back(superblock.get_variance());
        trerror  .push_back(superblock.truncation_error);
        length   .push_back(superblock.chain_length);
    }


    void push_back(int    iter_step_,
                   double time_step_,
                   double phys_time_,
                   double wall_time_)
    {
        iter_step.push_back(iter_step_);
        time_step.push_back(time_step_);
        phys_time.push_back(phys_time_);
        wall_time.push_back(wall_time_);
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
        chi  .push_back(chi_cur_);
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
