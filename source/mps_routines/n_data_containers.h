//
// Created by david on 8/17/17.
//

#ifndef DMRG_N_DATA_CONTAINERS_H
#define DMRG_N_DATA_CONTAINERS_H
#include <vector>
#include <general/class_hdf5_dataset_buffer.h>
#include <sim_parameters/n_model.h>
#include "class_observables.h"

class output_data_container{
private:
    class_hdf5 *hdf5_out; //Pointer is not owned! do not delete
    bool data_has_been_written_to_file = false;
public:
    class_hdf5_dataset_buffer<long>   chi;
    class_hdf5_dataset_buffer<long>   chi_max;
    class_hdf5_dataset_buffer<double> energy;
    class_hdf5_dataset_buffer<double> energy_ex; //Exact energy
    class_hdf5_dataset_buffer<double> entropy;
    class_hdf5_dataset_buffer<double> variance1;
    class_hdf5_dataset_buffer<double> variance2;
    class_hdf5_dataset_buffer<double> trerror;
    class_hdf5_dataset_buffer<double> e_error;
    class_hdf5_dataset_buffer<int>    iter_step;
    class_hdf5_dataset_buffer<double> time_step;
    class_hdf5_dataset_buffer<double> phys_time;
    class_hdf5_dataset_buffer<double> wall_time;
    class_hdf5_dataset_buffer<long>   length;

    output_data_container(const std::string &group_name,
                           const int iteration
                            ):
            output_data_container(nullptr,group_name, iteration)
    {
    }
    output_data_container(class_hdf5 * hdf5_out_ ,const std::string &group_name,
                          const int iteration
    ):      hdf5_out   (hdf5_out_),
            chi        (group_name, iteration, "chi"    ),
            chi_max    (group_name, iteration, "chi_max"    ),
            energy     (group_name, iteration, "energy"     ),
            energy_ex  (group_name, iteration, "energy_ex"  ),
            entropy    (group_name, iteration, "entropy"    ),
            variance1  (group_name, iteration, "variance1"   ),
            variance2  (group_name, iteration, "variance2"   ),
            trerror    (group_name, iteration, "truncation_error"    ),
            e_error    (group_name, iteration, "energy_error"    ),
            iter_step  (group_name, iteration, "iter_step"  ),
            time_step  (group_name, iteration, "time_step"  ),
            phys_time  (group_name, iteration, "phys_time"  ),
            wall_time  (group_name, iteration, "wall_time"  ),
            length     (group_name, iteration, "length"     )
    {
    }

    ~output_data_container(){
        if (hdf5_out != nullptr){
            write_data(*hdf5_out);
        }else if (!data_has_been_written_to_file){
            std::cerr << "Output data has not saved to file, yet it is being discarded!\n" << std::endl;
        }
    }

    void push_back(class_observables &observables){
        chi      .push_back(observables.get_chi());
        chi_max  .push_back(observables.get_chi_max());
        energy   .push_back(observables.get_energy());
        energy_ex.push_back(Model::energy_exact);
        entropy  .push_back(observables.get_entropy());
        variance1 .push_back(observables.get_variance1());
        variance2 .push_back(observables.get_variance2());
        trerror  .push_back(observables.get_truncation_error());
        length   .push_back(observables.get_chain_length());
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
                   double variance1_,
                   double variance2_,
                   double variance3_,
                   double trunc_error_,
                   double time_step_,
                   double phys_time_,
                   double wall_time_)
    {
        chi      .push_back(chi_cur_);
        chi_max  .push_back(chi_max_);
        energy   .push_back(energy_);
        entropy  .push_back(entropy_);
        variance1 .push_back(variance1_);
        variance2 .push_back(variance2_);
        trerror  .push_back(trunc_error_);
        time_step.push_back(time_step_);
        phys_time.push_back(phys_time_);
        wall_time.push_back(wall_time_);
    }

    void write_data(class_hdf5 & hdf5_out){
        chi      .write_buffer_to_file(hdf5_out);
        chi_max  .write_buffer_to_file(hdf5_out);
        energy   .write_buffer_to_file(hdf5_out);
        energy_ex.write_buffer_to_file(hdf5_out);
        entropy  .write_buffer_to_file(hdf5_out);
        variance1 .write_buffer_to_file(hdf5_out);
        variance2 .write_buffer_to_file(hdf5_out);
        trerror  .write_buffer_to_file(hdf5_out);
        e_error  .write_buffer_to_file(hdf5_out);
        iter_step.write_buffer_to_file(hdf5_out);
        time_step.write_buffer_to_file(hdf5_out);
        phys_time.write_buffer_to_file(hdf5_out);
        wall_time.write_buffer_to_file(hdf5_out);
        length   .write_buffer_to_file(hdf5_out);
        data_has_been_written_to_file = true;
    }

};


#endif //DMRG_N_DATA_CONTAINERS_H
