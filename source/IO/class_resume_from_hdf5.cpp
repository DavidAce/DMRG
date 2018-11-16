//
// Created by david on 2018-11-15.
//


#include "class_resume_from_hdf5.h"
#include <IO/class_hdf5_file.h>
#include <mps_routines/class_superblock.h>
#include <mps_routines/class_finite_chain_sweeper.h>
#include <memory>


class_resume_from_hdf5::class_resume_from_hdf5(
        std::shared_ptr<class_hdf5_file>                hdf5_,
        std::shared_ptr<class_superblock>               superblock_,
        std::shared_ptr<class_finite_chain_sweeper>     env_storage_,
        std::string sim_name_,
        SimulationType sim_type_
        )
        :
        hdf5           (std::move(hdf5_)),
        superblock     (std::move(superblock_)),
        env_storage    (std::move(env_storage_)),
        sim_name       (std::move(sim_name_)),
        sim_type       (std::move(sim_type_))
{
        if (settings::hdf5::resume_from_file){
                std::cout << "Attempting to resume" << std::endl;
                load_mps_chain();
        }

}


void class_resume_from_hdf5::load_mps_chain(){
        using Scalar = std::complex<double> ;
        Eigen::Tensor<Scalar,3> test;
        hdf5->read_dataset(test,"/xDMRG/chain/MPS/A_0");
        std::cout << "THE TENSOR \n" << test << std::endl;
}
