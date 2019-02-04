//
// Created by david on 2019-01-29.
//

#ifndef CLASS_FINITE_CHAIN_HDF5_H
#define CLASS_FINITE_CHAIN_HDF5_H
#include <memory>
#include <complex>
class class_finite_chain_state;
class class_hdf5_file;
//class class_superblock;


class class_finite_chain_hdf5{
private:
    using state_ptr = std::shared_ptr<class_finite_chain_state>;
    using hdf5_ptr  = std::shared_ptr<class_hdf5_file>;
    using Scalar    = std::complex<double>;
public:

    // Functions relating to HDF5 file storage
    static void write_all_to_hdf5(state_ptr state, hdf5_ptr hdf5, std::string sim_name);
    static void write_bond_matrices_to_hdf5(state_ptr state, hdf5_ptr hdf5, std::string sim_name);
    static void write_2site_mps_to_hdf5(state_ptr state, hdf5_ptr hdf5, std::string sim_name);
    static void write_2site_mpo_to_hdf5(state_ptr state, hdf5_ptr hdf5, std::string sim_name);
    static void write_2site_env_to_hdf5(state_ptr state, hdf5_ptr hdf5, std::string sim_name);
    static void write_2site_env2_to_hdf5(state_ptr state, hdf5_ptr hdf5, std::string sim_name);
    static void write_full_mps_to_hdf5(state_ptr state, hdf5_ptr hdf5, std::string sim_name);
    static void write_full_mpo_to_hdf5(state_ptr state, hdf5_ptr hdf5, std::string sim_name);
    static void write_hamiltonian_params_to_hdf5(state_ptr state, hdf5_ptr hdf5, std::string sim_name);
    static void write_entanglement_to_hdf5(state_ptr state, hdf5_ptr hdf5, std::string sim_name);

};


#endif //DMRG_CLASS_FINITE_CHAIN_HDF5_H

