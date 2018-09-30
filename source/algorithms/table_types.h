//
// Created by david on 2018-05-24.
//

#ifndef DMRG_TABLE_TYPES_H
#define DMRG_TABLE_TYPES_H

#include <vector>
#include <array>
#include <hdf5.h>
#include <hdf5_hl.h>



class class_table_dmrg {
private:
    struct data_struct {
        int     iteration;
        int     chain_length;
        int     position;
        long    chi;
        long    chi_max;
        double  energy_mpo; double  energy_ham; double  energy_mom;
        double  energy_min; double  energy_max; double  energy_tgt;
        double  variance_mpo; double  variance_ham; double  variance_mom;
        double  entanglement_entropy;
        double  truncation_error;
        double  wall_time;

        data_struct(
        int iteration_,
        int chain_length_,
        int position_,
        long chi_,
        long chi_max_,
        double energy_mpo_, double energy_ham_, double energy_mom_,
        double energy_min_, double energy_max_, double energy_tgt_,
        double variance_mpo_, double variance_ham_, double variance_mom_,
        double entropy_,
        double truncation_error_,
        double wall_time_) :
        iteration(iteration_),
        chain_length(chain_length_),
        position(position_),
        chi(chi_),
        chi_max(chi_max_),
        energy_mpo(energy_mpo_), energy_ham(energy_ham_), energy_mom(energy_mom_),
        energy_min(energy_min_), energy_max(energy_max_), energy_tgt(energy_tgt_),
        variance_mpo(variance_mpo_), variance_ham(variance_ham_), variance_mom(variance_mom_),
        entanglement_entropy(entropy_),
        truncation_error(truncation_error_),
        wall_time(wall_time_)
        {}
    };
    struct meta_struct {
        constexpr static hsize_t NFIELDS = 17;
        size_t dst_size = sizeof(data_struct);
        std::array<size_t, NFIELDS> dst_offsets = {HOFFSET(data_struct, iteration),
                                                   HOFFSET(data_struct, chain_length),
                                                   HOFFSET(data_struct, position),
                                                   HOFFSET(data_struct, chi),
                                                   HOFFSET(data_struct, chi_max),
                                                   HOFFSET(data_struct, energy_mpo), HOFFSET(data_struct, energy_ham), HOFFSET(data_struct, energy_mom),
                                                   HOFFSET(data_struct, energy_min), HOFFSET(data_struct, energy_max), HOFFSET(data_struct, energy_tgt),
                                                   HOFFSET(data_struct, variance_mpo), HOFFSET(data_struct, variance_ham), HOFFSET(data_struct, variance_mom),
                                                   HOFFSET(data_struct, entanglement_entropy),
                                                   HOFFSET(data_struct, truncation_error),
                                                   HOFFSET(data_struct, wall_time)
        };
        std::array<size_t, NFIELDS> dst_sizes = {
                sizeof(data_struct::iteration),
                sizeof(data_struct::chain_length),
                sizeof(data_struct::position),
                sizeof(data_struct::chi),
                sizeof(data_struct::chi_max),
                sizeof(data_struct::energy_mpo), sizeof(data_struct::energy_ham), sizeof(data_struct::energy_mom),
                sizeof(data_struct::energy_min), sizeof(data_struct::energy_max), sizeof(data_struct::energy_tgt),
                sizeof(data_struct::variance_mpo), sizeof(data_struct::variance_ham), sizeof(data_struct::variance_mom),
                sizeof(data_struct::entanglement_entropy), sizeof(data_struct::truncation_error), sizeof(data_struct::wall_time)
        };
        std::array<const char *, NFIELDS> field_names = {
                "iteration",
                "chain_length",
                "position",
                "chi",
                "chi_max",
                "energy_mpo","energy_ham","energy_mom",
                "energy_min","energy_max","energy_tgt",
                "variance_mpo","variance_ham","variance_mom",
                "entanglement_entropy",
                "truncation_error",
                "wall_time"
        };

        std::array<hid_t, NFIELDS> field_types = {H5T_NATIVE_INT,
                                                  H5T_NATIVE_INT,
                                                  H5T_NATIVE_INT,
                                                  H5T_NATIVE_LONG,
                                                  H5T_NATIVE_LONG,
                                                  H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                  H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                  H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                  H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE};
        hsize_t chunk_size = 1;
        void *fill_data = nullptr;
        int compress = 0;
    };
public:
    class_table_dmrg() = default;
    meta_struct meta;
    std::vector<data_struct> buffer;

};


class class_table_tebd{
private:
    struct data_struct {
        int     iteration;
        long    chi;
        long    chi_max;
        double  time_step;
        double  energy_mpo; double  energy_ham; double  energy_mom;
        double  variance_mpo; double  variance_ham; double  variance_mom;
        double  entanglement_entropy;
        double  truncation_error;
        double  phys_time;
        double  wall_time;

        data_struct(
                int    iteration_,
                long   chi_,
                long   chi_max_,
                double time_step_,
                double energy_mpo_, double energy_ham_, double energy_mom_,
                double variance_mpo_, double variance_ham_, double variance_mom_,
                double entropy_,
                double truncation_error_,
                double phys_time_,
                double wall_time_) :
                iteration(iteration_),
                chi(chi_),
                chi_max(chi_max_),
                time_step(time_step_),
                energy_mpo(energy_mpo_), energy_ham(energy_ham_), energy_mom(energy_mom_),
                variance_mpo(variance_mpo_), variance_ham(variance_ham_), variance_mom(variance_mom_),
                entanglement_entropy(entropy_),
                truncation_error(truncation_error_),
                phys_time(phys_time_),
                wall_time(wall_time_)
        {}
    };
    struct meta_struct {
        constexpr static hsize_t NFIELDS = 14;
        size_t dst_size = sizeof(data_struct);
        std::array<size_t, NFIELDS> dst_offsets = {HOFFSET(data_struct, iteration),
                                                   HOFFSET(data_struct, chi),
                                                   HOFFSET(data_struct, chi_max),
                                                   HOFFSET(data_struct, time_step),
                                                   HOFFSET(data_struct, energy_mpo), HOFFSET(data_struct, energy_ham), HOFFSET(data_struct, energy_mom),
                                                   HOFFSET(data_struct, variance_mpo), HOFFSET(data_struct, variance_ham), HOFFSET(data_struct, variance_mom),
                                                   HOFFSET(data_struct, entanglement_entropy),
                                                   HOFFSET(data_struct, truncation_error),
                                                   HOFFSET(data_struct, phys_time),
                                                   HOFFSET(data_struct, wall_time)
        };
        std::array<size_t, NFIELDS> dst_sizes = {
                sizeof(data_struct::iteration),
                sizeof(data_struct::chi),
                sizeof(data_struct::chi_max),
                sizeof(data_struct::time_step),
                sizeof(data_struct::energy_mpo), sizeof(data_struct::energy_ham), sizeof(data_struct::energy_mom),
                sizeof(data_struct::variance_mpo), sizeof(data_struct::variance_ham), sizeof(data_struct::variance_mom),
                sizeof(data_struct::entanglement_entropy),
                sizeof(data_struct::truncation_error),
                sizeof(data_struct::phys_time),
                sizeof(data_struct::wall_time)
        };
        std::array<const char *, NFIELDS> field_names = {"iteration",
                                                         "chi",
                                                         "chi_max",
                                                         "time_step",
                                                         "energy_mpo","energy_ham","energy_mom",
                                                         "variance_mpo","variance_ham","variance_mom",
                                                         "entanglement_entropy",
                                                         "truncation_error",
                                                         "phys_time",
                                                         "wall_time"
        };

        std::array<hid_t, NFIELDS> field_types = {H5T_NATIVE_INT,
                                                  H5T_NATIVE_LONG,
                                                  H5T_NATIVE_LONG,
                                                  H5T_NATIVE_DOUBLE,
                                                  H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                  H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                  H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                  H5T_NATIVE_DOUBLE};
        hsize_t chunk_size = 1;
        void *fill_data = nullptr;
        int compress = 0;
    };
public:
    class_table_tebd() = default;
    meta_struct meta;
    std::vector<data_struct> buffer;
};



class class_table_finite_chain{
private:
    struct data_struct {
        int     iteration;
        int     chain_length;
        int     position;
        long    chi;
        double  energy;
        double  entanglement_entropy;
        double  truncation_error;

        data_struct(   int iteration_,
                int chain_length_,
                int position_,
                long chi_,
                double energy_,
                double entanglement_entropy_,
                double truncation_error_):
                iteration(iteration_),
                chain_length(chain_length_),
                position(position_),
                chi(chi_),
                energy(energy_),
                entanglement_entropy(entanglement_entropy_),
                truncation_error(truncation_error_)
        {}
    };
    struct meta_struct {
        constexpr static hsize_t NFIELDS = 7;
        size_t dst_size = sizeof(data_struct);
        std::array<size_t, NFIELDS> dst_offsets = {HOFFSET(data_struct, iteration),
                                                   HOFFSET(data_struct, chain_length),
                                                   HOFFSET(data_struct, position),
                                                   HOFFSET(data_struct, chi),
                                                   HOFFSET(data_struct, energy),
                                                   HOFFSET(data_struct, entanglement_entropy),
                                                   HOFFSET(data_struct, truncation_error)
                                                    };
        std::array<size_t, NFIELDS> dst_sizes = {
                sizeof(data_struct::iteration),
                sizeof(data_struct::chain_length),
                sizeof(data_struct::position),
                sizeof(data_struct::chi),
                sizeof(data_struct::energy),
                sizeof(data_struct::entanglement_entropy),
                sizeof(data_struct::truncation_error)
        };
        std::array<const char *, NFIELDS> field_names = {"iteration",
                                                         "chain_length",
                                                         "position",
                                                         "chi",
                                                         "energy",
                                                         "entanglement_entropy",
                                                         "truncation_error"
        };

        std::array<hid_t, NFIELDS> field_types = {H5T_NATIVE_INT,
                                                  H5T_NATIVE_INT,
                                                  H5T_NATIVE_INT,
                                                  H5T_NATIVE_LONG,
                                                  H5T_NATIVE_DOUBLE,
                                                  H5T_NATIVE_DOUBLE,
                                                  H5T_NATIVE_DOUBLE
                                                };
        hsize_t chunk_size = 1;
        void *fill_data = nullptr;
        int compress = 0;
    };
public:
    class_table_finite_chain() = default;
    meta_struct meta;
    std::vector<data_struct> buffer;

};


//Profiling table definition
class class_table_profiling{
private:
    struct data_struct{
        int    iteration;
        double t_tot;
        double t_opt;
        double t_sim;
        double t_svd;
        double t_env;
        double t_evo;
        double t_udt;
        double t_sto;
        double t_ste;
        double t_prt;
        double t_obs;
        double t_mps;
        double t_chi;

        data_struct(int    iteration_,
             double t_tot_,  double t_opt_, double t_sim_,
             double t_svd_,  double t_env_, double t_evo_,
             double t_udt_,  double t_sto_, double t_ste_,
             double t_prt_,  double t_obs_, double t_mps_,
             double t_chi_)
                :iteration(iteration_),
                 t_tot(t_tot_), t_opt(t_opt_), t_sim(t_sim_),
                 t_svd(t_svd_), t_env(t_env_), t_evo(t_evo_),
                 t_udt(t_udt_), t_sto(t_sto_), t_ste(t_ste_),
                 t_prt(t_prt_), t_obs(t_obs_), t_mps(t_mps_),
                 t_chi(t_chi_)
        {}
    };
    struct meta_struct{
        constexpr static hsize_t                NFIELDS     = 14;
        size_t           dst_size                           = sizeof (data_struct);
        std::array       <size_t,NFIELDS>       dst_offsets = {HOFFSET(data_struct, iteration),
                                                               HOFFSET(data_struct, t_tot), HOFFSET(data_struct, t_opt), HOFFSET(data_struct, t_sim),
                                                               HOFFSET(data_struct, t_svd), HOFFSET(data_struct, t_env), HOFFSET(data_struct, t_evo),
                                                               HOFFSET(data_struct, t_udt), HOFFSET(data_struct, t_sto), HOFFSET(data_struct, t_ste),
                                                               HOFFSET(data_struct, t_prt), HOFFSET(data_struct, t_obs), HOFFSET(data_struct, t_mps),
                                                               HOFFSET(data_struct, t_chi)
        };
        std::array       <size_t,NFIELDS>       dst_sizes   = {sizeof(data_struct::iteration),
                                                               sizeof(data_struct::t_tot), sizeof(data_struct::t_opt), sizeof(data_struct::t_sim),
                                                               sizeof(data_struct::t_svd), sizeof(data_struct::t_env), sizeof(data_struct::t_evo),
                                                               sizeof(data_struct::t_udt), sizeof(data_struct::t_sto), sizeof(data_struct::t_ste),
                                                               sizeof(data_struct::t_prt), sizeof(data_struct::t_obs), sizeof(data_struct::t_mps),
                                                               sizeof(data_struct::t_chi)
        };
        std::array       <const char*,NFIELDS>  field_names = {"iteration",
                                                               "t_tot", "t_opt", "t_sim",
                                                               "t_svd", "t_env", "t_evo",
                                                               "t_udt", "t_sto", "t_ste",
                                                               "t_prt", "t_obs", "t_mps",
                                                               "t_con"
        };

        std::array       <hid_t,NFIELDS>        field_types = {H5T_NATIVE_INT,
                                                               H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                               H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                               H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                               H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                               H5T_NATIVE_DOUBLE
        };

        hsize_t          chunk_size                         = 1;
        void             *fill_data                         = nullptr;
        int              compress                           = 0;
    };
public:
    class_table_profiling() = default;
    meta_struct meta;
    std::vector<data_struct> buffer;
};


#endif //DMRG_TABLE_TYPES_H
