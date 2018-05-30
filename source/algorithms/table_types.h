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
    struct data {
        int     iteration;
        int     chain_length;
        int     position;
        int     chi;
        int     chi_max;
        double  energy_mpo; double  energy_ham; double  energy_mom;
        double  variance_mpo; double  variance_ham; double  variance_mom;
        double  entanglement_entropy;
        double  truncation_error;
        double  wall_time;

        data(
        int iteration_,
        int chain_length_,
        int position_,
        int chi_,
        int chi_max_,
        double energy_mpo_, double energy_ham_, double energy_mom_,
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
        variance_mpo(variance_mpo_), variance_ham(variance_ham_), variance_mom(variance_mom_),
        entanglement_entropy(entropy_),
        truncation_error(truncation_error_),
        wall_time(wall_time_)
        {}
    };
    struct meta_struct {
        constexpr static hsize_t NFIELDS = 14;
        size_t dst_size = sizeof(data);
        std::array<size_t, NFIELDS> dst_offsets = {HOFFSET(data, iteration),
                                                   HOFFSET(data, chain_length),
                                                   HOFFSET(data, position),
                                                   HOFFSET(data, chi),
                                                   HOFFSET(data, chi_max),
                                                   HOFFSET(data, energy_mpo), HOFFSET(data, energy_ham), HOFFSET(data, energy_mom),
                                                   HOFFSET(data, variance_mpo), HOFFSET(data, variance_ham), HOFFSET(data, variance_mom),
                                                   HOFFSET(data, entanglement_entropy),
                                                   HOFFSET(data, truncation_error),
                                                   HOFFSET(data, wall_time)
        };
        std::array<size_t, NFIELDS> dst_sizes = {
                sizeof(data::iteration),
                sizeof(data::chain_length),
                sizeof(data::position),
                sizeof(data::chi),
                sizeof(data::chi_max),
                sizeof(data::energy_mpo), sizeof(data::energy_ham), sizeof(data::energy_mom),
                sizeof(data::variance_mpo), sizeof(data::variance_ham), sizeof(data::variance_mom),
                sizeof(data::entanglement_entropy), sizeof(data::truncation_error), sizeof(data::wall_time)
        };
        std::array<const char *, NFIELDS> field_names = {"sweep",
                "chain_length",
                "iteration",
                "chi",
                "chi_max",
                "energy_mpo","energy_ham","energy_mom",
                "variance_mpo","variance_ham","variance_mom",
                "entanglement_entropy",
                "truncation_error",
                "wall_time"
        };

        std::array<hid_t, NFIELDS> field_types = {H5T_NATIVE_INT,
                                                  H5T_NATIVE_INT,
                                                  H5T_NATIVE_INT,
                                                  H5T_NATIVE_INT,
                                                  H5T_NATIVE_INT,
                                                  H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                  H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                  H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE};
        hsize_t chunk_size = 100;
        void *fill_data = nullptr;
        int compress = 0;
    };
public:
    class_table_dmrg() = default;
    meta_struct meta;
    std::vector<data> buffer;

};


class class_table_tebd{
private:
    struct data {
        int     iteration;
        int     chi;
        int     chi_max;
        double  time_step;
        double  energy_mpo; double  energy_ham; double  energy_mom;
        double  variance_mpo; double  variance_ham; double  variance_mom;
        double  entanglement_entropy;
        double  truncation_error;
        double  phys_time;
        double  wall_time;

        data(
                int iteration_,
                int chi_,
                int chi_max_,
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
        size_t dst_size = sizeof(data);
        std::array<size_t, NFIELDS> dst_offsets = {HOFFSET(data, iteration),
                                                   HOFFSET(data, chi),
                                                   HOFFSET(data, chi_max),
                                                   HOFFSET(data, time_step),
                                                   HOFFSET(data, energy_mpo), HOFFSET(data, energy_ham), HOFFSET(data, energy_mom),
                                                   HOFFSET(data, variance_mpo), HOFFSET(data, variance_ham), HOFFSET(data, variance_mom),
                                                   HOFFSET(data, entanglement_entropy),
                                                   HOFFSET(data, truncation_error),
                                                   HOFFSET(data, phys_time),
                                                   HOFFSET(data, wall_time)
        };
        std::array<size_t, NFIELDS> dst_sizes = {
                sizeof(data::iteration),
                sizeof(data::chi),
                sizeof(data::chi_max),
                sizeof(data::time_step),
                sizeof(data::energy_mpo), sizeof(data::energy_ham), sizeof(data::energy_mom),
                sizeof(data::variance_mpo), sizeof(data::variance_ham), sizeof(data::variance_mom),
                sizeof(data::entanglement_entropy),
                sizeof(data::truncation_error),
                sizeof(data::phys_time),
                sizeof(data::wall_time)
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
                                                  H5T_NATIVE_INT,
                                                  H5T_NATIVE_INT,
                                                  H5T_NATIVE_DOUBLE,
                                                  H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                  H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                  H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                                                  H5T_NATIVE_DOUBLE};
        hsize_t chunk_size = 100;
        void *fill_data = nullptr;
        int compress = 0;
    };
public:
    class_table_tebd() = default;
    meta_struct meta;
    std::vector<data> buffer;
};



class class_table_finite_chain{
private:
    struct data {
        int     sweeps;
        int     chain_length;
        int     position;
        double  energy;
        double  variance;

        data(   int sweeps_,
                int chain_length_,
                int position_,
                double energy_,
                double variance_):
                sweeps(sweeps_),
                chain_length(chain_length_),
                position(position_),
                energy(energy_),
                variance(variance_)
        {}
    };
    struct meta_struct {
        constexpr static hsize_t NFIELDS = 5;
        size_t dst_size = sizeof(data);
        std::array<size_t, NFIELDS> dst_offsets = {HOFFSET(data, sweeps),
                                                   HOFFSET(data, chain_length),
                                                   HOFFSET(data, position),
                                                   HOFFSET(data, energy),
                                                   HOFFSET(data, variance)};
        std::array<size_t, NFIELDS> dst_sizes = {
                sizeof(data::sweeps),
                sizeof(data::chain_length),
                sizeof(data::position),
                sizeof(data::energy),
                sizeof(data::variance)
        };
        std::array<const char *, NFIELDS> field_names = {"sweeps",
                                                         "chain_length",
                                                         "iteration",
                                                         "energy",
                                                         "variance"
        };

        std::array<hid_t, NFIELDS> field_types = {H5T_NATIVE_INT,
                                                  H5T_NATIVE_INT,
                                                  H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE};
        hsize_t chunk_size = 100;
        void *fill_data = nullptr;
        int compress = 0;
    };
public:
    class_table_finite_chain() = default;
    meta_struct meta;
    std::vector<data> buffer;

};


//Profiling table definition
class class_table_profiling{
private:
    struct data{
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

        data(int    iteration_,
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
        size_t           dst_size                           = sizeof (data);
        std::array       <size_t,NFIELDS>       dst_offsets = {HOFFSET(data, iteration),
                                                               HOFFSET(data, t_tot), HOFFSET(data, t_opt), HOFFSET(data, t_sim),
                                                               HOFFSET(data, t_svd), HOFFSET(data, t_env), HOFFSET(data, t_evo),
                                                               HOFFSET(data, t_udt), HOFFSET(data, t_sto), HOFFSET(data, t_ste),
                                                               HOFFSET(data, t_prt), HOFFSET(data, t_obs), HOFFSET(data, t_mps),
                                                               HOFFSET(data, t_chi)
        };
        std::array       <size_t,NFIELDS>       dst_sizes   = {sizeof(data::iteration),
                                                               sizeof(data::t_tot), sizeof(data::t_opt), sizeof(data::t_sim),
                                                               sizeof(data::t_svd), sizeof(data::t_env), sizeof(data::t_evo),
                                                               sizeof(data::t_udt), sizeof(data::t_sto), sizeof(data::t_ste),
                                                               sizeof(data::t_prt), sizeof(data::t_obs), sizeof(data::t_mps),
                                                               sizeof(data::t_chi)
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

        hsize_t          chunk_size                         = 100;
        void             *fill_data                         = nullptr;
        int              compress                           = 0;
    };
public:
    class_table_profiling() = default;
    meta_struct meta;
    std::vector<data> buffer;
};


#endif //DMRG_TABLE_TYPES_H
