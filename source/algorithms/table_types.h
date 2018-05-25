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
        double  energy1; double  energy2; double  energy3;
        double  energy4; double  energy5; double  energy6;
        double  variance1; double  variance2; double  variance3;
        double  variance4; double  variance5; double  variance6;
        double  entanglement_entropy;
        double  truncation_error;
        double  wall_time;

        data(
        int iteration_,
        int chain_length_,
        int position_,
        int chi_,
        int chi_max_,
        double energy1_, double energy2_, double energy3_,
        double energy4_, double energy5_, double energy6_,
        double variance1_, double variance2_, double variance3_,
        double variance4_, double variance5_, double variance6_,
        double entropy_,
        double truncation_error_,
        double wall_time_) :
        iteration(iteration_),
        chain_length(chain_length_),
        position(position_),
        chi(chi_),
        chi_max(chi_max_),
        energy1(energy1_), energy2(energy2_), energy3(energy3_),
        energy4(energy4_), energy5(energy5_), energy6(energy6_),
        variance1(variance1_), variance2(variance2_), variance3(variance3_),
        variance4(variance4_), variance5(variance5_), variance6(variance6_),
        entanglement_entropy(entropy_),
        truncation_error(truncation_error_),
        wall_time(wall_time_)
        {}
    };
    struct meta_struct {
        constexpr static hsize_t NFIELDS = 20;
        size_t dst_size = sizeof(data);
        std::array<size_t, NFIELDS> dst_offsets = {HOFFSET(data, iteration),
                                                   HOFFSET(data, chain_length),
                                                   HOFFSET(data, position),
                                                   HOFFSET(data, chi),
                                                   HOFFSET(data, chi_max),
                                                   HOFFSET(data, energy1), HOFFSET(data, energy2), HOFFSET(data, energy3),
                                                   HOFFSET(data, energy4), HOFFSET(data, energy5), HOFFSET(data, energy6),
                                                   HOFFSET(data, variance1), HOFFSET(data, variance2), HOFFSET(data, variance3),
                                                   HOFFSET(data, variance4), HOFFSET(data, variance5), HOFFSET(data, variance6),
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
                sizeof(data::energy1), sizeof(data::energy2), sizeof(data::energy3),
                sizeof(data::energy4), sizeof(data::energy5), sizeof(data::energy6),
                sizeof(data::variance1), sizeof(data::variance2), sizeof(data::variance3),
                sizeof(data::variance4), sizeof(data::variance5), sizeof(data::variance6),
                sizeof(data::entanglement_entropy), sizeof(data::truncation_error), sizeof(data::wall_time)
        };
        std::array<const char *, NFIELDS> field_names = {"sweep",
                "chain_length",
                "iteration",
                "chi;",
                "chi_max",
                "energy1","energy2","energy3",
                "energy4","energy5","energy6",
                "variance1","variance2","variance3",
                "variance4","variance5","variance6",
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
        double  energy1; double  energy2; double  energy3;
        double  energy4; double  energy5; double  energy6;
        double  variance1; double  variance2; double  variance3;
        double  variance4; double  variance5; double  variance6;
        double  entanglement_entropy;
        double  truncation_error;
        double  phys_time;
        double  wall_time;

        data(
                int iteration_,
                int chi_,
                int chi_max_,
                double time_step_,
                double energy1_, double energy2_, double energy3_,
                double energy4_, double energy5_, double energy6_,
                double variance1_, double variance2_, double variance3_,
                double variance4_, double variance5_, double variance6_,
                double entropy_,
                double truncation_error_,
                double phys_time_,
                double wall_time_) :
                iteration(iteration_),
                chi(chi_),
                chi_max(chi_max_),
                time_step(time_step_),
                energy1(energy1_), energy2(energy2_), energy3(energy3_),
                energy4(energy4_), energy5(energy5_), energy6(energy6_),
                variance1(variance1_), variance2(variance2_), variance3(variance3_),
                variance4(variance4_), variance5(variance5_), variance6(variance6_),
                entanglement_entropy(entropy_),
                truncation_error(truncation_error_),
                phys_time(phys_time_),
                wall_time(wall_time_)
        {}
    };
    struct meta_struct {
        constexpr static hsize_t NFIELDS = 20;
        size_t dst_size = sizeof(data);
        std::array<size_t, NFIELDS> dst_offsets = {HOFFSET(data, iteration),
                                                   HOFFSET(data, chi),
                                                   HOFFSET(data, chi_max),
                                                   HOFFSET(data, time_step),
                                                   HOFFSET(data, energy1), HOFFSET(data, energy2), HOFFSET(data, energy3),
                                                   HOFFSET(data, energy4), HOFFSET(data, energy5), HOFFSET(data, energy6),
                                                   HOFFSET(data, variance1), HOFFSET(data, variance2), HOFFSET(data, variance3),
                                                   HOFFSET(data, variance4), HOFFSET(data, variance5), HOFFSET(data, variance6),
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
                sizeof(data::energy1), sizeof(data::energy2), sizeof(data::energy3),
                sizeof(data::energy4), sizeof(data::energy5), sizeof(data::energy6),
                sizeof(data::variance1), sizeof(data::variance2), sizeof(data::variance3),
                sizeof(data::variance4), sizeof(data::variance5), sizeof(data::variance6),
                sizeof(data::entanglement_entropy),
                sizeof(data::truncation_error),
                sizeof(data::phys_time),
                sizeof(data::wall_time)
        };
        std::array<const char *, NFIELDS> field_names = {"sweep",
                                                         "chi",
                                                         "chi_max",
                                                         "time_step",
                                                         "energy1","energy2","energy3",
                                                         "energy4","energy5","energy6",
                                                         "variance1","variance2","variance3",
                                                         "variance4","variance5","variance6",
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
                                                               "t_chi"
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
