//
// Created by david on 2018-05-24.
//

#ifndef DMRG_LOG_TYPES_H
#define DMRG_LOG_TYPES_H

#include <vector>
#include <array>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <simulation/class_simulation_status.h>



class class_log_finite_dmrg_measurements {
public:
    struct data_struct {
        int     step;
        int     iteration;
        int     position;
        int     length;
        int     bond_dimension_midchain;
        int     bond_dimension_current;
        double  entanglement_entropy_midchain;
        double  entanglement_entropy_current;
        double  norm;
        double  energy;
        double  energy_per_site;
        double  energy_variance;
        double  energy_variance_per_site;
        double  spin_component_sx;
        double  spin_component_sy;
        double  spin_component_sz;
        double  truncation_error;
        double  wall_time;
    };

    struct meta_struct {
        constexpr static hsize_t NFIELDS = 18;
        size_t dst_size = sizeof(data_struct);
        std::array<size_t, NFIELDS> dst_offsets = {
                HOFFSET(data_struct, step),
                HOFFSET(data_struct, iteration),
                HOFFSET(data_struct, position),
                HOFFSET(data_struct, length),
                HOFFSET(data_struct, bond_dimension_midchain),
                HOFFSET(data_struct, bond_dimension_current),
                HOFFSET(data_struct, entanglement_entropy_midchain),
                HOFFSET(data_struct, entanglement_entropy_current),
                HOFFSET(data_struct, norm),
                HOFFSET(data_struct, energy),
                HOFFSET(data_struct, energy_per_site),
                HOFFSET(data_struct, energy_variance),
                HOFFSET(data_struct, energy_variance_per_site),
                HOFFSET(data_struct, spin_component_sx),
                HOFFSET(data_struct, spin_component_sy),
                HOFFSET(data_struct, spin_component_sz),
                HOFFSET(data_struct, truncation_error),
                HOFFSET(data_struct, wall_time)
        };
        std::array<size_t, NFIELDS> dst_sizes = {
                sizeof(data_struct::step),
                sizeof(data_struct::iteration),
                sizeof(data_struct::position),
                sizeof(data_struct::length),
                sizeof(data_struct::bond_dimension_midchain),
                sizeof(data_struct::bond_dimension_current),
                sizeof(data_struct::entanglement_entropy_midchain),
                sizeof(data_struct::entanglement_entropy_current),
                sizeof(data_struct::norm),
                sizeof(data_struct::energy),
                sizeof(data_struct::energy_per_site),
                sizeof(data_struct::energy_variance),
                sizeof(data_struct::energy_variance_per_site),
                sizeof(data_struct::spin_component_sx),
                sizeof(data_struct::spin_component_sy),
                sizeof(data_struct::spin_component_sz),
                sizeof(data_struct::truncation_error),
                sizeof(data_struct::wall_time),
        };
        std::array<const char *, NFIELDS> field_names = {
                "step",
                "iteration",
                "position",
                "length",
                "bond_dimension_midchain",
                "bond_dimension_current",
                "entanglement_entropy_midchain",
                "entanglement_entropy_current",
                "norm",
                "energy",
                "energy_per_site",
                "energy_variance",
                "energy_variance_per_site",
                "spin_component_sx",
                "spin_component_sy",
                "spin_component_sz",
                "truncation_error",
                "wall_time"
        };

        std::array<hid_t, NFIELDS> field_types = {
                H5T_NATIVE_INT,
                H5T_NATIVE_INT,
                H5T_NATIVE_INT,
                H5T_NATIVE_INT,
                H5T_NATIVE_INT,
                H5T_NATIVE_INT,
                H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,
                H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE};
        hsize_t chunk_size = 32;
        void *fill_data = nullptr;
        int compress = 0;
    };
public:
    class_log_finite_dmrg_measurements() = default;
    meta_struct meta;
    std::vector<data_struct> buffer;

};


class class_log_tebd{
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
                                                         "energy","energy_per_site_ham","energy_per_site_mom",
                                                         "variance_mpo","variance_ham","variance_mom",
                                                         "entanglement_entropy_midchain",
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
        hsize_t chunk_size = 32;
        void *fill_data = nullptr;
        int compress = 0;
    };
public:
    class_log_tebd() = default;
    meta_struct meta;
    std::vector<data_struct> buffer;
};



class class_log_profiling{
public:
    struct data_struct{
        int    step      = 0;
        int    iteration = 0;
        int    position  = 0;
        double t_tot = 0;
        double t_run = 0;
        double t_eig = 0;
        double t_svd = 0;
        double t_ene = 0;
        double t_var = 0;
        double t_ent = 0;
        double t_hdf = 0;
        double t_prj = 0;
        double t_opt = 0;
        double t_chk = 0;
        double t_ene_mpo = 0;
        double t_ene_ham = 0;
        double t_ene_mom = 0;
        double t_var_mpo = 0;
        double t_var_ham = 0;
        double t_var_mom = 0;
        double t_env = 0;
        double t_evo = 0;
        double t_udt = 0;
        double t_ste = 0;
        double t_prt = 0;
        double t_obs = 0;
        double t_mps = 0;
        double t_chi = 0;


    };
private:
    struct meta_struct{
        constexpr static hsize_t                NFIELDS     = 28;
        size_t           dst_size                           = sizeof (data_struct);
        std::array       <size_t,NFIELDS>       dst_offsets = {
                HOFFSET(data_struct, step),
                HOFFSET(data_struct, iteration),
                HOFFSET(data_struct, position),
                HOFFSET(data_struct, t_tot),
                HOFFSET(data_struct, t_run),
                HOFFSET(data_struct, t_eig),
                HOFFSET(data_struct, t_svd),
                HOFFSET(data_struct, t_ene),
                HOFFSET(data_struct, t_var),
                HOFFSET(data_struct, t_ent),
                HOFFSET(data_struct, t_hdf),
                HOFFSET(data_struct, t_prj),
                HOFFSET(data_struct, t_opt),
                HOFFSET(data_struct, t_chk),
                HOFFSET(data_struct, t_ene_mpo),
                HOFFSET(data_struct, t_ene_ham),
                HOFFSET(data_struct, t_ene_mom),
                HOFFSET(data_struct, t_var_mpo),
                HOFFSET(data_struct, t_var_ham),
                HOFFSET(data_struct, t_var_mom),
                HOFFSET(data_struct, t_env),
                HOFFSET(data_struct, t_evo),
                HOFFSET(data_struct, t_udt),
                HOFFSET(data_struct, t_ste),
                HOFFSET(data_struct, t_prt),
                HOFFSET(data_struct, t_obs),
                HOFFSET(data_struct, t_mps),
                HOFFSET(data_struct, t_chi),
        };
        std::array       <size_t,NFIELDS>       dst_sizes   = {
                sizeof(data_struct::step),
                sizeof(data_struct::iteration),
                sizeof(data_struct::position),
                sizeof(data_struct::t_tot),
                sizeof(data_struct::t_run),
                sizeof(data_struct::t_eig),
                sizeof(data_struct::t_svd),
                sizeof(data_struct::t_ene),
                sizeof(data_struct::t_var),
                sizeof(data_struct::t_ent),
                sizeof(data_struct::t_hdf),
                sizeof(data_struct::t_prj),
                sizeof(data_struct::t_opt),
                sizeof(data_struct::t_chk),
                sizeof(data_struct::t_ene_mpo),
                sizeof(data_struct::t_ene_ham),
                sizeof(data_struct::t_ene_mom),
                sizeof(data_struct::t_var_mpo),
                sizeof(data_struct::t_var_ham),
                sizeof(data_struct::t_var_mom),
                sizeof(data_struct::t_env),
                sizeof(data_struct::t_evo),
                sizeof(data_struct::t_udt),
                sizeof(data_struct::t_ste),
                sizeof(data_struct::t_prt),
                sizeof(data_struct::t_obs),
                sizeof(data_struct::t_mps),
                sizeof(data_struct::t_chi)
        };
        std::array       <const char*,NFIELDS>  field_names = {
                "step",
                "iteration",
                "position",
                "t_tot",
                "t_run",
                "t_eig",
                "t_svd",
                "t_ene",
                "t_var",
                "t_ent",
                "t_hdf",
                "t_prj",
                "t_opt",
                "t_chk",
                "t_ene_mpo",
                "t_ene_ham",
                "t_ene_mom",
                "t_var_mpo",
                "t_var_ham",
                "t_var_mom",
                "t_env",
                "t_evo",
                "t_udt",
                "t_ste",
                "t_prt",
                "t_obs",
                "t_mps",
                "t_chi"
        };

        std::array       <hid_t,NFIELDS>        field_types = {
                H5T_NATIVE_INT, H5T_NATIVE_INT, H5T_NATIVE_INT,
                H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,
                H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,
                H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,
                H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,
                H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE, H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE,H5T_NATIVE_DOUBLE
        };

        hsize_t          chunk_size                         = 32;
        void             *fill_data                         = nullptr;
        int              compress                           = 0;
    };
public:
    class_log_profiling() = default;
    meta_struct meta;
    std::vector<data_struct> buffer;
};


class class_log_simulation_status{
private:
    struct meta_struct{
        constexpr static hsize_t                NFIELDS     = 38;
        size_t           dst_size                           = sizeof (status_data);
        std::array       <size_t,NFIELDS>       dst_offsets =
                {
                    HOFFSET(status_data, iteration                     ),
                    HOFFSET(status_data, moves                         ),
                    HOFFSET(status_data, step                          ),
                    HOFFSET(status_data, position                      ),
                    HOFFSET(status_data, num_resets                    ),
                    HOFFSET(status_data, chi_temp                      ),
                    HOFFSET(status_data, chi_max                       ),
                    HOFFSET(status_data, min_sweeps                    ),
                    HOFFSET(status_data, energy_min                    ),
                    HOFFSET(status_data, energy_max                    ),
                    HOFFSET(status_data, energy_target                 ),
                    HOFFSET(status_data, energy_ubound                 ),
                    HOFFSET(status_data, energy_lbound                 ),
                    HOFFSET(status_data, energy_dens                   ),
                    HOFFSET(status_data, energy_dens_target            ),
                    HOFFSET(status_data, energy_dens_window            ),
                    HOFFSET(status_data, phys_time                     ),
                    HOFFSET(status_data, wall_time                     ),
                    HOFFSET(status_data, simu_time                     ),
                    HOFFSET(status_data, delta_t                       ),
                    HOFFSET(status_data, time_step_has_converged       ),
                    HOFFSET(status_data, simulation_has_converged      ),
                    HOFFSET(status_data, simulation_has_saturated      ),
                    HOFFSET(status_data, simulation_has_succeeded      ),
                    HOFFSET(status_data, simulation_has_to_stop        ),
                    HOFFSET(status_data, bond_dimension_has_reached_max),
                    HOFFSET(status_data, entanglement_has_converged    ),
                    HOFFSET(status_data, entanglement_has_saturated    ),
                    HOFFSET(status_data, variance_mpo_has_converged    ),
                    HOFFSET(status_data, variance_mpo_has_saturated    ),
                    HOFFSET(status_data, variance_ham_has_converged    ),
                    HOFFSET(status_data, variance_ham_has_saturated    ),
                    HOFFSET(status_data, variance_mom_has_converged    ),
                    HOFFSET(status_data, variance_mom_has_saturated    ),
                    HOFFSET(status_data, entanglement_saturated_for    ),
                    HOFFSET(status_data, variance_mpo_saturated_for    ),
                    HOFFSET(status_data, variance_ham_saturated_for    ),
                    HOFFSET(status_data, variance_mom_saturated_for    )

                };

        std::array       <size_t,NFIELDS>       dst_sizes   = {
                sizeof(status_data::iteration                     ),
                sizeof(status_data::moves                         ),
                sizeof(status_data::step                          ),
                sizeof(status_data::position                      ),
                sizeof(status_data::num_resets                    ),
                sizeof(status_data::chi_temp                      ),
                sizeof(status_data::chi_max                       ),
                sizeof(status_data::min_sweeps                    ),
                sizeof(status_data::energy_min                    ),
                sizeof(status_data::energy_max                    ),
                sizeof(status_data::energy_target                 ),
                sizeof(status_data::energy_ubound                 ),
                sizeof(status_data::energy_lbound                 ),
                sizeof(status_data::energy_dens                   ),
                sizeof(status_data::energy_dens_target            ),
                sizeof(status_data::energy_dens_window            ),
                sizeof(status_data::phys_time                     ),
                sizeof(status_data::wall_time                     ),
                sizeof(status_data::simu_time                     ),
                sizeof(status_data::delta_t                       ),
                sizeof(status_data::time_step_has_converged       ),
                sizeof(status_data::simulation_has_converged      ),
                sizeof(status_data::simulation_has_saturated      ),
                sizeof(status_data::simulation_has_succeeded      ),
                sizeof(status_data::simulation_has_to_stop        ),
                sizeof(status_data::bond_dimension_has_reached_max),
                sizeof(status_data::entanglement_has_converged    ),
                sizeof(status_data::entanglement_has_saturated    ),
                sizeof(status_data::variance_mpo_has_converged    ),
                sizeof(status_data::variance_mpo_has_saturated    ),
                sizeof(status_data::variance_ham_has_converged    ),
                sizeof(status_data::variance_ham_has_saturated    ),
                sizeof(status_data::variance_mom_has_converged    ),
                sizeof(status_data::variance_mom_has_saturated    ),
                sizeof(status_data::entanglement_saturated_for    ),
                sizeof(status_data::variance_mpo_saturated_for    ),
                sizeof(status_data::variance_ham_saturated_for    ),
                sizeof(status_data::variance_mom_saturated_for    )
        };

        std::array       <const char*,NFIELDS>  field_names =
                {
                    "iteration",
                    "moves",
                    "step",
                    "position",
                    "num_resets",
                    "chi_temp",
                    "chi_max",
                    "min_sweeps",
                    "energy_min",
                    "energy_max",
                    "energy_target",
                    "energy_ubound",
                    "energy_lbound",
                    "energy_dens",
                    "energy_dens_target",
                    "energy_dens_window",
                    "phys_time",
                    "wall_time",
                    "simu_time",
                    "delta_t",
                    "time_step_has_converged",
                    "simulation_has_converged",
                    "simulation_has_saturated",
                    "simulation_has_succeeded",
                    "simulation_has_to_stop",
                    "bond_dimension_has_reached_max",
                    "entanglement_has_converged",
                    "entanglement_has_saturated",
                    "variance_mpo_has_converged",
                    "variance_mpo_has_saturated",
                    "variance_ham_has_converged",
                    "variance_ham_has_saturated",
                    "variance_mom_has_converged",
                    "variance_mom_has_saturated",
                    "entanglement_saturated_for",
                    "variance_mpo_saturated_for",
                    "variance_ham_saturated_for",
                    "variance_mom_saturated_for"
                };

        std::array       <hid_t,NFIELDS>        field_types =
                {
                        H5T_NATIVE_UINT,
                        H5T_NATIVE_UINT,
                        H5T_NATIVE_UINT,
                        H5T_NATIVE_UINT,
                        H5T_NATIVE_UINT,
                        H5T_NATIVE_LONG,
                        H5T_NATIVE_LONG,
                        H5T_NATIVE_UINT,
                        H5T_NATIVE_DOUBLE,
                        H5T_NATIVE_DOUBLE,
                        H5T_NATIVE_DOUBLE,
                        H5T_NATIVE_DOUBLE,
                        H5T_NATIVE_DOUBLE,
                        H5T_NATIVE_DOUBLE,
                        H5T_NATIVE_DOUBLE,
                        H5T_NATIVE_DOUBLE,
                        H5T_NATIVE_DOUBLE,
                        H5T_NATIVE_DOUBLE,
                        H5T_NATIVE_DOUBLE,
                        H5T_NATIVE_DOUBLE,
                        H5T_NATIVE_HBOOL,
                        H5T_NATIVE_HBOOL,
                        H5T_NATIVE_HBOOL,
                        H5T_NATIVE_HBOOL,
                        H5T_NATIVE_HBOOL,
                        H5T_NATIVE_HBOOL,
                        H5T_NATIVE_HBOOL,
                        H5T_NATIVE_HBOOL,
                        H5T_NATIVE_HBOOL,
                        H5T_NATIVE_HBOOL,
                        H5T_NATIVE_HBOOL,
                        H5T_NATIVE_HBOOL,
                        H5T_NATIVE_HBOOL,
                        H5T_NATIVE_HBOOL,
                        H5T_NATIVE_UINT,
                        H5T_NATIVE_UINT,
                        H5T_NATIVE_UINT,
                        H5T_NATIVE_UINT
                };

        hsize_t          chunk_size                         = 32;
        void             *fill_data                         = nullptr;
        int              compress                           = 0;
    };
public:
    class_log_simulation_status() = default;
    meta_struct meta;
    std::vector<status_data> buffer;
};




#endif //DMRG_LOG_TYPES_H
