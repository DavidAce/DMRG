//
// Created by david on 2018-05-24.
//

#pragma once

#include <array>
#include <h5pp/details/h5ppHid.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <simulation/class_simulation_status.h>
#include <vector>

class h5pp_table_measurements_finite {
    public:
    static inline h5pp::hid::h5t h5_type;

    struct table {
        uint64_t              iter;
        uint64_t              step;
        uint64_t              position;
        uint64_t              length;
        long                  bond_dimension_midchain;
        long                  bond_dimension_current;
        long                  bond_dimension_limit;
        long                  bond_dimension_maximum;
        double                entanglement_entropy_midchain;
        double                entanglement_entropy_current;
        double                norm;
        double                energy;
        double                energy_per_site;
        double                energy_variance;
        double                energy_variance_per_site;
        double                energy_variance_lowest;
        double                energy_variance_per_site_lowest;
        std::array<double, 3> spin_components;
        double                truncation_error;
        double                wall_time;
    };

    h5pp_table_measurements_finite() { register_table_type(); }
    static void register_table_type() {
        if(h5_type.valid()) return;
        std::array<hsize_t, 1> spin_dims    = {3};
        h5pp::hid::h5t         H5_SPIN_TYPE = H5Tarray_create(H5T_NATIVE_DOUBLE, spin_dims.size(), spin_dims.data());

        h5_type = H5Tcreate(H5T_COMPOUND, sizeof(table));
        H5Tinsert(h5_type, "iter", HOFFSET(table, iter), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "step", HOFFSET(table, step), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "position", HOFFSET(table, position), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "length", HOFFSET(table, length), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "bond_dimension_midchain", HOFFSET(table, bond_dimension_midchain), H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "bond_dimension_current", HOFFSET(table, bond_dimension_current), H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "bond_dimension_limit", HOFFSET(table, bond_dimension_limit), H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "bond_dimension_maximum", HOFFSET(table, bond_dimension_maximum), H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "entanglement_entropy_midchain", HOFFSET(table, entanglement_entropy_midchain), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "entanglement_entropy_current", HOFFSET(table, entanglement_entropy_current), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "norm", HOFFSET(table, norm), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy", HOFFSET(table, energy), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_per_site", HOFFSET(table, energy_per_site), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_variance", HOFFSET(table, energy_variance), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_variance_per_site", HOFFSET(table, energy_variance_per_site), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_variance_lowest", HOFFSET(table, energy_variance_lowest), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_variance_per_site_lowest", HOFFSET(table, energy_variance_per_site_lowest), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "spin_components", HOFFSET(table, spin_components), H5_SPIN_TYPE);
        H5Tinsert(h5_type, "truncation_error", HOFFSET(table, truncation_error), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "wall_time", HOFFSET(table, wall_time), H5T_NATIVE_DOUBLE);
    }
};

class h5pp_table_measurements_infinite {
    public:
    static inline h5pp::hid::h5t h5_type;

    struct table {
        uint64_t iter;
        uint64_t step;
        uint64_t position;
        uint64_t length;
        long     bond_dimension;
        long     bond_dimension_limit;
        long     bond_dimension_maximum;
        double   entanglement_entropy;
        double   norm;
        double   energy_mpo;
        double   energy_per_site_mpo;
        double   energy_per_site_ham;
        double   energy_per_site_mom;
        double   energy_variance_mpo;
        double   energy_variance_per_site_mpo;
        double   energy_variance_per_site_ham;
        double   energy_variance_per_site_mom;
        double   truncation_error;
        double   wall_time;
        double   phys_time;
        double   time_step; // Only used in itebd
    };

    h5pp_table_measurements_infinite() { register_table_type(); }
    static void register_table_type() {
        if(h5_type.valid()) return;
        h5_type = H5Tcreate(H5T_COMPOUND, sizeof(table));
        H5Tinsert(h5_type, "iter", HOFFSET(table, iter), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "step", HOFFSET(table, step), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "position", HOFFSET(table, position), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "length", HOFFSET(table, length), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "bond_dimension", HOFFSET(table, bond_dimension), H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "bond_dimension_limit", HOFFSET(table, bond_dimension_limit), H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "bond_dimension_maximum", HOFFSET(table, bond_dimension_maximum), H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "entanglement_entropy", HOFFSET(table, entanglement_entropy), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "norm", HOFFSET(table, norm), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_mpo", HOFFSET(table, energy_mpo), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_per_site_mpo", HOFFSET(table, energy_per_site_mpo), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_per_site_ham", HOFFSET(table, energy_per_site_ham), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_per_site_mom", HOFFSET(table, energy_per_site_mom), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_variance_mpo", HOFFSET(table, energy_variance_mpo), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_variance_per_site_mpo", HOFFSET(table, energy_variance_per_site_mpo), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_variance_per_site_ham", HOFFSET(table, energy_variance_per_site_ham), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_variance_per_site_mom", HOFFSET(table, energy_variance_per_site_mom), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "truncation_error", HOFFSET(table, truncation_error), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "wall_time", HOFFSET(table, wall_time), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "phys_time", HOFFSET(table, phys_time), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "time_step", HOFFSET(table, time_step), H5T_NATIVE_DOUBLE);
    }
};

class h5pp_table_profiling {
    public:
    static inline h5pp::hid::h5t h5_type;
    struct table {
        uint64_t iter      = 0;
        uint64_t step      = 0;
        uint64_t position  = 0;
        double   t_tot     = 0;
        double   t_pre     = 0;
        double   t_pos     = 0;
        double   t_sim     = 0;
        double   t_con     = 0;
        double   t_eig     = 0;
        double   t_svd     = 0;
        double   t_ham     = 0;
        double   t_hsq     = 0;
        double   t_mps     = 0;
        double   t_mpo     = 0;
        double   t_opt     = 0;
        double   t_evo     = 0;
        double   t_env     = 0;
        double   t_ent     = 0;
        double   t_ene     = 0;
        double   t_var     = 0;
        double   t_prj     = 0;
        double   t_chk     = 0;
        double   t_hdf     = 0;
        double   t_ene_ham = 0;
        double   t_ene_mom = 0;
        double   t_var_ham = 0;
        double   t_var_mom = 0;
        double   t_vH2v    = 0;
        double   t_vHv     = 0;
        double   t_vH2     = 0;
        double   t_vH      = 0;
        double   t_op      = 0;
    };

    h5pp_table_profiling() { register_table_type(); }
    static void register_table_type() {
        if(h5_type.valid()) return;
        h5_type = H5Tcreate(H5T_COMPOUND, sizeof(table));
        H5Tinsert(h5_type, "iter", HOFFSET(table, iter), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "step", HOFFSET(table, step), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "position", HOFFSET(table, position), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "t_tot", HOFFSET(table, t_tot), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "t_pre", HOFFSET(table, t_pre), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "t_pos", HOFFSET(table, t_pos), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "t_sim", HOFFSET(table, t_sim), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "t_con", HOFFSET(table, t_con), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "t_eig", HOFFSET(table, t_eig), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "t_svd", HOFFSET(table, t_svd), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "t_ham", HOFFSET(table, t_ham), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "t_hsq", HOFFSET(table, t_hsq), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "t_mps", HOFFSET(table, t_mps), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "t_mpo", HOFFSET(table, t_mpo), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "t_opt", HOFFSET(table, t_opt), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "t_evo", HOFFSET(table, t_evo), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "t_env", HOFFSET(table, t_env), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "t_ent", HOFFSET(table, t_ent), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "t_ene", HOFFSET(table, t_ene), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "t_var", HOFFSET(table, t_var), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "t_prj", HOFFSET(table, t_prj), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "t_chk", HOFFSET(table, t_chk), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "t_hdf", HOFFSET(table, t_hdf), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "t_ene_ham", HOFFSET(table, t_ene_ham), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "t_ene_mom", HOFFSET(table, t_ene_mom), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "t_var_ham", HOFFSET(table, t_var_ham), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "t_var_mom", HOFFSET(table, t_var_mom), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "t_vH2v", HOFFSET(table, t_vH2v), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "t_vHv", HOFFSET(table, t_vHv), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "t_vH2", HOFFSET(table, t_vH2), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "t_vH", HOFFSET(table, t_vH), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "t_op", HOFFSET(table, t_op), H5T_NATIVE_DOUBLE);
    }
};

class h5pp_table_sim_status {
    public:
    static inline h5pp::hid::h5t h5_type;
    using table = class_simulation_status;

    h5pp_table_sim_status() { register_table_type(); }
    static void register_table_type() {
        if(h5_type.valid()) return;
        h5_type = H5Tcreate(H5T_COMPOUND, sizeof(table));
        H5Tinsert(h5_type, "iter", HOFFSET(table, iter), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "step", HOFFSET(table, step), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "position", HOFFSET(table, position), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "direction", HOFFSET(table, direction), H5T_NATIVE_INT);
        H5Tinsert(h5_type, "num_resets", HOFFSET(table, num_resets), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "min_iters", HOFFSET(table, min_iters), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "chi_max", HOFFSET(table, chi_max), H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "chi_lim", HOFFSET(table, chi_lim), H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "energy_min", HOFFSET(table, energy_min), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_max", HOFFSET(table, energy_max), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_target", HOFFSET(table, energy_target), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_ubound", HOFFSET(table, energy_ubound), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_lbound", HOFFSET(table, energy_lbound), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_dens", HOFFSET(table, energy_dens), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_dens_target", HOFFSET(table, energy_dens_target), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_dens_window", HOFFSET(table, energy_dens_window), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "phys_time", HOFFSET(table, phys_time), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "wall_time", HOFFSET(table, wall_time), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "simu_time", HOFFSET(table, simu_time), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "delta_t", HOFFSET(table, delta_t), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "simulation_has_stuck_for", HOFFSET(table, simulation_has_stuck_for), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "entanglement_saturated_for", HOFFSET(table, entanglement_saturated_for), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "variance_mpo_saturated_for", HOFFSET(table, variance_mpo_saturated_for), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "variance_ham_saturated_for", HOFFSET(table, variance_ham_saturated_for), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "variance_mom_saturated_for", HOFFSET(table, variance_mom_saturated_for), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "simulation_has_finished", HOFFSET(table, simulation_has_finished), H5T_NATIVE_UINT8);
        H5Tinsert(h5_type, "simulation_has_converged", HOFFSET(table, simulation_has_converged), H5T_NATIVE_UINT8);
        H5Tinsert(h5_type, "simulation_has_saturated", HOFFSET(table, simulation_has_saturated), H5T_NATIVE_UINT8);
        H5Tinsert(h5_type, "simulation_has_succeeded", HOFFSET(table, simulation_has_succeeded), H5T_NATIVE_UINT8);
        H5Tinsert(h5_type, "simulation_has_got_stuck", HOFFSET(table, simulation_has_got_stuck), H5T_NATIVE_UINT8);
        H5Tinsert(h5_type, "simulation_has_to_stop", HOFFSET(table, simulation_has_to_stop), H5T_NATIVE_UINT8);
        H5Tinsert(h5_type, "chi_lim_has_reached_chi_max", HOFFSET(table, chi_lim_has_reached_chi_max), H5T_NATIVE_UINT8);
        H5Tinsert(h5_type, "entanglement_has_converged", HOFFSET(table, entanglement_has_converged), H5T_NATIVE_UINT8);
        H5Tinsert(h5_type, "entanglement_has_saturated", HOFFSET(table, entanglement_has_saturated), H5T_NATIVE_UINT8);
        H5Tinsert(h5_type, "variance_mpo_has_converged", HOFFSET(table, variance_mpo_has_converged), H5T_NATIVE_UINT8);
        H5Tinsert(h5_type, "variance_mpo_has_saturated", HOFFSET(table, variance_mpo_has_saturated), H5T_NATIVE_UINT8);
        H5Tinsert(h5_type, "variance_ham_has_converged", HOFFSET(table, variance_ham_has_converged), H5T_NATIVE_UINT8);
        H5Tinsert(h5_type, "variance_ham_has_saturated", HOFFSET(table, variance_ham_has_saturated), H5T_NATIVE_UINT8);
        H5Tinsert(h5_type, "variance_mom_has_converged", HOFFSET(table, variance_mom_has_converged), H5T_NATIVE_UINT8);
        H5Tinsert(h5_type, "variance_mom_has_saturated", HOFFSET(table, variance_mom_has_saturated), H5T_NATIVE_UINT8);
        H5Tinsert(h5_type, "time_step_has_converged", HOFFSET(table, time_step_has_converged), H5T_NATIVE_UINT8);
    }
};

class h5pp_table_memory_usage {
    public:
    static inline h5pp::hid::h5t h5_type;

    struct table {
        uint64_t iter;
        uint64_t step;
        double   rss;
        double   hwm;
        double   vm;
    };

    h5pp_table_memory_usage() { register_table_type(); }
    static void register_table_type() {
        if(h5_type.valid()) return;
        h5_type = H5Tcreate(H5T_COMPOUND, sizeof(table));
        H5Tinsert(h5_type, "iter", HOFFSET(table, iter), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "step", HOFFSET(table, step), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "rss", HOFFSET(table, rss), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "hwm", HOFFSET(table, hwm), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "vm", HOFFSET(table, vm), H5T_NATIVE_DOUBLE);
    }
};
