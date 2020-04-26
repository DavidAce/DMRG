//
// Created by david on 2018-10-01.
//

#pragma once

#include <any>
#include <array>
#include <h5pp/details/h5ppHid.h>
#include <hdf5.h>
#include <hdf5_hl.h>

class h5tb_sdual_trf_ising {
    public:
    static inline h5pp::hid::h5t h5_type;

    struct table {
        double   J_rnd            = 0;         /*!< Randomly distributed nearest neighbour coupling */
        double   h_rnd            = 0;         /*!< Randomly distributed on-site field */
        double   J_ptb            = 0;         /*!< Perturbation to the coupling, std::pow(J_rnd + J_ptb,1-alpha) */
        double   h_ptb            = 0;         /*!< Perturbation to the coupling, std::pow(J_rnd + J_ptb,1-alpha) */
        double   J_avg            = 0;         /*!< Average of J_rnd between all sites*/
        double   h_avg            = 0;         /*!< Average of h_rnd on all sites */
        double   J_mean           = 0;         /*!< Mean for the distrbution of J_rnd */
        double   h_mean           = 0;         /*!< Mean for the distrbution of h_rnd */
        double   J_sigma          = 0;         /*!< Standard deviation for the distribution of J_rnd */
        double   h_sigma          = 0;         /*!< Standard deviation for the distribution of h_rnd */
        double   lambda           = 0;         /*!< Factor involved in next-nearest neighbor interaction */
        double   delta            = 0;         /*!< Difference J_mean - h_mean  */
        uint64_t spin_dim         = 0;         /*!< Spin dimension */
        char     distribution[16] = "uniform"; /*!< The random distribution of J_rnd and h_rnd. Choose between lognormal, normal or uniform */
    };

    h5tb_sdual_trf_ising() { register_table_type(); }
    static void register_table_type() {
        if(h5_type.valid()) return;

        // Create a type for the char array from the template H5T_C_S1
        // The template describes a string with a single char.
        // Set the size with H5Tset_size, or h5pp::hdf5::setStringSize(...)
        h5pp::hid::h5t h5t_custom_string = H5Tcopy(H5T_C_S1);
        H5Tset_size(h5t_custom_string, 16);

        // Optionally set the null terminator '\0'
        H5Tset_strpad(h5t_custom_string, H5T_STR_NULLTERM);

        h5_type = H5Tcreate(H5T_COMPOUND, sizeof(table));
        H5Tinsert(h5_type, "J_rnd", HOFFSET(table, J_rnd), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_rnd", HOFFSET(table, h_rnd), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J_ptb", HOFFSET(table, J_ptb), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_ptb", HOFFSET(table, h_ptb), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J_avg", HOFFSET(table, J_avg), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_avg", HOFFSET(table, h_avg), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J_mean", HOFFSET(table, J_mean), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_mean", HOFFSET(table, h_mean), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J_sigma", HOFFSET(table, J_sigma), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_sigma", HOFFSET(table, h_sigma), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "lambda", HOFFSET(table, lambda), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "delta", HOFFSET(table, delta), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "spin_dim", HOFFSET(table, spin_dim), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "distribution", HOFFSET(table, distribution), h5t_custom_string);
    }
};

class h5tb_tf_ising {
    public:
    static inline h5pp::hid::h5t h5_type;

    struct table {
        double   J_nn             = 0;         /*!< Nearest neighbor coupling */
        double   J_nnn            = 0;         /*!< Next-nearest neighbor coupling */
        double   h_field          = 0;         /*!< On-site magnetic field */
        double   h_rnd            = 0;         /*!< Random field value */
        double   h_ptb            = 0;         /*!< Perturbation */
        double   h_mean           = 0;         /*!< Mean of the distribution for the random field */
        double   h_sigma          = 0;         /*!< Randomness strength. In distribution this is N(r_mean,r_sigma) or U(r_mean-r_sigma,r_mean+r_sigma) */
        uint64_t spin_dim         = 0;         /*!< Spin dimension */
        char     distribution[16] = "uniform"; /*!< The random distribution of J_rnd and h_rnd. Choose between lognormal, normal or uniform */
    };

    h5tb_tf_ising() { register_table_type(); }

    static void register_table_type() {
        if(h5_type.valid()) return;
        // Create a type for the char array from the template H5T_C_S1
        // The template describes a string with a single char.
        // Set the size with H5Tset_size, or h5pp::hdf5::setStringSize(...)
        h5pp::hid::h5t h5t_custom_string = H5Tcopy(H5T_C_S1);
        H5Tset_size(h5t_custom_string, 16);
        // Optionally set the null terminator '\0'
        H5Tset_strpad(h5t_custom_string, H5T_STR_NULLTERM);
        h5_type = H5Tcreate(H5T_COMPOUND, sizeof(table));
        H5Tinsert(h5_type, "J_nn", HOFFSET(table, J_nn), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J_nnn", HOFFSET(table, J_nnn), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_field", HOFFSET(table, h_field), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_rnd", HOFFSET(table, h_rnd), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_ptb", HOFFSET(table, h_ptb), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_mean", HOFFSET(table, h_mean), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_sigma", HOFFSET(table, h_sigma), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "spin_dim", HOFFSET(table, spin_dim), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "distribution", HOFFSET(table, distribution), h5t_custom_string);
    }
};
