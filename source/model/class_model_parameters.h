//
// Created by david on 2018-10-01.
//

#pragma once

#include <any>
#include <array>
#include <h5pp/details/h5ppHid.h>
#include <hdf5.h>
#include <hdf5_hl.h>

class h5tb_ising_selfdual_tf_rf_nn {
    public:
    static inline h5pp::hid::h5t h5_type;

    struct table {
        double   J_rand           = 0;         /*!< Randomly distributed nearest neighbour coupling */
        double   h_rand           = 0;         /*!< Randomly distributed on-site field */
        double   J_pert           = 0;         /*!< Perturbation to the coupling, std::pow(J_rnd + J_ptb,1-alpha) */
        double   h_pert           = 0;         /*!< Perturbation to the coupling, std::pow(J_rnd + J_ptb,1-alpha) */
        double   J_avrg           = 0;         /*!< Average of J_rnd between all sites*/
        double   h_avrg           = 0;         /*!< Average of h_rnd on all sites */
        double   J_mean           = 0;         /*!< Mean for the distrbution of J_rnd */
        double   h_mean           = 0;         /*!< Mean for the distrbution of h_rnd */
        double   J_stdv           = 0;         /*!< Standard deviation for the distribution of J_rnd */
        double   h_stdv           = 0;         /*!< Standard deviation for the distribution of h_rnd */
        double   lambda           = 0;         /*!< Factor involved in next-nearest neighbor interaction */
        double   delta            = 0;         /*!< Difference J_mean - h_mean  */
        uint64_t spin_dim         = 2;         /*!< Spin dimension */
        char     distribution[16] = "uniform"; /*!< The random distribution of J_rnd and h_rnd. Choose between lognormal, normal or uniform */
    };

    h5tb_ising_selfdual_tf_rf_nn() { register_table_type(); }
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
        H5Tinsert(h5_type, "J_rand", HOFFSET(table, J_rand), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_rand", HOFFSET(table, h_rand), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J_pert", HOFFSET(table, J_pert), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_pert", HOFFSET(table, h_pert), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J_avrg", HOFFSET(table, J_avrg), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_avrg", HOFFSET(table, h_avrg), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J_mean", HOFFSET(table, J_mean), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_mean", HOFFSET(table, h_mean), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J_stdv", HOFFSET(table, J_stdv), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_stdv", HOFFSET(table, h_stdv), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "lambda", HOFFSET(table, lambda), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "delta", HOFFSET(table, delta), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "spin_dim", HOFFSET(table, spin_dim), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "distribution", HOFFSET(table, distribution), h5t_custom_string);
    }
};

class h5tb_ising_tf_rf_nn {
    public:
    static inline h5pp::hid::h5t h5_type;

    struct table {
        double   J1               = 0;         /*!< Nearest neighbor coupling */
        double   J2               = 0;         /*!< Next-nearest neighbor coupling */
        double   h_tran           = 0;         /*!< Transverse field strength */
        double   h_mean           = 0;         /*!< Random field mean of distribution */
        double   h_stdv           = 0;         /*!< Random field standard deviation. In distribution this is N(h_mean,h_stdv) or U(h_mean-h_stdv,h_mean+h_stdv) */
        double   h_rand           = 0;         /*!< Random field value */
        double   h_pert           = 0;         /*!< Perturbation */
        uint64_t spin_dim         = 2;         /*!< Spin dimension */
        char     distribution[16] = "uniform"; /*!< The random distribution of J_rnd and h_rnd. Choose between lognormal, normal or uniform */
    };

    h5tb_ising_tf_rf_nn() { register_table_type(); }

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
        H5Tinsert(h5_type, "J1", HOFFSET(table, J1), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J2", HOFFSET(table, J2), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_tran", HOFFSET(table, h_tran), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_mean", HOFFSET(table, h_mean), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_stdv", HOFFSET(table, h_stdv), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_rand", HOFFSET(table, h_rand), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_pert", HOFFSET(table, h_pert), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "spin_dim", HOFFSET(table, spin_dim), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "distribution", HOFFSET(table, distribution), h5t_custom_string);
    }
};
