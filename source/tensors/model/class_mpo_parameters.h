//
// Created by david on 2018-10-01.
//

#pragma once

#include <any>
#include <array>
#include <stdexcept>

#include <h5pp/details/h5ppHid.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <tools/common/log.h>

template<typename SrcType, typename TgtType, size_t size>
void copy_c_str(const SrcType &src, TgtType (&tgt)[size])
// Use to copy the distribution char array from string
{
    tgt[src.copy(tgt, size - 1)] = 0; // Null terminates
}

class h5tb_ising_sdual {
    public:
    struct table {
        double J_mean           = 0;         /*!< Mean for the distrbution of J_rnd */
        double J_stdv           = 0;         /*!< Standard deviation for the distribution of J_rnd */
        double J_rand           = 0;         /*!< Randomly distributed nearest neighbour coupling */
        double J_avrg           = 0;         /*!< Average of J_rnd between all sites*/
        double J_pert           = 0;         /*!< Perturbation to the coupling, std::pow(J_rnd + J_ptb,1-alpha) */
        double h_mean           = 0;         /*!< Mean for the distrbution of h_rnd */
        double h_stdv           = 0;         /*!< Standard deviation for the distribution of h_rnd */
        double h_rand           = 0;         /*!< Randomly distributed on-site field */
        double h_avrg           = 0;         /*!< Average of h_rnd on all sites */
        double h_pert           = 0;         /*!< Perturbation to the coupling, std::pow(J_rnd + J_ptb,1-alpha) */
        double lambda           = 0;         /*!< Factor involved in next-nearest neighbor interaction */
        double delta            = 0;         /*!< Difference J_mean - h_mean  */
        long   spin_dim         = 2;         /*!< Spin dimension */
        char   distribution[16] = "uniform"; /*!< The random distribution of J_rnd and h_rnd. Choose between lognormal, normal or uniform */
    };
    static inline h5pp::hid::h5t h5_type;
    table                        param;

    h5tb_ising_sdual() { register_table_type(); }
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
        H5Tinsert(h5_type, "J_mean", HOFFSET(table, J_mean), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J_stdv", HOFFSET(table, J_stdv), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J_rand", HOFFSET(table, J_rand), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J_avrg", HOFFSET(table, J_avrg), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J_pert", HOFFSET(table, J_pert), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_mean", HOFFSET(table, h_mean), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_stdv", HOFFSET(table, h_stdv), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_rand", HOFFSET(table, h_rand), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_avrg", HOFFSET(table, h_avrg), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_pert", HOFFSET(table, h_pert), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "lambda", HOFFSET(table, lambda), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "delta", HOFFSET(table, delta), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "spin_dim", HOFFSET(table, spin_dim), H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "distribution", HOFFSET(table, distribution), h5t_custom_string);
    }

    static std::vector<std::string> get_parameter_names() {
        return {"J_mean", "J_stdv", "J_rand", "J_avrg", "J_pert", "h_mean",   "h_stdv",
                "h_rand", "h_avrg", "h_pert", "lambda", "delta",  "spin_dim", "distribution"};
    }
    static void print_parameter_names() {
        auto name = get_parameter_names();
        tools::log->info("{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12}", name[0], name[1], name[2],
                         name[3], name[4], name[5], name[6], name[7], name[8], name[9], name[10], name[11], name[12], name[13]);
    }
    void print_parameter_values() const {
        tools::log->info(
            "{:<+12.8f} {:<+12.8f} {:<+12.8f} {:<+12.8f} {:<+12.8f} {:<+12.8f} {:<+12.8f} {:<+12.8f} {:<+12.8f} {:<+12.8f} {:<+12.8f} {:<+12} {:<12}",
            param.J_mean, param.J_stdv, param.J_rand, param.J_avrg, param.J_pert, param.h_mean, param.h_stdv, param.h_rand, param.h_avrg, param.h_pert,
            param.lambda, param.delta, param.spin_dim, param.distribution);
    }
};

class h5tb_ising_tf_rf {
    public:
    struct table {
        double J1               = 0;         /*!< Nearest neighbor coupling */
        double J2               = 0;         /*!< Next-nearest neighbor coupling */
        double h_tran           = 0;         /*!< Transverse field strength */
        double h_mean           = 0;         /*!< Random field mean of distribution */
        double h_stdv           = 0;         /*!< Random field standard deviation. In distribution this is N(h_mean,h_stdv) or U(h_mean-h_stdv,h_mean+h_stdv) */
        double h_rand           = 0;         /*!< Random field value */
        double h_pert           = 0;         /*!< Perturbation */
        long   spin_dim         = 2;         /*!< Spin dimension */
        char   distribution[16] = "uniform"; /*!< The random distribution of J_rnd and h_rnd. Choose between lognormal, normal or uniform */
    };

    static inline h5pp::hid::h5t h5_type;
    table                        param;

    h5tb_ising_tf_rf() { register_table_type(); }

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
        H5Tinsert(h5_type, "spin_dim", HOFFSET(table, spin_dim), H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "distribution", HOFFSET(table, distribution), h5t_custom_string);
    }

    static std::vector<std::string> get_parameter_names() { return {"J1", "J2", "h_tran", "h_mean", "h_stdv", "h_rand", "h_pert", "spin_dim", "distribution"}; }
    static void                     print_parameter_names() {
        auto name = get_parameter_names();
        tools::log->info("{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12}", name[0], name[1], name[2], name[3], name[4], name[5], name[6],
                         name[7], name[8]);
    }
    void print_parameter_values() const {
        tools::log->info("{:<+12.8f} {:<+12.8f} {:<+12.8f} {:<+12.8f} {:<+12.8f} {:<+12.8f} {:<+12.8f} {:<+12} {:<12}", param.J1, param.J2, param.h_tran,
                         param.h_mean, param.h_stdv, param.h_rand, param.h_pert, param.spin_dim, param.distribution);
    }
};

class h5tb_lbit {
    public:
    struct table {
        double J1               = 0;         /*!< On-site interaction */
        double J2               = 0;         /*!< Two-body interaction */
        double J3               = 0;         /*!< Three-body interaction */
        double w1               = 0;         /*!< Width of the uniform box distribution U(-w1,w1) */
        double w2               = 0;         /*!< Width of the uniform box distribution U(-w2,w2) */
        double w3               = 0;         /*!< Width of the uniform box distribution U(-w3,w3) */
        double J1_pert          = 0;         /*!< On-site perturbation */
        double J2_pert          = 0;         /*!< Two-body perturbation */
        double J3_pert          = 0;         /*!< Three-body perturbation */
        long   spin_dim         = 2;         /*!< Spin dimension */
        char   distribution[16] = "uniform"; /*!< The random distribution of J_rnd and h_rnd. Choose between lognormal, normal or uniform */
    };

    static inline h5pp::hid::h5t h5_type;
    table                        param;

    h5tb_lbit() { register_table_type(); }

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
        H5Tinsert(h5_type, "J3", HOFFSET(table, J3), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "w1", HOFFSET(table, w1), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "w2", HOFFSET(table, w2), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "w3", HOFFSET(table, w3), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J1_pert", HOFFSET(table, J1_pert), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J2_pert", HOFFSET(table, J2_pert), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J3_pert", HOFFSET(table, J3_pert), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "spin_dim", HOFFSET(table, spin_dim), H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "distribution", HOFFSET(table, distribution), h5t_custom_string);
    }

    static std::vector<std::string> get_parameter_names() {
        return {"J1", "J2", "J3", "w1", "w2", "w3", "J1_pert", "J2_pert", "J3_pert", "spin_dim", "distribution"};
    }
    static void print_parameter_names() {
        auto name = get_parameter_names();
        tools::log->info("{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12}", name[0], name[1], name[2], name[3], name[4], name[5],
                         name[6], name[7], name[8], name[9], name[10]);
    }
    void print_parameter_values() const {
        tools::log->info("{:<+12.8f} {:<+12.8f} {:<+12.8f} {:<+12.8f} {:<+12.8f} {:<+12.8f} {:<+12.8f} {:<+12.8f} {:<+12.8f} {:<+12} {:<12}", param.J1,
                         param.J2, param.J3, param.w1, param.w2, param.w3, param.J1_pert, param.J2_pert, param.J3_pert, param.spin_dim, param.distribution);
    }
};
