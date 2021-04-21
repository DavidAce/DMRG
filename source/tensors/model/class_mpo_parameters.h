//
// Created by david on 2018-10-01.
//

#pragma once

#include <any>
#include <array>
#include <h5pp/details/h5ppHid.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdexcept>
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
        tools::log->info("{:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8}", name[0], name[1], name[2], name[3], name[4],
                         name[5], name[6], name[7], name[8], name[9], name[10], name[11], name[12], name[13]);
    }
    void print_parameter_values() const {
        tools::log->info("{:<+8.4f} {:<+8.4f} {:<+8.4f} {:<+8.4f} {:<+8.4f} {:<+8.4f} {:<+8.4f} {:<+8.4f} {:<+8.4f} {:<+8.4f} {:<+8.4f} {:<+8.4f} {:<8} {:<8}",
                         param.J_mean, param.J_stdv, param.J_rand, param.J_avrg, param.J_pert, param.h_mean, param.h_stdv, param.h_rand, param.h_avrg,
                         param.h_pert, param.lambda, param.delta, param.spin_dim, param.distribution);
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
        H5Tinsert(h5_type, "J1_rand", HOFFSET(table, J1), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J2_rand", HOFFSET(table, J2), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_tran", HOFFSET(table, h_tran), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_mean", HOFFSET(table, h_mean), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_stdv", HOFFSET(table, h_stdv), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_rand", HOFFSET(table, h_rand), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_pert", HOFFSET(table, h_pert), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "spin_dim", HOFFSET(table, spin_dim), H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "distribution", HOFFSET(table, distribution), h5t_custom_string);
    }

    static std::vector<std::string> get_parameter_names() {
        return {"J1_rand", "J2_rand", "h_tran", "h_mean", "h_stdv", "h_rand", "h_pert", "spin_dim", "distribution"};
    }
    static void print_parameter_names() {
        auto name = get_parameter_names();
        tools::log->info("{:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8}", name[0], name[1], name[2], name[3], name[4], name[5], name[6], name[7],
                         name[8]);
    }
    void print_parameter_values() const {
        tools::log->info("{:<+8.4f} {:<+8.4f} {:<+8.4f} {:<+8.4f} {:<+8.4f} {:<+8.4f} {:<+8.4f} {:<+12} {:<8}", param.J1, param.J2, param.h_tran, param.h_mean,
                         param.h_stdv, param.h_rand, param.h_pert, param.spin_dim, param.distribution);
    }
};

class h5tb_lbit {
    private:
    static constexpr size_t                 J2_size = 9;
    static constexpr std::array<hsize_t, 1> J2_dims = {J2_size};

    public:
    using J2Type = std::array<double, J2_size>;

    struct table {
        double   J1_rand          = 0;         /*!< On-site interaction */
        J2Type   J2_rand          = {};        /*!< Two-body interaction */
        double   J3_rand          = 0;         /*!< Three-body interaction */
        double   J1_mean          = 0;         /*!< Constant offset for on-site */
        double   J2_mean          = 0;         /*!< Constant offset for two-body interaction */
        double   J3_mean          = 0;         /*!< Constant offset for three-body interaction */
        double   J1_wdth          = 0;         /*!< Width of the uniform box distribution U(-w1,w1) */
        double   J2_wdth          = 0;         /*!< Width of the uniform box distribution U(-J2_wdth,J2_wdth) */
        double   J3_wdth          = 0;         /*!< Width of the uniform box distribution U(-J3_wdth,J3_wdth) */
        double   J2_base          = 0;         /*!< Base for power-decay of two-body interactions: J2_rand*J2_base^-|i-j| */
        size_t   J2_span          = 0;         /*!< Maximum allowed range for pairwise interactions, |i-j| <= J2_span. Note that J2_span + 1 MPOs are used */
        double   J1_pert          = 0;         /*!< On-site perturbation */
        double   J2_pert          = 0;         /*!< Two-body perturbation */
        double   J3_pert          = 0;         /*!< Three-body perturbation */
        double   f_mixer          = 0;         /*!< Mixing factor for unitary transformation to real-space */
        uint64_t u_layer          = 0;         /*!< Number of unitary 2-site layers which transform lbit <-> real spaces */
        long     spin_dim         = 2;         /*!< Spin dimension */
        char     distribution[16] = "uniform"; /*!< The random distribution of J_rnd and h_rnd. Choose between lognormal, normal or uniform */
    };

    static inline h5pp::hid::h5t h5_type;
    table                        param;

    h5tb_lbit() { register_table_type(); }

    static void register_table_type() {
        if(h5_type.valid()) return;
        h5pp::hid::h5t H5T_JARRAY_DOUBLE = H5Tarray_create(H5T_NATIVE_DOUBLE, J2_dims.size(), J2_dims.data());

        // Create a type for the char array from the template H5T_C_S1
        // The template describes a string with a single char.
        // Set the size with H5Tset_size, or h5pp::hdf5::setStringSize(...)
        h5pp::hid::h5t h5t_custom_string = H5Tcopy(H5T_C_S1);
        H5Tset_size(h5t_custom_string, 16);
        // Optionally set the null terminator '\0'
        H5Tset_strpad(h5t_custom_string, H5T_STR_NULLTERM);
        h5_type = H5Tcreate(H5T_COMPOUND, sizeof(table));
        H5Tinsert(h5_type, "J1_rand", HOFFSET(table, J1_rand), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J2_rand", HOFFSET(table, J2_rand), H5T_JARRAY_DOUBLE);
        H5Tinsert(h5_type, "J3_rand", HOFFSET(table, J3_rand), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J1_mean", HOFFSET(table, J1_mean), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J2_mean", HOFFSET(table, J2_mean), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J3_mean", HOFFSET(table, J3_mean), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J1_wdth", HOFFSET(table, J1_wdth), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J2_wdth", HOFFSET(table, J2_wdth), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J3_wdth", HOFFSET(table, J3_wdth), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J2_base", HOFFSET(table, J2_base), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J2_span", HOFFSET(table, J2_span), H5T_NATIVE_ULONG);
        H5Tinsert(h5_type, "J1_pert", HOFFSET(table, J1_pert), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J2_pert", HOFFSET(table, J2_pert), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J3_pert", HOFFSET(table, J3_pert), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "f_mixer", HOFFSET(table, f_mixer), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "u_layer", HOFFSET(table, u_layer), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "spin_dim", HOFFSET(table, spin_dim), H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "distribution", HOFFSET(table, distribution), h5t_custom_string);
    }

    [[nodiscard]] std::string J2_str() const { return fmt::format("[{:<+8.4f}]", fmt::join(param.J2_rand, " ")); }

    static std::vector<std::string> get_parameter_names() {
        return {"J1_rand", "J2_rand", "J3_rand", "J1_mean", "J2_mean", "J3_mean", "J1_wdth", "J2_wdth",  "J3_wdth",
                "J2_base", "J2_span", "J1_pert", "J2_pert", "J3_pert", "f_mixer", "u_layer", "spin_dim", "distribution"};
    }

    static void print_parameter_names() {
        auto name = get_parameter_names();
        tools::log->info("{:<8} {:<82} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8} {:<8}", name[0], name[1],
                         name[2], name[3], name[4], name[5], name[6], name[7], name[8], name[9], name[10], name[11], name[12], name[13], name[14], name[15],
                         name[16], name[17]);
    }

    void print_parameter_values() const {
        tools::log->info("{:<+8.4f} {:<} {:<+8.4f} {:<+8.4f} {:<+8.4f} {:<+8.4f} {:<+8.4f} {:<+8.4f} {:<+8.4f} {:<+8.4f} {:<8} "
                         "{:<+8.4f} {:<+8.4f} {:<+8.4f} {:<+8.4f} {:<8} {:<8} {:<8}",
                         param.J1_rand, J2_str(), param.J3_rand, param.J1_mean, param.J2_mean, param.J3_mean, param.J1_wdth, param.J2_wdth, param.J3_wdth,
                         param.J2_base, param.J2_span, param.J1_pert, param.J2_pert, param.J3_pert, param.f_mixer, param.u_layer, param.spin_dim,
                         param.distribution);
    }
};
