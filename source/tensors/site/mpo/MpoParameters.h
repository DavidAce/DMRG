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
        double J_pert           = 0;         /*!< Perturbation to the coupling, */
        double h_mean           = 0;         /*!< Mean for the distrbution of h_rnd */
        double h_stdv           = 0;         /*!< Standard deviation for the distribution of h_rnd */
        double h_rand           = 0;         /*!< Randomly distributed on-site field */
        double h_avrg           = 0;         /*!< Average of h_rnd on all sites */
        double h_pert           = 0;         /*!< Perturbation to the field */
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

    [[nodiscard]] std::string fmt_value(std::string_view p) const {
        /* clang-format off */
        if(p == "J_mean")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.J_mean);
        if(p == "J_stdv")           return fmt::format(FMT_STRING("{:<9.2e}") , param.J_stdv);
        if(p == "J_rand")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.J_rand);
        if(p == "J_avrg")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.J_avrg);
        if(p == "J_pert")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.J_pert);
        if(p == "h_mean")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.h_mean);
        if(p == "h_stdv")           return fmt::format(FMT_STRING("{:<9.2e}") , param.h_stdv);
        if(p == "h_rand")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.h_rand);
        if(p == "h_avrg")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.h_avrg);
        if(p == "h_pert")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.h_pert);
        if(p == "lambda")           return fmt::format(FMT_STRING("{:<7.4f}") , param.lambda);
        if(p == "delta")            return fmt::format(FMT_STRING("{:<+7.4f}"), param.delta);
        if(p == "spin_dim")         return fmt::format(FMT_STRING("{:>8}")    , param.spin_dim);
        if(p == "distribution")     return fmt::format(FMT_STRING("{:<12}")   , param.distribution);
        /* clang-format on */
        throw std::runtime_error(fmt::format("Unrecognized parameter: {}", p));
    }

    static std::vector<std::string> get_parameter_names() {
        return {"J_mean", "J_stdv", "J_rand", "J_avrg", "J_pert", "h_mean",   "h_stdv",
                "h_rand", "h_avrg", "h_pert", "lambda", "delta",  "spin_dim", "distribution"};
    }

    void print_parameter_names() const {
        std::string name_line;
        for(const auto &name : get_parameter_names()) { name_line.append(fmt::format(FMT_STRING("{:<{}} "), name, fmt_value(name).size())); }
        tools::log->info(name_line);
    }

    void print_parameter_values() const {
        std::string value_line;
        for(const auto &name : get_parameter_names()) { value_line.append(fmt::format(FMT_STRING("{} "), fmt_value(name))); }
        tools::log->info(value_line);
    }
};

class h5tb_ising_majorana {
    public:
    struct table {
        double J_mean           = 0; /*!< Mean for the distrbution of J_rnd */
        double J_boxw           = 0; /*!< Standard deviation for the distribution of J_rnd */
        double J_rand           = 0; /*!< Randomly distributed nearest neighbour coupling */
        double J_avrg           = 0; /*!< Average of J_rnd between all sites*/
        double J_pert           = 0; /*!< Perturbation to the coupling, */
        double h_mean           = 0; /*!< Mean for the distrbution of h_rnd */
        double h_boxw           = 0; /*!< Standard deviation for the distribution of h_rnd */
        double h_rand           = 0; /*!< Randomly distributed on-site field */
        double h_avrg           = 0; /*!< Average of h_rnd on all sites */
        double h_pert           = 0; /*!< Perturbation to the field */
        double g                = 0; /*!< Interaction parameter for nearest ZZ and next-nearest XX neighbor coupling */
        double delta            = 0; /*!< Delta defined as log(J_mean) - log(h_mean). We get J_mean and h_mean by fixing delta = 2lnW, W = J_boxw = 1/h_boxw */
        long   spin_dim         = 2; /*!< Spin dimension */
        char   distribution[16] = "uniform"; /*!< The random distribution of J_rnd and h_rnd. Choose between lognormal, normal or uniform */
    };
    static inline h5pp::hid::h5t h5_type;
    table                        param;

    h5tb_ising_majorana() { register_table_type(); }
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
        H5Tinsert(h5_type, "J_boxw", HOFFSET(table, J_boxw), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J_rand", HOFFSET(table, J_rand), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J_avrg", HOFFSET(table, J_avrg), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J_pert", HOFFSET(table, J_pert), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_mean", HOFFSET(table, h_mean), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_boxw", HOFFSET(table, h_boxw), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_rand", HOFFSET(table, h_rand), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_avrg", HOFFSET(table, h_avrg), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_pert", HOFFSET(table, h_pert), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "g", HOFFSET(table, g), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "delta", HOFFSET(table, delta), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "spin_dim", HOFFSET(table, spin_dim), H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "distribution", HOFFSET(table, distribution), h5t_custom_string);
    }

    [[nodiscard]] std::string fmt_value(std::string_view p) const {
        /* clang-format off */
        if(p == "J_mean")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.J_mean);
        if(p == "J_boxw")           return fmt::format(FMT_STRING("{:<9.2e}") , param.J_boxw);
        if(p == "J_rand")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.J_rand);
        if(p == "J_avrg")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.J_avrg);
        if(p == "J_pert")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.J_pert);
        if(p == "h_mean")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.h_mean);
        if(p == "h_boxw")           return fmt::format(FMT_STRING("{:<9.2e}") , param.h_boxw);
        if(p == "h_rand")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.h_rand);
        if(p == "h_avrg")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.h_avrg);
        if(p == "h_pert")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.h_pert);
        if(p == "g")                return fmt::format(FMT_STRING("{:<7.4f}") , param.g);
        if(p == "delta")            return fmt::format(FMT_STRING("{:<+7.4f}"), param.delta);
        if(p == "spin_dim")         return fmt::format(FMT_STRING("{:>8}")    , param.spin_dim);
        if(p == "distribution")     return fmt::format(FMT_STRING("{:<12}")   , param.distribution);
        /* clang-format on */
        throw std::runtime_error(fmt::format("Unrecognized parameter: {}", p));
    }

    static std::vector<std::string> get_parameter_names() {
        return {"J_mean", "J_boxw", "J_rand", "J_avrg", "J_pert", "h_mean", "h_boxw", "h_rand", "h_avrg", "h_pert", "g", "delta", "spin_dim", "distribution"};
    }

    void print_parameter_names() const {
        std::string name_line;
        for(const auto &name : get_parameter_names()) { name_line.append(fmt::format(FMT_STRING("{:<{}} "), name, fmt_value(name).size())); }
        tools::log->info(name_line);
    }

    void print_parameter_values() const {
        std::string value_line;
        for(const auto &name : get_parameter_names()) { value_line.append(fmt::format(FMT_STRING("{} "), fmt_value(name))); }
        tools::log->info(value_line);
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

    [[nodiscard]] std::string fmt_value(std::string_view p) const {
        /* clang-format off */
        if(p == "J1")               return fmt::format(FMT_STRING("{:<+9.2e}"), param.J1);
        if(p == "J2")               return fmt::format(FMT_STRING("{:<+9.2e}"), param.J2);
        if(p == "h_tran")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.h_tran);
        if(p == "h_mean")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.h_mean);
        if(p == "h_stdv")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.h_stdv);
        if(p == "h_rand")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.h_rand);
        if(p == "h_pert")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.h_pert);
        if(p == "spin_dim")         return fmt::format(FMT_STRING("{:>8}")    , param.spin_dim);
        if(p == "distribution")     return fmt::format(FMT_STRING("{:<12}")   , param.distribution);
        /* clang-format on */
        throw std::runtime_error(fmt::format("Unrecognized parameter: {}", p));
    }

    static std::vector<std::string> get_parameter_names() { return {"J1", "J2", "h_tran", "h_mean", "h_stdv", "h_rand", "h_pert", "spin_dim", "distribution"}; }

    void print_parameter_names() const {
        std::string name_line;
        for(const auto &name : get_parameter_names()) { name_line.append(fmt::format(FMT_STRING("{:<{}} "), name, fmt_value(name).size())); }
        tools::log->info(name_line);
    }

    void print_parameter_values() const {
        std::string value_line;
        for(const auto &name : get_parameter_names()) { value_line.append(fmt::format(FMT_STRING("{} "), fmt_value(name))); }
        tools::log->info(value_line);
    }
};

class h5tb_lbit {
    private:
    static constexpr size_t                 J2_size = 64;
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
        double   J2_xcls          = 0;         /*!< Exp. decay rate of two-body interactions: exp(-|i-j|/J2_xcls) * J2_rand */
        size_t   J2_span          = 0;         /*!< Maximum range for pairwise interactions, |i-j| <= J2_span. */
        size_t   J2_ctof          = 0;         /*!< Effective range for pairwise interactions, |i-j| <= std::min(J2_span, model_size-1). */
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
        H5Tinsert(h5_type, "J2_xcls", HOFFSET(table, J2_xcls), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J2_span", HOFFSET(table, J2_span), H5T_NATIVE_ULONG);
        H5Tinsert(h5_type, "J2_ctof", HOFFSET(table, J2_ctof), H5T_NATIVE_ULONG);
        H5Tinsert(h5_type, "J1_pert", HOFFSET(table, J1_pert), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J2_pert", HOFFSET(table, J2_pert), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J3_pert", HOFFSET(table, J3_pert), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "f_mixer", HOFFSET(table, f_mixer), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "u_layer", HOFFSET(table, u_layer), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "spin_dim", HOFFSET(table, spin_dim), H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "distribution", HOFFSET(table, distribution), h5t_custom_string);
    }

    [[nodiscard]] std::string J2_str() const {
        return fmt::format(FMT_STRING("[{:<+9.2e}]"), fmt::join(param.J2_rand.begin(), param.J2_rand.begin() + param.J2_ctof + 1, ","));
    }
    [[nodiscard]] std::string fmt_value(std::string_view p) const {
        /* clang-format off */
        if(p == "J1_rand")     return fmt::format(FMT_STRING("{:<+9.2e}"),param.J1_rand);
        if(p == "J2_rand")     return fmt::format(FMT_STRING("[{:<+9.2e}]"), fmt::join(param.J2_rand.begin(), param.J2_rand.begin() + param.J2_ctof + 1, ","));;
        if(p == "J3_rand")     return fmt::format(FMT_STRING("{:<+9.2e}"), param.J3_rand);
        if(p == "J1_mean")     return fmt::format(FMT_STRING("{:<+9.2e}"), param.J1_mean);
        if(p == "J2_mean")     return fmt::format(FMT_STRING("{:<+9.2e}"), param.J2_mean);
        if(p == "J3_mean")     return fmt::format(FMT_STRING("{:<+9.2e}"), param.J3_mean);
        if(p == "J1_wdth")     return fmt::format(FMT_STRING("{:<7.4f}"),  param.J1_wdth);
        if(p == "J2_wdth")     return fmt::format(FMT_STRING("{:<7.4f}"),  param.J2_wdth);
        if(p == "J3_wdth")     return fmt::format(FMT_STRING("{:<7.4f}"),  param.J3_wdth);
        if(p == "J2_xcls")     return fmt::format(FMT_STRING("{:<7.4f}"),  param.J2_xcls);
        if(p == "J2_span")     return fmt::format(FMT_STRING("{:>7}"),     param.J2_span);
        if(p == "J2_ctof")     return fmt::format(FMT_STRING("{:>7}"),     param.J2_ctof);
        if(p == "J1_pert")     return fmt::format(FMT_STRING("{:<+7.3f}"), param.J1_pert);
        if(p == "J2_pert")     return fmt::format(FMT_STRING("{:<+7.3f}"), param.J2_pert);
        if(p == "J3_pert")     return fmt::format(FMT_STRING("{:<+7.3f}"), param.J3_pert);
        if(p == "f_mixer")     return fmt::format(FMT_STRING("{:<7.4f}"),  param.f_mixer);
        if(p == "u_layer")     return fmt::format(FMT_STRING("{:>7}"),     param.u_layer);
        if(p == "spin_dim")    return fmt::format(FMT_STRING("{:>8}"),     param.spin_dim);
        if(p == "distribution")return fmt::format(FMT_STRING("{:<12}")    ,param.distribution);
        /* clang-format on */
        throw std::runtime_error(fmt::format("Unrecognized parameter: {}", p));
    }

    static std::vector<std::string> get_parameter_names() {
        return {"J1_rand", "J2_rand", "J3_rand", "J1_mean", "J2_mean", "J3_mean", "J1_wdth", "J2_wdth",  "J3_wdth",     "J2_xcls",
                "J2_span", "J2_ctof", "J1_pert", "J2_pert", "J3_pert", "f_mixer", "u_layer", "spin_dim", "distribution"};
    }

    void print_parameter_names() const {
        std::string name_line;
        for(const auto &name : get_parameter_names()) { name_line.append(fmt::format(FMT_STRING("{:<{}} "), name, fmt_value(name).size())); }
        tools::log->info(name_line);
    }

    void print_parameter_values() const {
        std::string value_line;
        for(const auto &name : get_parameter_names()) { value_line.append(fmt::format(FMT_STRING("{} "), fmt_value(name))); }
        tools::log->info(value_line);
    }
};
