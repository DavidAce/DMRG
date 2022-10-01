#pragma once

#include "debug/exceptions.h"
#include "tools/common/log.h"
#include <any>
#include <array>
#include <h5pp/details/h5ppHid.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdexcept>

class h5tb_ising_selfdual {
    public:
    struct table {
        double J_mean       = 0;       /*!< Mean for the distrbution of J_rand */
        double J_wdth       = 0;       /*!< Width for the distrbution of J_rand */
        double J_rand       = 0;       /*!< Randomly distributed nearest neighbour coupling */
        double h_mean       = 0;       /*!< Mean for the distrbution of h_rand */
        double h_wdth       = 0;       /*!< Width for the distrbution of h_rand */
        double h_rand       = 0;       /*!< Randomly distributed on-site field */
        double lambda       = 0;       /*!< Factor involved in next-nearest neighbor interaction */
        double delta        = 0;       /*!< Difference log(J_mean) - log(h_mean) */
        long   spin_dim     = 2;       /*!< Spin dimension */
        char  *distribution = nullptr; /*!< The random distribution of J_rnd and h_rnd. Choose between lognormal, normal or uniform */
        ~table() noexcept {
            if(distribution != nullptr) { free(distribution); }
        }
        std::string_view get_distribution() const {
            if(distribution == nullptr) throw except::runtime_error("distribution == nullptr");
            return distribution;
        }
        void set_distribution(std::string_view dist) {
            if(distribution != nullptr) free(distribution);
            distribution = static_cast<char *>(malloc((dist.size() + 1) * sizeof(char))); // Add +1 for null terminator
            strcpy(distribution, dist.data());
        }
    };
    static inline h5pp::hid::h5t h5_type;
    table                        param;

    h5tb_ising_selfdual() { register_table_type(); }
    static void register_table_type() {
        if(h5_type.valid()) return;

        // Create a type for the char array from the template H5T_C_S1
        h5pp::hid::h5t h5t_custom_string = H5Tcopy(H5T_C_S1);
        H5Tset_size(h5t_custom_string, H5T_VARIABLE);
        H5Tset_strpad(h5t_custom_string, H5T_STR_NULLTERM);

        h5_type = H5Tcreate(H5T_COMPOUND, sizeof(table));
        H5Tinsert(h5_type, "J_mean", HOFFSET(table, J_mean), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J_wdth", HOFFSET(table, J_wdth), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J_rand", HOFFSET(table, J_rand), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_mean", HOFFSET(table, h_mean), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_wdth", HOFFSET(table, h_wdth), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_rand", HOFFSET(table, h_rand), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "lambda", HOFFSET(table, lambda), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "delta", HOFFSET(table, delta), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "spin_dim", HOFFSET(table, spin_dim), H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "distribution", HOFFSET(table, distribution), h5t_custom_string);
    }

    [[nodiscard]] std::string fmt_value(std::string_view p) const {
        /* clang-format off */
        if(p == "J_mean")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.J_mean);
        if(p == "J_wdth")           return fmt::format(FMT_STRING("{:<9.2e}") , param.J_wdth);
        if(p == "J_rand")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.J_rand);
        if(p == "h_mean")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.h_mean);
        if(p == "h_wdth")           return fmt::format(FMT_STRING("{:<9.2e}") , param.h_wdth);
        if(p == "h_rand")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.h_rand);
        if(p == "lambda")           return fmt::format(FMT_STRING("{:<7.4f}") , param.lambda);
        if(p == "delta")            return fmt::format(FMT_STRING("{:<+7.4f}"), param.delta);
        if(p == "spin_dim")         return fmt::format(FMT_STRING("{:>8}")    , param.spin_dim);
        if(p == "distribution")     return fmt::format(FMT_STRING("{:<12}")   , param.get_distribution());
        /* clang-format on */
        throw except::runtime_error("Unrecognized parameter: {}", p);
    }

    static std::vector<std::string> get_parameter_names() {
        return {"J_mean", "J_wdth", "J_rand", "h_mean", "h_wdth", "h_rand", "lambda", "delta", "spin_dim", "distribution"};
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
        double J_mean       = 0; /*!< Mean for the distrbution of J_rand */
        double J_wdth       = 0; /*!< Width for the distrbution of J_rand */
        double J_rand       = 0; /*!< Randomly distributed nearest neighbour coupling */
        double h_mean       = 0; /*!< Mean for the distrbution of h_rand */
        double h_wdth       = 0; /*!< Width for the distrbution of h_rand */
        double h_rand       = 0; /*!< Randomly distributed on-site field */
        double g            = 0; /*!< Interaction parameter for nearest ZZ and next-nearest XX neighbor coupling */
        double delta        = 0; /*!< Delta defined as log(J_mean) - log(h_mean). We get J_mean and h_mean by fixing delta = 2lnW, W = J_wdth = 1/h_wdth */
        long   spin_dim     = 2; /*!< Spin dimension */
        char  *distribution = nullptr; /*!< The random distribution of J_rnd and h_rnd. Choose between lognormal, normal or uniform */
        ~table() noexcept {
            if(distribution != nullptr) { free(distribution); }
        }
        std::string_view get_distribution() const {
            if(distribution == nullptr) throw except::runtime_error("distribution == nullptr");
            return distribution;
        }
        void set_distribution(std::string_view dist) {
            if(distribution != nullptr) free(distribution);
            distribution = static_cast<char *>(malloc((dist.size() + 1) * sizeof(char))); // Add +1 for null terminator
            strcpy(distribution, dist.data());
        }
    };
    static inline h5pp::hid::h5t h5_type;
    table                        param;

    h5tb_ising_majorana() { register_table_type(); }
    static void register_table_type() {
        if(h5_type.valid()) return;

        // Create a type for the char array from the template H5T_C_S1
        h5pp::hid::h5t h5t_custom_string = H5Tcopy(H5T_C_S1);
        H5Tset_size(h5t_custom_string, H5T_VARIABLE);
        H5Tset_strpad(h5t_custom_string, H5T_STR_NULLTERM);

        h5_type = H5Tcreate(H5T_COMPOUND, sizeof(table));
        H5Tinsert(h5_type, "J_mean", HOFFSET(table, J_mean), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J_wdth", HOFFSET(table, J_wdth), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J_rand", HOFFSET(table, J_rand), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_mean", HOFFSET(table, h_mean), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_wdth", HOFFSET(table, h_wdth), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_rand", HOFFSET(table, h_rand), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "g", HOFFSET(table, g), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "delta", HOFFSET(table, delta), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "spin_dim", HOFFSET(table, spin_dim), H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "distribution", HOFFSET(table, distribution), h5t_custom_string);
    }

    [[nodiscard]] std::string fmt_value(std::string_view p) const {
        /* clang-format off */
        if(p == "J_mean")           return fmt::format(FMT_STRING("{:<8.2e}") , param.J_mean);
        if(p == "J_wdth")           return fmt::format(FMT_STRING("{:<8.2e}") , param.J_wdth);
        if(p == "J_rand")           return fmt::format(FMT_STRING("{:<8.2e}") , param.J_rand);
        if(p == "h_mean")           return fmt::format(FMT_STRING("{:<8.2e}") , param.h_mean);
        if(p == "h_wdth")           return fmt::format(FMT_STRING("{:<8.2e}") , param.h_wdth);
        if(p == "h_rand")           return fmt::format(FMT_STRING("{:<8.2e}") , param.h_rand);
        if(p == "g")                return fmt::format(FMT_STRING("{:<7.4f}") , param.g);
        if(p == "delta")            return fmt::format(FMT_STRING("{:<+7.4f}"), param.delta);
        if(p == "spin_dim")         return fmt::format(FMT_STRING("{:>8}")    , param.spin_dim);
        if(p == "distribution")     return fmt::format(FMT_STRING("{:<12}")   , param.get_distribution());
        /* clang-format on */
        throw except::runtime_error("Unrecognized parameter: {}", p);
    }

    static std::vector<std::string> get_parameter_names() {
        return {"J_mean", "J_wdth", "J_rand", "h_mean", "h_wdth", "h_rand", "g", "delta", "spin_dim", "distribution"};
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
        double J1           = 0;       /*!< Nearest neighbor coupling */
        double J2           = 0;       /*!< Next-nearest neighbor coupling */
        double h_tran       = 0;       /*!< Transverse field strength */
        double h_mean       = 0;       /*!< Random field mean of distribution */
        double h_wdth       = 0;       /*!< Random field width of distribution. */
        double h_rand       = 0;       /*!< Random field value */
        long   spin_dim     = 2;       /*!< Spin dimension */
        char  *distribution = nullptr; /*!< The random distribution of J_rnd and h_rnd. Choose between lognormal, normal or uniform */
        ~table() noexcept {
            if(distribution != nullptr) { free(distribution); }
        }
        std::string_view get_distribution() const {
            if(distribution == nullptr) throw except::runtime_error("distribution == nullptr");
            return distribution;
        }
        void set_distribution(std::string_view dist) {
            if(distribution != nullptr) free(distribution);
            distribution = static_cast<char *>(malloc((dist.size() + 1) * sizeof(char))); // Add +1 for null terminator
            strcpy(distribution, dist.data());
        }
    };

    static inline h5pp::hid::h5t h5_type;
    table                        param;

    h5tb_ising_tf_rf() { register_table_type(); }

    static void register_table_type() {
        if(h5_type.valid()) return;
        // Create a type for the char array from the template H5T_C_S1
        h5pp::hid::h5t h5t_custom_string = H5Tcopy(H5T_C_S1);
        H5Tset_size(h5t_custom_string, H5T_VARIABLE);
        H5Tset_strpad(h5t_custom_string, H5T_STR_NULLTERM);

        h5_type = H5Tcreate(H5T_COMPOUND, sizeof(table));
        H5Tinsert(h5_type, "J1", HOFFSET(table, J1), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J2", HOFFSET(table, J2), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_tran", HOFFSET(table, h_tran), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_mean", HOFFSET(table, h_mean), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_wdth", HOFFSET(table, h_wdth), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "h_rand", HOFFSET(table, h_rand), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "spin_dim", HOFFSET(table, spin_dim), H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "distribution", HOFFSET(table, distribution), h5t_custom_string);
    }

    [[nodiscard]] std::string fmt_value(std::string_view p) const {
        /* clang-format off */
        if(p == "J1")               return fmt::format(FMT_STRING("{:<+9.2e}"), param.J1);
        if(p == "J2")               return fmt::format(FMT_STRING("{:<+9.2e}"), param.J2);
        if(p == "h_tran")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.h_tran);
        if(p == "h_mean")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.h_mean);
        if(p == "h_wdth")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.h_wdth);
        if(p == "h_rand")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.h_rand);
        if(p == "spin_dim")         return fmt::format(FMT_STRING("{:>8}")    , param.spin_dim);
        if(p == "distribution")     return fmt::format(FMT_STRING("{:<12}")   , param.get_distribution());
        /* clang-format on */
        throw except::runtime_error("Unrecognized parameter: {}", p);
    }

    static std::vector<std::string> get_parameter_names() { return {"J1", "J2", "h_tran", "h_mean", "h_wdth", "h_rand", "spin_dim", "distribution"}; }

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
    public:
    struct table {
        double   J1_rand      = 0;       /*!< On-site interaction */
        hvl_t    J2_rand      = {};      /*!< Two-body interaction */
        double   J3_rand      = 0;       /*!< Three-body interaction */
        double   J1_mean      = 0;       /*!< Constant offset for on-site */
        double   J2_mean      = 0;       /*!< Constant offset for two-body interaction */
        double   J3_mean      = 0;       /*!< Constant offset for three-body interaction */
        double   J1_wdth      = 0;       /*!< Width of the distribution J1 */
        double   J2_wdth      = 0;       /*!< Width of the distribution J2 */
        double   J3_wdth      = 0;       /*!< Width of the distribution J3 */
        double   J2_xcls      = 0;       /*!< Exp. decay rate of two-body interactions: exp(-|i-j|/J2_xcls) * J2_rand */
        size_t   J2_span      = 0;       /*!< Maximum range for pairwise interactions, |i-j| <= J2_span. */
        size_t   J2_ctof      = 0;       /*!< Effective range for pairwise interactions, |i-j| <= std::min(J2_span, model_size-1). */
        double   f_mixer      = 0;       /*!< Mixing factor for unitary transformation to real-space */
        uint64_t u_layer      = 0;       /*!< Number of unitary 2-site layers which transform lbit <-> real spaces */
        long     spin_dim     = 2;       /*!< Spin dimension */
        char    *distribution = nullptr; /*!< The random distribution of J_rnd and h_rnd. Choose between lognormal, normal or uniform */
        table() noexcept : J2_rand(hvl_t{0, nullptr}) {}
        ~table() noexcept {
            free(distribution);
            free(J2_rand.p);
        }
        std::string_view get_distribution() const {
            if(distribution == nullptr) throw except::runtime_error("distribution == nullptr");
            return distribution;
        }
        void set_distribution(std::string_view dist) {
            free(distribution);
            distribution = static_cast<char *>(malloc((dist.size() + 1) * sizeof(char))); // Add +1 for null terminator
            strcpy(distribution, dist.data());
        }
        std::vector<double> get_J2_rand() const {
            if(J2_rand.p == nullptr) throw except::runtime_error("J2_rand.p == nullptr");
            auto bgn = static_cast<double *>(J2_rand.p);
            auto end = static_cast<double *>(J2_rand.p) + J2_rand.len;
            return std::vector<double>(bgn, end);
        }
        size_t get_J2_rand_size() const { return J2_rand.len; }
        void   set_J2_rand(const std::vector<double> &J2) {
            free(J2_rand.p);
            J2_rand.p   = static_cast<double *>(malloc(J2.size() * sizeof(double)));
            J2_rand.len = J2.size();
            std::copy_n(J2.data(), J2.size(), static_cast<double *>(J2_rand.p));
        }
    };

    static inline h5pp::hid::h5t h5_type;
    table                        param;

    h5tb_lbit() { register_table_type(); }

    static void register_table_type() {
        if(h5_type.valid()) return;
        h5pp::hid::h5t H5T_VARRAY_DOUBLE = H5Tvlen_create(H5T_NATIVE_DOUBLE);

        // Create a type for the char array from the template H5T_C_S1
        h5pp::hid::h5t h5t_custom_string = H5Tcopy(H5T_C_S1);
        H5Tset_size(h5t_custom_string, H5T_VARIABLE);
        // Optionally set the null terminator '\0'
        H5Tset_strpad(h5t_custom_string, H5T_STR_NULLTERM);
        h5_type = H5Tcreate(H5T_COMPOUND, sizeof(table));
        H5Tinsert(h5_type, "J1_rand", HOFFSET(table, J1_rand), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "J2_rand", HOFFSET(table, J2_rand), H5T_VARRAY_DOUBLE);
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
        H5Tinsert(h5_type, "f_mixer", HOFFSET(table, f_mixer), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "u_layer", HOFFSET(table, u_layer), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "spin_dim", HOFFSET(table, spin_dim), H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "distribution", HOFFSET(table, distribution), h5t_custom_string);
    }

    [[nodiscard]] std::string J2_str() const { return fmt::format(FMT_STRING("[{:<+9.2e}]"), fmt::join(param.get_J2_rand(), ",")); }
    [[nodiscard]] std::string fmt_value(std::string_view p) const {
        // J2_rand is special since it varies in length for each mpo. Let's just pad with nan to make it pretty
        if(p == "J2_rand") {
            auto J2_rand = param.get_J2_rand();
            J2_rand.reserve(param.J2_ctof + 1);
            for(size_t i = J2_rand.size(); i < param.J2_ctof + 1; i++) J2_rand.emplace_back(std::numeric_limits<double>::quiet_NaN());
            return fmt::format(FMT_STRING("[{:<+9.2e}]"), fmt::join(J2_rand, ","));
        }

        /* clang-format off */
        if(p == "J1_rand")     return fmt::format(FMT_STRING("{:<+9.2e}"),param.J1_rand);
        if(p == "J3_rand")     return fmt::format(FMT_STRING("{:<+9.2e}"), param.J3_rand);
        if(p == "J1_mean")     return fmt::format(FMT_STRING("{:<+9.2e}"), param.J1_mean);
        if(p == "J2_mean")     return fmt::format(FMT_STRING("{:<+9.2e}"), param.J2_mean);
        if(p == "J3_mean")     return fmt::format(FMT_STRING("{:<+9.2e}"), param.J3_mean);
        if(p == "J1_wdth")     return fmt::format(FMT_STRING("{:<7.4f}"),  param.J1_wdth);
        if(p == "J2_wdth")     return fmt::format(FMT_STRING("{:<7.4f}"),  param.J2_wdth);
        if(p == "J3_wdth")     return fmt::format(FMT_STRING("{:<7.4f}"),  param.J3_wdth);
        if(p == "J2_xcls")     return fmt::format(FMT_STRING("{:<7.4f}"),  param.J2_xcls);
        if(p == "J2_span")     return fmt::format(FMT_STRING("{:>7}"),     param.J2_span == -1ul ? -1l : static_cast<long>(param.J2_span));
        if(p == "J2_ctof")     return fmt::format(FMT_STRING("{:>7}"),     param.J2_ctof);
        if(p == "f_mixer")     return fmt::format(FMT_STRING("{:<7.4f}"),  param.f_mixer);
        if(p == "u_layer")     return fmt::format(FMT_STRING("{:>7}"),     param.u_layer);
        if(p == "spin_dim")    return fmt::format(FMT_STRING("{:>8}"),     param.spin_dim);
        if(p == "distribution")return fmt::format(FMT_STRING("{:<12}"),    param.get_distribution());
        /* clang-format on */
        throw except::runtime_error("Unrecognized parameter: {}", p);
    }

    static std::vector<std::string> get_parameter_names() {
        return {"J1_rand", "J2_rand", "J3_rand", "J1_mean", "J2_mean", "J3_mean", "J1_wdth",  "J2_wdth",
                "J3_wdth", "J2_xcls", "J2_span", "J2_ctof", "f_mixer", "u_layer", "spin_dim", "distribution"};
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
