#include "h5tb.h"
#include "debug/exceptions.h"
#include "tools/common/log.h"
#include <hdf5.h>

const h5pp::hid::h5t &h5tb_base::get_h5_type() const {
    register_table_type();
    return h5_type;
}
void h5tb_base::print_parameter_names() const noexcept {
    std::string name_line;
    for(const auto &name : get_parameter_names()) name_line.append(fmt::format(FMT_STRING("{:<{}} "), name, fmt_value(name).size()));
    tools::log->info(name_line);
}

void h5tb_base::print_parameter_values() const noexcept {
    std::string value_line;
    for(const auto &name : get_parameter_names()) value_line.append(fmt::format(FMT_STRING("{} "), fmt_value(name)));
    tools::log->info(value_line);
}

/*
 *
 *
 *  ISING SELFDUAL
 *
 *
 */

void h5tb_ising_selfdual::register_table_type() const {
    if(h5_type.valid()) return;

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
    H5Tinsert(h5_type, "distribution", HOFFSET(table, distribution), h5pp::vstr_t::get_h5type());
}

std::string h5tb_ising_selfdual::fmt_value(std::string_view p) const {
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
        if(p == "distribution")     return fmt::format(FMT_STRING("{:<12}")   , param.distribution);
    /* clang-format on */
    throw except::runtime_error("Unrecognized parameter: {}", p);
}

std::vector<std::string_view> h5tb_ising_selfdual::get_parameter_names() const noexcept {
    return {"J_mean", "J_wdth", "J_rand", "h_mean", "h_wdth", "h_rand", "lambda", "delta", "spin_dim", "distribution"};
}

/*
 *
 *
 *  ISING MAJORANA
 *
 *
 */

void h5tb_ising_majorana::register_table_type() const {
    if(h5_type.valid()) return;

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
    H5Tinsert(h5_type, "distribution", HOFFSET(table, distribution), h5pp::vstr_t::get_h5type());
}

std::string h5tb_ising_majorana::fmt_value(std::string_view p) const {
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
        if(p == "distribution")     return fmt::format(FMT_STRING("{:<12}")   , param.distribution);
    /* clang-format on */
    throw except::runtime_error("Unrecognized parameter: {}", p);
}

std::vector<std::string_view> h5tb_ising_majorana::get_parameter_names() const noexcept {
    return {"J_mean", "J_wdth", "J_rand", "h_mean", "h_wdth", "h_rand", "g", "delta", "spin_dim", "distribution"};
}

/*
 *
 *
 *  ISING TRANSVERSE RANDOM FIELD
 *
 *
 */

void h5tb_ising_tf_rf::register_table_type() const {
    if(h5_type.valid()) return;

    h5_type = H5Tcreate(H5T_COMPOUND, sizeof(table));
    H5Tinsert(h5_type, "J1", HOFFSET(table, J1), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5_type, "J2", HOFFSET(table, J2), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5_type, "h_tran", HOFFSET(table, h_tran), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5_type, "h_mean", HOFFSET(table, h_mean), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5_type, "h_wdth", HOFFSET(table, h_wdth), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5_type, "h_rand", HOFFSET(table, h_rand), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5_type, "spin_dim", HOFFSET(table, spin_dim), H5T_NATIVE_LONG);
    H5Tinsert(h5_type, "distribution", HOFFSET(table, distribution), h5pp::vstr_t::get_h5type());
}

std::string h5tb_ising_tf_rf::fmt_value(std::string_view p) const {
    /* clang-format off */
        if(p == "J1")               return fmt::format(FMT_STRING("{:<+9.2e}"), param.J1);
        if(p == "J2")               return fmt::format(FMT_STRING("{:<+9.2e}"), param.J2);
        if(p == "h_tran")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.h_tran);
        if(p == "h_mean")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.h_mean);
        if(p == "h_wdth")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.h_wdth);
        if(p == "h_rand")           return fmt::format(FMT_STRING("{:<+9.2e}"), param.h_rand);
        if(p == "spin_dim")         return fmt::format(FMT_STRING("{:>8}")    , param.spin_dim);
        if(p == "distribution")     return fmt::format(FMT_STRING("{:<12}")   , param.distribution);
    /* clang-format on */
    throw except::runtime_error("Unrecognized parameter: {}", p);
}

std::vector<std::string_view> h5tb_ising_tf_rf::get_parameter_names() const noexcept {
    return {"J1", "J2", "h_tran", "h_mean", "h_wdth", "h_rand", "spin_dim", "distribution"};
}

/*
 *
 *
 *  LBIT
 *
 *
 */

h5pp::hid::h5t &h5tb_lbit::get_h5t_enum_w8() {
    create_enum_w8();
    return enum_w8;
}

void h5tb_lbit::create_enum_w8() {
    if(enum_w8.valid()) return;
    enum_w8 = H5Tenum_create(H5T_NATIVE_INT);
    int val;
    H5Tenum_insert(enum_w8, "IDENTITY", (val = static_cast<int>(UnitaryGateWeight::IDENTITY), &val));
    H5Tenum_insert(enum_w8, "EXPDECAY", (val = static_cast<int>(UnitaryGateWeight::EXPDECAY), &val));
}

void h5tb_lbit::commit_enum_w8(const h5pp::hid::h5f &file_id) {
    if(H5Tcommitted(get_h5t_enum_w8()) > 0) return;
    herr_t err = H5Tcommit(file_id, "UnitaryGateWeight", get_h5t_enum_w8(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(err < 0) throw except::runtime_error("Failed to commit StorageEvent to file");
}

void h5tb_lbit::register_table_type() const {
    if(h5_type.valid()) return;
    h5_type = H5Tcreate(H5T_COMPOUND, sizeof(table));
    H5Tinsert(h5_type, "J1_rand", HOFFSET(table, J1_rand), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5_type, "J2_rand", HOFFSET(table, J2_rand), h5pp::varr_t<double>::get_h5type());
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
    H5Tinsert(h5_type, "u_depth", HOFFSET(table, u_depth), H5T_NATIVE_UINT64);
    H5Tinsert(h5_type, "u_fmix", HOFFSET(table, u_fmix), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5_type, "u_tstd", HOFFSET(table, u_tstd), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5_type, "u_cstd", HOFFSET(table, u_cstd), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5_type, "u_tgw8", HOFFSET(table, u_tgw8), get_h5t_enum_w8());
    H5Tinsert(h5_type, "u_cgw8", HOFFSET(table, u_cgw8), get_h5t_enum_w8());
    H5Tinsert(h5_type, "spin_dim", HOFFSET(table, spin_dim), H5T_NATIVE_LONG);
    H5Tinsert(h5_type, "distribution", HOFFSET(table, distribution), h5pp::vstr_t::get_h5type());
}

std::string h5tb_lbit::fmt_value(std::string_view p) const {
    // J2_rand is special since it varies in length for each mpo. Let's just pad with nan to make it pretty
    if(p == "J2_rand") {
        std::vector<double> J2_rand = param.J2_rand;
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
        if(p == "u_depth")     return fmt::format(FMT_STRING("{:>7}"),     param.u_depth);
        if(p == "u_fmix")      return fmt::format(FMT_STRING("{:<7.4f}"),  param.u_fmix);
        if(p == "u_tstd")      return fmt::format(FMT_STRING("{:<7.4f}"),  param.u_tstd);
        if(p == "u_cstd")      return fmt::format(FMT_STRING("{:<7.4f}"),  param.u_cstd);
        if(p == "u_tgw8")      return fmt::format(FMT_STRING("{}"),        enum2sv(param.u_tgw8));
        if(p == "u_cgw8")      return fmt::format(FMT_STRING("{}"),        enum2sv(param.u_cgw8));
        if(p == "spin_dim")    return fmt::format(FMT_STRING("{:>8}"),     param.spin_dim);
        if(p == "distribution")return fmt::format(FMT_STRING("{:<12}"),    param.distribution);
        /* clang-format on */
        throw except::runtime_error("Unrecognized parameter: {}", p);
}

std::vector<std::string_view> h5tb_lbit::get_parameter_names() const noexcept {
    return {"J1_rand", "J2_rand", "J3_rand", "J1_mean", "J2_mean", "J3_mean", "J1_wdth",  "J2_wdth",
            "J3_wdth", "J2_xcls", "J2_span", "J2_ctof", "u_fmix",  "u_depth", "spin_dim", "distribution"};
}
