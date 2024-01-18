#include "../lbit.h"
#include "config/enums.h"
#include "debug/exceptions.h"
#include "tools/common/log.h"
#include <h5pp/details/h5ppType.h>
#include <hdf5.h>

const h5pp::hid::h5t qm::lbit::UnitaryGateParameters::get_h5_type() {
    static h5pp::hid::h5t h5_type;
    if(h5_type.valid()) return h5_type;
    h5_type                               = H5Tcreate(H5T_COMPOUND, sizeof(qm::lbit::UnitaryGateParameters));
    std::array<hsize_t, 1> array2_dims    = {2};
    std::array<hsize_t, 1> array4_dims    = {4};
    h5pp::hid::h5t         H5_ARRAY2_TYPE = H5Tarray_create(H5T_NATIVE_ULONG, array2_dims.size(), array2_dims.data());
    h5pp::hid::h5t         H5_ARRAY4_TYPE = H5Tarray_create(H5T_NATIVE_DOUBLE, array4_dims.size(), array4_dims.data());

    H5Tinsert(h5_type, "layer", HOFFSET(qm::lbit::UnitaryGateParameters, layer), H5T_NATIVE_ULONG);
    H5Tinsert(h5_type, "sites", HOFFSET(qm::lbit::UnitaryGateParameters, sites), H5_ARRAY2_TYPE);
    H5Tinsert(h5_type, "f", HOFFSET(qm::lbit::UnitaryGateParameters, f), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5_type, "w", HOFFSET(qm::lbit::UnitaryGateParameters, w), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5_type, "theta", HOFFSET(qm::lbit::UnitaryGateParameters, theta), H5_ARRAY4_TYPE);
    H5Tinsert(h5_type, "c", HOFFSET(qm::lbit::UnitaryGateParameters, c), h5pp::type::getH5Type<decltype(qm::lbit::UnitaryGateParameters::c)>());
    H5Tinsert(h5_type, "g8w8", HOFFSET(qm::lbit::UnitaryGateParameters, g8w8), get_h5t_enum_uw());
    H5Tinsert(h5_type, "type", HOFFSET(qm::lbit::UnitaryGateParameters, type), get_h5t_enum_ut());
    return h5_type;
}

void qm::lbit::UnitaryGateParameters::print_parameter_names() const noexcept {
    std::string name_line;
    for(const auto &name : get_parameter_names()) name_line.append(fmt::format("{:<{}} ", name, fmt_value(name).size()));
    tools::log->info(name_line);
}

void qm::lbit::UnitaryGateParameters::print_parameter_values() const noexcept {
    std::string value_line;
    for(const auto &name : get_parameter_names()) value_line.append(fmt::format("{} ", fmt_value(name)));
    tools::log->info(value_line);
}

const h5pp::hid::h5t qm::lbit::UnitaryGateParameters::get_h5t_enum_ut() {
    static h5pp::hid::h5t enum_ut;
    if(enum_ut.valid()) return enum_ut;
    enum_ut = H5Tenum_create(H5T_NATIVE_INT);
    int val;
    H5Tenum_insert(enum_ut, "ANDERSON", (val = static_cast<int>(UnitaryGateType::ANDERSON), &val));
    H5Tenum_insert(enum_ut, "MBL", (val = static_cast<int>(UnitaryGateType::MBL), &val));
    return enum_ut;
}

const h5pp::hid::h5t qm::lbit::UnitaryGateParameters::get_h5t_enum_uw() {
    static h5pp::hid::h5t enum_uw;
    if(enum_uw.valid()) return enum_uw;
    enum_uw = H5Tenum_create(H5T_NATIVE_INT);
    int val;
    H5Tenum_insert(enum_uw, "IDENTITY", (val = static_cast<int>(UnitaryGateWeight::IDENTITY), &val));
    H5Tenum_insert(enum_uw, "EXPDECAY", (val = static_cast<int>(UnitaryGateWeight::EXPDECAY), &val));
    return enum_uw;
}

std::string qm::lbit::UnitaryGateParameters::fmt_value(std::string_view p) const {
    /* clang-format off */
        if(p == "layer")  return fmt::format("{:>7}"     , layer);
        if(p == "sites")  return fmt::format("{::4}"     , sites);
        if(p == "f")      return fmt::format("{:<7.4f}"  , f);
        if(p == "w")      return fmt::format("{:<9.2e}"  , w);
        if(p == "theta")  return fmt::format("{::<+9.2e}", theta);
        if(p == "c")      return fmt::format("{}"        , c);
        if(p == "g8w8")   return fmt::format("{}"        , enum2sv(g8w8));
        if(p == "type")   return fmt::format("{}"        , enum2sv(type));
    /* clang-format on */
    throw except::runtime_error("Unrecognized parameter: {}", p);
}

constexpr std::array<std::string_view, 8> qm::lbit::UnitaryGateParameters::get_parameter_names() noexcept {
    return {"layer", "sites", "f", "w", "theta", "c", "g8w8",  "type"};
}