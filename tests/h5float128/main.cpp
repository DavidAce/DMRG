
#if defined(USE_QUADMATH)
    #include "debug/exceptions.h"
    #include "io/fmt.h"
    #include "io/fmt_f128_t.h"
    #include "math/f128.h"
    #include <bitset>
    #include <fmt/format.h>
    #include <h5pp/h5pp.h>
template<typename T>
std::string get_binary(T val) {
    auto s = std::string();
    //    s.resize(sizeof(val));
    //    memcpy(s.data(), &val, sizeof(val));
    //    std::string b;
    //    for(const auto c : s) b += fmt::format("{:08b}", c);
    //    //    for(const auto c: s) b += fmt::format("{:08b}", static_cast<int>(c));
    //    std::reverse(b.begin(), b.end());
    //
    //    return b;

    auto p = reinterpret_cast<unsigned char *>(&val);
    for(size_t i = 0; i < sizeof(T); ++i) {
        auto bs = fmt::format("{:08b}", p[i]);
        std::reverse(bs.begin(), bs.end());
        //        fmt::print("p[{}] = {} -> {}\n", i, p[i], bs);
        s += bs;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

template<typename T>
T get_decimal(std::string s) {
    //    std::string b;
    std::reverse(s.begin(), s.end());
    std::vector<uint8_t> b;
    for(size_t i = 0; i < sizeof(T); ++i) {
        auto bs = s.substr(i * 8, 8);
        //        fmt::print("b[{}] = {}\n", i, bs);
        //        b += bs;
        std::reverse(bs.begin(), bs.end());

        uint8_t n = static_cast<uint8_t>(std::stoi(bs, 0, 2));
        b.emplace_back(static_cast<uint8_t>(n));
    }
    //    fmt::print("as bytes: {}\n", b);
    T val;
    memcpy(&val, b.data(), sizeof(T));
    return val;
}

int main() {
    auto tval = float(5.0);
    auto tbin = get_binary(tval);
    auto wval = __float128(0.1);
    auto wbin = get_binary(wval);
    fmt::print("tval = {} -> {} -> {}\n", tval, tbin, get_decimal<float>(tbin));
    fmt::print("wval = {} -> {} -> {}\n", f128_t(wval), wbin, f128_t(get_decimal<__float128>(wbin)));

    auto f = h5pp::File("h5float128.h5", h5pp::FilePermission::READWRITE, 2);

    //    auto h5_float128_t = h5pp::hid::h5t();
    hid_t h5_float128_t;
    #if __BYTE_ORDER == LITTLE_ENDIAN
    h5_float128_t = H5Tcopy(H5T_IEEE_F64LE);
    #else
    h5_float128_t = H5Tcopy(H5T_IEEE_F64BE);
    #endif
    auto sserr = H5Tset_size(h5_float128_t, 16);
    auto sperr = H5Tset_precision(h5_float128_t, 128);
    auto soerr = H5Tset_offset(h5_float128_t, 0);
    #if __BYTE_ORDER == LITTLE_ENDIAN
    auto sferr = H5Tset_fields(h5_float128_t, 127, 112, 15, 0, 112);
    #else
    auto sferr = H5Tset_fields(h5_real_t, 0, 1, 15, 16, 112);
    #endif
    auto sberr = H5Tset_ebias(h5_float128_t, 127);
    auto snerr = H5Tset_norm(h5_float128_t, H5T_norm_t::H5T_NORM_MSBSET);
    auto sierr = H5Tset_inpad(h5_float128_t, H5T_pad_t::H5T_PAD_ZERO);
    auto sparr = H5Tset_pad(h5_float128_t, H5T_pad_t::H5T_PAD_ZERO, H5T_pad_t::H5T_PAD_ZERO);
    if(sserr < 0) throw except::runtime_error("H5Tset_size returned {}", sserr);
    if(sperr < 0) throw except::runtime_error("H5Tset_precision returned {}", sperr);
    if(soerr < 0) throw except::runtime_error("H5Tset_offset returned {}", soerr);
    if(sferr < 0) throw except::runtime_error("H5Tset_fields returned {}", sferr);
    if(sberr < 0) throw except::runtime_error("H5Tset_ebias returned {}", sberr);
    if(snerr < 0) throw except::runtime_error("H5Tset_norm returned {}", snerr);
    if(sierr < 0) throw except::runtime_error("H5Tset_inpad returned {}", sierr);
    if(sparr < 0) throw except::runtime_error("H5Tset_pad returned {}", sparr);
    //    H5Tcommit(f.openFileHandle(), "H5T_IEEE_F128LE", h5_float128_t, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //    f.writeDataset(val, "float128", h5_float128_t);
    auto o     = h5pp::Options();
    o.linkPath = "float128";
    o.h5Type   = h5_float128_t;
    o.dataDims = {};

    auto dataInfo = h5pp::scan::scanDataInfo(wval, o);
    auto dsetInfo = h5pp::scan::makeDsetInfo(f.openFileHandle(), o, f.plists);
    h5pp::hdf5::createDataset(dsetInfo, f.plists);
    h5pp::hdf5::writeDataset(wval, dataInfo, dsetInfo, f.plists);

    auto rval = f.readDataset<__float128>("float128", std::nullopt, h5_float128_t);
    auto rbin = get_binary(rval);
    fmt::print("rval = {} -> {} -> {}\n", f128_t(rval), rbin, f128_t(get_decimal<__float128>(rbin)));
    if(rval != wval) throw std::runtime_error("rval != wval");
    return 0;
}
#else
int main() { return 0; }
#endif