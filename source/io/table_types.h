#pragma once

#include "algorithms/AlgorithmStatus.h"
#include "io/fmt.h"
#include "tid/enums.h"
#include <array>
#include <cstdio>
#include <h5pp/h5pp.h>
#include <memory>
#include <unordered_map>
#include <vector>

class h5pp_table_measurements_finite {
    public:
    static inline h5pp::hid::h5t h5_type;

    struct table {
        uint64_t              iter                          = 0;
        uint64_t              step                          = 0;
        long                  position                      = 0;
        uint64_t              length                        = 0;
        long                  bond_mid                      = 0;
        long                  bond_lim                      = 0;
        long                  bond_max                      = 0;
        double                entanglement_entropy_midchain = 0;
        double                number_entropy_midchain       = 0;
        double                norm                          = 0;
        double                energy                        = 0;
        double                energy_variance               = 0;
        double                energy_variance_lowest        = 0;
        std::array<double, 3> spin_components               = {0};
        double                truncation_error              = 0;
        double                total_time                    = 0;
        double                algorithm_time                = 0;
        double                physical_time                 = 0;
    };

    h5pp_table_measurements_finite() { register_table_type(); }
    static void register_table_type() {
        if(h5_type.valid()) return;
        std::array<hsize_t, 1> spin_dims    = {3};
        h5pp::hid::h5t         H5_SPIN_TYPE = H5Tarray_create(H5T_NATIVE_DOUBLE, spin_dims.size(), spin_dims.data());

        h5_type = H5Tcreate(H5T_COMPOUND, sizeof(table));
        H5Tinsert(h5_type, "iter", HOFFSET(table, iter), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "step", HOFFSET(table, step), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "position", HOFFSET(table, position), H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "length", HOFFSET(table, length), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "bond_mid", HOFFSET(table, bond_mid), H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "bond_lim", HOFFSET(table, bond_lim), H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "bond_max", HOFFSET(table, bond_max), H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "entanglement_entropy_midchain", HOFFSET(table, entanglement_entropy_midchain), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "number_entropy_midchain", HOFFSET(table, number_entropy_midchain), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "norm", HOFFSET(table, norm), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy", HOFFSET(table, energy), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_variance", HOFFSET(table, energy_variance), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_variance_lowest", HOFFSET(table, energy_variance_lowest), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "spin_components", HOFFSET(table, spin_components), H5_SPIN_TYPE);
        H5Tinsert(h5_type, "truncation_error", HOFFSET(table, truncation_error), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "total_time", HOFFSET(table, total_time), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "algorithm_time", HOFFSET(table, algorithm_time), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "physical_time", HOFFSET(table, physical_time), H5T_NATIVE_DOUBLE);
    }
};

class h5pp_table_measurements_infinite {
    public:
    static inline h5pp::hid::h5t h5_type;

    struct table {
        uint64_t             iter                         = 0;
        uint64_t             step                         = 0;
        long                 position                     = 0;
        uint64_t             length                       = 0;
        long                 bond_dim                     = 0;
        long                 bond_lim                     = 0;
        long                 bond_max                     = 0;
        double               entanglement_entropy         = 0;
        double               norm                         = 0;
        double               energy_mpo                   = 0;
        double               energy_per_site_mpo          = 0;
        double               energy_per_site_ham          = 0;
        double               energy_per_site_mom          = 0;
        double               energy_variance_mpo          = 0;
        double               energy_variance_per_site_mpo = 0;
        double               energy_variance_per_site_ham = 0;
        double               energy_variance_per_site_mom = 0;
        double               truncation_error             = 0;
        double               wall_time                    = 0;
        double               phys_time                    = 0;
        std::complex<double> time_step;
    };

    h5pp_table_measurements_infinite() { register_table_type(); }
    static void register_table_type() {
        if(h5_type.valid()) return;
        h5_type = H5Tcreate(H5T_COMPOUND, sizeof(table));
        H5Tinsert(h5_type, "iter", HOFFSET(table, iter), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "step", HOFFSET(table, step), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "position", HOFFSET(table, position), H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "length", HOFFSET(table, length), H5T_NATIVE_UINT64);
        H5Tinsert(h5_type, "bond_dim", HOFFSET(table, bond_dim), H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "bond_lim", HOFFSET(table, bond_lim), H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "bond_max", HOFFSET(table, bond_max), H5T_NATIVE_LONG);
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
        H5Tinsert(h5_type, "time_step", HOFFSET(table, time_step), h5pp::type::compound::H5T_COMPLEX<double>::h5type());
    }
};

class h5pp_table_algorithm_status {
    public:
    static inline h5pp::hid::h5t h5_type;
    static inline h5pp::hid::h5t h5_algo_type;
    static inline h5pp::hid::h5t h5_algo_stop;
    using table = AlgorithmStatus;

    h5pp_table_algorithm_status() { register_table_type(); }
    static void register_table_type() {
        if(h5_type.valid()) return;
        h5_type = H5Tcreate(H5T_COMPOUND, sizeof(table));

        if(not h5_algo_type.valid()) {
            h5_algo_type = H5Tcreate(H5T_ENUM, sizeof(AlgorithmType));
            int val;
            H5Tenum_insert(h5_algo_type, "iDMRG", (val = 0, &val));
            H5Tenum_insert(h5_algo_type, "fDMRG", (val = 1, &val));
            H5Tenum_insert(h5_algo_type, "xDMRG", (val = 2, &val));
            H5Tenum_insert(h5_algo_type, "iTEBD", (val = 3, &val));
            H5Tenum_insert(h5_algo_type, "flBIT", (val = 4, &val));
            H5Tenum_insert(h5_algo_type, "ANY", (val = 5, &val));
        }
        if(not h5_algo_stop.valid()) {
            h5_algo_stop = H5Tcreate(H5T_ENUM, sizeof(AlgorithmStop));
            int val;
            H5Tenum_insert(h5_algo_stop, "SUCCESS", (val = 0, &val));
            H5Tenum_insert(h5_algo_stop, "SATURATED", (val = 1, &val));
            H5Tenum_insert(h5_algo_stop, "MAX_ITERS", (val = 2, &val));
            H5Tenum_insert(h5_algo_stop, "MAX_RESET", (val = 3, &val));
            H5Tenum_insert(h5_algo_stop, "RANDOMIZE", (val = 4, &val));
            H5Tenum_insert(h5_algo_stop, "NONE", (val = 5, &val));
        }

        /* clang-format off */
        H5Tinsert(h5_type, "iter",                        HOFFSET(table, iter),                       H5T_NATIVE_ULONG);
        H5Tinsert(h5_type, "step",                        HOFFSET(table, step),                       H5T_NATIVE_ULONG);
        H5Tinsert(h5_type, "position",                    HOFFSET(table, position),                   H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "direction",                   HOFFSET(table, direction),                  H5T_NATIVE_INT);
        H5Tinsert(h5_type, "num_resets",                  HOFFSET(table, num_resets),                 H5T_NATIVE_ULONG);
        H5Tinsert(h5_type, "min_iters",                   HOFFSET(table, min_iters),                  H5T_NATIVE_ULONG);
        H5Tinsert(h5_type, "bond_lim",                    HOFFSET(table, bond_lim),                   H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "bond_max",                    HOFFSET(table, bond_max),                   H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "bond_init",                   HOFFSET(table, bond_init),                  H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "energy_min",                  HOFFSET(table, energy_min),                 H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_max",                  HOFFSET(table, energy_max),                 H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_tgt",                  HOFFSET(table, energy_tgt),                 H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_ulim",                 HOFFSET(table, energy_ulim),                H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_llim",                 HOFFSET(table, energy_llim),                H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_dens",                 HOFFSET(table, energy_dens),                H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_dens_target",          HOFFSET(table, energy_dens_target),         H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_dens_window",          HOFFSET(table, energy_dens_window),         H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_variance_lowest",      HOFFSET(table, energy_variance_lowest),     H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "energy_variance_max_digits",  HOFFSET(table, energy_variance_max_digits), H5T_NATIVE_ULONG);
        H5Tinsert(h5_type, "energy_variance_prec_limit",  HOFFSET(table, energy_variance_prec_limit), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "env_expansion_alpha",         HOFFSET(table, env_expansion_alpha),        H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "env_expansion_variance",      HOFFSET(table, env_expansion_variance),     H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "env_expansion_step",          HOFFSET(table, env_expansion_step),         H5T_NATIVE_ULONG);
        H5Tinsert(h5_type, "phys_time",                   HOFFSET(table, phys_time),                  H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "wall_time",                   HOFFSET(table, wall_time),                  H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "algo_time",                   HOFFSET(table, algo_time),                  H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "delta_t",                     HOFFSET(table, delta_t),                    h5pp::type::compound::H5T_COMPLEX<double>::h5type());
        H5Tinsert(h5_type, "algo_type",                   HOFFSET(table, algo_type),                  h5_algo_type);
        H5Tinsert(h5_type, "algo_stop",                   HOFFSET(table, algo_stop),                  h5_algo_stop);
        H5Tinsert(h5_type, "algorithm_has_finished",      HOFFSET(table, algorithm_has_finished),     H5T_NATIVE_UINT8);
        H5Tinsert(h5_type, "algorithm_has_succeeded",     HOFFSET(table, algorithm_has_succeeded),    H5T_NATIVE_UINT8);
        H5Tinsert(h5_type, "algorithm_has_to_stop",       HOFFSET(table, algorithm_has_to_stop),      H5T_NATIVE_UINT8);
        H5Tinsert(h5_type, "algorithm_has_stuck_for",     HOFFSET(table, algorithm_has_stuck_for),    H5T_NATIVE_ULONG);
        H5Tinsert(h5_type, "algorithm_saturated_for",     HOFFSET(table, algorithm_saturated_for),    H5T_NATIVE_ULONG);
        H5Tinsert(h5_type, "algorithm_converged_for",     HOFFSET(table, algorithm_converged_for),    H5T_NATIVE_ULONG);
        H5Tinsert(h5_type, "entanglement_converged_for",  HOFFSET(table, entanglement_converged_for), H5T_NATIVE_ULONG);
        H5Tinsert(h5_type, "entanglement_saturated_for",  HOFFSET(table, entanglement_saturated_for), H5T_NATIVE_ULONG);
        H5Tinsert(h5_type, "variance_mpo_converged_for",  HOFFSET(table, variance_mpo_converged_for), H5T_NATIVE_ULONG);
        H5Tinsert(h5_type, "variance_mpo_saturated_for",  HOFFSET(table, variance_mpo_saturated_for), H5T_NATIVE_ULONG);
        H5Tinsert(h5_type, "variance_ham_converged_for",  HOFFSET(table, variance_ham_converged_for), H5T_NATIVE_ULONG);
        H5Tinsert(h5_type, "variance_ham_saturated_for",  HOFFSET(table, variance_ham_saturated_for), H5T_NATIVE_ULONG);
        H5Tinsert(h5_type, "variance_mom_converged_for",  HOFFSET(table, variance_mom_converged_for), H5T_NATIVE_ULONG);
        H5Tinsert(h5_type, "variance_mom_saturated_for",  HOFFSET(table, variance_mom_saturated_for), H5T_NATIVE_ULONG);
        H5Tinsert(h5_type, "bond_limit_has_reached_max",  HOFFSET(table, bond_limit_has_reached_max), H5T_NATIVE_UINT8);
        H5Tinsert(h5_type, "spin_parity_has_converged",   HOFFSET(table, spin_parity_has_converged),  H5T_NATIVE_UINT8);
        H5Tinsert(h5_type, "time_step_has_converged",     HOFFSET(table, time_step_has_converged),    H5T_NATIVE_UINT8);
        H5Tinsert(h5_type, "fes_is_running",              HOFFSET(table, fes_is_running),             H5T_NATIVE_UINT8);
        /* clang-format on */
    }
};

class h5pp_table_memory_usage {
    public:
    static inline h5pp::hid::h5t h5_type;

    struct table {
        uint64_t iter;
        uint64_t step;
        int64_t  bond_lim;
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
        H5Tinsert(h5_type, "bond_lim", HOFFSET(table, bond_lim), H5T_NATIVE_LONG);
        H5Tinsert(h5_type, "rss", HOFFSET(table, rss), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "hwm", HOFFSET(table, hwm), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "vm", HOFFSET(table, vm), H5T_NATIVE_DOUBLE);
    }
};

template<typename T>
class h5pp_table_data {
    private:
    struct meta {
        size_t      type_size;
        size_t      data_size;
        std::string fieldname;
        //        bool operator == (const meta & rhs) const = default;
        bool operator==(const meta &rhs) const { return type_size == rhs.type_size and data_size == rhs.data_size and fieldname == rhs.fieldname; }
    };
    struct metaHasher {
        auto operator()(const meta &m) const { return std::hash<std::string>{}(fmt::format("{}|{}|{}", m.type_size, m.data_size, m.fieldname)); }
    };

    public:
    static inline std::unordered_map<meta, h5pp::hid::h5t, metaHasher> h5_types;

    static std::vector<std::byte> make_entry(uint64_t iter, uint64_t step, int64_t bond_lim, const T *const data, size_t data_size) {
        size_t                 total_size = 2 * sizeof(uint64_t) + 1 * sizeof(int64_t) + data_size * sizeof(T);
        std::vector<std::byte> entry(total_size);
        std::memcpy(entry.data() + 0 * sizeof(uint64_t), &iter, sizeof(uint64_t));
        std::memcpy(entry.data() + 1 * sizeof(uint64_t), &step, sizeof(uint64_t));
        std::memcpy(entry.data() + 2 * sizeof(int64_t), &bond_lim, sizeof(int64_t));
        std::memcpy(entry.data() + 3 * sizeof(uint64_t), data, data_size * sizeof(T));
        return entry;
    }

    [[nodiscard]] static h5pp::hid::h5t &register_table_type(size_t data_size, std::string_view fieldname) {
        size_t data_offset = 2 * sizeof(uint64_t) + 1 * sizeof(int64_t);
        size_t total_size  = data_offset + data_size * sizeof(T);
        meta   m           = {sizeof(T), data_size, std::string(fieldname)};
        if(h5_types.find(m) == h5_types.end()) {
            h5pp::hid::h5t h5_type = H5Tcreate(H5T_COMPOUND, total_size);
            H5Tinsert(h5_type, "iter", 0 * sizeof(uint64_t), H5T_NATIVE_UINT64);
            H5Tinsert(h5_type, "step", 1 * sizeof(uint64_t), H5T_NATIVE_UINT64);
            H5Tinsert(h5_type, "bond_lim", 2 * sizeof(uint64_t), H5T_NATIVE_INT64);
            auto h5type = h5pp::util::getH5Type<T>();
            if(data_size == 1)
                H5Tinsert(h5_type, std::string(fieldname).c_str(), data_offset, h5type);
            else
                for(size_t elem = 0; elem < data_size; elem++) {
                    H5Tinsert(h5_type, fmt::format("{}{}", fieldname, elem).c_str(), data_offset + elem * sizeof(T), h5type);
                }
            h5_types[m] = h5_type;
        }
        return h5_types[m];
    }
};

class h5pp_ur {
    public:
    static inline h5pp::hid::h5t h5_type;
    static inline h5pp::hid::h5t h5_level_type;
    static constexpr size_t      nlen = 255;
    struct item {
        char  *name;
        double time;
        double sum;
        double pcnt;
        double avg;
        int    level;
        size_t count;
        item() = default;
        item(const item &it) : time(it.time), sum(it.sum), pcnt(it.pcnt), avg(it.avg), level(it.level), count(it.count) { copy_name(it.name); };
        item(std::string_view name_, double time, double sum, double pcnt, double avg, int level, size_t count)
            : time(time), sum(sum), pcnt(pcnt), avg(avg), level(level), count(count) {
            copy_name(name_);
        }
        // Automatically deallocate when no longer needed
        ~item() { free(name); }

        private:
        void copy_name(std::string_view name_) {
            // Note that we must use C-style malloc/free here rather than C++-style new/delete, since that is what HDF5 uses internally.
            // This is particularly important when we read data from file, and let HDF5 allocate the vlen buffer
            name = static_cast<char *>(malloc((name_.size() + 1) * sizeof(char))); // Add +1 for null terminator
            strncpy(name, name_.data(), name_.size());
            name[name_.size()] = '\0'; // Make sure name is null-terminated!
        }
    };

    h5pp_ur() { register_table_type(); }
    static void register_table_type() {
        if(h5_type.valid()) return;
        if(not h5_level_type.valid()) {
            h5_level_type = H5Tcreate(H5T_ENUM, sizeof(tid::level));
            int val;
            H5Tenum_insert(h5_level_type, "parent", (val = tid::level::parent, &val));
            H5Tenum_insert(h5_level_type, "normal", (val = tid::level::normal, &val));
            H5Tenum_insert(h5_level_type, "extra", (val = tid::level::extra, &val));
            H5Tenum_insert(h5_level_type, "detailed", (val = tid::level::detailed, &val));
        }

        h5_type = H5Tcreate(H5T_COMPOUND, sizeof(item));

        // Create a type for the char array from the template H5T_C_S1
        // The template describes a string with a single char.
        // Set the size with H5Tset_size, or h5pp::hdf5::setStringSize(...)
        h5pp::hid::h5t h5t_custom_string = H5Tcopy(H5T_C_S1);
        H5Tset_size(h5t_custom_string, H5T_VARIABLE);

        // Optionally set the null terminator or pad with '\0'
        H5Tset_strpad(h5t_custom_string, H5T_STR_NULLTERM);
        H5Tset_cset(h5t_custom_string, H5T_CSET_UTF8);

        H5Tinsert(h5_type, "name", HOFFSET(item, name), h5t_custom_string);
        H5Tinsert(h5_type, "time", HOFFSET(item, time), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "sum", HOFFSET(item, sum), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "pcnt", HOFFSET(item, pcnt), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "avg", HOFFSET(item, avg), H5T_NATIVE_DOUBLE);
        H5Tinsert(h5_type, "level", HOFFSET(item, level), h5_level_type);
        H5Tinsert(h5_type, "count", HOFFSET(item, count), H5T_NATIVE_UINT64);
    }
};