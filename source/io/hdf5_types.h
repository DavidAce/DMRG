#pragma once

#include "algorithms/AlgorithmStatus.h"
#include "debug/exceptions.h"
#include "io/fmt.h"
#include "tid/enums.h"
#include <array>
#include <cstdio>
#include <h5pp/h5pp.h>
#include <memory>
#include <unordered_map>
#include <vector>

class h5_enum_storage_event {
    private:
    static inline h5pp::hid::h5t h5_storage_event;

    public:
    static h5pp::hid::h5t &get_h5t();
    static void            create();
    static void            commit(const h5pp::hid::h5f &file_id);
};

class h5_enum_storage_level {
    private:
    static inline h5pp::hid::h5t h5_storage_level;

    public:
    static h5pp::hid::h5t &get_h5t();
    static void            create();
    static void            commit(const h5pp::hid::h5f &file_id);
};

class h5pp_table_measurements_finite {
    public:
    static inline h5pp::hid::h5t h5_type;

    struct table {
        uint64_t              iter                          = 0;
        uint64_t              step                          = 0;
        long                  position                      = 0;
        StorageEvent          event                         = StorageEvent::NONE;
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

    h5pp_table_measurements_finite();
    static void register_table_type();
};

class h5pp_table_measurements_infinite {
    public:
    static inline h5pp::hid::h5t h5_type;

    struct table {
        uint64_t             iter                         = 0;
        uint64_t             step                         = 0;
        long                 position                     = 0;
        StorageEvent         event                        = StorageEvent::NONE;
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

    h5pp_table_measurements_infinite();
    static void register_table_type();
};

class h5pp_table_algorithm_status {
    public:
    static inline h5pp::hid::h5t h5_type;
    static inline h5pp::hid::h5t h5_algo_type;
    static inline h5pp::hid::h5t h5_algo_stop;
    using table = AlgorithmStatus;

    h5pp_table_algorithm_status();
    static void register_table_type();
};

class h5pp_table_memory_usage {
    public:
    static inline h5pp::hid::h5t h5_type;

    struct table {
        uint64_t     iter     = 0;
        uint64_t     step     = 0;
        long         position = 0;
        StorageEvent event    = StorageEvent::NONE;
        int64_t      bond_lim = 0;
        double       rss      = 0;
        double       hwm      = 0;
        double       vm       = 0;
    };

    h5pp_table_memory_usage();
    static void register_table_type();
};

template<typename T>
class h5pp_table_data {
    private:
    struct meta {
        size_t                  type_size;
        size_t                  data_size;
        std::string             fieldname;
        static constexpr size_t nfieldsFixed = 6;
        bool operator==(const meta &rhs) const { return type_size == rhs.type_size and data_size == rhs.data_size and fieldname == rhs.fieldname; }
        static std::array<size_t, nfieldsFixed> get_extent(size_t data_size) {
            return {sizeof(uint64_t), sizeof(uint64_t), sizeof(int64_t), sizeof(StorageEvent), sizeof(int64_t), data_size * sizeof(T)};
        }

        static std::array<size_t, nfieldsFixed> get_offset(size_t data_size) {
            std::array<size_t, nfieldsFixed> ext = get_extent(data_size);
            std::array<size_t, nfieldsFixed> off = {0};
            std::partial_sum(ext.begin(), ext.end(), off.begin() + 1, std::plus<size_t>());
            return off;
        }
        static size_t get_entry_size(size_t data_size) {
            auto ext = get_extent(data_size);
            return std::accumulate(ext.begin(), ext.end(), 0ul);
        }
    };
    struct metaHasher {
        auto operator()(const meta &m) const { return std::hash<std::string>{}(fmt::format("{}|{}|{}", m.type_size, m.data_size, m.fieldname)); }
    };

    public:
    static inline std::unordered_map<meta, h5pp::hid::h5t, metaHasher> h5_types;

    static std::vector<std::byte> make_entry(uint64_t iter, uint64_t step, int64_t position, StorageEvent storage_event, int64_t bond_lim, const T *const data,
                                             size_t data_size) {
        auto                   extent     = meta::get_extent(data_size);
        auto                   offset     = meta::get_offset(data_size);
        size_t                 entry_size = meta::get_entry_size(data_size);
        std::vector<std::byte> entry(entry_size);
        std::memcpy(entry.data() + offset[0], &iter, extent[0]);
        std::memcpy(entry.data() + offset[1], &step, extent[1]);
        std::memcpy(entry.data() + offset[2], &position, extent[2]);
        std::memcpy(entry.data() + offset[3], &storage_event, extent[3]);
        std::memcpy(entry.data() + offset[4], &bond_lim, extent[4]);
        std::memcpy(entry.data() + offset[5], data, extent[5]);
        return entry;
    }

    [[nodiscard]] static h5pp::hid::h5t &register_table_type(const h5pp::hid::h5t &h5elem_t, size_t data_size, std::string_view fieldname) {
        meta m = {sizeof(T), data_size, std::string(fieldname)};
        if(h5_types.find(m) == h5_types.end()) {
            h5_enum_storage_event::create();
            auto           offset     = meta::get_offset(data_size);
            size_t         entry_size = meta::get_entry_size(data_size);
            h5pp::hid::h5t h5comp_t   = H5Tcreate(H5T_COMPOUND, entry_size);
            H5Tinsert(h5comp_t, "iter", offset[0], H5T_NATIVE_UINT64);
            H5Tinsert(h5comp_t, "step", offset[1], H5T_NATIVE_UINT64);
            H5Tinsert(h5comp_t, "position", offset[2], H5T_NATIVE_INT64);
            H5Tinsert(h5comp_t, "event", offset[3], h5_enum_storage_event::get_h5t());
            H5Tinsert(h5comp_t, "bond_lim", offset[4], H5T_NATIVE_INT64);

            h5pp::hid::h5t h5data_t;
            if constexpr(std::is_same_v<T, hvl_t>)
                h5data_t = h5pp::hid::h5t(H5Tvlen_create(h5elem_t));
            else if constexpr(h5pp::type::sfinae::is_std_array_v<T>) {
                auto dims = std::array<hsize_t, 1>{std::tuple_size<T>::value}; // Gets N in std::array<S,N>
                h5data_t  = H5Tarray_create(h5elem_t, dims.size(), dims.data());
            } else
                h5data_t = h5elem_t;
            if(data_size == 1)
                H5Tinsert(h5comp_t, std::string(fieldname).c_str(), offset.back(), h5data_t);
            else
                for(size_t elem = 0; elem < data_size; elem++) {
                    H5Tinsert(h5comp_t, fmt::format("{}{}", fieldname, elem).c_str(), offset.back() + elem * sizeof(T), h5data_t);
                }
            h5_types[m] = h5comp_t;
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
        item(const item &it);
        item(std::string_view name_, double time, double sum, double pcnt, double avg, int level, size_t count);
        ~item() noexcept; // Automatically deallocate when no longer needed

        private:
        void copy_name(std::string_view name_);
    };

    h5pp_ur();
    static void register_table_type();
};