#pragma once
#include "config/enums.h"
#include "general/text.h"
#include <deque>
#include <h5pp/details/h5ppFormat.h>
#include <h5pp/details/h5ppHid.h>
#include <h5pp/details/h5ppInfo.h>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

struct BufferedTableInfo {
    private:
    struct ContiguousBuffer {
        hsize_t                offset = 0; // In table units
        hsize_t                extent = 0; // In table units
        std::vector<std::byte> rawdata;
    };
    hsize_t count = 0;

    public:
    h5pp::TableInfo              *info = nullptr;
    std::vector<ContiguousBuffer> recordBuffer;
    size_t                        maxRecords = 10000;
    BufferedTableInfo();
    BufferedTableInfo(h5pp::TableInfo *info_);
    BufferedTableInfo &operator=(h5pp::TableInfo *info_);
    ~BufferedTableInfo();

    void insert(const std::vector<std::byte> &entry, hsize_t index /* index in units of table entries */);
    void insert(const std::vector<std::byte>::const_iterator begin, const std::vector<std::byte>::const_iterator end, hsize_t index);

    hsize_t get_count(); /* Number of inserts */
    void    flush();
};

struct FileStats {
    size_t               files = 0;
    size_t               count = 0;
    uintmax_t            bytes = 0;
    double               elaps = 0.0;
    [[nodiscard]] double get_speed() const {
        if(elaps == 0)
            return 0.;
        else
            return static_cast<double>(count) / elaps;
    }
};

struct FileId {
    long   seed      = -1;
    char   path[512] = {};
    size_t hash      = {};
    FileId()         = default;
    FileId(long seed_, std::string_view path_, size_t hash_);
    [[nodiscard]] std::string string() const;
};

struct lbit {
    double                   J1_mean, J2_mean, J3_mean;
    double                   J1_wdth, J2_wdth, J3_wdth;
    size_t                   J2_span;
    double                   xi_Jcls;
    size_t                   u_depth;
    double                   u_fmix;
    double                   u_tstd, u_cstd;
    UnitaryGateWeight        u_g8w8;
    UnitaryGateType          u_type;
    long                     u_bond = -1;
    std::vector<std::string> fields = {"J1_mean", "J2_mean", "J3_mean", "J1_wdth", "J2_wdth", "J3_wdth", "J2_span", "xi_Jcls",
                                       "u_fmix",  "u_depth", "u_tstd",  "u_cstd",  "u_g8w8",  "u_type",  "u_bond"};
};

struct sdual {
    double                   J_mean;
    double                   J_wdth;
    double                   h_mean;
    double                   h_wdth;
    double                   lambda;
    double                   delta;
    std::vector<std::string> fields = {"J_mean", "J_wdth", "h_mean", "h_wdth", "lambda", "delta"};
};

struct majorana {
    double                   g;
    double                   delta;
    std::vector<std::string> fields = {"g", "delta"};
};

template<typename Param>
struct ModelId {
    Param       p;
    size_t      model_size;
    std::string model_type;
    std::string distribution;
    std::string algorithm;
    std::string key;
    std::string path;
    std::string basepath;
};

struct PathId {
    public:
    std::string src_path, tgt_path;
    std::string base, algo, state;
    PathId(std::string_view base_, std::string_view algo_, std::string_view state_);
    [[nodiscard]] bool        match(std::string_view algo_pattern, std::string_view state_pattern) const;
    [[nodiscard]] std::string dset_path(std::string_view dsetname) const;
    [[nodiscard]] std::string table_path(std::string_view tablename) const;
    [[nodiscard]] std::string crono_path(std::string_view tablename, size_t iter) const;
    [[nodiscard]] std::string scale_path(std::string_view tablename, size_t bond) const;
    [[nodiscard]] std::string fesle_path(std::string_view tablename, size_t bond) const;
    template<typename KeyT>
    [[nodiscard]] std::string create_path(std::string_view tablename, size_t idx) const;
};

template<typename InfoType>
struct InfoId {
    private:
    bool                              modified = false;
    std::unordered_map<long, hsize_t> db;

    public:
    InfoType info = InfoType();
    InfoId()      = default;
    InfoId(long seed_, hsize_t index_);
    InfoId(const InfoType &info_);
    bool    db_modified() const { return modified; }
    bool    has_index(long seed) const { return db.find(seed) != db.end(); }
    hsize_t get_index(long seed) const {
        auto res = db.find(seed);
        if(res != db.end())
            return res->second;
        else
            return db.size();
    }
    void insert(long seed, hsize_t index) {
        auto res = db.insert({seed, index});
        if(res.second) modified = true;
    }
    [[nodiscard]] const std::unordered_map<long, hsize_t> &get_db() const { return db; }
};

template<>
struct InfoId<BufferedTableInfo> {
    private:
    bool                              modified = false;
    std::unordered_map<long, hsize_t> db;

    public:
    h5pp::TableInfo   info = h5pp::TableInfo();
    BufferedTableInfo buff = BufferedTableInfo();
    InfoId()               = default;
    InfoId(long seed_, hsize_t index_);
    InfoId(const h5pp::TableInfo &info_);
    InfoId &operator=(const h5pp::TableInfo &info_);
    bool    db_modified() const { return modified; }
    bool    has_index(long seed) const { return db.find(seed) != db.end(); }
    hsize_t get_index(long seed) const {
        auto res = db.find(seed);
        if(res != db.end())
            return res->second;
        else
            return db.size();
    }
    void insert(long seed, hsize_t index) {
        auto res = db.insert({seed, index});
        if(res.second) modified = true;
    }

    [[nodiscard]] const std::unordered_map<long, hsize_t> &get_db() const { return db; }
};

struct SeedId {
    long    seed  = -1;
    hsize_t index = std::numeric_limits<hsize_t>::max();
    SeedId()      = default;
    SeedId(long seed_, hsize_t index_) : seed(seed_), index(index_) {}
    [[nodiscard]] std::string string() const { return h5pp::format("seed {} | index {}", index, seed); }
};

class H5T_FileId {
    public:
    static inline h5pp::hid::h5t h5_type;
    H5T_FileId();
    static void register_table_type();
};

class H5T_SeedId {
    public:
    static inline h5pp::hid::h5t h5_type;
    H5T_SeedId();
    static void register_table_type();
};

class H5T_profiling {
    public:
    static inline h5pp::hid::h5t h5_type;

    struct item {
        double time;
        double avg;
        size_t count;
    };

    H5T_profiling();
    static void register_table_type();
};