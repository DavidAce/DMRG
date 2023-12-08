#include "id.h"
#include "general/sfinae.h"
#include "logger.h"
#include "meta.h"
#include "tid/tid.h"
#include <h5pp/details/h5ppHdf5.h>
#include <h5pp/details/h5ppInfo.h>

BufferedTableInfo::BufferedTableInfo() = default;
BufferedTableInfo::BufferedTableInfo(h5pp::TableInfo *info_) : info(info_) {}
BufferedTableInfo &BufferedTableInfo::operator=(h5pp::TableInfo *info_) {
    // Guard self assignment
    if(info == info_) return *this;
    *this = BufferedTableInfo(info_);
    return *this;
}
BufferedTableInfo::~BufferedTableInfo() { flush(); }

void BufferedTableInfo::insert(const std::vector<std::byte>::const_iterator begin, const std::vector<std::byte>::const_iterator end, hsize_t index) {
    if(info == nullptr) throw std::runtime_error("insert: info is nullptr");
    if(std::distance(begin, end) != static_cast<long>(info->recordBytes.value())) throw std::runtime_error("insert: record and entry size mismatch");
    if(recordBuffer.size() >= maxRecords) flush();

    // We need to find out if there is already a contigous buffer where we can append this entry. If not, we start a new contiguous buffer
    for(auto &r : recordBuffer) {
        if(r.offset + r.extent == index) {
            r.rawdata.insert(r.rawdata.end(), begin, end);
            r.extent += 1;
            //            h5pp::print("Inserting 1 records at index {} into offset {} extent {}\n", index, r.offset, r.extent );
            count++;
            return;
        }
    }
    // None was found, so we make a new one
    recordBuffer.emplace_back(ContiguousBuffer{index, 1ul, std::vector<std::byte>(begin, end)});
    count++;
    //    auto & r = recordBuffer.back();
    //    h5pp::print("Inserting 1 records at index {} into offset {} extent {}\n",  index, r.offset, r.extent );
}

void BufferedTableInfo::insert(const std::vector<std::byte> &entry, hsize_t index) {
    if(info == nullptr) throw std::runtime_error("insert: info is nullptr");
    if(entry.size() != info->recordBytes.value()) throw std::runtime_error("insert: record and entry size mismatch");
    if(recordBuffer.size() >= maxRecords) flush();

    // We need to find out if there is already a contigous buffer where we can append this entry. If not, we start a new contiguous buffer
    for(auto &r : recordBuffer) {
        if(r.offset + r.extent == index) {
            r.rawdata.insert(r.rawdata.end(), entry.begin(), entry.end());
            r.extent += 1;
            //            h5pp::print("Inserting 1 records at index {} into offset {} extent {}\n", index, r.offset, r.extent );
            count++;
            return;
        }
    }
    // None was found, so we make a new one
    recordBuffer.emplace_back(ContiguousBuffer{index, 1ul, entry});
    count++;
    //    auto & r = recordBuffer.back();
    //    h5pp::print("Inserting 1 records at index {} into offset {} extent {}\n",  index, r.offset, r.extent );
}

void BufferedTableInfo::flush() {
    if(recordBuffer.empty()) return;
    if(info == nullptr) throw std::runtime_error("flush: info is nullptr");
    for(const auto &r : recordBuffer) { h5pp::hdf5::writeTableRecords(r.rawdata, *info, r.offset, r.extent); }

    recordBuffer.clear();
}

hsize_t BufferedTableInfo::get_count() { return count; }

FileId::FileId(long seed_, std::string_view path_, size_t hash_) : seed(seed_), hash(hash_) {
    strncpy(path, path_.data(), sizeof(path) - 1);
    path[sizeof(path) - 1] = '\0'; /* This line is important to make sure we don't get UB */
}
std::string FileId::string() const { return h5pp::format("path [{}] | seed {} | hash {}", path, seed, hash); }

template<typename InfoType>
InfoId<InfoType>::InfoId(long seed_, hsize_t index_) {
    db[seed_] = index_;
}
template<typename InfoType>
InfoId<InfoType>::InfoId(const InfoType &info_) : info(info_) {}
template struct InfoId<h5pp::DsetInfo>;
template struct InfoId<h5pp::TableInfo>;

InfoId<BufferedTableInfo>::InfoId(long seed_, hsize_t index_) { db[seed_] = index_; }
InfoId<BufferedTableInfo>::InfoId(const h5pp::TableInfo &info_) : info(info_), buff(BufferedTableInfo(&info)) {
    //    h5pp::print("Copy constructor for buffered table {}\n", info.tablePath.value());
}
InfoId<BufferedTableInfo> &InfoId<BufferedTableInfo>::operator=(const h5pp::TableInfo &info_) {
    //        h5pp::print("Assignment operator for buffered table {}\n", info_.tablePath.value());
    if(&info == &info_) return *this;
    info = info_;
    buff = BufferedTableInfo(&info);
    return *this;
}

PathId::PathId(std::string_view base_, std::string_view algo_, std::string_view state_) : base(base_), algo(algo_), state(state_) {
    src_path = fmt::format("{}/{}", algo, state);
    tgt_path = fmt::format("{}/{}/{}", base, algo, state);
}

bool PathId::match(std::string_view algo_pattern, std::string_view state_pattern) const {
    return text::match_pattern(algo, algo_pattern) and text::match_pattern(state, state_pattern);
}

[[nodiscard]] std::string PathId::dset_path(std::string_view dsetname) const { return h5pp::format("{}/{}/{}/dsets/{}", base, algo, state, dsetname); }
[[nodiscard]] std::string PathId::table_path(std::string_view tablename) const { return h5pp::format("{}/{}/{}/tables/{}", base, algo, state, tablename); }
[[nodiscard]] std::string PathId::crono_path(std::string_view tablename, size_t iter) const {
    /*
     * When collecting a "crono" kind of table:
     *      - the source path is <base>/<algo>/<state>/tables/<tablename>
     *      - we find entries for all iterations in <tablename>
     *      - we collect the contribution from each realization to each iteration separately
     *      - the target path <base>/<algo>/<state>/cronos/iter_<iter>/<tablename>, collects all the realizations
     */
    return h5pp::format("{}/{}/{}/cronos/iter_{}/{}", base, algo, state, iter, tablename);
}
[[nodiscard]] std::string PathId::scale_path(std::string_view tablename, size_t bond) const {
    /*
     * When collecting a "scale" kind of table:
     *      - the source path is <base>/<algo>/<state>/tablename>
     *      - we want all entries of event type BOND_STATE in each <tablename>
     *      - we collect the contribution from each realization to each bond_# separately
     *      - the target path <base>/<algo>/<state>/fes/bond_<#>/<tablename>, collects all the contributions
     */
    return h5pp::format("{}/{}/{}/fes-inc/bond_{}/{}", base, algo, state, bond, tablename);
}

[[nodiscard]] std::string PathId::fesle_path(std::string_view tablename, size_t bond) const {
    /*
     * When collecting a "FES" kind of table:
     *      - the source paths are at <base>/<algo>/<state>/<tablename>
     *      - we want all entries of event type FES_STATE in each <tablename>
     *      - we collect the contribution from each realization to each bond_# separately
     *      - the target path <base>/<algo>/<state>/fes/bond_<#>/<tablename>, collects all the contributions
     */
    return h5pp::format("{}/{}/{}/fes-dec/bond_{}/{}", base, algo, state, bond, tablename);
}

template<typename KeyT>
[[nodiscard]] std::string PathId::create_path(std::string_view tablename, size_t idx) const {
    static_assert(sfinae::is_any_v<KeyT, FesDnKey, FesUpKey, CronoKey>);
    if constexpr(std::is_same_v<KeyT, FesDnKey>) return fesle_path(tablename, idx);
    if constexpr(std::is_same_v<KeyT, FesUpKey>) return scale_path(tablename, idx);
    if constexpr(std::is_same_v<KeyT, CronoKey>) return crono_path(tablename, idx);
}
template std::string PathId::create_path<FesDnKey>(std::string_view tablename, size_t idx) const;
template std::string PathId::create_path<FesUpKey>(std::string_view tablename, size_t idx) const;
template std::string PathId::create_path<CronoKey>(std::string_view tablename, size_t idx) const;

H5T_FileId::H5T_FileId() { register_table_type(); }
void H5T_FileId::register_table_type() {
    if(h5_type.valid()) return;
    auto           t_scope  = tid::tic_scope(__FUNCTION__);
    h5pp::hid::h5t H5T_PATH = H5Tcopy(H5T_C_S1);
    H5Tset_size(H5T_PATH, 512);
    // Optionally set the null terminator '\0'
    H5Tset_strpad(H5T_PATH, H5T_STR_NULLTERM);
    h5_type = H5Tcreate(H5T_COMPOUND, sizeof(FileId));
    H5Tinsert(h5_type, "seed", HOFFSET(FileId, seed), H5T_NATIVE_LONG);
    H5Tinsert(h5_type, "path", HOFFSET(FileId, path), H5T_PATH);
    H5Tinsert(h5_type, "hash", HOFFSET(FileId, hash), H5T_NATIVE_ULONG);
}

H5T_SeedId::H5T_SeedId() { register_table_type(); }
void H5T_SeedId::register_table_type() {
    auto t_scope = tid::tic_scope(__FUNCTION__);
    if(h5_type.valid()) return;
    h5_type = H5Tcreate(H5T_COMPOUND, sizeof(SeedId));
    H5Tinsert(h5_type, "seed", HOFFSET(SeedId, seed), H5T_NATIVE_LONG);
    H5Tinsert(h5_type, "index", HOFFSET(SeedId, index), H5T_NATIVE_HSIZE);
}

H5T_profiling::H5T_profiling() { register_table_type(); }
void H5T_profiling::register_table_type() {
    if(h5_type.valid()) return;
    h5_type = H5Tcreate(H5T_COMPOUND, sizeof(item));
    H5Tinsert(h5_type, "time", HOFFSET(item, time), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5_type, "avg", HOFFSET(item, avg), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5_type, "count", HOFFSET(item, count), H5T_NATIVE_UINT64);
}