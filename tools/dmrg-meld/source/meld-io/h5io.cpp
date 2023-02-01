#include "h5io.h"
#include "debug/exceptions.h"
#include "general/prof.h"
#include "general/text.h"
#include "h5xf.h"
#include "meld-io/h5dbg.h"
#include "meld-io/id.h"
#include "meld-io/logger.h"
#include "meld-io/meta.h"
#include "meld-io/parse.h"
#include "meld-io/type.h"
#include "mpi/mpi-tools.h"
#include "tid/tid.h"
#include <cstdlib>
#include <functional>
#include <h5pp/h5pp.h>
#include <h5tb.h>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace tools::h5io {

    namespace internal {
        struct SearchResult {
            std::string                      root;
            std::string                      key;
            long                             hits;
            long                             depth;
            mutable std::vector<std::string> result;
            bool                             operator==(const SearchResult &rhs) const {
                return this->root == rhs.root and this->key == rhs.key and this->hits == rhs.hits and this->depth == rhs.depth;
            }
        };

        struct SearchHasher {
            auto operator()(const SearchResult &r) const { return std::hash<std::string>{}(fmt::format("{}|{}|{}|{}", r.root, r.key, r.hits, r.depth)); }
        };
    }

    std::string get_tmp_dirname(std::string_view exename) { return fmt::format("{}.{}", h5pp::fs::path(exename).filename(), getenv("USER")); }

    size_t get_num_decimals(double val) {
        auto s       = fmt::format("{}", val);
        auto pos_dot = s.find('.');
        if(pos_dot == std::string::npos) return 0;
        return s.substr(pos_dot + 1).size();
    }
    size_t get_max_decimals(std::initializer_list<double> l) {
        if(l.size() == 0) throw std::runtime_error("l is empty");
        std::vector<size_t> decs(l.size());
        std::transform(l.begin(), l.end(), decs.begin(), get_num_decimals);
        return *std::max_element(decs.begin(), decs.end());
    }

    template<typename T>
    std::string get_standardized_base(const ModelId<T> &H) {
        auto t_scope = tid::tic_scope(__FUNCTION__);
        if constexpr(std::is_same_v<T, sdual>) {
            size_t decimals = 4;
            return h5pp::format("L{1}/l{2:.{0}f}/d{3:+.{0}f}", decimals, H.model_size, H.p.lambda, H.p.delta);
        }
        if constexpr(std::is_same_v<T, majorana>) {
            size_t decimals = 4;
            return h5pp::format("L{1}/g{2:.{0}f}/d{3:+.{0}f}", decimals, H.model_size, H.p.g, H.p.delta);
        }
        if constexpr(std::is_same_v<T, lbit>) {
            size_t J_dec = get_max_decimals({H.p.J1_mean, H.p.J1_wdth, H.p.J2_mean, H.p.J2_wdth, H.p.J3_mean, H.p.J3_wdth});
            size_t x_dec = get_max_decimals({H.p.xi_Jcls});
            size_t u_dec = get_max_decimals({H.p.u_fmix, H.p.u_tstd, H.p.u_cstd});
            auto   L_str = fmt::format(FMT_COMPILE("L{}"), H.model_size);
            auto   J_str = fmt::format(FMT_COMPILE("J[{1:+.{0}f}±{2:.{0}f}_{3:+.{0}f}±{4:.{0}f}_{5:+.{0}f}±{6:.{0}f}]"), J_dec, H.p.J1_mean, H.p.J1_wdth,
                                       H.p.J2_mean, H.p.J2_wdth, H.p.J3_mean, H.p.J3_wdth);
            auto   x_str = fmt::format(FMT_COMPILE("x{1:.{0}f}"), x_dec, H.p.xi_Jcls);
            // J2_span is special since it can be -1ul, meaning long range. We prefer putting L in the path rather than 18446744073709551615
            auto r_str = H.p.J2_span == -1ul ? fmt::format("rL") : fmt::format("r{}", H.p.J2_span);
            auto u_str = fmt::format(FMT_COMPILE("u[d{1}_f{2:.{0}f}_tw{3:.{0}f}{4:.2}_cw{5:.{0}f}{6:.2}]"), u_dec, H.p.u_depth, H.p.u_fmix, H.p.u_tstd,
                                     enum2sv(H.p.u_tgw8), H.p.u_cstd, enum2sv(H.p.u_cgw8));
            auto base  = fmt::format(FMT_COMPILE("{0}/{1}/{2}/{3}/{4}"), L_str, J_str, x_str, r_str, u_str);
            tools::logger::log->info("creating base: {}", base);
            return base;
        }
    }
    template std::string get_standardized_base(const ModelId<sdual> &H);
    template std::string get_standardized_base(const ModelId<majorana> &H);
    template std::string get_standardized_base(const ModelId<lbit> &H);

    std::vector<std::string> findKeys(const h5pp::File &h5_src, const std::string &root, const std::vector<std::string> &expectedKeys, long hits, long depth,
                                      bool usecache) {
        auto                                                                      t_scope = tid::tic_scope(__FUNCTION__);
        static std::unordered_set<internal::SearchResult, internal::SearchHasher> cache;
        std::vector<std::string>                                                  result;
        for(const auto &key : expectedKeys) {
            internal::SearchResult searchQuery = {root, key, hits, depth, {}};
            std::string            cacheMsg;
            auto                   cachedIt = cache.find(searchQuery);
            bool                   cacheHit = cachedIt != cache.end() and                         // The search query has been processed earlier
                            ((hits > 0 and static_cast<long>(cachedIt->result.size()) >= hits) or // Asked for a up to "hits" reslts
                             (hits <= 0 and static_cast<long>(cachedIt->result.size()) >= 1)      // Asked for any number of results
                            );
            if(cacheHit and usecache) {
                for(auto &item : cachedIt->result) {
                    if(result.size() > 1 and std::find(result.begin(), result.end(), item) != result.end()) continue;
                    result.emplace_back(item);
                }
                cacheMsg = " | cache hit";
            } else {
                std::vector<std::string> found;
                auto                     key_split = text::split(key, "/");
                if(key_split.size() > 1) {
                    // We try to find the last key given the current root at infinite depth. Then we keep
                    // groups if all the key_split elements are present in the path
                    auto key_last = key_split.back();
                    found         = h5_src.findGroups(key_last, root, -1, -1);
                    // Prune found items if they do not have the expected path depth
                    text::erase_if(found, [&key_split](auto &item) -> bool {
                        auto rel_depth = std::count(item.begin(), item.end(), '/');
                        return rel_depth + 1 != static_cast<long>(key_split.size());
                    });

                    // Prune found items if they do not have the same groups as key_split
                    text::erase_if(found, [&key_split](auto &item) -> bool {
                        // Check that all parts of key_split exist in item
                        for(const auto &s : key_split) {
                            auto fuzz_pos    = s.find_first_of('*', 0);
                            auto key_in_item = item.find(s.substr(0, fuzz_pos)) != std::string::npos;
                            if(not key_in_item) return true;
                        }
                        return false;
                    });

                } else {
                    if(key.empty())
                        found.emplace_back(key);
                    else if(key.back() == '*') {
                        std::string_view key_match = std::string_view(key).substr(0, key.size() - 2);
                        for(auto &item : h5_src.findGroups(key_match, root, hits, depth)) {
                            if(not text::startsWith(item, key_match)) continue;
                            if(found.size() > 1 and std::find(found.begin(), found.end(), item) != found.end()) continue;
                            found.emplace_back(item);
                        }
                    } else {
                        for(auto &item : h5_src.findGroups(key, root, hits, depth)) {
                            tools::logger::log->info("Found item [{}]", item);
                            if(not text::endsWith(item, key)) continue;
                            if(found.size() > 1 and std::find(found.begin(), found.end(), item) != found.end()) continue;
                            found.emplace_back(item);
                        }
                    }
                }
                auto [it, ok] = cache.insert(searchQuery);
                it->result    = found;
                // Now we have built a cache hit. Add it to results
                for(auto &item : it->result) {
                    if(result.size() > 1 and std::find(result.begin(), result.end(), item) != result.end()) continue;
                    result.emplace_back(item);
                }
            }
            // Sort
            std::sort(result.begin(), result.end(), text::natcomp);
            tools::logger::log->trace("Search: root {} | key [{}] | result {}{}", root, key, result, cacheMsg);
        }
        return result;
    }
    template<typename T>
    std::vector<ModelKey> loadModel(const h5pp::File &h5_src, std::unordered_map<std::string, ModelId<T>> &srcModelDb, const std::vector<ModelKey> &srcKeys) {
        auto                  t_scope = tid::tic_scope(__FUNCTION__);
        std::vector<ModelKey> keys;
        for(const auto &srcKey : srcKeys) {
            auto path = fmt::format("{}/{}/{}", srcKey.algo, srcKey.model, srcKey.name);
            auto key  = fmt::format("{}|{}", h5pp::fs::path(h5_src.getFilePath()).parent_path(), path);
            if(srcModelDb.find(key) == srcModelDb.end() and h5_src.linkExists(path)) {
                // Copy the model from the attributes in h5_src to a struct ModelId
                srcModelDb[key]   = ModelId<T>();
                auto &srcModelId  = srcModelDb[key];
                auto &hamiltonian = srcModelId.p;
                if constexpr(std::is_same_v<T, sdual>) {
                    auto h5tb_hamiltonian   = h5_src.readTableRecords<h5tb_ising_selfdual::table>(path, h5pp::TableSelection::FIRST);
                    hamiltonian.J_mean      = h5tb_hamiltonian.J_mean;
                    hamiltonian.J_wdth      = h5tb_hamiltonian.J_wdth;
                    hamiltonian.h_mean      = h5tb_hamiltonian.h_mean;
                    hamiltonian.h_wdth      = h5tb_hamiltonian.h_wdth;
                    hamiltonian.lambda      = h5tb_hamiltonian.lambda;
                    hamiltonian.delta       = h5tb_hamiltonian.delta;
                    srcModelId.distribution = h5tb_hamiltonian.distribution;
                }
                if constexpr(std::is_same_v<T, majorana>) {
                    auto h5tb_hamiltonian   = h5_src.readTableRecords<h5tb_ising_majorana::table>(path, h5pp::TableSelection::FIRST);
                    hamiltonian.g           = h5tb_hamiltonian.g;
                    hamiltonian.delta       = h5tb_hamiltonian.delta;
                    srcModelId.distribution = h5tb_hamiltonian.distribution;
                }
                if constexpr(std::is_same_v<T, lbit>) {
                    auto h5tb_hamiltonian   = h5_src.readTableRecords<h5tb_lbit::table>(path, h5pp::TableSelection::FIRST);
                    hamiltonian.J1_mean     = h5tb_hamiltonian.J1_mean;
                    hamiltonian.J2_mean     = h5tb_hamiltonian.J2_mean;
                    hamiltonian.J3_mean     = h5tb_hamiltonian.J3_mean;
                    hamiltonian.J1_wdth     = h5tb_hamiltonian.J1_wdth;
                    hamiltonian.J2_wdth     = h5tb_hamiltonian.J2_wdth;
                    hamiltonian.J3_wdth     = h5tb_hamiltonian.J3_wdth;
                    hamiltonian.xi_Jcls     = h5tb_hamiltonian.xi_Jcls;
                    hamiltonian.J2_span     = h5tb_hamiltonian.J2_span;
                    hamiltonian.u_depth     = h5tb_hamiltonian.u_depth;
                    hamiltonian.u_fmix      = h5tb_hamiltonian.u_fmix;
                    hamiltonian.u_tstd      = h5tb_hamiltonian.u_tstd;
                    hamiltonian.u_cstd      = h5tb_hamiltonian.u_cstd;
                    hamiltonian.u_tgw8      = h5tb_hamiltonian.u_tgw8;
                    hamiltonian.u_cgw8      = h5tb_hamiltonian.u_cgw8;
                    srcModelId.distribution = h5tb_hamiltonian.distribution;
                    tools::logger::log->info("{}: u_tgw8 {} u_cgw8 {}", path, enum2sv(hamiltonian.u_tgw8), enum2sv(hamiltonian.u_cgw8));
                }
                srcModelId.model_size = h5_src.readAttribute<size_t>(path, "model_size");
                srcModelId.model_type = h5_src.readAttribute<std::string>(path, "model_name");
                srcModelId.algorithm  = srcKey.algo;
                srcModelId.key        = key;
                srcModelId.path       = path;
                srcModelId.basepath   = get_standardized_base(srcModelId);
                tools::logger::log->info("srcModelId.basepath = {} | key {}", srcModelId.basepath, key);
            }
            keys.emplace_back(srcKey);
            keys.back().key = key;
        }
        return keys;
    }
    template std::vector<ModelKey> loadModel(const h5pp::File &h5_src, std::unordered_map<std::string, ModelId<sdual>> &srcModelDb,
                                             const std::vector<ModelKey> &srcKeys);
    template std::vector<ModelKey> loadModel(const h5pp::File &h5_src, std::unordered_map<std::string, ModelId<majorana>> &srcModelDb,
                                             const std::vector<ModelKey> &srcKeys);
    template std::vector<ModelKey> loadModel(const h5pp::File &h5_src, std::unordered_map<std::string, ModelId<lbit>> &srcModelDb,
                                             const std::vector<ModelKey> &srcKeys);

    template<typename T>
    void saveModel([[maybe_unused]] const h5pp::File &h5_src, h5pp::File &h5_tgt, std::unordered_map<std::string, InfoId<h5pp::TableInfo>> &tgtModelDb,
                   const ModelId<T> &modelId, const FileId &fileId) {
        auto              t_scope      = tid::tic_scope(__FUNCTION__);
        const std::string tgtModelPath = fmt::format("{}/{}", modelId.basepath, modelId.path);

        if(tgtModelDb.find(tgtModelPath) == tgtModelDb.end()) {
            tools::logger::log->debug("Copying model to tgtPath {}", tgtModelPath);
            tools::logger::log->trace("modelId key    {}", modelId.key);
            tools::logger::log->trace("modelId path   {}", modelId.path);
            tools::logger::log->trace("modelId bpath  {}", modelId.basepath);
            tools::logger::log->trace("modelId fields {}", modelId.p.fields);
            tgtModelDb[tgtModelPath] = h5_tgt.getTableInfo(tgtModelPath);

            auto &tgtId   = tgtModelDb[tgtModelPath];
            auto &tgtInfo = tgtModelDb[tgtModelPath].info;

            // It doesn't make sense to copy the whole hamiltonian table here:
            // It is specific to a single realization, but here we collect the fields common to all realizations.
            const auto &srcModelInfo = h5_src.getTableInfo(modelId.path); // This will return a TableInfo to a an existing table in the source file
            auto        modelPath    = fmt::format("{}/{}/model", modelId.basepath, modelId.algorithm); // Generate a new path for the target group
            auto        tablePath    = fmt::format("{}/hamiltonian", modelPath);                        // Generate a new path for the target table

            // Update an entry of the hamiltonian table with the relevant fields
            tools::logger::log->trace("Copying model {}", modelId.basepath);
            auto h5t_model = h5pp::util::getFieldTypeId(srcModelInfo, modelId.p.fields); // Generate a h5t with the relevant fields
            auto modelData =
                h5_src.readTableField<std::vector<std::byte>>(srcModelInfo, h5t_model, h5pp::TableSelection::LAST); // Read those fields into a buffer
            tgtInfo = h5_tgt.createTable(h5t_model, tablePath, h5pp::format("{} Hamiltonian", modelId.algorithm));
            h5_tgt.writeTableRecords(modelData, tablePath);
            // Update the database
            tgtId.insert(fileId.seed, 0);

            // Now copy some helpful scalar datasets. This data is available in the attributes of the table above but this is also handy
            h5_tgt.writeDataset(modelId.model_size, fmt::format("{}/{}/model/model_size", modelId.basepath, modelId.algorithm));
            h5_tgt.writeDataset(modelId.model_type, fmt::format("{}/{}/model/model_type", modelId.basepath, modelId.algorithm));
            h5_tgt.writeDataset(modelId.distribution, fmt::format("{}/{}/model/distribution", modelId.basepath, modelId.algorithm));
            if constexpr(std::is_same_v<T, sdual>) {
                h5_tgt.writeDataset(modelId.p.J_mean, fmt::format("{}/{}/model/J_mean", modelId.basepath, modelId.algorithm));
                h5_tgt.writeDataset(modelId.p.J_wdth, fmt::format("{}/{}/model/J_wdth", modelId.basepath, modelId.algorithm));
                h5_tgt.writeDataset(modelId.p.h_mean, fmt::format("{}/{}/model/h_mean", modelId.basepath, modelId.algorithm));
                h5_tgt.writeDataset(modelId.p.h_wdth, fmt::format("{}/{}/model/h_wdth", modelId.basepath, modelId.algorithm));
                h5_tgt.writeDataset(modelId.p.lambda, fmt::format("{}/{}/model/lambda", modelId.basepath, modelId.algorithm));
                h5_tgt.writeDataset(modelId.p.delta, fmt::format("{}/{}/model/delta", modelId.basepath, modelId.algorithm));
            }
            if constexpr(std::is_same_v<T, majorana>) {
                h5_tgt.writeDataset(modelId.p.g, fmt::format("{}/{}/model/g", modelId.basepath, modelId.algorithm));
                h5_tgt.writeDataset(modelId.p.delta, fmt::format("{}/{}/model/delta", modelId.basepath, modelId.algorithm));
            }
            if constexpr(std::is_same_v<T, lbit>) {
                h5_tgt.writeDataset(modelId.p.J1_mean, fmt::format("{}/J1_mean", modelPath));
                h5_tgt.writeDataset(modelId.p.J2_mean, fmt::format("{}/J2_mean", modelPath));
                h5_tgt.writeDataset(modelId.p.J3_mean, fmt::format("{}/J3_mean", modelPath));
                h5_tgt.writeDataset(modelId.p.J1_wdth, fmt::format("{}/J1_wdth", modelPath));
                h5_tgt.writeDataset(modelId.p.J2_wdth, fmt::format("{}/J2_wdth", modelPath));
                h5_tgt.writeDataset(modelId.p.J3_wdth, fmt::format("{}/J3_wdth", modelPath));
                h5_tgt.writeDataset(modelId.p.J2_span, fmt::format("{}/J2_span", modelPath));
                h5_tgt.writeDataset(modelId.p.xi_Jcls, fmt::format("{}/xi_Jcls", modelPath));
                h5_tgt.writeDataset(modelId.p.u_depth, fmt::format("{}/u_depth", modelPath));
                h5_tgt.writeDataset(modelId.p.u_fmix, fmt::format("{}/u_fmix", modelPath));
                h5_tgt.writeDataset(modelId.p.u_tstd, fmt::format("{}/u_tstd", modelPath));
                h5_tgt.writeDataset(modelId.p.u_cstd, fmt::format("{}/u_cstd", modelPath));
                h5_tgt.writeDataset(modelId.p.u_tgw8, fmt::format("{}/u_tgw8", modelPath));
                h5_tgt.writeDataset(modelId.p.u_cgw8, fmt::format("{}/u_cgw8", modelPath));
            }
        }
    }
    template void saveModel(const h5pp::File &h5_src, h5pp::File &h5_tgt, std::unordered_map<std::string, InfoId<h5pp::TableInfo>> &tgtModelDb,
                            const ModelId<sdual> &modelId, const FileId &fileId);
    template void saveModel(const h5pp::File &h5_src, h5pp::File &h5_tgt, std::unordered_map<std::string, InfoId<h5pp::TableInfo>> &tgtModelDb,
                            const ModelId<majorana> &modelId, const FileId &fileId);
    template void saveModel(const h5pp::File &h5_src, h5pp::File &h5_tgt, std::unordered_map<std::string, InfoId<h5pp::TableInfo>> &tgtModelDb,
                            const ModelId<lbit> &modelId, const FileId &fileId);

    std::vector<DsetKey> gatherDsetKeys(const h5pp::File &h5_src, std::unordered_map<std::string, h5pp::DsetInfo> &srcDsetDb, const PathId &pathid,
                                        const std::vector<DsetKey> &srcKeys) {
        auto                 t_scope = tid::tic_scope(__FUNCTION__);
        std::vector<DsetKey> keys;

        std::string   srcParentPath = h5pp::fs::path(h5_src.getFilePath()).parent_path();
        h5pp::Options options;
        for(auto &srcKey : srcKeys) {
            // Make sure to only collect keys for objects that we have asked for.
            if(not pathid.match(srcKey.algo, srcKey.state)) {
                tools::logger::log->debug("Failed to match pathid {} to algo: [{}] and state: [{}]", pathid.src_path, srcKey.algo, srcKey.state);
                continue;
            }

            // A "srcKey" is a struct describing a dataset that we want, with name, type, size, etc
            // Here we look for a matching dataset in the given group path, and generate a key for it.
            // The srcKeys that match an existing dataset are returned in "keys", with an additional parameter ".key" which
            // is an unique identifier to find the DsetInfo object in the map srcDsetDb
            // The database srcDsetDb  is used to keep track the DsetInfo structs for found datasets.
            auto path = fmt::format("{}/{}", pathid.src_path, srcKey.name);
            auto key  = fmt::format("{}|{}", srcParentPath, path);
            if(srcDsetDb.find(key) == srcDsetDb.end()) {
                srcDsetDb[key] = h5_src.getDatasetInfo(path);
                if(srcDsetDb[key].dsetExists.value())
                    tools::logger::log->debug("Detected new source dataset {}", key);
                else
                    tools::logger::log->trace("Missing source dataset {}", key);
            } else {
                auto t_read = tid::tic_scope("readDsetInfo");

                auto &srcInfo = srcDsetDb[key];
                // We reuse the struct srcDsetDb[dsetKey] for every source file,
                // but each time have to renew the following fields
                srcInfo.h5File     = h5_src.openFileHandle();
                srcInfo.h5Dset     = std::nullopt;
                srcInfo.h5Space    = std::nullopt;
                srcInfo.dsetExists = std::nullopt;
                srcInfo.dsetSize   = std::nullopt;
                srcInfo.dsetDims   = std::nullopt;
                srcInfo.dsetByte   = std::nullopt;
                srcInfo.dsetPath   = path;
                h5pp::scan::readDsetInfo(srcInfo, srcInfo.h5File.value(), options, h5_src.plists);
            }
            auto &srcInfo = srcDsetDb[key];
            if(srcInfo.dsetExists and srcInfo.dsetExists.value()) {
                keys.emplace_back(srcKey);
                keys.back().key = key;
            } else {
                tools::logger::log->debug("Missing dataset [{}] in file [{}]", path, h5_src.getFilePath());
            }
        }
        return keys;
    }

    template<StrictTableSize strictTableSize, typename KeyT>
    std::vector<KeyT> gatherTableKeys(const h5pp::File &h5_src, std::unordered_map<std::string, h5pp::TableInfo> &srcTableDb, const PathId &pathid,
                                      const std::vector<KeyT> &srcKeys) {
        auto              t_scope = tid::tic_scope(__FUNCTION__);
        std::vector<KeyT> keys;
        h5pp::Options     options;
        std::string       srcParentPath = h5pp::fs::path(h5_src.getFilePath()).parent_path();
        for(auto &srcKey : srcKeys) {
            // Make sure to only collect keys for objects that we have asked for.
            if(not pathid.match(srcKey.algo, srcKey.state)) continue;

            // A "srcKey" is a struct describing a table that we want, with name, type, size, etc
            // Here we look for a matching table in the given group path, and generate a key for it.
            // The srcKeys that match an existing table are returned in "keys", with an additional parameter ".key" which
            // is an unique identifier to find the TableInfo object in the map srcTableDb
            // The database srcTableDb  is used to keep track the TableInfo structs for found datasets.
            auto                   path = fmt::format("{}/{}", pathid.src_path, srcKey.name);
            auto                   key  = fmt::format("{}|{}", srcParentPath, path);
            std::optional<hsize_t> numRecords_old;
            if(srcTableDb.find(key) == srcTableDb.end()) {
                srcTableDb[key] = h5_src.getTableInfo(path);
                if(srcTableDb[key].tableExists.value()) tools::logger::log->debug("Detected new source {}: {}", srcKey.classtag, key);
            } else {
                auto t_read = tid::tic_scope("readTableInfo");

                auto &srcInfo = srcTableDb[key];
                // We reuse the struct srcDsetDb[key] for every source file,
                // but each time have to renew the following fields
                numRecords_old      = srcInfo.numRecords;
                srcInfo.h5File      = h5_src.openFileHandle();
                srcInfo.h5Dset      = std::nullopt;
                srcInfo.numRecords  = std::nullopt;
                srcInfo.tableExists = std::nullopt;
                srcInfo.tablePath   = path;
                h5pp::scan::readTableInfo(srcInfo, srcInfo.h5File.value(), options, h5_src.plists);
            }

            auto &srcInfo = srcTableDb[key];
            if(srcInfo.tableExists and srcInfo.tableExists.value()) {
                // Check that the tables are the same size before and after
                if constexpr(strictTableSize == StrictTableSize::TRUE) {
                    if(srcKey.expected_size != -1ul and srcInfo.numRecords.value() != srcKey.expected_size)
                        throw except::runtime_error("Table size mismatch: num records {} | expected {}: {}", srcInfo.numRecords.value(), srcKey.expected_size,
                                                    key);
                    //                    if(numRecords_old and numRecords_old.value() != srcInfo.numRecords.value())
                    //                        throw except::runtime_error("Table size mismatch: num records {} | expected {}: {}", srcInfo.numRecords.value(),
                    //                        numRecords_old.value(), key);
                }

                keys.emplace_back(srcKey);
                keys.back().key = key;
            } else {
                srcInfo.h5File = std::nullopt;
                tools::logger::log->debug("Missing {} [{}] in file [{}]", srcKey.classtag, path, h5_src.getFilePath());
            }
        }
        return keys;
    }

    template<typename ModelType>
    void merge(h5pp::File &h5_tgt, const h5pp::File &h5_src, const FileId &fileId, const tools::h5db::Keys &keys, tools::h5db::TgtDb &tgtdb) {
        auto t_scope = tid::tic_scope(__FUNCTION__);
        // Define reusable source Info
        static tools::h5db::SrcDb<ModelId<ModelType>> srcdb;
        h5pp::fs::path                                parent_path = h5pp::fs::path(h5_src.getFilePath()).parent_path();
        if(srcdb.parent_path != parent_path) {
            // Clear when moving to another set of seeds (new point on the phase diagram)
            srcdb.clear();
            srcdb.parent_path = parent_path;
        }

        // Start finding the required components in the source
        auto groups = tools::h5io::findKeys(h5_src, "/", keys.get_algos(), -1, 0);
        for(const auto &algo : groups) {
            // Start by extracting the model
            auto modelKeys = tools::h5io::loadModel(h5_src, srcdb.model, keys.models);
            if(modelKeys.size() != 1) throw std::runtime_error("Exactly 1 model has to be loaded into keys");
            auto &modelId = srcdb.model[modelKeys.back().key];
            // Save the model to file if it hasn't
            tools::h5io::saveModel(h5_src, h5_tgt, tgtdb.model, modelId, fileId);
            auto tgt_base = modelId.basepath;
            // Next search for tables and datasets in the source file
            // and transfer them to the target file
            auto state_groups = tools::h5io::findKeys(h5_src, algo, keys.get_states(), -1, 0);
            for(const auto &state : state_groups) {
                auto pathid = PathId(tgt_base, algo, state);
                // Try gathering all the tables
                try {
                    auto t_dsets  = tid::tic_scope("dsets");
                    auto dsetKeys = tools::h5io::gatherDsetKeys(h5_src, srcdb.dset, pathid, keys.dsets);
                    tools::h5xf::transferDatasets(h5_tgt, tgtdb.dset, h5_src, srcdb.dset, pathid, dsetKeys, fileId);
                } catch(const std::runtime_error &ex) { tools::logger::log->warn("Dset transfer failed in [{}]: {}", pathid.src_path, ex.what()); }

                try {
                    auto t_table   = tid::tic_scope("table");
                    auto tableKeys = tools::h5io::gatherTableKeys<StrictTableSize::FALSE>(h5_src, srcdb.table, pathid, keys.tables);
                    tools::h5xf::transferTables(h5_tgt, tgtdb.table, srcdb.table, pathid, tableKeys, fileId);
                } catch(const std::runtime_error &ex) { tools::logger::log->error("Table transfer failed in [{}]: {}", pathid.src_path, ex.what()); }
                try {
                    auto t_fesup   = tid::tic_scope("fesup");
                    auto fesupKeys = tools::h5io::gatherTableKeys<StrictTableSize::FALSE>(h5_src, srcdb.fesup, pathid, keys.fesups);
                    tools::h5xf::transferSeries(h5_tgt, tgtdb.fesup, srcdb.fesup, pathid, fesupKeys, fileId);
                } catch(const std::runtime_error &ex) { tools::logger::log->error("Fesup transfer failed in[{}]: {}", pathid.src_path, ex.what()); }
                try {
                    auto t_fesdn   = tid::tic_scope("fesdn");
                    auto fesdnKeys = tools::h5io::gatherTableKeys<StrictTableSize::TRUE>(h5_src, srcdb.fesdn, pathid, keys.fesdns);
                    tools::h5xf::transferSeries(h5_tgt, tgtdb.fesdn, srcdb.fesdn, pathid, fesdnKeys, fileId);
                } catch(const std::runtime_error &ex) { tools::logger::log->error("Fesdn transfer failed in[{}]: {}", pathid.src_path, ex.what()); }
                try {
                    auto t_crono   = tid::tic_scope("crono");
                    auto cronoKeys = tools::h5io::gatherTableKeys<StrictTableSize::TRUE>(h5_src, srcdb.crono, pathid, keys.cronos);
                    tools::h5xf::transferSeries(h5_tgt, tgtdb.crono, srcdb.crono, pathid, cronoKeys, fileId);
                } catch(const std::runtime_error &ex) { tools::logger::log->error("Crono transfer failed in[{}]: {}", pathid.src_path, ex.what()); }
            }
        }

        auto t_close = tid::tic_scope("close");

        // Check that there are no errors hiding in the HDF5 error-stack
        auto num_errors = H5Eget_num(H5E_DEFAULT);
        if(num_errors > 0) {
            H5Eprint(H5E_DEFAULT, stderr);
            throw std::runtime_error(fmt::format("Error when treating file [{}]", h5_src.getFilePath()));
        }
    }
    template void merge<sdual>(h5pp::File &h5_tgt, const h5pp::File &h5_src, const FileId &fileId, const tools::h5db::Keys &keys, tools::h5db::TgtDb &tgtdb);
    template void merge<majorana>(h5pp::File &h5_tgt, const h5pp::File &h5_src, const FileId &fileId, const tools::h5db::Keys &keys, tools::h5db::TgtDb &tgtdb);
    template void merge<lbit>(h5pp::File &h5_tgt, const h5pp::File &h5_src, const FileId &fileId, const tools::h5db::Keys &keys, tools::h5db::TgtDb &tgtdb);

    void writeProfiling(h5pp::File &h5_tgt) {
        auto t_scope = tid::tic_scope(__FUNCTION__);
        H5T_profiling::register_table_type();
        for(const auto &t : tid::get_tree("", tid::level::normal)) {
            auto tablepath = h5pp::format(".db/prof_{}/{}", mpi::world.id, t->get_label());
            if(not h5_tgt.linkExists(tablepath)) h5_tgt.createTable(H5T_profiling::h5_type, tablepath, "H5MBL Profiling", {100});
            H5T_profiling::item entry{t->get_time(), t->get_time_avg(), t->get_tic_count()};
            h5_tgt.writeTableRecords(entry, tablepath, 0);
        }
    }
}