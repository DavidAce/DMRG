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
#include <regex>
#include <set>
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

    std::string get_tmp_dirname(std::string_view exename) { return fmt::format("{}.{}", h5pp::fs::path(exename).filename().string(), getenv("USER")); }

    template<typename T>
    std::string get_standardized_base(const ModelId<T> &H, int decimals) {
        auto t_scope = tid::tic_scope(__FUNCTION__);
        if constexpr(std::is_same_v<T, sdual>) return h5pp::format("L_{1}/l_{2:.{0}f}/d_{3:+.{0}f}", decimals, H.model_size, H.p.lambda, H.p.delta);
        if constexpr(std::is_same_v<T, majorana>) return h5pp::format("L_{1}/g_{2:.{0}f}/d_{3:+.{0}f}", decimals, H.model_size, H.p.g, H.p.delta);
        if constexpr(std::is_same_v<T, lbit>) {
            std::string J_mean_str = fmt::format("J[{1:+.{0}f}_{2:+.{0}f}_{3:+.{0}f}]", decimals, H.p.J1_mean, H.p.J2_mean, H.p.J3_mean);
            std::string J_wdth_str = fmt::format("w[{1:+.{0}f}_{2:+.{0}f}_{3:+.{0}f}]", decimals, H.p.J1_wdth, H.p.J2_wdth, H.p.J3_wdth);
            std::string x_str      = fmt::format("x_{1:.{0}f}", decimals, H.p.J2_xcls);
            std::string f_str      = fmt::format("f_{1:.{0}f}", decimals, H.p.f_mixer);
            std::string u_str      = fmt::format("u_{}", H.p.u_layer);
            std::string r_str;
            // J2_span is special since it can be -1ul. We prefer putting -1 in the path rather than 18446744073709551615
            if(H.p.J2_span == -1ul)
                r_str = fmt::format("r_L");
            else
                r_str = fmt::format("r_{}", H.p.J2_span);
            auto base = h5pp::format("L_{}/{}/{}/{}/{}/{}/{}", H.model_size, J_mean_str, J_wdth_str, x_str, f_str, u_str, r_str);
            tools::logger::log->info("creating base with {} decimals: {}", decimals, base);
            return base;
        }
    }
    template std::string get_standardized_base(const ModelId<sdual> &H, int decimals);
    template std::string get_standardized_base(const ModelId<majorana> &H, int decimals);
    template std::string get_standardized_base(const ModelId<lbit> &H, int decimals);

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
                if(key.empty())
                    found.emplace_back(key);
                else if(key.back() == '*') {
                    std::string_view key_match = std::string_view(key).substr(0, key.size() - 2);
                    tools::logger::log->info("key_match: ", key_match);
                    for(auto &item : h5_src.findGroups(key_match, root, hits, depth)) {
                        if(not text::startsWith(item, key_match)) continue;
                        if(found.size() > 1 and std::find(found.begin(), found.end(), item) != found.end()) continue;
                        found.emplace_back(item);
                    }
                } else {
                    for(auto &item : h5_src.findGroups(key, root, hits, depth)) {
                        if(not text::endsWith(item, key)) continue;
                        if(found.size() > 1 and std::find(found.begin(), found.end(), item) != found.end()) continue;
                        found.emplace_back(item);
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
            tools::logger::log->trace("Search: key [{}] | result {}{}", key, result, cacheMsg);
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
                    hamiltonian.J2_xcls     = h5tb_hamiltonian.J2_xcls;
                    hamiltonian.J2_span     = h5tb_hamiltonian.J2_span;
                    hamiltonian.f_mixer     = h5tb_hamiltonian.f_mixer;
                    hamiltonian.u_layer     = h5tb_hamiltonian.u_layer;
                    srcModelId.distribution = h5tb_hamiltonian.distribution;
                }
                srcModelId.model_size = h5_src.readAttribute<size_t>(path, "model_size");
                srcModelId.model_type = h5_src.readAttribute<std::string>(path, "model_name");
                srcModelId.algorithm  = srcKey.algo;
                srcModelId.key        = key;
                srcModelId.path       = path;
                srcModelId.basepath   = get_standardized_base(srcModelId);
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
                h5_tgt.writeDataset(modelId.p.J2_xcls, fmt::format("{}/J2_xcls", modelPath));
                h5_tgt.writeDataset(modelId.p.J2_span, fmt::format("{}/J2_span", modelPath));
                h5_tgt.writeDataset(modelId.p.f_mixer, fmt::format("{}/f_mixer", modelPath));
                h5_tgt.writeDataset(modelId.p.u_layer, fmt::format("{}/u_layer", modelPath));
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
            if(not pathid.match(srcKey.algo, srcKey.state, srcKey.point)) continue;

            // A "srcKey" is a struct describing a dataset that we want, with name, type, size, etc
            // Here we look for a matching dataset in the given group path, and generate a key for it.
            // The srcKeys that match an existing dataset are returned in "keys", with an additional parameter ".key" which
            // is an unique identifier to find the DsetInfo object in the map srcDsetDb
            // The database srcDsetDb  is used to keep track the DsetInfo structs for found datasets.
            auto path = fmt::format("{}/{}", pathid.src_path, srcKey.name);
            auto key  = fmt::format("{}|{}", srcParentPath, path);
            if(srcDsetDb.find(key) == srcDsetDb.end()) {
                srcDsetDb[key] = h5_src.getDatasetInfo(path);
                if(srcDsetDb[key].dsetExists.value()) tools::logger::log->debug("Detected new source dataset {}", key);
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

    std::vector<ScaleKey> gatherScaleKeys(const h5pp::File &h5_src, std::unordered_map<std::string, h5pp::TableInfo> &srcTableDb, const PathId &pathid,
                                          const std::vector<ScaleKey> &srcKeys) {
        auto                  t_scope = tid::tic_scope(__FUNCTION__);
        std::vector<ScaleKey> keys;
        h5pp::Options         options;
        std::string           srcParentPath = h5pp::fs::path(h5_src.getFilePath()).parent_path();
        for(const auto &srcKey : srcKeys) {
            // Make sure to only collect keys for objects that we have asked for.
            if(not pathid.match(srcKey.algo, srcKey.state, srcKey.point)) continue;

            // A "srcKey" is a struct describing a table that we want, with name, type, size, etc
            // Here we look for a matching table in the given group path, and generate a key for it.
            // The srcKeys that match an existing table are returned in "keys", with an additional parameter ".key" which
            // is an unique identifier to find the TableInfo object in the map srcTableDb
            // The database srcTableDb  is used to keep track the TableInfo structs for found datasets.

            // scaleKeys should be a vector [ "bond_8", "bond_16", "bond_24" ...] with tables for each scaling measurement
            // Find the largest existing dim dimension. We need this to fill missing keys
            std::string              scaleKeyMax;
            std::vector<std::string> scaleKeys;

            if(srcKey.allKeys.empty()) {
                scaleKeys   = findKeys(h5_src, pathid.src_path, {srcKey.scale}, -1, 1);
                scaleKeyMax = scaleKeys.back();
            } else {
                auto foundKeys = findKeys(h5_src, pathid.src_path, {srcKey.scale}, -1, 1, false);
                if(not foundKeys.empty()) scaleKeyMax = foundKeys.back();
                scaleKeys = std::vector<std::string>(srcKey.allKeys.begin(), srcKey.allKeys.end());
            }

            for(const auto &scaleKey : scaleKeys) {
                auto path = fmt::format("{}/{}/{}", pathid.src_path, scaleKey, srcKey.name);
                auto key  = fmt::format("{}|{}", srcParentPath, path);
                auto dim  = tools::parse::extract_parameter_from_path<size_t>(path, "bond_");
                if(srcTableDb.find(key) == srcTableDb.end()) {
                    srcTableDb[key] = h5_src.getTableInfo(path);
                    if(srcTableDb[key].tableExists.value()) tools::logger::log->debug("Detected new source scale {}", key);
                } else {
                    auto t_read = tid::tic_scope("readTableInfo");

                    auto &srcInfo = srcTableDb[key];
                    // We reuse the struct srcDsetDb[dsetKey] for every source file,
                    // but each time have to renew the following fields

                    srcInfo.h5File      = h5_src.openFileHandle();
                    srcInfo.h5Dset      = std::nullopt;
                    srcInfo.numRecords  = std::nullopt;
                    srcInfo.tableExists = std::nullopt;
                    srcInfo.tablePath   = path;
                    h5pp::scan::readTableInfo(srcInfo, srcInfo.h5File.value(), options, h5_src.plists);
                }

                auto &srcInfo = srcTableDb[key];
                if(srcInfo.tableExists and srcInfo.tableExists.value()) {
                    keys.emplace_back(srcKey);
                    keys.back().key = key;
                    keys.back().dim = dim; // Here we save the dim value to use later for generating a target path
                } else {
                    tools::logger::log->debug("Missing scale [{}] in file [{}].", path, h5_src.getFilePath());
                    if(not scaleKeyMax.empty()) {
                        tools::logger::log->debug("Extrapolating with [{}]", scaleKeyMax);
                        auto pathMax = fmt::format("{}/{}/{}", pathid.src_path, scaleKeyMax, srcKey.name);
                        auto keyMax  = fmt::format("{}|{}", srcParentPath, pathMax);
                        //                        auto dimMax     = tools::parse::extract_parameter_from_path<size_t>(pathMax, "bond_");
                        srcTableDb[key] = srcTableDb[keyMax]; // This now points to the largest bond_# group that actually exists in the source file.

                        // Here we mock a larger dim dimension by reusing the largest one
                        keys.emplace_back(srcKey);
                        keys.back().key = keyMax;
                        keys.back().dim = dim; // Here we save the missing dim dimension to use later for generating a target path
                    }
                }
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
            if(not pathid.match(srcKey.algo, srcKey.state, srcKey.point)) continue;

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
                auto point_groups = tools::h5io::findKeys(h5_src, fmt::format("{}/{}", algo, state), keys.get_points(), -1, 1);
                for(const auto &point : point_groups) {
                    auto pathid = PathId(tgt_base, algo, state, point);
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
                        auto t_scale   = tid::tic_scope("scale");
                        auto scaleKeys = tools::h5io::gatherScaleKeys(h5_src, srcdb.scale, pathid, keys.scales);
                        tools::h5xf::transferScales(h5_tgt, tgtdb.scale, srcdb.scale, pathid, scaleKeys, fileId);
                    } catch(const std::runtime_error &ex) { tools::logger::log->error("Scale transfer failed in[{}]: {}", pathid.src_path, ex.what()); }
                    try {
                        auto t_bondd   = tid::tic_scope("bondd");
                        auto bonddKeys = tools::h5io::gatherTableKeys<StrictTableSize::TRUE>(h5_src, srcdb.bondd, pathid, keys.bondds);
                        tools::h5xf::transferSeries(h5_tgt, tgtdb.bondd, srcdb.bondd, pathid, bonddKeys, fileId);
                    } catch(const std::runtime_error &ex) { tools::logger::log->error("Bondd transfer failed in[{}]: {}", pathid.src_path, ex.what()); }
                    try {
                        auto t_crono   = tid::tic_scope("crono");
                        auto cronoKeys = tools::h5io::gatherTableKeys<StrictTableSize::TRUE>(h5_src, srcdb.crono, pathid, keys.cronos);
                        tools::h5xf::transferSeries(h5_tgt, tgtdb.crono, srcdb.crono, pathid, cronoKeys, fileId);
                    } catch(const std::runtime_error &ex) { tools::logger::log->error("Crono transfer failed in[{}]: {}", pathid.src_path, ex.what()); }
                }
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