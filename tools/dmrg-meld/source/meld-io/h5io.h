#pragma once
#include "meld-io/h5db.h"
#include "meld-io/id.h"
#include "meld-io/meta.h"
#include <deque>
#include <h5pp/details/h5ppInfo.h>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>
namespace h5pp {
    class File;
}

namespace tools::h5io {

    inline std::string tmp_path;
    inline std::string tgt_path;

    std::string get_tmp_dirname(std::string_view exename);

    template<typename T>
    std::string get_standardized_base(const ModelId<T> &H, int decimals = 4);

    std::vector<std::string> findKeys(const h5pp::File &h5_src, const std::string &root, const std::vector<std::string> &expectedKeys, long hits = -1,
                                      long depth = 0, bool usecache = true);

    template<typename T>
    std::vector<ModelKey> loadModel(const h5pp::File &h5_src, std::unordered_map<std::string, ModelId<T>> &srcModelDb, const std::vector<ModelKey> &srcKeys);

    template<typename T>
    void saveModel(const h5pp::File &h5_src, h5pp::File &h5_tgt, std::unordered_map<std::string, InfoId<h5pp::TableInfo>> &tgtTableDb,
                   const ModelId<T> &modelId, const FileId &fileId);

    std::vector<DsetKey>  gatherDsetKeys(const h5pp::File &h5_src, std::unordered_map<std::string, h5pp::DsetInfo> &srcDsetDb, const PathId &pathid,
                                         const std::vector<DsetKey> &srcKeys);
    std::vector<TableKey> gatherTableKeys(const h5pp::File &h5_src, std::unordered_map<std::string, h5pp::TableInfo> &srcTableDb, const PathId &pathid,
                                          const std::vector<TableKey> &tables);

    template<typename ModelType>
    void merge(h5pp::File &h5_tgt, const h5pp::File &h5_src, const FileId &fileId, const tools::h5db::Keys &keys, tools::h5db::TgtDb &tgtdb);

    void writeProfiling(h5pp::File &h5_tgt);

}