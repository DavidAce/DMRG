#pragma once
#include "general/enums.h"
#include "meld-io/id.h"
#include "meld-io/meta.h"
#include <h5pp/details/h5ppFilesystem.h>
#include <h5pp/details/h5ppInfo.h>
#include <string>
#include <unordered_map>
#include <vector>
namespace h5pp {
    class File;
}

namespace tools::h5db {

    struct Keys {
        std::vector<DsetKey>  dsets;
        std::vector<TableKey> tables;
        std::vector<CronoKey> cronos;
        std::vector<ScaleKey> scales;
        std::vector<BonddKey> bondds;
        std::vector<ModelKey> models;

        [[nodiscard]] std::vector<std::string> get_algos() const {
            std::vector<std::string> v;
            for(auto &d : dsets) v.emplace_back(d.algo);
            for(auto &d : tables) v.emplace_back(d.algo);
            for(auto &d : cronos) v.emplace_back(d.algo);
            for(auto &d : scales) v.emplace_back(d.algo);
            for(auto &d : bondds) v.emplace_back(d.algo);
            for(auto &d : models) v.emplace_back(d.algo);
            std::sort(v.begin(), v.end());
            v.erase(std::unique(v.begin(), v.end()), v.end());
            return v;
        }
        [[nodiscard]] std::vector<std::string> get_states() const {
            std::vector<std::string> v;
            for(auto &d : dsets) v.emplace_back(d.state);
            for(auto &d : tables) v.emplace_back(d.state);
            for(auto &d : cronos) v.emplace_back(d.state);
            for(auto &d : scales) v.emplace_back(d.state);
            for(auto &d : bondds) v.emplace_back(d.state);
            std::sort(v.begin(), v.end());
            v.erase(std::unique(v.begin(), v.end()), v.end());
            return v;
        }
        [[nodiscard]] std::vector<std::string> get_points() const {
            std::vector<std::string> v;
            for(auto &d : dsets) v.emplace_back(d.point);
            for(auto &d : tables) v.emplace_back(d.point);
            for(auto &d : cronos) v.emplace_back(d.point);
            for(auto &d : scales) v.emplace_back(d.point);
            for(auto &d : bondds) v.emplace_back(d.point);
            std::sort(v.begin(), v.end());
            v.erase(std::unique(v.begin(), v.end()), v.end());
            return v;
        }
    };

    template<typename ModelType>
    struct SrcDb {
        h5pp::fs::path                                   parent_path;
        std::unordered_map<std::string, h5pp::DsetInfo>  dset;
        std::unordered_map<std::string, h5pp::TableInfo> table;
        std::unordered_map<std::string, h5pp::TableInfo> crono;
        std::unordered_map<std::string, h5pp::TableInfo> scale;
        std::unordered_map<std::string, h5pp::TableInfo> bondd;
        std::unordered_map<std::string, ModelType>       model;
        void                                             clear() {
            dset.clear();
            table.clear();
            crono.clear();
            scale.clear();
            bondd.clear();
            model.clear();
        }
    };

    struct TgtDb {
        std::unordered_map<std::string, FileId>                    file;
        std::unordered_map<std::string, InfoId<h5pp::DsetInfo>>    dset;
        std::unordered_map<std::string, InfoId<h5pp::TableInfo>>   table;
        std::unordered_map<std::string, InfoId<BufferedTableInfo>> crono;
        std::unordered_map<std::string, InfoId<BufferedTableInfo>> scale;
        std::unordered_map<std::string, InfoId<BufferedTableInfo>> bondd;
        std::unordered_map<std::string, InfoId<h5pp::TableInfo>>   model;
        void                                                       clear() {
            file.clear();
            dset.clear();
            table.clear();
            crono.clear();
            scale.clear();
            bondd.clear();
            model.clear();
        }
    };

    std::unordered_map<std::string, FileId> loadFileDatabase(const h5pp::File &h5_tgt);

    template<typename InfoType, typename KeyType>
    std::unordered_map<std::string, InfoId<InfoType>> loadDatabase(const h5pp::File &h5_tgt, const std::vector<KeyType> &keys);

    void saveDatabase(h5pp::File &h5_tgt, const std::unordered_map<std::string, FileId> &fileIdDb);

    template<typename InfoType>
    void saveDatabase(h5pp::File &h5_tgt, std::unordered_map<std::string, InfoId<InfoType>> &infoDb);

    FileIdStatus getFileIdStatus(const std::unordered_map<std::string, FileId> &fileIdDb, const FileId &newFileId);

}
