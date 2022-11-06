#include "h5xf.h"
#include "config/enums.h"
#include "debug/exceptions.h"
#include "meld-io/logger.h"
#include "tid/tid.h"
#include <h5pp/h5pp.h>
namespace tools::h5xf {

    namespace internal {
        void copy_dset(h5pp::File &h5_tgt, const h5pp::File &h5_src, h5pp::DsetInfo &tgtInfo, h5pp::DsetInfo &srcInfo, hsize_t index, size_t axis) {
            auto           t_scope = tid::tic_scope(__FUNCTION__);
            auto           data    = h5_src.readDataset<std::vector<std::byte>>(srcInfo); // Read the data into a generic buffer.
            h5pp::DataInfo dataInfo;

            // The axis parameter must be  < srcInfo.rank+1
            if(axis >= srcInfo.dsetDims->size()) {
                dataInfo.dataByte = srcInfo.dsetByte;
                dataInfo.dataSize = srcInfo.dsetSize;
                dataInfo.dataRank = tgtInfo.dsetRank; // E.g. if stacking rank 2 matrices then axis == 2 and so data rank must be 3.
                dataInfo.dataDims = std::vector<hsize_t>(axis + 1, 1ull);
                std::copy_n(srcInfo.dsetDims->begin(), srcInfo.dsetDims->size(), dataInfo.dataDims->begin());
                dataInfo.h5Space = h5pp::util::getMemSpace(dataInfo.dataSize.value(), dataInfo.dataDims.value());
            } else {
                dataInfo.dataDims = srcInfo.dsetDims;
                dataInfo.dataSize = srcInfo.dsetSize;
                dataInfo.dataByte = srcInfo.dsetByte;
                dataInfo.dataRank = srcInfo.dsetRank;
                dataInfo.h5Space  = srcInfo.h5Space;
            }
            tgtInfo.resizePolicy                   = h5pp::ResizePolicy::GROW;
            tgtInfo.dsetSlab                       = h5pp::Hyperslab();
            tgtInfo.dsetSlab->extent               = dataInfo.dataDims;
            tgtInfo.dsetSlab->offset               = std::vector<hsize_t>(tgtInfo.dsetDims->size(), 0);
            tgtInfo.dsetSlab->offset.value()[axis] = index;

            if(tgtInfo.dsetSlab->extent->size() != tgtInfo.dsetSlab->offset->size()) {
                throw except::logic_error("dsetSlab has mismatching ranks: \n "
                                          "target: offset {}:{} != extent {}:{} | dset {}\n"
                                          "source: dsetDimms {}:{} | dset {} \n"
                                          "axis: {}",
                                          tgtInfo.dsetSlab->offset.value(), tgtInfo.dsetSlab->offset->size(), tgtInfo.dsetSlab->extent.value(),
                                          tgtInfo.dsetSlab->extent->size(), tgtInfo.dsetPath.value(), srcInfo.dsetDims.value(), srcInfo.dsetDims->size(),
                                          srcInfo.dsetPath.value(), axis

                );
            }

            h5_tgt.appendToDataset(data, dataInfo, tgtInfo, static_cast<size_t>(axis));
            // Restore previous settings
            tgtInfo.dsetSlab = std::nullopt;
        }

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

    void transferDatasets(h5pp::File &h5_tgt, std::unordered_map<std::string, InfoId<h5pp::DsetInfo>> &tgtDsetDb, const h5pp::File &h5_src,
                          std::unordered_map<std::string, h5pp::DsetInfo> &srcDsetDb, const PathId &pathid, const std::vector<DsetKey> &srcDsetKeys,
                          const FileId &fileId) {
        auto t_scope = tid::tic_scope(__FUNCTION__);
        for(const auto &srcKey : srcDsetKeys) {
            if(srcDsetDb.find(srcKey.key) == srcDsetDb.end()) throw std::logic_error(h5pp::format("Key [{}] was not found in source map", srcKey.key));
            auto &srcInfo = srcDsetDb[srcKey.key];
            if(not srcInfo.dsetExists or not srcInfo.dsetExists.value()) continue;
            auto tgtName = h5pp::fs::path(srcInfo.dsetPath.value()).filename().string();
            auto tgtPath = h5pp::format("{}/{}", pathid.tgt_path, tgtName);
            if(tgtDsetDb.find(tgtPath) == tgtDsetDb.end()) {
                auto t_create = tid::tic_scope("createDataset");
                auto tgtDims  = srcInfo.dsetDims.value();
                if(tgtDims.empty()) tgtDims = {0}; // In case src is a scalar.
                while(tgtDims.size() < srcKey.axis + 1) tgtDims.push_back(1);
                tgtDims[srcKey.axis] = 0; // Create with 0 extent in the new axis direction, so that the dataset starts empty (zero volume)

                // Determine a good chunksize between 10 and 500 elements
                auto tgtChunk         = tgtDims;
                auto chunkSize        = std::clamp(5e5 / static_cast<double>(srcInfo.dsetByte.value()), 10., 1000.);
                tgtChunk[srcKey.axis] = static_cast<hsize_t>(chunkSize); // number of elements in the axis that we append into.

                if(srcKey.size == Size::VAR) {
                    auto        srcGroupPath    = h5pp::fs::path(srcInfo.dsetPath.value()).parent_path().string();
                    std::string statusTablePath = fmt::format("{}/status", srcGroupPath);
                    long        bond_max        = h5_src.readTableField<long>(statusTablePath, "bond_max", h5pp::TableSelection::FIRST);
                    tgtDims[0]                  = static_cast<hsize_t>(bond_max);
                }

                tools::logger::log->debug("Adding target dset {} | dims {} | chnk {}", tgtPath, tgtDims, tgtChunk);
                tgtDsetDb[tgtPath] = h5_tgt.createDataset(tgtPath, srcInfo.h5Type.value(), H5D_CHUNKED, tgtDims, tgtChunk);
            }
            auto &tgtId   = tgtDsetDb[tgtPath];
            auto &tgtInfo = tgtId.info;
            // Determine the target index where to copy this record
            hsize_t index = tgtId.get_index(fileId.seed); // Positive if it has been read already
            index         = index != std::numeric_limits<hsize_t>::max() ? index : tgtInfo.dsetDims.value().at(srcKey.axis);
            internal::copy_dset(h5_tgt, h5_src, tgtInfo, srcInfo, index, srcKey.axis);

            // Update the database
            tgtId.insert(fileId.seed, index);
        }
    }

    void transferTables(h5pp::File &h5_tgt, std::unordered_map<std::string, InfoId<h5pp::TableInfo>> &tgtTableDb,
                        std::unordered_map<std::string, h5pp::TableInfo> &srcTableDb, const PathId &pathid, const std::vector<TableKey> &srcTableKeys,
                        const FileId &fileId) {
        auto                      t_scope = tid::tic_scope(__FUNCTION__);
        std::vector<StorageEvent> srcEvents; // This stores the "event" row in the table. We need to make sure to only transfer StorageEvent::ITER_STATE
        for(const auto &srcKey : srcTableKeys) {
            if(srcTableDb.find(srcKey.key) == srcTableDb.end()) throw std::runtime_error(h5pp::format("Key [{}] was not found in source map", srcKey.key));
            auto &srcInfo = srcTableDb[srcKey.key];
            if(not srcInfo.tableExists or not srcInfo.tableExists.value()) continue;
            auto tgtName = h5pp::fs::path(srcInfo.tablePath.value()).filename().string();
            auto tgtPath = pathid.table_path(tgtName);
            if(tgtTableDb.find(tgtPath) == tgtTableDb.end()) {
                auto t_create = tid::tic_scope("createTable");
                tools::logger::log->debug("Adding target table {}", tgtPath);
                h5pp::TableInfo tableInfo = h5_tgt.getTableInfo(tgtPath);
                if(not tableInfo.tableExists.value())
                    tableInfo = h5_tgt.createTable(srcInfo.h5Type.value(), tgtPath, srcInfo.tableTitle.value(), std::nullopt, true);
                tgtTableDb[tgtPath] = tableInfo;
            }
            auto &tgtId   = tgtTableDb[tgtPath];
            auto &tgtInfo = tgtId.info;

            // Determine the target index where to copy this record
            hsize_t index = tgtId.get_index(fileId.seed); // Positive if it has been read already, numeric_limits::max if it's new
            tools::logger::log->debug("Transferring table {} -> {}", tgtPath, srcInfo.numRecords.value(), index);
            auto t_copy = tid::tic_scope("copyTableRecords");
            h5_tgt.copyTableRecords(srcInfo, tgtInfo, h5pp::TableSelection::LAST, index);
            // Update the database
            tgtId.insert(fileId.seed, index);
        }
    }

    template<typename KeyT>
    void transferSeries(h5pp::File &h5_tgt, std::unordered_map<std::string, InfoId<BufferedTableInfo>> &tgtTableDb,
                        std::unordered_map<std::string, h5pp::TableInfo> &srcTableDb, const PathId &pathid, const std::vector<KeyT> &srcKeys,
                        const FileId &fileId) {
        // In this function we take series data from each srcTable and create multiple tables tgtTable, one for each
        // index (e.g. iter, bond dim, etc). Each entry in tgtTable corresponds to the same index point on different realizations.
        auto                      t_scope = tid::tic_scope(__FUNCTION__);
        std::vector<std::byte>    srcReadBuffer;
        std::vector<size_t>       srcIndices; // We can assume all tables have the same indexing numbers. Only update on mismatch
        std::vector<StorageEvent> srcEvents;  // This stores the "event" row in the table. We need to make sure to only transfer StorageEvent::ITER_STATE
        for(const auto &srcKey : srcKeys) {
            if(srcTableDb.find(srcKey.key) == srcTableDb.end()) throw except::logic_error("Key [{}] was not found in source map", srcKey.key);
            auto &srcInfo = srcTableDb[srcKey.key];
            tools::logger::log->trace("Transferring {} table {} | records {}", srcKey.classtag, srcInfo.tablePath.value(), srcInfo.numRecords.value());

            // Iterate over all table elements. These should be a time series measured at every iteration
            // Note that there could in principle exist duplicate entries, which is why we can't trust the
            // "rec" iterator but have to get the iteration number from the table directly.
            // Try getting the iteration number, which is more accurate.

            try {
                if(srcIndices.size() != srcInfo.numRecords.value()) { // Update iteration numbers if it's not the same that we have already.
                    auto t_read     = tid::tic_scope("readTableField");
                    auto indexfield = text::match(srcKey.index, srcInfo.fieldNames.value());
                    auto eventfield = text::match({"event"}, srcInfo.fieldNames.value());
                    if(not indexfield)
                        throw std::runtime_error(fmt::format("Index field was not found in table [{}]\n Expected fields {}\n Existing fields {}",
                                                             srcInfo.tablePath.value(), srcKey.index, srcInfo.fieldNames.value()));
                    if(not eventfield)
                        throw std::runtime_error(
                            fmt::format("Event field was not found in table [{}]\n Existing fields {}", srcInfo.tablePath.value(), srcInfo.fieldNames.value()));
                    h5pp::hdf5::readTableField(srcIndices, srcInfo, {indexfield.value()});
                    h5pp::hdf5::readTableField(srcEvents, srcInfo, {eventfield.value()});
                }
            } catch(const std::exception &ex) { throw std::logic_error(fmt::format("Failed to get iteration numbers: {}", ex.what())); }
            for(size_t rec = 0; rec < srcInfo.numRecords.value(); rec++) {
                size_t srcIndex = rec;
                if(not srcEvents.empty() and srcEvents[rec] != StorageEvent::ITER_STATE) continue; // Only take cronos from iteration events
                if(not srcIndices.empty()) srcIndex = srcIndices[rec];                             // Get the actual index number from the table
                auto tgtName = h5pp::fs::path(srcInfo.tablePath.value()).filename().string();
                auto tgtPath = pathid.create_path<KeyT>(tgtName, srcIndex);

                if(tgtTableDb.find(tgtPath) == tgtTableDb.end()) {
                    auto t_create = tid::tic_scope("createTable");
                    tools::logger::log->debug("Adding target {} {}", srcKey.classtag, tgtPath);
                    h5pp::TableInfo tableInfo = h5_tgt.getTableInfo(tgtPath);
                    // Disabling compression is supposed to give a nice speedup. Read here:
                    // https://support.hdfgroup.org/HDF5/doc1.8/Advanced/DirectChunkWrite/UsingDirectChunkWrite.pdf
                    if(not tableInfo.tableExists.value())
                        tableInfo = h5_tgt.createTable(srcInfo.h5Type.value(), tgtPath, srcInfo.tableTitle.value(), std::nullopt, true);

                    tgtTableDb[tgtPath] = tableInfo;
                }
                auto &tgtId   = tgtTableDb[tgtPath];
                auto &tgtBuff = tgtId.buff;

                // Determine the target index where to copy this record
                // tgtInfo.numRecords is the total number of realizations registered until now
                // tgtId.get_index(fileId.seed) Gets the index number if it already exists, else
                // the size of the table (so we can append)
                hsize_t tgtIndex = tgtId.get_index(fileId.seed);

                // read a source record at "srcIndex" into the "tgtIndex" position in the buffer
                auto t_buffer = tid::tic_scope(fmt::format("bufferRecords-{}", srcKey.classtag));
                srcReadBuffer.resize(srcInfo.recordBytes.value());
                h5pp::hdf5::readTableRecords(srcReadBuffer, srcInfo, rec, 1);
                tgtBuff.insert(srcReadBuffer, tgtIndex);

                // Update the database
                tgtId.insert(fileId.seed, tgtIndex);
            }
        }
    }
    template void transferSeries(h5pp::File &h5_tgt, std::unordered_map<std::string, InfoId<BufferedTableInfo>> &tgtTableDb,
                                 std::unordered_map<std::string, h5pp::TableInfo> &srcTableDb, const PathId &pathid, const std::vector<CronoKey> &srcKeys,
                                 const FileId &fileId);
    template void transferSeries(h5pp::File &h5_tgt, std::unordered_map<std::string, InfoId<BufferedTableInfo>> &tgtTableDb,
                                 std::unordered_map<std::string, h5pp::TableInfo> &srcTableDb, const PathId &pathid, const std::vector<BonddKey> &srcKeys,
                                 const FileId &fileId);

    void transferScales(h5pp::File &h5_tgt, std::unordered_map<std::string, InfoId<BufferedTableInfo>> &tgtTableDb,
                        std::unordered_map<std::string, h5pp::TableInfo> &srcTableDb, const PathId &pathid, const std::vector<ScaleKey> &srcScaleKeys,
                        const FileId &fileId) {
        // In this function we take time series data from each srcTable and create multiple tables tgtTable, one for each
        // time point (iteration). Each entry in tgtTable corresponds to the same time point on different realizations.
        auto                   t_scope = tid::tic_scope(__FUNCTION__);
        std::vector<std::byte> srcReadBuffer;
        for(const auto &srcKey : srcScaleKeys) {
            if(srcTableDb.find(srcKey.key) == srcTableDb.end()) throw std::logic_error(h5pp::format("Key [{}] was not found in source map", srcKey.key));
            auto &srcInfo    = srcTableDb[srcKey.key];
            auto  srcRecords = srcInfo.numRecords.value();
            tools::logger::log->trace("Transferring scale table {}", srcInfo.tablePath.value());
            auto tgtName = h5pp::fs::path(srcInfo.tablePath.value()).filename().string();
            auto tgtPath = pathid.scale_path(tgtName, srcKey.dim);

            if(tgtTableDb.find(tgtPath) == tgtTableDb.end()) {
                auto t_create = tid::tic_scope("createTable");
                tools::logger::log->debug("Adding target scale {}", tgtPath);
                h5pp::TableInfo tableInfo = h5_tgt.getTableInfo(tgtPath);
                // Disabling compression is supposed to give a nice speedup. Read here:
                // https://support.hdfgroup.org/HDF5/doc1.8/Advanced/DirectChunkWrite/UsingDirectChunkWrite.pdf
                if(not tableInfo.tableExists.value())
                    tableInfo = h5_tgt.createTable(srcInfo.h5Type.value(), tgtPath, srcInfo.tableTitle.value(), std::nullopt, true);

                tgtTableDb[tgtPath] = tableInfo;
                //                    tgtTableDb[tgtPath].info.assertWriteReady();
            }
            auto &tgtId   = tgtTableDb[tgtPath];
            auto &tgtBuff = tgtId.buff;
            //                auto &tgtInfo = tgtTableDb[tgtPath].info;

            // Determine the target index where to copy this record
            // Under normal circumstances, the "index" counts the number of realizations, or simulation seeds.
            // tgtInfo.numRecords is the total number of realizations registered until now

            hsize_t index = tgtId.get_index(fileId.seed); // Positive if it has been read already
                                                          //            if(index != std::numeric_limits<hsize_t>::max()) {
                                                          //                // The table entry for the current realization may already have been added
            //                // This happens for instance if we make extra entries in "finished" after adding them to "checkpoint" and/or "savepoint"
            //                // However, these entries should be identical, so no need to copy them again.
            // #pragma message "Skipping here may not be wise?"
            //                tools::logger::log->info("Skip copying existing scale entry: {} | index {}", tgtPath, index);
            //                continue;
            //            }

            //            index = index == std::numeric_limits<hsize_t>::max() ? tgtBuff.get_count() : index;

            // read a source record at "iter" into the "index" position in the buffer
            auto t_buffer = tid::tic_scope("bufferScaleRecords");
            srcReadBuffer.resize(srcInfo.recordBytes.value());
            h5pp::hdf5::readTableRecords(srcReadBuffer, srcInfo, srcRecords - 1, 1); // Read the last entry
            tgtBuff.insert(srcReadBuffer, index);

            // copy/append a source record at "iter" into the "index" position on the table.
            //                tools::logger::log->trace("Copying scale index {} -> {}: {}", rec, index, tgtPath);
            //                auto t_copy = tid::tic_scope("copyTableRecords");
            //                h5_tgt.copyTableRecords(srcInfo, rec, 1, tgtInfo, static_cast<hsize_t>(index));
            // Update the database
            tgtId.insert(fileId.seed, index);
        }
    }

}
