#include "h5xf.h"
#include "config/enums.h"
#include "debug/exceptions.h"
#include "general/seq.h"
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
                // In this case we are collecting datasets onto a superdataset with dimension +1, enumerating them by that axis
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
                tgtDims[srcKey.axis] = 0;          // Create with 0 extent in the new axis direction, so that the dataset starts empty (zero volume)

                // Determine a good chunksize between 10 and 500 elements
                auto tgtChunk           = tgtDims;
                auto tgtMaxDims         = tgtDims;
                auto chunkSize          = std::clamp(5e4 / static_cast<double>(srcInfo.dsetByte.value()), 10., 10000.);
                tgtChunk[srcKey.axis]   = static_cast<hsize_t>(chunkSize); // number of elements in the axis that we append into.
                tgtMaxDims[srcKey.axis] = H5S_UNLIMITED;
                if(srcKey.size == Size::VAR) {
                    auto        srcGroupPath    = h5pp::fs::path(srcInfo.dsetPath.value()).parent_path().string();
                    std::string statusTablePath = fmt::format("{}/status", srcGroupPath);
                    long        bond_max        = h5_src.readTableField<long>(statusTablePath, "bond_max", h5pp::TableSelection::FIRST);
                    tgtDims[0]                  = static_cast<hsize_t>(bond_max);
                }

                tools::logger::log->info("Adding target dset {} | dims {} | chnk {}", tgtPath, tgtDims, tgtChunk);
                tgtDsetDb[tgtPath] = h5_tgt.createDataset(tgtPath, srcInfo.h5Type.value(), H5D_CHUNKED, tgtDims, tgtChunk, tgtMaxDims, 3);
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
            if(srcTableDb.find(srcKey.key) == srcTableDb.end()) throw except::range_error("{}: Key [{}] was not found in source map", __FUNCTION__, srcKey.key);
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
            tools::logger::log->debug("Transferring table {} | src idx {} -> tgt idx {}", tgtPath, srcInfo.numRecords.value(), index);
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
        std::vector<size_t>       srcIndices; // We can assume all tables have the same indexing numbers. Only update on mismatch
        std::vector<StorageEvent> srcEvents;  // This stores the "event" row in the table. We need to make sure to only transfer StorageEvent::ITER_STATE
        std::vector<std::byte>    srcReadBuffer;
        for(const auto &srcKey : srcKeys) {
            if(srcTableDb.find(srcKey.key) == srcTableDb.end()) throw except::logic_error("Key [{}] was not found in source map", srcKey.key);
            auto &srcInfo = srcTableDb[srcKey.key];

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
                        throw except::range_error("{}: Index field was not found in table [{}]\n Expected fields {}\n Existing fields {}", __FUNCTION__,
                                                  srcInfo.tablePath.value(), srcKey.index, srcInfo.fieldNames.value());
                    if(not eventfield)
                        throw except::range_error("{}: Event field was not found in table [{}]\n Existing fields {}", __FUNCTION__, srcInfo.tablePath.value(),
                                                  srcInfo.fieldNames.value());
                    h5pp::hdf5::readTableField(srcIndices, srcInfo, {indexfield.value()});
                    h5pp::hdf5::readTableField(srcEvents, srcInfo, {eventfield.value()});
                }
            } catch(const std::exception &ex) { throw std::logic_error(fmt::format("Failed to get iteration numbers: {}", ex.what())); }
            for(size_t rec = 0; rec < srcInfo.numRecords.value(); rec++) {
                size_t srcIndex = rec;
                if(not srcIndices.empty()) srcIndex = srcIndices[rec]; // Get the actual index number from the table
                if(not srcEvents.empty()) {
                    if(not seq::has(srcKey.event, srcEvents[rec])) continue;
                }
                auto tgtName = h5pp::fs::path(srcInfo.tablePath.value()).filename().string();
                auto tgtPath = pathid.create_path<KeyT>(tgtName, srcIndex);

                if(tgtTableDb.find(tgtPath) == tgtTableDb.end()) {
                    auto t_create = tid::tic_scope("createTable");
                    tools::logger::log->debug("Adding target {} {}", srcKey.classtag, tgtPath);
                    h5pp::TableInfo tableInfo = h5_tgt.getTableInfo(tgtPath);
                    // Disabling compression is supposed to give a nice speedup. Read here:
                    // https://support.hdfgroup.org/HDF5/doc1.8/Advanced/DirectChunkWrite/UsingDirectChunkWrite.pdf
                    if(not tableInfo.tableExists.value()) {
                        auto recordBytes = H5Tget_size(srcInfo.h5Type.value());
                        auto targetBytes = 500 * 1024; // 500 kB
                        auto chunkSize   = static_cast<hsize_t>(std::ceil(targetBytes / static_cast<double>(recordBytes)));
                        //                        tools::logger::log->info("Chunk size: {} | {}", chunkSize, tgtPath);
                        tableInfo = h5_tgt.createTable(srcInfo.h5Type.value(), tgtPath, srcInfo.tableTitle.value(), chunkSize, 6);
                    }
                    tgtTableDb[tgtPath] = tableInfo;
                }
                auto &tgtId   = tgtTableDb[tgtPath];
                auto &tgtBuff = tgtId.buff;

                // Determine the target index where to copy this record
                // tgtInfo.numRecords is the total number of realizations registered until now
                // tgtId.get_index(fileId.seed) Gets the index number if it already exists, else
                // the size of the table (so we can append)
                hsize_t tgtIndex = tgtId.get_index(fileId.seed);

                tools::logger::log->trace("Transferring {} table {} ({} records) | {} {} @ idx {} -> tgt {}", srcKey.classtag, srcInfo.tablePath.value(),
                                          srcInfo.numRecords.value(), srcKey.index, srcIndex, rec, tgtIndex);

                // read a source record at "srcIndex" into the "tgtIndex" position in the buffer
                auto t_buffer  = tid::tic_scope(fmt::format("bufferRecords-{}", srcKey.classtag));
                srcInfo.h5Type = H5Dget_type(srcInfo.h5Dset.value());
                h5pp::hdf5::readTableRecords(srcReadBuffer, srcInfo, rec, 1, h5_tgt.plists);
                tgtBuff.insert(srcReadBuffer, tgtIndex);
                if(srcInfo.reclaimInfo.has_value()) {
                    // One or more table fields are vlen arrays.
                    // This means that HDF5 has allocated memory for those vlen arrays on a pointer somewhere in srcReadBuffer.
                    // This memory needs to be de-allocated to avoid a memory leak, which unfortunately means that we have to flush
                    // the buffer immediately before srcReadBuffer is overwritten in the next round.
                    // When calling tgtBuff.insert we copy the pointer to the data, but not the data itself.
                    // If we then overwrite srcReadBuffer, we lose track of the vlen memory pointed to previously.
                    tgtBuff.flush();
                    srcInfo.reclaim();
                }
                // Update the database
                tgtId.insert(fileId.seed, tgtIndex);
            }
        }
    }
    template void transferSeries(h5pp::File &h5_tgt, std::unordered_map<std::string, InfoId<BufferedTableInfo>> &tgtTableDb,
                                 std::unordered_map<std::string, h5pp::TableInfo> &srcTableDb, const PathId &pathid, const std::vector<CronoKey> &srcKeys,
                                 const FileId &fileId);
    template void transferSeries(h5pp::File &h5_tgt, std::unordered_map<std::string, InfoId<BufferedTableInfo>> &tgtTableDb,
                                 std::unordered_map<std::string, h5pp::TableInfo> &srcTableDb, const PathId &pathid, const std::vector<FesDnKey> &srcKeys,
                                 const FileId &fileId);
    template void transferSeries(h5pp::File &h5_tgt, std::unordered_map<std::string, InfoId<BufferedTableInfo>> &tgtTableDb,
                                 std::unordered_map<std::string, h5pp::TableInfo> &srcTableDb, const PathId &pathid, const std::vector<FesUpKey> &srcKeys,
                                 const FileId &fileId);

}
