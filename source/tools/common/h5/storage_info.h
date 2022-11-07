#pragma once
#include "config/enums.h"
#include <string_view>
class AlgorithmStatus;
struct StorageInfo {
    public:
    size_t        iter      = -1ul;
    size_t        step      = -1ul;
    long          position  = -1;
    long          direction = 0;
    long          bond_lim  = -1;
    long          bond_max  = -1;
    double        trnc_lim  = -1;
    AlgorithmType algo_type = AlgorithmType::ANY;
    std::string   algo_name;
    StorageEvent &storage_event; // Reference to the status.event
    StorageLevel  storage_level = StorageLevel::NONE;
    std::string   state_name;
    void          assert_well_defined() const;
    std::string   get_state_prefix() const;
    std::string   get_mps_prefix() const;
    StorageInfo(const AlgorithmStatus &status, std::string_view state_name, StorageEvent event = StorageEvent::NONE);
    ~StorageInfo() noexcept;
};

// Metadata needed to uniquely identify a save point for a given HDF5 link
struct StorageAttrs {
    size_t       iter          = -1ul;
    size_t       step          = -1ul;
    long         bond_lim      = -1;
    long         bond_max      = -1;
    double       trnc_lim      = -1;
    StorageEvent storage_event = StorageEvent::NONE;
    StorageLevel storage_level = StorageLevel::NONE;
    mutable bool link_exists   = false;
    StorageAttrs()             = default;
    StorageAttrs(const StorageInfo &sinfo);
    bool operator==(const StorageAttrs &sinfo);
    bool operator==(const StorageInfo &sinfo);
};
