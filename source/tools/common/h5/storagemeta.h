#pragma once
#include <config/enums.h>
#include <string>
#include <vector>

struct StorageMeta {
    StorageLevel             storage_level;
    StorageReason            storage_reason;
    std::string              algo_name;
    std::string              state_name;
    std::string              state_prefix;
    std::string              model_prefix;
    std::string              timer_prefix;
    std::vector<std::string> table_prefxs;
};
