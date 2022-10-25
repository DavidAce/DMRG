#pragma once
#include <string>

// Define the allowed items
enum class Type { INT, LONG, DOUBLE, COMPLEX, TID };
enum class Size { FIX, VAR };
// struct DsetKey {
//     Type        type;
//     Size        size = Size::FIX;
//     std::string name;
//     std::string key;
// };

struct Key {
    std::string algo, state, point;
    std::string name;
    std::string key;
    size_t      expected_size = -1ul;
    Key()                     = default;
    Key(std::string_view algo_, std::string_view state_, std::string_view point_, std::string_view name_, size_t expected_size_ = -1ul)
        : algo(algo_), state(state_), point(point_), name(name_), expected_size(expected_size_) {}
};

struct DsetKey : public Key {
    static constexpr std::string_view classtag = "dset";
    Size                              size     = Size::FIX;
    size_t                            axis;
    DsetKey(std::string_view algo_, std::string_view state_, std::string_view point_, std::string_view name_, Size size_, size_t axis_)
        : Key(algo_, state_, point_, name_), size(size_), axis(axis_) {}
};

struct TableKey : public Key {
    using Key::Key;
    static constexpr std::string_view classtag = "table";
};

struct BonddKey : public Key {
    using Key::Key;
    static constexpr std::string_view                classtag = "bondd"; // Name of this type of key
    static constexpr std::array<std::string_view, 3> index    = {"bond_lim", "bond_limit",
                                                                 "bond_dimension_limit"}; // Index this data according to this table field (first match)
};

struct CronoKey : public Key {
    using Key::Key;
    static constexpr std::string_view                classtag = "crono";  // Name of this type of key
    static constexpr std::array<std::string_view, 1> index    = {"iter"}; // Index this data according to this table field (first match)
};

struct ScaleKey : public Key {
    static constexpr std::string_view classtag = "scale";
    std::vector<std::string>          allKeys;
    std::string                       scale; // The group pattern, usually something like "bond_*"
    size_t                            dim;   // The actual value of dim found
    ScaleKey(std::string_view algo_, std::string_view state_, std::string_view point_, std::string_view scale_, std::string_view name_,
             std::vector<int> bondlims_ = {})
        : Key(algo_, state_, point_, name_), scale(scale_), dim(-1ul) {
        allKeys.reserve(bondlims_.size());
        for(const auto &b : bondlims_) allKeys.emplace_back(h5pp::format("bond_{}", b));
    }
};

// struct BonddKey : public Key {
//     std::string group; // The group pattern, usually something like "bond_*"
//     size_t dim; // The actual value of the dim dimension found
//     std::vector<size_t> dims = {}; // All the dim dimensions that are expected. This has to be entered manually, and extrapolation occurs when dim isn't
//     found BonddKey(std::string_view algo_, std::string_view state_, std::string_view point_, std::string_view group_, std::string_view name_, const
//     std::vector<size_t> & dims_)
//         : Key(algo_, state_, point_, name_), group(group_), dim(-1ul), dims(dims_){}
// };

struct ModelKey {
    static constexpr std::string_view tag = "model";
    std::string                       algo, model;
    std::string                       name;
    std::string                       key;
    ModelKey(std::string_view algo_, std::string_view model_, std::string_view name_) : algo(algo_), model(model_), name(name_) {}
};
