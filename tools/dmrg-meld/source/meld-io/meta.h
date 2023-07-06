#pragma once
#include "config/enums.h"
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
    std::string algo, state;
    std::string name;
    std::string key;
    mutable size_t      expected_size = -1ul;
    Key()                     = default;
    Key(std::string_view algo_, std::string_view state_, std::string_view name_, size_t expected_size_ = -1ul)
        : algo(algo_), state(state_), name(name_), expected_size(expected_size_) {}
};

struct DsetKey : public Key {
    static constexpr std::string_view classtag = "dset";
    Size                              size     = Size::FIX;
    size_t                            axis;
    DsetKey(std::string_view algo_, std::string_view state_, std::string_view name_, Size size_, size_t axis_)
        : Key(algo_, state_, name_), size(size_), axis(axis_) {}
};

struct TableKey : public Key {
    using Key::Key;
    static constexpr std::string_view            classtag = "table";
    static constexpr std::array<StorageEvent, 1> event    = {StorageEvent::LAST_STATE}; // Grab table rows with these event types
};

struct FesUpKey : public Key {
    using Key::Key;
    static constexpr std::string_view                classtag = "fesup";                       // Name of this type of key
    static constexpr std::array<StorageEvent, 1>     event    = {StorageEvent::BOND_INCREASE}; // Grab table rows with these event types
    static constexpr std::array<std::string_view, 1> index    = {"bond_lim"};                  // Index this data according to this table field (first match)
};

struct FesDnKey : public Key {
    using Key::Key;
    static constexpr std::string_view                classtag = "fesdn";                   // Name of this type of key
    static constexpr std::array<StorageEvent, 1>     event    = {StorageEvent::FES_STATE}; // Grab table rows with these event types
    static constexpr std::array<std::string_view, 1> index    = {"bond_lim"};              // Index this data according to this table field (first match)
};

struct CronoKey : public Key {
    using Key::Key;
    static constexpr std::string_view                classtag = "crono";                    // Name of this type of key
    static constexpr std::array<StorageEvent, 1>     event    = {StorageEvent::ITER_STATE}; // Grab table rows with these event types
    static constexpr std::array<std::string_view, 1> index    = {"iter"};                   // Index this data according to this table field (first match)
};

struct ModelKey {
    static constexpr std::string_view classtag = "model";
    std::string                       algo, model;
    std::string                       name;
    std::string                       key;
    ModelKey(std::string_view algo_, std::string_view model_, std::string_view name_) : algo(algo_), model(model_), name(name_) {}
};
