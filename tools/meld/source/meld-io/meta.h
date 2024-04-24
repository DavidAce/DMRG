#pragma once
#include "config/enums.h"
#include "general/enums.h"
#include <string>

// Define the allowed items
enum class Type { INT, LONG, DOUBLE, COMPLEX, TID };
enum class Size { FIX, VAR };

struct Key {
    std::string    algo, state, subgroup;
    std::string    name;
    std::string    key;
    mutable size_t expected_size = -1ul;
                   Key()         = default;
                   Key(std::string_view algo_, std::string_view state_, std::string_view subgroup_, std::string_view name_, size_t expected_size_ = -1ul)
        : algo(algo_), state(state_), subgroup(subgroup_), name(name_), expected_size(expected_size_) {}
    virtual std::string create_source_path() const {
        std::string srcpath;
        if(!algo.empty()) srcpath += fmt::format("/{}", algo);
        if(!state.empty()) srcpath += fmt::format("/{}", state);
        if(!subgroup.empty()) srcpath += fmt::format("/{}", subgroup);
        if(!name.empty()) srcpath += fmt::format("/{}", name);
        return srcpath;
    };
};

struct DsetKey : public Key {
    static constexpr std::string_view classtag = "dset";
    Size                              size     = Size::FIX;
    size_t                            axis;
    SlabSelect                        ssel = SlabSelect::FULL;
    DsetKey(std::string_view algo_, std::string_view state_, std::string_view subgroup, std::string_view name_, Size size_, size_t axis_,
            SlabSelect ssel_ = SlabSelect::FULL)
        : Key(algo_, state_, subgroup, name_), size(size_), axis(axis_), ssel(ssel_) {}
};

struct TableKey : public Key {
    using Key::Key;
    static constexpr std::string_view classtag = "table";
    static constexpr std::array       eventkey = {StorageEvent::FINISHED}; // Grab table rows with these event types
};

struct FesUpKey : public Key {
    using Key::Key;
    static constexpr std::string_view classtag = "fesup";                     // Name of this type of key
    static constexpr std::array       eventkey = {StorageEvent::BOND_UPDATE}; // Grab table rows with these event types
    static constexpr std::array       indexkey = {"bond_lim"};                // Index this data according to this table field (first match)
    static constexpr std::string_view indexpfx = "bond_lim_";                 // Prefix for the target subgroup (e.g. /xDMRG/state_emid/fes/bond_2/<name>)
};

struct RbdsKey : public Key {
    using Key::Key;
    static constexpr std::string_view classtag = "rbds";                    // Name of this type of key
    static constexpr std::array       eventkey = {StorageEvent::RBDS_STEP}; // Grab table rows with these event types
    static constexpr std::array       indexkey = {"bond_lim"};              // Index according to this column name (first match or uses the record number)
    static constexpr std::string_view indexpfx = "bond_lim_";               // Prefix for the target subgroup (e.g. /xDMRG/state_emid/rbds/bond_lim_2/<name>)
};

struct RtesKey : public Key {
    using Key::Key;
    static constexpr std::string_view classtag = "rtes";                    // Name of this type of key
    static constexpr std::array       eventkey = {StorageEvent::RTES_STEP}; // Grab table rows with these event types
    static constexpr std::array       indexkey = {""};                      // Index according to this column name (first match or uses the record number)
    static constexpr std::string_view indexpfx = "trnc_lim_";               // Prefix for the target subgroup (e.g. /xDMRG/state_emid/rtes/trnc_lim_2/<name>)
};

struct CronoKey : public Key {
    using Key::Key;
    static constexpr std::string_view classtag = "crono";                   // Name of this type of key
    static constexpr std::array       eventkey = {StorageEvent::ITERATION}; // Grab table rows with these event types
    static constexpr std::array       indexkey = {"iter"};                  // Index this data according to this table field (first match)
    static constexpr std::string_view indexpfx = "iter_";                   // Prefix for the target subgroup (e.g. /xDMRG/state_emid/crono/iter_2/<name>)
};

struct ModelKey {
    static constexpr std::string_view classtag = "model";
    std::string                       algo, model;
    std::string                       name;
    std::string                       key;
    ModelKey(std::string_view algo_, std::string_view model_, std::string_view name_) : algo(algo_), model(model_), name(name_) {}
};
