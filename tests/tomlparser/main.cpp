#include "config/enums.h"
#include "debug/exceptions.h"
#include "general/sfinae.h"
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <iostream>
#include <toml.hpp>

template<typename T>
T get_value(const toml::table &tbl, std::string_view path) {
    auto node = tbl.at_path(path);
    if(!node) throw except::runtime_error("node [{}] is null", path);
    if constexpr(std::is_enum_v<T>) {
        auto value = node.value<std::string_view>();
        if(value.has_value())
            return sv2enum<T>(value.value());
        else
            throw std::runtime_error("Missing value for node");
    } else if constexpr(sfinae::is_std_vector_v<T>) {
        using value_type = typename T::value_type;
        auto *array      = node.as_array();
        T     res;
        if(array) {
            for(auto &&item : *array) {
                if(item.value<value_type>().has_value()) res.emplace_back(item.value<value_type>().value());
            }
            return res;
        }

    } else if constexpr(sfinae::is_std_array_v<T>) {
        if(!node.is_array()) throw except::runtime_error("toml: node is not an array");
        using value_type = typename T::value_type;
        T     res;
        auto *array = node.as_array();
        if(array == nullptr) throw except::runtime_error("array is nullptr");
        // auto val = array->as<double>()->get();
        if(res.size() != array->size())
            throw except::runtime_error("tomlplusplus: expected {}, but parsed array of size = {}", sfinae::type_name<T>(), array->size());
        // std::copy(array->begin(), array->end(), res.begin());
        std::transform(array->begin(), array->end(), res.begin(), [](auto &&elem) -> value_type {
            if(!elem.is_value()) throw except::runtime_error("tomlplusplus: elem is not a value");
            // if(!elem.has_value()) throw except::runtime_error("tomlplusplus:: value has no value");
            return elem.template as<value_type>()->get();
            // auto ptr = elem.template as<value_type>();
            // if(!ptr)  throw except::runtime_error("tomlplusplus:: &elem.as<value_type>() does not point to a value");
            // return ptr->value();
        });
        return res;
    }
}

// void print_node(const toml::node &node) {
//     fmt::print("{} \n", std::string_view(node));
// }

template<typename T>
struct tomlitem {
    std::string key;
    T value;
    std::string comment;
};


void print_table(const toml::table &tbl, std::string_view path) {
    size_t depth = std::count(path.begin(), path.end(), '.');
    if(!path.empty()) fmt::print("[{}]\n", path);
    for(const auto &[key, node] : tbl) {
        auto ksv = std::string_view{key};
        if(node.is_table())
            if(path.empty()) {
                print_table(*node.as_table(), fmt::format("{}", ksv));
            } else {
                print_table(*node.as_table(), fmt::format("{}.{}", path, ksv));
            }
        else { fmt::print(fmt::runtime(" {} {:>{}} {}\n"),   ksv, "=", std::max(50ul, ksv.size()) - ksv.size(), "value"); }

        // if (node.is_table()) print_node(node);
        // fmt::print("{} | is_table {} | source {}\n", std::string_view(key), node.is_table(), std::string_view(node.source().begin(), node.source().end()));
        // fmt::print("{}\n", std::string_view(key));
        // if(node.is_table()) {
        //     for (auto & [key, snode] : *node.as_table()) {
        //         fmt::print(" {}\n", );
        //
        //     }
        // }
    }
}

int main(int argc, char **argv) {
    toml::table tbl;
    try {
        tbl = toml::parse_file(TEST_TOML_FILE);
        //        auto policy = get_value<StoragePolicy>(tbl["storage"]["dataset"]["lbit_analysis"]["policy"]);
        auto policy = get_value<StoragePolicy>(tbl, "storage.dataset.lbit_analysis.policy");
        fmt::print("policy:{}\n", enum2sv(policy));
        fmt::print("testarr:{}\n", get_value<std::vector<double>>(tbl, "model.lbit.testarr"));
        fmt::print("testarr:{}\n", get_value<std::vector<int>>(tbl, "model.lbit.testarr"));
        fmt::print("testarr:{}\n", get_value<std::array<double, 3>>(tbl, "model.lbit.testarr"));
        // fmt::print("testarr:{}\n", get_value<std::array<double, 2>>(tbl, "model.lbit.testarr"));
        // fmt::print("testarr:{}\n", get_value<std::vector<double>>(tbl, "model.lbit.testarr2"));

    } catch(const toml::parse_error &err) {
        std::cerr << "Parsing failed:\n" << err << "\n";
        return 1;
    }

    // print_table(tbl, "");
    // tbl.insert("xdmrg.testitem", 4.32);
    tbl["xdmrg"].as_table()->insert("testitem", 4.32);
    // tbl.insert("udmrg");

    print_table(tbl, "");

    std::cout << tbl << "\n";

    return 0;
}
