#include "config/enums.h"
#include "debug/exceptions.h"
#include "general/sfinae.h"
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <iostream>
#include <toml.hpp>

// template<typename T>
// T get_value(const toml::table &tbl, std::string_view path) {
//     auto node = tbl.at_path(path);
//     if(!node) throw except::runtime_error("node [{}] is null", path);
//     if constexpr(std::is_enum_v<T>) {
//         auto value = node.value<std::string_view>();
//         if(value.has_value())
//             return sv2enum<T>(value.value());
//         else
//             throw std::runtime_error("Missing value for node");
//     } else if constexpr(sfinae::is_std_vector_v<T>) {
//         using value_type = typename T::value_type;
//         auto *array      = node.as_array();
//         T     res;
//         if(array) {
//             for(auto &&item : *array) {
//                 if(item.value<value_type>().has_value()) res.emplace_back(item.value<value_type>().value());
//             }
//             return res;
//         }
//
//     } else if constexpr(sfinae::is_std_array_v<T>) {
//         if(!node.is_array()) throw except::runtime_error("toml: node is not an array");
//         using value_type = typename T::value_type;
//         T     res;
//         auto *array = node.as_array();
//         if(array == nullptr) throw except::runtime_error("array is nullptr");
//         // auto val = array->as<double>()->get();
//         if(res.size() != array->size())
//             throw except::runtime_error("tomlplusplus: expected {}, but parsed array of size = {}", sfinae::type_name<T>(), array->size());
//         // std::copy(array->begin(), array->end(), res.begin());
//         std::transform(array->begin(), array->end(), res.begin(), [](auto &&elem) -> value_type {
//             if(!elem.is_value()) throw except::runtime_error("tomlplusplus: elem is not a value");
//             // if(!elem.has_value()) throw except::runtime_error("tomlplusplus:: value has no value");
//             return elem.template as<value_type>()->get();
//             // auto ptr = elem.template as<value_type>();
//             // if(!ptr)  throw except::runtime_error("tomlplusplus:: &elem.as<value_type>() does not point to a value");
//             // return ptr->value();
//         });
//         return res;
//     }
// }
//

// void print_node(const toml::node &node) {
//     fmt::print("{} \n", std::string_view(node));
// }
//
// void print_table(const toml::table &tbl) {
//     for (const auto &[key, node]: tbl) {
//         // if (node.is_table()) print_node(node);
//         fmt::print("{} | is_table {} | source {}\n", std::string_view(key), node.is_table(), std::string_view(node.source().begin(), node.source().end()));
//         // fmt::print("{} | is_table {}\n", std::string_view(key), node.is_table());
//     }
// }

// auto find_recurse(auto data, std::string key) {
//     auto subkey = key.substr(0, key.find_first_of('.'));
//     fmt::print("{}\n", subkey);
//     if(subkey.empty()) {
//         return key;
//     } else
//         return find_recurse(data, subkey);
// }

template<typename T>
auto find_recurse(auto data, std::string key) {
    while(true) {
        auto dotpos = key.find_first_of('.');
        if(dotpos == std::string::npos) return toml::find<T>(data, key);
        // advance
        fmt::print("could not find {} in data\n", key);
        data = toml::find(data, key.substr(0, dotpos));
        key  = key.substr(dotpos + 1);
    }
}

size_t get_maxkeywidth(const toml::ordered_value &table) {
    // fmt::print("{}\n");
    size_t maxkeywidth = 0;
    if(table.is_table()) {
        auto maxkeyit = std::max_element(table.as_table().begin(), table.as_table().end(), [](auto a, auto b) -> bool { return a.first.size() < b.first.size(); });
        if(maxkeyit != table.as_table().end()) {maxkeywidth = std::max(maxkeywidth, maxkeyit->first.size());}
        for(const auto &[key, elem] : table.as_table()) {
            if(elem.is_table())
                maxkeywidth = std::max(maxkeywidth, get_maxkeywidth(elem));
            else {
                auto keyval = fmt::format("{}{:>{}} {}", key, "=", std::max<size_t>(key.size()+1, maxkeywidth+1) - key.size(), toml::format(elem));
                fmt::print("{}{:>{}} {}\n", keyval, "", std::max(30ul, keyval.size()) - keyval.size(), elem.comments().empty() ? "#" : elem.comments().front());
            }
        }
    }
}

// template<typename T>
void print_table(const toml::ordered_value &table) {
    // fmt::print("{}\n");
    if(table.is_table()) {
        size_t maxkeywidth = 0;
        auto maxkeyit = std::max_element(table.as_table().begin(), table.as_table().end(), [](auto a, auto b) -> bool { return a.first.size() < b.first.size(); });
        if(maxkeyit != table.as_table().end()) {maxkeywidth = maxkeyit->first.size();}
        for(const auto &[key, elem] : table.as_table()) {
            if(elem.is_table())
                print_table(elem);
            else {
                auto keyval = fmt::format("{}{:>{}} {}", key, "=", std::max<size_t>(key.size()+1, maxkeywidth+1) - key.size(), toml::format(elem));
                fmt::print("{}{:>{}} {}\n", keyval, "", std::max(30ul, keyval.size()) - keyval.size(), elem.comments().empty() ? "#" : elem.comments().front());
            }
        }
    }
}

int main(int argc, char **argv) {
    auto data = toml::parse(TEST_TOML_FILE);
    fmt::print("{}\n", toml::find<std::array<double, 3>>(data.at("model").at("lbit"), "testarr"));
    data.at("model").at("lbit").at("testarr") = std::array<double, 3>({2., 3., 4.});
    fmt::print("{} | {}\n", toml::find<std::array<double, 3>>(data.at("model").at("lbit"), "testarr"), data.at("model").at("lbit").at("testarr").comments());
    data.at("model");

    auto tbl = toml::ordered_value(
        toml::ordered_table{
            {"a", 42},
            {"b", "foo"},
            {"c", "f0000000000oo"},
        },
        toml::table_format_info{.fmt = toml::table_format::multiline});
    tbl.comments().emplace_back("# This is a main comment");
    tbl["a"].comments().emplace_back("# This is a comment");
    tbl["a"].comments().emplace_back("# This is an inline comment");
    tbl["test"]        = toml::ordered_value(toml::ordered_table{}, toml::table_format_info{.fmt = toml::table_format::dotted});
    tbl["test"]["val"] = 1.4;
    print_table(tbl);
    // tbl.emplace_back(std::make_pair("test", toml::table()));

    // print_table(tbl);
    // std::string formatted_table = toml::format(tbl, toml::spec::v(1, 0, 1));
    // fmt::print("{}\n", formatted_table);

    // std::cout << toml::format(tbl) << std::endl;

    // tbl.emplace_back(toml::basic_value<toml::table>());
    // auto model = toml::find(data, "model");
    // auto lbit  = toml::find(model, "lbit");
    // fmt::print("testarr:{}\n", toml::find<std::vector<double>>(lbit, "testarr"));
    // fmt::print("recurse:{}\n", find_recurse<std::vector<double>>(data, "model.lbit.testarr"));

    // try {

    //     //        auto policy = get_value<StoragePolicy>(tbl["storage"]["dataset"]["lbit_analysis"]["policy"]);
    //     auto policy = get_value<StoragePolicy>(tbl, "storage.dataset.lbit_analysis.policy");
    //     fmt::print("policy:{}\n", enum2sv(policy));
    //     fmt::print("testarr:{}\n", get_value<std::vector<double>>(tbl, "model.lbit.testarr"));
    //     fmt::print("testarr:{}\n", get_value<std::vector<int>>(tbl, "model.lbit.testarr"));
    //     fmt::print("testarr:{}\n", get_value<std::array<double, 3>>(tbl, "model.lbit.testarr"));
    //     // fmt::print("testarr:{}\n", get_value<std::array<double, 2>>(tbl, "model.lbit.testarr"));
    //     // fmt::print("testarr:{}\n", get_value<std::vector<double>>(tbl, "model.lbit.testarr2"));
    //
    // } catch(const toml::parse_error &err) {
    //     std::cerr << "Parsing failed:\n" << err << "\n";
    //     return 1;
    // }
    //
    // print_table(tbl);
    return 0;
}
