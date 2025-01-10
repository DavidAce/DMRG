#include "config/enums.h"
#include "general/sfinae.h"
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <iostream>
#include <toml.hpp>

template<typename T>
T get_value(toml::node_view<toml::node> node) {
    if constexpr(std::is_enum_v<T>) {
        auto value = node.value<std::string_view>();
        if(value.has_value())
            return sv2enum<T>(value.value());
        else
            throw std::runtime_error("Missing value for node");
    } else if constexpr(sfinae::is_std_vector_v<T>) {
        using value_type = typename T::value_type;
        T            res;
        toml::array *array = node.as_array();
        if(array) {
            for(auto &&item : *array) {
                if(item.value<value_type>().has_value()) res.emplace_back(item.value<value_type>().value());
            }
            return res;
        }

    } else if constexpr(sfinae::is_std_array_v<T>) {
        using value_type = typename T::value_type;
        T            res;
        toml::array *array = node.as_array();
        if(array and array->size() == res.size()) {
            for(size_t i = 0; i < res.size(); ++i) {
                if(array->at(i).value<value_type>().has_value()) res[i] = array->at(i).value<value_type>().value();
            }
            return res;
        }
    }
}

int main(int argc, char **argv) {
    toml::table tbl;
    try {
        tbl = toml::parse_file(TEST_TOML_FILE);
        //        auto policy = get_value<StoragePolicy>(tbl["storage"]["dataset"]["lbit_analysis"]["policy"]);
        auto policy = get_value<StoragePolicy>(tbl.at_path("storage.dataset.lbit_analysis.policy"));
        fmt::print("policy:{}\n", enum2sv(policy));
        fmt::print("testarr:{}\n", get_value<std::vector<int>>(tbl.at_path("model.lbit.testarr")));
        fmt::print("testarr:{}\n", get_value<std::array<int, 3>>(tbl.at_path("model.lbit.testarr")));
    } catch(const toml::parse_error &err) {
        std::cerr << "Parsing failed:\n" << err << "\n";
        return 1;
    }
    return 0;
}
