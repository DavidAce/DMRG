#include "tid.h"
#include <fmt/format.h>
#include <regex>
#include <string>
namespace tid {
    namespace internal {

        const ur *ur_ref_t::operator->() const { return &ref.get(); }

        [[nodiscard]] std::string ur_ref_t::str() const {
            return fmt::format(FMT_STRING("{0:<{1}} {2:>8.3f} s | sum {3:>8.3f} s | {4:>6.2f} % | avg {5:>8.2e} s | level {6} | count {7}"), key,
                               tree_max_key_size, ref.get().get_time(), sum, 100 * frac, ref.get().get_time_avg(), level2sv(ref.get().get_level()),
                               ref.get().get_tic_count());
        }

        template<typename T>
        T split(std::string_view strv, std::string_view delims) {
            T output;
            //            output.reserve(strv.length() / 4);
            auto first = strv.begin();
            while(first != strv.end()) {
                const auto second = std::find_first_of(first, std::cend(strv), std::cbegin(delims), std::cend(delims));
                const auto pos    = static_cast<std::string_view::size_type>(std::distance(strv.begin(), first));
                const auto len    = static_cast<std::string_view::size_type>(std::distance(first, second));
                if(first != second) output.emplace_back(strv.substr(pos, len));
                if(second == strv.end()) break;
                first = std::next(second);
            }
            //            output.shrink_to_fit();
            return output;
        }
        template std::vector<std::string_view> split(std::string_view strv, std::string_view delims);
        template std::deque<std::string_view>  split(std::string_view strv, std::string_view delims);
    }

    ur &get(std::deque<std::string_view> &keys, level l, ur &u) {
        if(keys.empty()) return u;
        auto &u_sub = u.insert(keys[0], l);
        keys.pop_front();
        return get(keys, l, u_sub);
    }
    ur &get(std::string_view key, level l, std::optional<std::string_view> prefix) {
        std::string parsed_key(key);
        std::string prefix_key;
        if(prefix)
            prefix_key = std::string(prefix.value());
        else
            prefix_key = internal::ur_prefix;

        // Use prepended '.' to go to parent scope
        while(true) {
            auto pos_dd = parsed_key.rfind('.', 0);
            if(pos_dd != std::string::npos) {
                parsed_key.erase(pos_dd, 1);
                auto pos_sc = prefix_key.rfind('.');
                if(pos_sc != std::string::npos) prefix_key.erase(pos_sc, prefix_key.size() - pos_sc);
            } else
                break;
        }

        if(prefix_key.empty()) {
            prefix_key = std::string(parsed_key);
        } else
            prefix_key = fmt::format("{}.{}", prefix_key, parsed_key);
        // If the element does not exist we insert it here
        if(prefix_key.empty()) throw std::runtime_error(fmt::format("Invalid key: {}", prefix_key));
        auto  sp       = tid::internal::split<std::deque<std::string_view>>(prefix_key, ".");
        auto  result   = tid::internal::tid_db.insert(std::make_pair(sp[0], tid::ur(sp[0])));
        auto &ur_found = result.first->second;
        if(result.second and l != level::parent) ur_found.set_level(l);
        sp.pop_front();
        return get(sp, l, ur_found);
    }

    ur &get_unscoped(std::string_view key, level l) {
        return get(key, l, "");

        //        if(key.empty()) throw std::runtime_error(fmt::format("Invalid key: {}", key));
        //        auto result = internal::tid_db.insert(std::make_pair(key, tid::ur(key)));
        //        auto & ur = result.first->second;
        //        if(result.second and l != level::parent)  ur.set_level(l); // Set the level on creation only
        //        return ur;
    }

    token tic_token(std::string_view key, level l) { return tid::get(key, l).tic_token(); }

    token tic_scope(std::string_view key, level l) { return tid::get(key, l).tic_token(key); }

    void tic(std::string_view key, level l) { get(key, l).tic(); }

    void toc(std::string_view key, level l) { get(key, l).toc(); }

    void reset(std::string_view expr) {
        for(auto &[key, ur] : tid::internal::tid_db) {
            if(std::regex_match(std::string(key), std::regex(std::string(expr)))) ur.reset();
        }
    }

    void reset(const std::vector<std::string> &excl) {
        for(auto &[key, ur] : tid::internal::tid_db) {
            for(const auto &e : excl) {
                if(key == e) continue;
            }
            ur.reset();
        }
    }

    void set_prefix(std::string_view prefix) { internal::ur_prefix = std::string(prefix); }

    std::vector<internal::ur_ref_t> get_tree(const tid::ur &u, std::string_view prefix, level l) {
        std::string key;
        if(prefix.empty())
            key = u.get_label();
        else
            key = fmt::format("{}.{}", prefix, u.get_label());

        std::vector<internal::ur_ref_t> tree = {internal::ur_ref_t{key, u, 0.0, 1.0}};
        for(const auto &un : u.ur_under) {
            tree.front().sum += un.second->get_time(); // Add up times under
            for(const auto &t : get_tree(*un.second, key, l)) {
                //                if(un.second->get_time() == 0) {
                //                    // If the intermediate node did not measure time, add the times under it instead
                //                    tree.front().sum += t.sum;
                //                }
                tree.emplace_back(t);
            }
        }

        // Sort the tree
        if(tree.size() > 1)
            std::sort(std::next(tree.begin()), tree.end(), [](const internal::ur_ref_t &t1, const internal::ur_ref_t &t2) -> bool { return t1.key < t2.key; });

        // Prune the tree based on level
        tree.erase(std::remove_if(tree.begin(), tree.end(), [&l](auto &t) -> bool { return t->get_level() > l; }), tree.end());

        // Find the longest key in the tree
        auto   max_it       = std::max_element(tree.begin(), tree.end(), [](const auto &a, const auto &b) -> bool { return a.key.size() < b.key.size(); });
        size_t max_key_size = 0;
        if(max_it != tree.end()) max_key_size = max_it->key.size();

        // Calculate the fractions and set the maximum key size
        for(auto &t : tree) {
            t.tree_max_key_size = max_key_size;
            if(tree.front()->get_time() == 0) break;
            if(&t == &tree.front()) continue;
            auto t_parent = tree.front()->get_time() == 0.0 ? tree.front().sum : tree.front()->get_time();
            t.frac        = t->get_time() / t_parent;
        }

        return tree;
    }

    std::vector<internal::ur_ref_t> get_tree(std::string_view prefix, level l) {
        std::vector<internal::ur_ref_t> tree;
        for(const auto &[key, u] : tid::internal::tid_db) {
            if(key == prefix or prefix.empty()) {
                auto t = get_tree(u, "", l);
                tree.insert(tree.end(), t.begin(), t.end());
            }
        }
        // Find the longest key in the tree
        auto max_it = std::max_element(tree.begin(), tree.end(), [](const auto &a, const auto &b) -> bool { return a.key.size() < b.key.size(); });
        if(max_it != tree.end()) {
            // Set the max key size
            for(auto &t : tree) t.tree_max_key_size = max_it->key.size();
        }

        return tree;
    }

    std::vector<internal::ur_ref_t> search(const tid::ur &u, std::string_view match) {
        std::vector<internal::ur_ref_t> matches;
        for(const auto &t : get_tree(u, "", level::pedant)) {
            if(t.key.find(match) != std::string_view::npos) matches.push_back(t);
        }
        return matches;
    }

    std::vector<internal::ur_ref_t> search(std::string_view match) {
        std::vector<internal::ur_ref_t> matches;
        for(const auto &t : get_tree("", level::pedant)) {
            if(t.key.find(match) != std::string_view::npos) matches.push_back(t);
        }
        return matches;
    }

    void print_tree(const tid::ur &u, std::string_view prefix, level l) {
        for(const auto &t : tid::get_tree(u, prefix)) {
            if(t->get_level() <= l) fmt::print("{}\n", t.str());
        }
    }
    void print_tree(std::string_view prefix, level l) {
        for(const auto &t : tid::get_tree(prefix)) {
            if(t->get_level() <= l) fmt::print("{}\n", t.str());
        }
    }
}