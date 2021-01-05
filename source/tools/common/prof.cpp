//
// Created by david on 2019-06-08.
//

#include <fstream>
#include <sstream>

#include <config/enums.h>
#include <config/nmspc_settings.h>
#include <tools/common/log.h>
#include <tools/common/prof.h>
#include <general/class_tic_toc.h>

namespace tools::common::profile::internal {
    // Implement a custom ordered map
    // In the map, we want keep the order in which they were appended. This behavior does not exist in stl ordered maps.
    template<typename KeyT, typename ValT> using iterator = typename insert_ordered_map<KeyT,ValT>::iterator;
    template<typename KeyT, typename ValT> using const_iterator = typename insert_ordered_map<KeyT,ValT>::const_iterator;

    template<typename KeyT, typename ValT>
    ValT & insert_ordered_map<KeyT,ValT>::operator[](const KeyT &key) {
        auto it = find(key);
        if(it == data.end()){
            if constexpr(std::is_convertible_v<KeyT,std::string>)
                throw std::runtime_error(fmt::format("Invalid key: {}",key));
            else throw std::runtime_error("Invalid key");
        }
        return it->second;
    }

    template<typename KeyT, typename ValT>
    void insert_ordered_map<KeyT,ValT>::append(const KeyT & key, ValT val){
        auto it = find(key);
        if(it == data.end())
            data.emplace_back(std::make_pair(key, std::move(val)));
    }

    template<typename KeyT, typename ValT> typename insert_ordered_map<KeyT,ValT>::iterator insert_ordered_map<KeyT,ValT>::begin() { return data.begin(); }
    template<typename KeyT, typename ValT> typename insert_ordered_map<KeyT,ValT>::iterator insert_ordered_map<KeyT,ValT>::end() { return data.end(); }
    template<typename KeyT, typename ValT> typename insert_ordered_map<KeyT,ValT>::const_iterator insert_ordered_map<KeyT,ValT>::begin() const { return data.begin(); }
    template<typename KeyT, typename ValT> typename insert_ordered_map<KeyT,ValT>::const_iterator insert_ordered_map<KeyT,ValT>::end() const { return data.end(); }
    template<typename KeyT, typename ValT> typename insert_ordered_map<KeyT,ValT>::iterator insert_ordered_map<KeyT,ValT>::find(const KeyT & key){
        return std::find_if(data.begin(), data.end(), [&key](const auto &element) { return element.first == key; });
    }
    template<typename KeyT, typename ValT> typename insert_ordered_map<KeyT,ValT>::const_iterator insert_ordered_map<KeyT,ValT>::find(const KeyT & key) const {
        return std::find_if(data.begin(), data.end(), [&key](const auto &element) { return element.first == key; });
    }

    template class insert_ordered_map<std::string, std::unique_ptr<class_tic_toc>>;
    template class insert_ordered_map<AlgorithmType, MapTicTocUnique>;
}




tools::common::profile::internal::MapTicTocUnique &tools::common::profile::get_default_prof() {
    if(not internal::default_algo_type) throw std::runtime_error("No default algorithm type has been set for profiling");
    return prof[internal::default_algo_type.value()];
}
void          tools::common::profile::set_default_prof(AlgorithmType algo_type) { internal::default_algo_type = algo_type; }
AlgorithmType tools::common::profile::get_current_algo_type() {
    if(not internal::default_algo_type) throw std::runtime_error("No default algorithm type has been set");
    return internal::default_algo_type.value();
}

void tools::common::profile::print_profiling_all() { tools::common::profile::print_profiling(std::nullopt); }

void tools::common::profile::print_profiling(std::optional<AlgorithmType> algo_type) {
    if(t_tot == nullptr) return;
    if(settings::profiling::on) {
        static double last_print_time = 0;
        if(std::abs(t_tot->get_measured_time() - last_print_time) < 5.0) return; // Do not print if there's already been a print within the last 5 seconds
        auto t_tot_percent = 100.0 / std::max(1.0, t_tot->get_measured_time());
        //        auto t_sim_percent = 100.0 / std::max(1.0, t_sim->get_measured_time());

        /* clang-format off */
        tools::log->info("{:<6}{:<30}{:>10.3f} s"                       ,"DMRG ", t_tot->get_name(), t_tot->get_measured_time());

        if(algo_type)
            for(const auto &[t_key, t_val] : prof[algo_type.value()]) {
                if(t_val->get_measured_time() == 0) continue;
                tools::log->info("{:<8}  {:<30}{:>10.3f} s ({:<3.2f} % of total)", enum2str(algo_type.value()), t_val->get_name(), t_val->get_measured_time(),
                                 t_val->get_measured_time() * t_tot_percent);
            }
        else
            for(const auto &[a_key, a_val] : prof)
                for(const auto &[t_key, t_val] : a_val) {
                    if(t_val->get_measured_time() == 0) continue;
                    tools::log->info("{:<8}  {:<30}{:>10.3f} s ({:<3.2f} % of total)", enum2str(a_key), t_val->get_name(), t_val->get_measured_time(),
                                     t_val->get_measured_time() * t_tot_percent);
                }

        last_print_time = t_tot->get_measured_time();
    }
}

void tools::common::profile::print_profiling_laps(std::optional<AlgorithmType> algo_type) {
    if(t_tot == nullptr) return;
    if(settings::profiling::on) {
        auto t_tot_percent = 100.0 / std::max(1.0, t_tot->get_lap());
        //        auto t_sim_percent = 100.0 / std::max(1.0, t_sim->get_measured_time());

        /* clang-format off */
        tools::log->info("{:<6}{:<30}{:>10.3f} s"                       ,"DMRG ", t_tot->get_name(), t_tot->get_lap());
        t_tot->start_lap();
        if(algo_type)
            for(const auto &[t_key, t_val] : prof[algo_type.value()]) {
                if(t_val->get_lap() == 0) continue;
                tools::log->info("{:<6}{:<30}{:>10.3f} s ({:<3.2f} % of total)", enum2str(algo_type.value()), t_val->get_name(), t_val->get_lap(),
                                 t_val->get_lap() * t_tot_percent);
                t_val->start_lap();
            }
        else
            for(const auto &[a_key, a_val] : prof)
                for(const auto &[t_key, t_val] : a_val) {
                    if(t_val->get_lap() == 0) continue;
                    tools::log->info("{:<6}{:<30}{:>10.3f} s ({:<3.2f} % of total)", enum2str(a_key), t_val->get_name(), t_val->get_lap(),
                                     t_val->get_lap() * t_tot_percent);
                    t_val->start_lap();
                }
    }
}

void tools::common::profile::init_profiling() {
    if(t_tot != nullptr) return;
    /* clang-format off */
    t_tot             = std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "+ Total Time              ");

    prof.append(AlgorithmType::xDMRG, internal::MapTicTocUnique());
    prof.append(AlgorithmType::fDMRG, internal::MapTicTocUnique());
    prof.append(AlgorithmType::iDMRG, internal::MapTicTocUnique());
    prof.append(AlgorithmType::iTEBD, internal::MapTicTocUnique());
    prof.append(AlgorithmType::fLBIT, internal::MapTicTocUnique());
    prof.append(AlgorithmType::ANY, internal::MapTicTocUnique());

    prof[AlgorithmType::xDMRG].append("t_pre",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Preprocess             "));
    prof[AlgorithmType::xDMRG].append("t_rnd",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Randomize state        "));
    prof[AlgorithmType::xDMRG].append("t_pos",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Postprocess            "));
    prof[AlgorithmType::xDMRG].append("t_sim",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|+ Simulation             "));
    prof[AlgorithmType::xDMRG].append("t_con",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Convergence checks    "));
    prof[AlgorithmType::xDMRG].append("t_eig",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Eig. decomp.          "));
    prof[AlgorithmType::xDMRG].append("t_svd",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Svd. decomp.          "));
    prof[AlgorithmType::xDMRG].append("t_env",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Environment upd.      "));
    prof[AlgorithmType::xDMRG].append("t_ent",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Entanglement entropy  "));
    prof[AlgorithmType::xDMRG].append("t_spn",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Spin components       "));
    prof[AlgorithmType::xDMRG].append("t_ene",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Energy                "));
    prof[AlgorithmType::xDMRG].append("t_var",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Variance              "));
    prof[AlgorithmType::xDMRG].append("t_prj",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Projections           "));
    prof[AlgorithmType::xDMRG].append("t_dbg",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Debugging checks      "));
    prof[AlgorithmType::xDMRG].append("t_out",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Print status          "));
    prof[AlgorithmType::xDMRG].append("t_hdf",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- h5pp storage          "));
    prof[AlgorithmType::xDMRG].append("t_mps",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Multisite-MPS         "));
    prof[AlgorithmType::xDMRG].append("t_mpo",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Multisite-MPO         "));
    prof[AlgorithmType::xDMRG].append("t_opt",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |+ Optimization          "));
    prof[AlgorithmType::xDMRG].append("t_opt_dir",       std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " ||+ Direct               "));
    prof[AlgorithmType::xDMRG].append("t_opt_dir_bfgs",  std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |||+ L-BFGS              "));
    prof[AlgorithmType::xDMRG].append("t_opt_dir_step",  std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " || |- step               "));
    prof[AlgorithmType::xDMRG].append("t_opt_dir_vH2",   std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " || |- vH2                "));
    prof[AlgorithmType::xDMRG].append("t_opt_dir_vH2v",  std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " || |- vH2v               "));
    prof[AlgorithmType::xDMRG].append("t_opt_dir_vH",    std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " || |- vH                 "));
    prof[AlgorithmType::xDMRG].append("t_opt_dir_vHv",   std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " || |- vHv                "));
    prof[AlgorithmType::xDMRG].append("t_opt_sub",       std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " ||+ Subspace             "));
    prof[AlgorithmType::xDMRG].append("t_opt_sub_ham",   std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " | |- Hamiltonian Matrix  "));
    prof[AlgorithmType::xDMRG].append("t_opt_sub_hsq",   std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " | |- Hamiltonian Matrix² "));
    prof[AlgorithmType::xDMRG].append("t_opt_sub_lu",    std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " | |- LU decomposition    "));
    prof[AlgorithmType::xDMRG].append("t_opt_sub_eig",   std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " | |- Eigenvalue decomp   "));
    prof[AlgorithmType::xDMRG].append("t_opt_sub_bfgs",  std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " | |+ L-BFGS              "));
    prof[AlgorithmType::xDMRG].append("t_opt_sub_step",  std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |  |- step               "));
    prof[AlgorithmType::xDMRG].append("t_opt_sub_vH2",   std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |  |- vH2                "));
    prof[AlgorithmType::xDMRG].append("t_opt_sub_vH2v",  std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |  |- vH2v               "));
    prof[AlgorithmType::xDMRG].append("t_opt_sub_vH",    std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |  |- vH                 "));
    prof[AlgorithmType::xDMRG].append("t_opt_sub_vHv",   std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |  |- vHv                "));

    prof[AlgorithmType::fDMRG].append("t_pre",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Preprocess             "));
    prof[AlgorithmType::fDMRG].append("t_rnd",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Randomize state        "));
    prof[AlgorithmType::fDMRG].append("t_pos",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Postprocess            "));
    prof[AlgorithmType::fDMRG].append("t_sim",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|+ Simulation             "));
    prof[AlgorithmType::fDMRG].append("t_con",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Convergence checks    "));
    prof[AlgorithmType::fDMRG].append("t_eig",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Eig. decomp.          "));
    prof[AlgorithmType::fDMRG].append("t_svd",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Svd. decomp.          "));
    prof[AlgorithmType::fDMRG].append("t_env",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Environment upd.      "));
    prof[AlgorithmType::fDMRG].append("t_ent",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Entanglement entropy  "));
    prof[AlgorithmType::fDMRG].append("t_spn",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Spin components       "));
    prof[AlgorithmType::fDMRG].append("t_ene",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Energy                "));
    prof[AlgorithmType::fDMRG].append("t_var",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Variance              "));
    prof[AlgorithmType::fDMRG].append("t_prj",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Projections           "));
    prof[AlgorithmType::fDMRG].append("t_dbg",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Debugging checks      "));
    prof[AlgorithmType::fDMRG].append("t_out",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Print status          "));
    prof[AlgorithmType::fDMRG].append("t_hdf",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- h5pp storage          "));
    prof[AlgorithmType::fDMRG].append("t_mps",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Multisite-MPS         "));
    prof[AlgorithmType::fDMRG].append("t_mpo",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Multisite-MPO         "));

    prof[AlgorithmType::fLBIT].append("t_pre",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Preprocess             "));
    prof[AlgorithmType::fLBIT].append("t_rnd",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Randomize state        "));
    prof[AlgorithmType::fLBIT].append("t_pos",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Postprocess            "));
    prof[AlgorithmType::fLBIT].append("t_sim",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|+ Simulation             "));
    prof[AlgorithmType::fLBIT].append("t_con",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Convergence checks    "));
    prof[AlgorithmType::fLBIT].append("t_eig",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Eig. decomp.          "));
    prof[AlgorithmType::fLBIT].append("t_evo",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Time evolution        "));
    prof[AlgorithmType::fLBIT].append("t_map",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Map lbit<-->real      "));
    prof[AlgorithmType::fLBIT].append("t_svd",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Svd. decomp.          "));
    prof[AlgorithmType::fLBIT].append("t_env",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Environment upd.      "));
    prof[AlgorithmType::fLBIT].append("t_ent",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Entanglement entropy  "));
    prof[AlgorithmType::fLBIT].append("t_spn",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Spin components       "));
    prof[AlgorithmType::fLBIT].append("t_ene",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Energy                "));
    prof[AlgorithmType::fLBIT].append("t_var",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Variance              "));
    prof[AlgorithmType::fLBIT].append("t_prj",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Projections           "));
    prof[AlgorithmType::fLBIT].append("t_dbg",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Debugging checks      "));
    prof[AlgorithmType::fLBIT].append("t_out",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Print status          "));
    prof[AlgorithmType::fLBIT].append("t_hdf",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- h5pp storage          "));
    prof[AlgorithmType::fLBIT].append("t_mps",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Multisite-MPS         "));
    prof[AlgorithmType::fLBIT].append("t_mpo",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Multisite-MPO         "));

    prof[AlgorithmType::iDMRG].append("t_pre",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Preprocess             "));
    prof[AlgorithmType::iDMRG].append("t_rnd",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Randomize state        "));
    prof[AlgorithmType::iDMRG].append("t_pos",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Postprocess            "));
    prof[AlgorithmType::iDMRG].append("t_sim",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|+ Simulation             "));
    prof[AlgorithmType::iDMRG].append("t_con",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Convergence checks    "));
    prof[AlgorithmType::iDMRG].append("t_eig",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Eig. decomp.          "));
    prof[AlgorithmType::iDMRG].append("t_svd",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Svd. decomp.          "));
    prof[AlgorithmType::iDMRG].append("t_env",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Environment upd.      "));
    prof[AlgorithmType::iDMRG].append("t_ent",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Entanglement entropy  "));
    prof[AlgorithmType::iDMRG].append("t_ene",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Energy                "));
    prof[AlgorithmType::iDMRG].append("t_var",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Variance              "));
    prof[AlgorithmType::iDMRG].append("t_dbg",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Debugging checks      "));
    prof[AlgorithmType::iDMRG].append("t_out",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Print status          "));
    prof[AlgorithmType::iDMRG].append("t_hdf",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- h5pp storage          "));
    prof[AlgorithmType::iDMRG].append("t_ene_ham",       std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Energy (HAM)          "));
    prof[AlgorithmType::iDMRG].append("t_ene_mom",       std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Energy (MOM)          "));
    prof[AlgorithmType::iDMRG].append("t_var_ham",       std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Variance (HAM)        "));
    prof[AlgorithmType::iDMRG].append("t_var_mom",       std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Variance (MOM)        "));

    prof[AlgorithmType::iTEBD].append("t_pre",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Preprocessing          "));
    prof[AlgorithmType::iTEBD].append("t_pos",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|- Postprocessing         "));
    prof[AlgorithmType::iTEBD].append("t_sim",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, "|+ Simulation             "));
    prof[AlgorithmType::iTEBD].append("t_con",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Convergence checks    "));
    prof[AlgorithmType::iTEBD].append("t_svd",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Svd. decomp.          "));
    prof[AlgorithmType::iTEBD].append("t_evo",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Time evolution        "));
    prof[AlgorithmType::iTEBD].append("t_ent",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Entanglement entropy  "));
    prof[AlgorithmType::iTEBD].append("t_dbg",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Debugging checks      "));
    prof[AlgorithmType::iTEBD].append("t_out",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Print status          "));
    prof[AlgorithmType::iTEBD].append("t_hdf",           std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- h5pp storage          "));
    prof[AlgorithmType::iTEBD].append("t_ene_ham",       std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Energy (HAM)          "));
    prof[AlgorithmType::iTEBD].append("t_ene_mom",       std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Energy (MOM)          "));
    prof[AlgorithmType::iTEBD].append("t_var_ham",       std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Variance (HAM)        "));
    prof[AlgorithmType::iTEBD].append("t_var_mom",       std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Variance (MOM)        "));

    prof[AlgorithmType::ANY].append("t_gate_move",       std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Gate: Move            "));
    prof[AlgorithmType::ANY].append("t_gate_apply",      std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Gate: Apply           "));
    prof[AlgorithmType::ANY].append("t_gate_merge",      std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Gate: Merge           "));
    prof[AlgorithmType::ANY].append("t_gate_return",     std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Gate: Return          "));
    prof[AlgorithmType::ANY].append("t_merge_split",     std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Merge: Split          "));
    prof[AlgorithmType::ANY].append("t_merge_merge",     std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Merge: Merge          "));
    prof[AlgorithmType::ANY].append("t_split_svdm",      std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Split: SVD main       "));
    prof[AlgorithmType::ANY].append("t_split_svda",      std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Split: SVD left       "));
    prof[AlgorithmType::ANY].append("t_split_svdb",      std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Split: SVD right      "));
    prof[AlgorithmType::ANY].append("t_splitA_svd",      std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- SplitA:SVD            "));
    prof[AlgorithmType::ANY].append("t_splitB_svd",      std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- SplitB:SVD            "));
    prof[AlgorithmType::ANY].append("t_write_h5pp",      std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Write h5pp            "));
    prof[AlgorithmType::ANY].append("t_map_norm"  ,      std::make_unique<class_tic_toc>(settings::profiling::on, settings::profiling::precision, " |- Map normalize         "));

    /* clang-format on */


}


void tools::common::profile::reset_profiling(std::optional<AlgorithmType> algo_type, const std::vector<std::string> &excl) {
    if(t_tot == nullptr) throw std::runtime_error("Profiling timers have not been initialized");
    // Do not reset total time!
    tools::log->trace("Resetting profiling timers");
    if(algo_type)
        for(auto &[t_key, t_val] : prof[algo_type.value()]) {
            for(const auto &e : excl)
                if(t_key == e) continue;
            t_val->reset();
        }
    else
        for(auto &[a_key, a_val] : prof)
            for(auto &[t_key, t_val] : a_val) {
                for(const auto &e : excl)
                    if(t_key == e) continue;
                t_val->reset();
            }
}

//void tools::common::profile::reset_for_run_algorithm(std::optional<AlgorithmType> algo_type, const std::vector<std::string> &excl) {
//    if(t_tot == nullptr) throw std::runtime_error("Profiling timers have not been initialized");
//    if(not algo_type) algo_type = tools::common::profile::get_current_algo_type();
//    // Do not reset total time!
//    tools::log->trace("Resetting profiling timers for new {} run", enum2str(algo_type.value()));
//    if(algo_type)
//        for(auto &[t_key, t_val] : prof[algo_type.value()]) {
//            for(const auto &e : excl)
//                if(t_key == e) continue;
//            t_val->reset();
//        }
//    else
//        for(auto &[a_key, a_val] : prof)
//            for(auto &[t_key, t_val] : a_val) {
//                for(const auto &e : excl)
//                    if(t_key == e) continue;
//                t_val->reset();
//            }
//}

double tools::common::profile::mem_usage_in_mb(std::string_view name) {
    std::ifstream filestream("/proc/self/status");
    std::string   line;
    while(std::getline(filestream, line)) {
        std::istringstream is_line(line);
        std::string        key;
        if(std::getline(is_line, key, ':')) {
            if(key == name) {
                std::string value_str;
                if(std::getline(is_line, value_str)) {
                    // Filter non-digit characters
                    value_str.erase(std::remove_if(value_str.begin(), value_str.end(), [](auto const &c) -> bool { return not std::isdigit(c); }),
                                    value_str.end());
                    // Extract the number
                    long long value = 0;
                    try {
                        std::string::size_type sz; // alias of size_t
                        value = std::stoll(value_str, &sz);
                    } catch(const std::exception &ex) {
                        tools::log->error("Could not read mem usage from /proc/self/status: Failed to parse string [{}]: {}", value_str, ex.what());
                    }
                    // Now we have the value in kb
                    return static_cast<double>(value) / 1024.0;
                }
            }
        }
    }
    return -1.0;
}

double tools::common::profile::mem_rss_in_mb() { return mem_usage_in_mb("VmRSS"); }
double tools::common::profile::mem_hwm_in_mb() { return mem_usage_in_mb("VmHWM"); }
double tools::common::profile::mem_vm_in_mb() { return mem_usage_in_mb("VmPeak"); }

void tools::common::profile::print_mem_usage() {
    tools::log->info("{:<30}{:>10.1f} MB", "Memory RSS", mem_rss_in_mb());
    tools::log->info("{:<30}{:>10.1f} MB", "Memory Peak", mem_hwm_in_mb());
    tools::log->info("{:<30}{:>10.1f} MB", "Memory Vm", mem_vm_in_mb());
}
