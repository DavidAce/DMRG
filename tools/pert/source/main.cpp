#include <config/settings.h>
#include <config/threading.h>
#include <h5pp/h5pp.h>
#include <math/tenx/threads.h>
#include <omp.h>
#include <tensors/model/ModelFinite.h>
#include <tensors/site/mpo/MpoSite.h>
#include <tensors/state/StateFinite.h>
#include <tools/common/h5.h>
#include <tools/common/log.h>
#include <tools/finite/h5.h>
#include <tools/finite/measure.h>

struct PertData {
    long   seed            = 0;
    long   model_size      = 0;
    double delta           = 0.0;
    double g               = 0;
    double g_pert          = 0;
    double energy          = 0;
    double energy_pert     = 0;
    double energy_change   = 0;
    double variance        = 0;
    double variance_pert   = 0;
    double variance_change = 0;
};
struct spins {
    double x,y,z;
};

int main() {
    tools::log                       = tools::Logger::setLogger("pert", 2, true);
    settings::threading::num_threads = 8; // omp_get_max_threads();
    settings::configure_threads();
    double         g_perturb = 1e-5;
    h5pp::fs::path h5dir     = "/mnt/WDB-AN1500/mbl_transition/xdmrg3-letsgo/output/L12/g0.000";
    auto           h5tgt     = h5pp::File("../../../../output/pert-L12-g1e-5.h5", h5pp::FileAccess::REPLACE);

    h5pp::hid::h5t h5_type = H5Tcreate(H5T_COMPOUND, sizeof(PertData));
    H5Tinsert(h5_type, "seed", HOFFSET(PertData, seed), H5T_NATIVE_LONG);
    H5Tinsert(h5_type, "model_size", HOFFSET(PertData, model_size), H5T_NATIVE_LONG);
    H5Tinsert(h5_type, "delta", HOFFSET(PertData, delta), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5_type, "g", HOFFSET(PertData, g), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5_type, "g_pert", HOFFSET(PertData, g_pert), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5_type, "energy", HOFFSET(PertData, energy), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5_type, "energy_pert", HOFFSET(PertData, energy_pert), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5_type, "energy_change", HOFFSET(PertData, energy_change), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5_type, "variance", HOFFSET(PertData, variance), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5_type, "variance_pert", HOFFSET(PertData, variance_pert), H5T_NATIVE_DOUBLE);
    H5Tinsert(h5_type, "variance_change", HOFFSET(PertData, variance_change), H5T_NATIVE_DOUBLE);

    // Collect and sort all the files in h5dir
    using h5iter = h5pp::fs::recursive_directory_iterator;
    long counter = 0;
    for(const auto &file : h5iter(h5dir)) {
        if(not file.is_regular_file()) continue;
        if(file.path().extension() != ".h5") continue;

        // tools::log->info("file: {}\n", file.path().string());
        auto h5src = h5pp::File(file, h5pp::FileAccess::READONLY);

        auto var_old = h5src.readTableField<double>("xDMRG/state_emid/measurements", "energy_variance", h5pp::TableSelection::LAST);
        if(var_old > 1e-12) continue;
        auto ene_old = h5src.readTableField<double>("xDMRG/state_emid/measurements", "energy", h5pp::TableSelection::LAST);
        auto spinZ   = h5src.readTableField<spins>("xDMRG/state_emid/measurements", "spin_global", h5pp::TableSelection::LAST).z;

        auto seed        = h5src.readDataset<long>("common/seed");
        auto model_size  = h5src.readAttribute<long>("xDMRG/model/hamiltonian", "model_size");
        auto model_type  = h5src.readAttribute<ModelType>("xDMRG/model/hamiltonian", "model_type");
        auto position    = h5src.readAttribute<long>("xDMRG/state_emid/mps", "position");
        auto model_g     = h5src.readTableField<double>("xDMRG/model/hamiltonian", "g", h5pp::TableSelection::FIRST);
        auto model_delta = h5src.readTableField<double>("xDMRG/model/hamiltonian", "delta", h5pp::TableSelection::FIRST);

        settings::model::model_size            = model_size;
        settings::model::model_type            = model_type;
        settings::model::ising_majorana::g     = model_g;
        settings::model::ising_majorana::delta = model_delta;

        tools::log->set_level(spdlog::level::warn);
        auto model = ModelFinite(ModelType::ising_majorana, model_size);
        tools::finite::h5::load::model(h5src, AlgorithmType::xDMRG, model);

        auto mpsinfo = tools::finite::h5::load::MpsInfo();
        auto state   = StateFinite(AlgorithmType::xDMRG, model_size, position);
        tools::finite::h5::load::state(h5src, "xDMRG/state_emid", state, mpsinfo);
        tools::log->set_level(spdlog::level::info);

        // auto expH_old    = tools::finite::measure::expectation_value(state, state, mposEne_old);

        // auto mposEne_old = model.get_all_mpo_tensors(MposWithEdges::ON);
        // auto mposVar_old = model.get_compressed_mpos_squared(MposWithEdges::ON);

        // auto expH2_old   = tools::finite::measure::expectation_value(state, state, mposVar_old);
        // auto var_old     = std::abs(std::abs(expH2_old) - std::abs(expH_old * expH_old));
        // auto var_file    = h5src.readTableField<double>("xDMRG/state_emid/measurements", "energy_variance", h5pp::TableSelection::LAST);

        // tools::log->info("Old H {:<20.16f} H² {:.3e} | var {:.3e}", std::real(ene_old), std::real(var_old), std::real(var_old) - std::real(ene_old *
        // ene_old));
        for(auto &mpo : model.MPO) { mpo->set_parameter("g", g_perturb); }
        model.build_mpo();
        auto mposEne_pert = model.get_all_mpo_tensors(MposWithEdges::ON);
        auto ene_pert     = tools::finite::measure::expectation_value(state, state, mposEne_pert);

        for(const auto &mpo : model.MPO) mpo->set_energy_shift_mpo(cplx(std::real(ene_pert) / static_cast<double>(model_size), 0.0));
        for(auto &mpo : model.MPO) { mpo->set_parity_shift_mpo_squared(spinZ > 0 ? 1 : -1, "z"); }

        model.clear_cache();
        model.build_mpo();
        model.build_mpo_squared();

        auto mposVar_pert = model.get_compressed_mpos_squared(MposWithEdges::ON);
        auto var_pert     = std::abs(tools::finite::measure::expectation_value(state, state, mposVar_pert)); // Has reduced energy so this is immediately the variance

        auto ene_chn = std::real(ene_pert - ene_old);
        auto var_chn = std::abs(std::abs(var_pert) - std::abs(var_old));

        tools::log->info("seed {} | change energy {:.3e} | variance {:.3e} -> {:.3e} change {:.3e}", seed, ene_chn, var_old, var_pert, var_chn);

        auto tablename = fmt::format("L{}-d{:.2f}-gpert{:.1e}", model_size, model_delta, g_perturb);
        if(not h5tgt.linkExists(tablename)) { h5tgt.createTable(h5_type, tablename, "PertData", 100, 3); }
        auto pert            = PertData();
        pert.seed            = seed;
        pert.model_size      = model_size;
        pert.delta           = model_delta;
        pert.g               = model_g;   // The perturbation
        pert.g_pert          = g_perturb; // The perturbation
        pert.energy          = std::real(ene_old);
        pert.energy_pert     = std::real(ene_pert);
        pert.energy_change   = ene_chn;
        pert.variance        = std::real(var_old);
        pert.variance_pert   = std::real(var_pert);
        pert.variance_change = var_chn;
        h5tgt.appendTableRecords(pert, tablename);
        //
        counter++;
        // if(counter > 100) exit(0);
    }

    return 0;
}