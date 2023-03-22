//
// Created by david on 2021-10-12.
//

#include "parse.h"
#include "enums.h"
#include "loader.h"
#include "settings.h"
#include "tools/common/log.h"
#include <CLI/CLI.hpp>
#include <h5pp/h5pp.h>

template<typename T>
std::string filename_append_number(const std::string &filename, const T number) {
    if constexpr(std::is_signed_v<T>)
        if(number < 0) return filename;
    if constexpr(std::is_unsigned_v<T>)
        if(number == std::numeric_limits<T>::max()) return filename;
    // Append the seed_model to the output filename
    h5pp::fs::path oldFileName = filename;
    h5pp::fs::path newFileName = filename;
    if(oldFileName.stem().string().find(std::to_string(number)) != std::string::npos) return filename;
    newFileName.replace_filename(fmt::format("{}_{}{}", oldFileName.stem().string(), number, oldFileName.extension().string()));
    tools::log->info("Appended number [{}] to filename: [{}]", number, newFileName.string());
    return newFileName.string();
}

template<>
constexpr auto sv2enum<h5pp::LogLevel>(std::string_view item) {
    if(item == "trace") return h5pp::LogLevel::trace;
    if(item == "debug")
        return h5pp::LogLevel::debug;
    else
        return h5pp::LogLevel::info;
}

template<>
constexpr auto sv2enum<spdlog::level::level_enum>(std::string_view item) {
    if(item == "trace") return spdlog::level::level_enum::trace;
    if(item == "debug")
        return spdlog::level::level_enum::debug;
    else
        return spdlog::level::level_enum::info;
}

// MWE: https://godbolt.org/z/jddxod53d
int settings::parse(int argc, char **argv) {
    using namespace settings;
    using namespace h5pp;
    using namespace spdlog;

    auto s2e_log     = mapStr2Enum<spdlog::level::level_enum>("trace", "debug", "info");
    auto s2e_logh5pp = mapStr2Enum<h5pp::LogLevel>("trace", "debug", "info");
    auto s2e_model   = ModelType_s2e;
    int  dummy       = 0;

    auto preload = [&argc, &argv, &s2e_log]() -> int {
        CLI::App pre;
        pre.get_formatter()->column_width(90);
        pre.option_defaults()->always_capture_default();
        pre.allow_extras(true);
        pre.set_help_flag("--help-preload", "Help for preloading configuration");
        /* clang-format off */
        pre.add_option("-c,--config"                       , input::config_filename , "Path to a .cfg or .h5 file from a previous simulation");
        pre.add_option("-v,--log,--verbosity,--loglevel"   , console::loglevel      , "Log level of DMRG++")->transform(CLI::CheckedTransformer(s2e_log, CLI::ignore_case))->type_name("ENUM");
        pre.add_option("--timestamp"                       , console::timestamp     , "Log timestamp");
        /* clang-format on */
        pre.parse(argc, argv);
        tools::log = tools::Logger::setLogger("DMRG++ config", settings::console::loglevel, settings::console::timestamp);
        tools::log->info("Preloading {}", input::config_filename);
        //  Try loading the given config file.
        //  Note that there is a default "input/input.config" if none was given
        Loader dmrg_config(settings::input::config_filename);
        if(dmrg_config.file_exists) {
            dmrg_config.load();
            settings::load(dmrg_config);
        } else if(pre.get_option("--config")->empty()) {
            tools::log->warn("The default config file does not exist: {}", input::config_filename);
        } else
            throw except::runtime_error("Could not find config file: {}", settings::input::config_filename); // Invalid file
        return 0;
    };
    preload();

    CLI::App app;
    app.description("DMRG++: An MPS-based algorithm to calculate 1D quantum-states");
    app.get_formatter()->column_width(90);
    app.option_defaults()->always_capture_default();
    app.allow_extras(false);
    /* clang-format off */
    app.add_flag("--help-preload"                      , "Print help related to preloading configuration");
    app.add_option("-c,--config"                       , input::config_filename         , "Path to a .cfg or .h5 file from a previous simulation");
    app.add_option("-m,--model"                        , model::model_type              , "Select the Hamiltonian")->transform(CLI::CheckedTransformer(s2e_model, CLI::ignore_case));
    app.add_option("-b,--bitfield"                     , input::bitfield                , "Integer whose bitfield sets the initial product state. Negative is unused");
    app.add_option("-o,--outfile"                      , storage::output_filepath       , "Path to the output file. The seed number gets appended by default (see -x)");
    app.add_option("-s,--seed"                         , input::seed                    , "Positive number seeds the random number generator");
    app.add_option("-t,--threads"                      , threading::num_threads         , "Total number of threads (omp + std threads). Use env OMP_NUM_THREADS to control omp.");
    app.add_flag("--show-threads"                      , threading::show_threads        , "Show information about threading and exit immediately");
    app.add_flag  ("--append-seed, !--no-append-seed"  , storage::output_append_seed    , "Append seed to the output filename")->default_val(true);
    app.add_option("-z,--compression"                  , storage::compression_level     , "Compression level of h5pp")->check(CLI::Range(0,9));
    app.add_flag  ("-r,--resume"                                                        , "Resume simulation from last iteration");
    app.add_option("--resume-iter"                     , storage::file_resume_iter      , "Resume from iteration");
    app.add_option("--resume-name"                     , storage::file_resume_name      , "Resume from state matching this name");
    app.add_option("-v,--log,--verbosity,--loglevel"   , console::loglevel              , "Log level of DMRG++")->transform(CLI::CheckedTransformer(s2e_log, CLI::ignore_case))->type_name("ENUM");
    app.add_option("-V,--logh5pp"                      , console::logh5pp               , "Log level of h5pp")->transform(CLI::CheckedTransformer(s2e_logh5pp, CLI::ignore_case))->type_name("ENUM");
    app.add_option("--timestamp"                       , console::timestamp             , "Log timestamp");
    app.add_option("--dummyrange"                      , dummy                          , "Dummy")->check(CLI::Range(0,3));
    /* clang-format on */

    app.parse(argc, argv);

    if(app.count("--resume") > 0 or app.count("--resume-iter") > 0 or app.count("--resume-name") > 0) {
        tools::log->info("Resuming from iter {}", storage::file_resume_iter);
        settings::storage::file_collision_policy = FileCollisionPolicy::RESUME;
    }

    //    for(const auto &res : app.get_options()) fmt::print("{:<32} = {}\n", res->get_name(), res->results());

    // Generate the correct output filename based on given seeds
    if(storage::output_append_seed) {
        settings::storage::output_filepath = filename_append_number(settings::storage::output_filepath, settings::input::seed);
        settings::storage::output_filepath = filename_append_number(settings::storage::output_filepath, settings::input::bitfield);
    }

    return 0;
}