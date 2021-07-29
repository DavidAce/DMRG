#include <algorithm>
#include <algorithms/AlgorithmStatus.h>
#include <chrono>
#include <config/settings.h>
#include <cstdlib>
#include <fstream>
#include <h5pp/h5pp.h>
#include <io/filesystem.h>
#include <io/fmt.h>
#include <regex>
#include <tid/tid.h>
#include <tools/common/h5.h>
#include <tools/common/log.h>

std::string get_dirname() { return "DMRG." + std::string(getenv("USER")); }

std::string replace_substr(std::string text, const std::string &search, const std::string &replace) {
    size_t pos = 0;
    while((pos = text.find(search, pos)) != std::string::npos) {
        text.replace(pos, search.length(), replace);
        pos += replace.length();
    }
    return text;
}

std::string tools::common::h5::tmp::internal::get_tmp_dir() {
    if(fs::exists(settings::output::temp_dir))
        return settings::output::temp_dir;
    else if(fs::exists("/dev/shm"))
        return "/dev/shm";
    else if(fs::exists("/scratch/local"))
        return "/scratch/local";
    else
        return fs::temp_directory_path();
}

const tools::common::h5::tmp::internal::pathpair &tools::common::h5::tmp::internal::register_paths(const std::string &filepath) {
    if(not fs::path(filepath).has_filename()) throw std::runtime_error(fmt::format("Given output file path has no filename: [{}]", filepath));
    std::string filename = fs::path(filepath).filename();
    if(internal::file_register.find(filename) != internal::file_register.end()) return internal::file_register[filename];
    std::string original_path         = fs::absolute(filepath);
    std::string temporary_path        = fs::path(get_tmp_dir()) / fs::path(get_dirname()) / filename;
    internal::file_register[filename] = {original_path, temporary_path};
    return internal::file_register[filename];
}

const tools::common::h5::tmp::internal::pathpair &tools::common::h5::tmp::internal::get_paths(const std::string &filepath) {
    if(not fs::path(filepath).has_filename()) throw std::runtime_error(fmt::format("Given output file path has no filename: [{}]", filepath));
    std::string filename = fs::path(filepath).filename();
    return internal::file_register[filename];
}

const std::string &tools::common::h5::tmp::get_temporary_filepath(const std::string &filepath) { return internal::get_paths(filepath).temporary_path; }
const std::string &tools::common::h5::tmp::get_original_filepath(const std::string &filepath) { return internal::get_paths(filepath).original_path; }

void tools::common::h5::tmp::register_new_file(const std::string &filepath) { internal::register_paths(filepath); }

// std::string tools::common::h5::tmp::set_tmp_prefix(const std::string &filepath) {
//    return init::register_paths(filepath).temporary_path;

//
//
//    fs::path temp_path = fs::path(get_tmp_dir()) / fs::path(get_dirname());
//    if (output_filepath.find(temp_path.string()) != std::string::npos) {
//        tools::log->debug("Already temporary path [{}]", output_filepath);
//        return output_filepath;
//    } else {
//        fs::path concat_filepath;
//        if(fs::path(output_filepath).is_absolute()){
//            concat_filepath = temp_path.string() + output_filepath;
//        }else{
//            concat_filepath = temp_path / output_filepath;
//        }
//        // The filename may be relative or full, just make sure to replace .. with __
//        fs::path new_filepath = replace_substr(concat_filepath, "..", "__");
//        tools::log->debug("Temp path       [{}]", temp_path.string());
//        tools::log->debug("output_filepath [{}]", output_filepath);
//        tools::log->debug("concat_filepath [{}]", concat_filepath.string());
//        tools::log->debug("new_filepath    [{}]", new_filepath.string());
//        tools::log->debug("Set temporary path [{}] -> [{}]", output_filepath, new_filepath.string());
//        return new_filepath;
//    }
//}

// std::string tools::common::h5::tmp::unset_tmp_prefix(const std::string &filepath) {
//    return init::register_paths(filepath).original_path;
//    fs::path temp_path = fs::path(get_tmp_dir()) / fs::path(get_dirname());
//    std::string::size_type pos = output_filepath.find(temp_path.string());
//    if (pos != std::string::npos){
//        std::string new_filepath = output_filepath;
//        new_filepath.erase(pos, temp_path.string().length());
//        new_filepath = fs::current_path() / std::regex_replace(new_filepath, std::regex("__"), "..");
//        tools::log->debug("Unset temporary path [{}] -> [{}]", output_filepath, new_filepath);
//        return new_filepath;
//    }else{
//        tools::log->debug("Already final path [{}]", output_filepath);
//        return output_filepath;
//    }

//}

void tools::common::h5::tmp::create_directory(const std::string &path) {
    if(path.empty()) return;
    if(tools::log == nullptr) {
        if(spdlog::get("DMRG") == nullptr)
            tools::log = spdlog::default_logger();
        else
            tools::log = spdlog::get("DMRG");
    }

    fs::path dir = fs::absolute(path);
    if(dir.has_filename() and dir.has_extension()) dir = dir.parent_path();

    try {
        if(fs::create_directories(dir)) {
            tools::log->info("Created directory: {}", dir.string());
        } else {
            tools::log->info("Directory already exists: {}", dir.string());
        }
    } catch(std::exception &ex) { throw std::runtime_error("Failed to create directory: " + std::string(ex.what())); }
}

void tools::common::h5::tmp::copy_from_tmp(const AlgorithmStatus &status, const h5pp::File &h5ppfile, StorageReason storage_reason,
                                           std::optional<CopyPolicy> copy_policy) {
    if(not settings::output::use_temp_dir) return;
    if(not copy_policy) return copy_from_tmp(status, h5ppfile, storage_reason, CopyPolicy::TRY);
    if(copy_policy == CopyPolicy::OFF) return;

    // Check if we already copied the file this iteration and step
    static std::unordered_map<std::string, std::pair<uint64_t, uint64_t>> save_log;
    auto                                                                  save_point = std::make_pair(status.iter, status.step);

    if(copy_policy == CopyPolicy::TRY) {
        if(save_log[h5ppfile.getFilePath()] == save_point) return;
        switch(storage_reason) {
            case StorageReason::SAVEPOINT:
            case StorageReason::CHECKPOINT:
                if(status.iter % settings::output::copy_from_temp_freq != 0) return; // Check that we write according to the frequency given
            case StorageReason::FINISHED:
            case StorageReason::CHI_UPDATE:
            case StorageReason::PROJ_STATE:
            case StorageReason::INIT_STATE:
            case StorageReason::EMIN_STATE:
            case StorageReason::EMAX_STATE:
            case StorageReason::MODEL: break;
        }
        tools::common::h5::tmp::copy_from_tmp(h5ppfile.getFilePath());
    } else if(copy_policy == CopyPolicy::FORCE)
        tools::common::h5::tmp::copy_from_tmp(h5ppfile.getFilePath());

    save_log[h5ppfile.getFilePath()] = save_point;
}

void tools::common::h5::tmp::copy_from_tmp(const std::string &filepath) {
    if(filepath.empty()) return;
    if(not fs::exists(get_temporary_filepath(filepath)) and internal::file_register.find(filepath) != internal::file_register.end()) {
        tools::log->debug("Temporary file is already deleted: [{}]", get_temporary_filepath(filepath));
        return;
    }

    copy_file(get_temporary_filepath(filepath), get_original_filepath(filepath));
}

void tools::common::h5::tmp::copy_into_tmp(const std::string &filepath) {
    if(filepath.empty()) return;
    copy_file(get_original_filepath(filepath), get_temporary_filepath(filepath));
}

void tools::common::h5::tmp::copy_file(const std::string &src, const std::string &tgt) {
    if(src == tgt) {
        tools::log->info("Skipping file copy. Identical paths: [{}] == [{}]", src, tgt);
        return;
    }
    fs::path target_path = tgt;
    fs::path source_path = src;
    if(not fs::exists(source_path)) {
        tools::log->warn("Could not copy: source file does not exist [{}]", src);
        return;
    }
    auto t_copy = tid::tic_token("copy");
    if(not fs::exists(target_path.parent_path())) { tools::common::h5::tmp::create_directory(target_path); }
    if(fs::exists(target_path)) {
        std::ifstream                          target_stream(target_path.string(), std::ios_base::binary);
        std::ifstream                          source_stream(source_path.string(), std::ios_base::binary);
        typedef std::istreambuf_iterator<char> isbuf_it;
        if(std::equal(isbuf_it(target_stream.rdbuf()), isbuf_it(), isbuf_it(source_stream.rdbuf()), isbuf_it())) {
            tools::log->debug("Source [{}] identical to target [{}] ... Skipping copy", src, tgt);
            return;
        }
    }
    tools::log->debug("Copying file: {} -> {}", src, tgt);
    h5pp::hdf5::copyFile(source_path, target_path, h5pp::FilePermission::REPLACE);
    //    fs::copy(source_path, target_path, fs::copy_options::overwrite_existing);
}

void tools::common::h5::tmp::remove_from_tmp(const std::string &filepath) {
    if(filepath.empty()) {
        tools::log->trace("Nothing to delete\n");
        return;
    }
    const auto &[orig, temp] = internal::get_paths(filepath);
    if(temp == orig) {
        tools::log->debug("Final path and temporary paths are identical. Skipping removal: [{}] == [{}]", orig, temp);
        return;
    }
    if(temp.empty()) {
        tools::log->debug("No temp file path defined - nothing to delete");
    } else {
        if(fs::exists(temp)) {
            tools::log->debug("Deleting temporary file: {}", temp);
            fs::remove(temp);
        } else {
            tools::log->debug("File does not exist [{}]. Nothing to delete", temp);
        }
    }
}
