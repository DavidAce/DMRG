#include "h5dbg.h"
#include "tid/tid.h"

std::string tools::h5dbg::get_hid_string_details(const hid_t &id) {
    auto        t_scope  = tid::tic_scope(__FUNCTION__);
    auto        hid_type = H5Iget_type(id);
    std::string msg      = fmt::format("id {}: ", id);
    switch(hid_type) {
        case H5I_type_t::H5I_DATASET: {
            msg.append("[DATASET]: ");
            break;
        }
        case H5I_type_t::H5I_ATTR: {
            msg.append("[ATTRIBUTE]: ");
            break;
        }
        case H5I_type_t::H5I_DATATYPE: {
            msg.append("[DATATYPE]: ");
            break;
        }
        case H5I_type_t::H5I_DATASPACE: {
            h5pp::hid::h5s hid = id;
            msg.append("[DATASPACE]: ");
            break;
        }
        case H5I_type_t::H5I_GROUP: {
            msg.append("[GROUP]: ");
            break;
        }
        case H5I_type_t::H5I_FILE: {
            msg.append("[FILE]: ");
            break;
        }
        case H5I_type_t::H5I_BADID: {
            msg.append("[BADID] ");
            break;
        }
        default: msg.append("[UNKNOWN]");
    }
    msg.append(h5pp::format("ref count {}: ", H5Iget_ref(id)));
    msg.append(h5pp::hdf5::getName(id));

    return msg;
}

void tools::h5dbg::print_dangling_ids(const h5pp::File &file, std::string_view func, int line, unsigned obj) {
    // Check that there are no open HDF5 handles
    h5pp::hid::h5f handle       = file.openFileHandle();
    auto           ssize_objids = H5Fget_obj_count(handle, obj);
    if(ssize_objids > 1) {
        tools::logger::log->info("File [{}] has {} open ids at {} line {}", file.getFileName(), ssize_objids, func, line);
        auto               size_objids = static_cast<size_t>(ssize_objids);
        std::vector<hid_t> objids(size_objids);
        H5Fget_obj_ids(handle, obj, size_objids, objids.data());
        for(const auto &id : objids) {
            //            auto rc = H5Iget_ref(id);
            //            if( rc > 1){
            tools::logger::log->info(" -- {}", tools::h5dbg::get_hid_string_details(id));
            //            }
        }
    } else if(ssize_objids < 0) {
        tools::logger::log->info("File [{}] failed to count ids at {} line {}: {}", file.getFileName(), func, line, ssize_objids);
    }
    // Check that there are no errors lurking in the HDF5 error-stack
    auto num_errors = H5Eget_num(H5E_DEFAULT);
    if(num_errors > 0) {
        H5Eprint(H5E_DEFAULT, stderr);
        throw std::runtime_error(fmt::format("Error when treating file [{}]", file.getFileName()));
    }
}
void tools::h5dbg::assert_no_dangling_ids(const h5pp::File &file, std::string_view func, int line, unsigned obj) {
    // Check that there are no open HDF5 handles
    h5pp::hid::h5f handle       = file.openFileHandle();
    auto           ssize_objids = H5Fget_obj_count(handle, obj);
    if(ssize_objids > 1) {
        tools::logger::log->info("File [{}] has {} open ids at {} line {}", file.getFileName(), ssize_objids, func, line);
        auto               size_objids = static_cast<size_t>(ssize_objids);
        std::vector<hid_t> objids(size_objids);
        H5Fget_obj_ids(handle, obj, size_objids, objids.data());
        for(auto &id : objids) tools::logger::log->error(" -- {}", tools::h5dbg::get_hid_string_details(id));
        if(ssize_objids > 1) throw std::runtime_error("Found dangling object id's (ref count > 1)");
    } else if(ssize_objids < 0) {
        tools::logger::log->info("File [{}] failed to count ids at {} line {}: {}", file.getFileName(), func, line, ssize_objids);
    }
    // Check that there are no errors lurking in the HDF5 error-stack
    auto num_errors = H5Eget_num(H5E_DEFAULT);
    if(num_errors > 0) {
        H5Eprint(H5E_DEFAULT, stderr);
        throw std::runtime_error(fmt::format("Error when treating file [{}]", file.getFileName()));
    }
}