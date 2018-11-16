//
// Created by david on 2018-05-24.
//
#include <iostream>
#include <iterator>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <IO/class_hdf5_file.h>
#include <IO/class_hdf5_table_buffer2.h>
#include <algorithms/table_types.h>


/*! \brief Prints the content of a vector nicely */
template<typename T, size_t N>
std::ostream &operator<<(std::ostream &out, const std::array<T,N> &v) {
    if (!v.empty()) {
        out << "[ ";
        std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, " "));
        out << "]";
    }
    return out;
}

template<typename table_type>
class_hdf5_table<table_type>::class_hdf5_table(std::shared_ptr<class_hdf5_file> hdf5_out_,
                 std::string group_name_,
                 std::string table_name_,
                 bool mpi_on_):
        hdf5_file(std::move(hdf5_out_)),
        group_name(group_name_),
        table_name(table_name_),
        mpi_on(mpi_on_)
{
    table_entries      = std::make_unique<table_type>();
    table_path = group_name + "/" + table_name;
    initialize_table();
}


template<typename table_type>
class_hdf5_table<table_type>::~class_hdf5_table(){
    if (hdf5_file){
        write_buffer_to_file();
    }else if (!buffer_is_empty){
        std::cerr << "Warning: Output data_struct has not been saved to file, yet it is being discarded!\n" << std::endl;
    }
}


template<typename table_type>
void class_hdf5_table<table_type>::initialize_table(){
    if (table_entries->buffer.empty() and !table_is_ready) {
        hsize_t NRECORDS = table_entries->buffer.size();
        hdf5_file->create_group_link(group_name);
        H5TBmake_table(table_name.c_str(), hdf5_file->file, table_path.c_str(), table_entries->meta.NFIELDS, NRECORDS,
                       table_entries->meta.dst_size, table_entries->meta.field_names.data(), table_entries->meta.dst_offsets.data(), table_entries->meta.field_types.data(),
                       table_entries->meta.chunk_size, table_entries->meta.fill_data, table_entries->meta.compress, table_entries->buffer.data());
        table_is_ready = true;
        H5Fflush(hdf5_file->file,H5F_SCOPE_GLOBAL);

    }
}

template<typename table_type>
void class_hdf5_table<table_type>:: write_buffer_to_file() {
    if (!table_entries->buffer.empty() and table_is_ready) {
        hsize_t NRECORDS = table_entries->buffer.size();
        hdf5_file->create_group_link(group_name);
        H5TBappend_records(hdf5_file->file, table_path.c_str(), NRECORDS, table_entries->meta.dst_size,
                           table_entries->meta.dst_offsets.data(), table_entries->meta.dst_sizes.data(), table_entries->buffer.data());

        table_entries->buffer.clear();
        recorded_elements += NRECORDS;
        H5Fflush(hdf5_file->file,H5F_SCOPE_GLOBAL);
    }
    buffer_is_empty = true;
}



//Explicit instantiations

template class class_hdf5_table<class_table_dmrg>;
template class class_hdf5_table<class_table_tebd>;
template class class_hdf5_table<class_table_finite_chain>;
template class class_hdf5_table<class_table_profiling>;
