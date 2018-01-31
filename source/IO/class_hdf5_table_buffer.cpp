//
// Created by david on 2018-01-26.
//

#include <IO/class_hdf5_table_buffer.h>
#include <IO/class_hdf5_file.h>
#include <H5Cpp.h>
#include <hdf5_hl.h>



class_hdf5_table_buffer::class_hdf5_table_buffer(std::shared_ptr<class_hdf5_file> hdf5_out_,
                                                 const std::string group_name_,
                                                 const std::string table_name_)
        :
        hdf5_out(std::move(hdf5_out_)),
        group_name(group_name_),
        table_name(table_name_)
{
    this->reserve(max_elements);
};

class_hdf5_table_buffer::class_hdf5_table_buffer(std::nullptr_t nullp,
                                                 const std::string group_name_,
                                                 const std::string table_name_)
        :
        group_name(group_name_),
        table_name(table_name_)
{
    this->reserve(max_elements);
};


class_hdf5_table_buffer::class_hdf5_table_buffer(const std::string group_name_,
                                                 const std::string table_name_)
        : class_hdf5_table_buffer(nullptr, group_name_,table_name_)
{

}




//template<typename TableType>
void class_hdf5_table_buffer::write_buffer_to_file(class_hdf5_file &hdf5_out) {
    if (!this->empty()) {
        hsize_t NRECORDS = this->size();
        if(table_is_empty) {

            H5TBmake_table("Table Title", hdf5_out.file, table_name.c_str(), meta.NFIELDS, NRECORDS,
                           meta.dst_size, meta.field_names, meta.dst_offset, meta.field_type,
                           meta.chunk_size, meta.fill_data, meta.compress, this->data());
            table_is_empty = false;
        }
        else{
            H5TBappend_records(hdf5_out.file, table_name.c_str(), NRECORDS, meta.dst_size, meta.dst_offset, meta.dst_sizes,
                               this->data());
        }
    }
    this->clear();
    buffer_is_empty = true;
}
