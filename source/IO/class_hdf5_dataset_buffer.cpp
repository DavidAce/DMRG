//
// Created by david on 2017-12-02.
//
#include "class_hdf5_file.h"
#include "class_hdf5_dataset_buffer.h"

template<typename DataType, typename AttrType,typename IterType>
std::string class_hdf5_dataset_buffer<DataType, AttrType, IterType>::left_pad(const std::string &temp){
    std::ostringstream ostr;
    ostr << std::setfill('0') << std::setw(pad_to_width) << temp;
    return ostr.str();
}


template<typename DataType, typename AttrType,typename IterType>
class_hdf5_dataset_buffer<DataType, AttrType, IterType>::class_hdf5_dataset_buffer(std::shared_ptr<class_hdf5_file> hdf5_out_,
                                                                         const std::string &group_name_,
                                                                         const IterType &iteration_,
                                                                         const std::string &dataset_name_)
: hdf5_out(std::move(hdf5_out_)),
group_name     (group_name_),
iteration      (iteration_),
dataset_name   (dataset_name_)
{
    this->reserve(max_elements);
}

template<typename DataType, typename AttrType,typename IterType>
class_hdf5_dataset_buffer<DataType, AttrType, IterType>::class_hdf5_dataset_buffer(std::shared_ptr<class_hdf5_file> hdf5_out_,
                                                                                   const std::string &group_name_,
                                                                                   const IterType &iteration_,
                                                                                   const std::string &dataset_name_,
                                                                                   const AttrType &attribute_,
                                                                                   const std::string &attribute_name_)
        :
        hdf5_out(std::move(hdf5_out_)),
        group_name(group_name_),
        iteration(iteration_),
        dataset_name(dataset_name_),
        attribute(attribute_),
        attribute_name(attribute_name_),
        attribute_set(true)
{
    this->reserve(max_elements);
}

template<typename DataType, typename AttrType,typename IterType>
class_hdf5_dataset_buffer<DataType, AttrType, IterType>::class_hdf5_dataset_buffer(const std::string &group_name_,
                                                                         const IterType &iteration_,
                                                                         const std::string &dataset_name_)
:class_hdf5_dataset_buffer(nullptr, group_name_, iteration_, dataset_name_)
{}
template<typename DataType, typename AttrType,typename IterType>
class_hdf5_dataset_buffer<DataType, AttrType, IterType>::class_hdf5_dataset_buffer(const std::string &group_name_,
                                                                         const IterType &iteration_,
                                                                         const std::string &dataset_name_,
                                                                         const AttrType &attribute_,
                                                                         const std::string &attribute_name_)
        :
        class_hdf5_dataset_buffer(nullptr, group_name_, iteration_, dataset_name_, attribute_, attribute_name_)
{}


template<typename DataType, typename AttrType,typename IterType>
void class_hdf5_dataset_buffer<DataType, AttrType, IterType>::hdf5_set_group_name    (const std::string &group_name_)   {
    group_name     = group_name_;
}


template<typename DataType, typename AttrType,typename IterType>
void class_hdf5_dataset_buffer<DataType, AttrType, IterType>::hdf5_set_dataset_name  (const std::string &dataset_name_) {
    dataset_name   = dataset_name_;
}

template<typename DataType, typename AttrType,typename IterType>
void class_hdf5_dataset_buffer<DataType, AttrType, IterType>::hdf5_set_iteration     (const IterType &iteration_)             {
    iteration      = iteration_;
}

template<typename DataType, typename AttrType,typename IterType>
void class_hdf5_dataset_buffer<DataType, AttrType, IterType>::hdf5_set_width     (const int &pad_to_width_)             {
    pad_to_width      = pad_to_width_;
}

template<typename DataType, typename AttrType,typename IterType>
void class_hdf5_dataset_buffer<DataType, AttrType, IterType>::hdf5_set_all(const std::string &group_name_, const IterType iteration_, const std::string &dataset_name_){
    group_name     = group_name_;
    iteration      = iteration_;
    dataset_name   = dataset_name_;
}

template<typename DataType, typename AttrType,typename IterType>
void class_hdf5_dataset_buffer<DataType, AttrType, IterType>::hdf5_set_attribute(const AttrType &attribute_, const std::string &attribute_name_){
    attribute_name  = attribute_name_;
    attribute       = attribute_;
    attribute_set   = true;
}

template<typename DataType, typename AttrType,typename IterType>
void class_hdf5_dataset_buffer<DataType, AttrType, IterType>::hdf5_refresh_relative_name(){
    dataset_relative_name = "/" + group_name + "_" + left_pad(std::to_string(iteration)) + "/" + dataset_name ;
}


template<typename DataType, typename AttrType,typename IterType>
void class_hdf5_dataset_buffer<DataType, AttrType, IterType>::write_buffer_to_file() {
    if (!this->empty()) {
        hdf5_refresh_relative_name();
        hdf5_out->write_dataset(static_cast<std::vector<DataType>>(*this), dataset_relative_name);

        if (attribute_set) {
            hdf5_out->write_attribute_to_dataset(dataset_relative_name, attribute, attribute_name);
        }
        this->clear();
    }
    data_has_been_written_to_file = true;
    H5Fflush(hdf5_out->file,H5F_SCOPE_GLOBAL);
}


template class class_hdf5_dataset_buffer<double>;
template class class_hdf5_dataset_buffer<int>;
template class class_hdf5_dataset_buffer<long>;